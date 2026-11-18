import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from variant_cohort.cli import run
from variant_cohort.phase import stitch

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'src'
BASE = 'ACGT' * 60
HEADER = '##fileformat=VCFv4.2\n'


def fasta(seq=BASE, name='chr1'):
    return f'>{name}\n{seq}\n'


def data(*fragments, seq=BASE):
    out = []
    for fragment_id, samples, rows in fragments:
        lines = [HEADER.rstrip('\n'),
                 '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(samples)]
        for pos, ref, alts, cells in rows:
            lines.append('\t'.join(['chr1', str(pos), '.', ref, ','.join(alts),
                                    '.', 'PASS', '.', 'GT:PS'] + list(cells)))
        out.append({'id': fragment_id, 'vcf': '\n'.join(lines) + '\n'})
    return {'fasta': fasta(seq), 'fragments': out}


class ChainAndOrientation(unittest.TestCase):
    def test_three_block_chain_with_one_swap(self):
        a = ('a', ['S'], [(1, 'A', ['C'], ['1|0:1']), (5, 'A', ['G'], ['1|0:1'])])
        b = ('b', ['S'], [(5, 'A', ['G'], ['0|1:1']), (9, 'A', ['T'], ['1|0:1'])])
        c = ('c', ['S'], [(9, 'A', ['T'], ['1|0:1']), (13, 'A', ['C'], ['0|1:1'])])
        result = stitch(data(a, b, c))
        self.assertEqual(len(result['blocks']), 1)
        block = result['blocks'][0]
        self.assertEqual(block['anchor'], {'fragment_id': 'a', 'phase_set': 1})
        members = {m['fragment_id']: m['flipped'] for m in block['members']}
        self.assertEqual(members, {'a': False, 'b': True, 'c': True})
        copies = {(allele['position'], allele['alt']): allele['copies'] for allele in block['alleles']}
        self.assertEqual(copies, {(1, 'C'): [1, 0], (5, 'G'): [1, 0],
                                  (9, 'T'): [0, 1], (13, 'C'): [1, 0]})

    def test_two_blocks_keep_direction(self):
        a = ('a', ['S'], [(1, 'A', ['C'], ['1|0:1']), (5, 'A', ['G'], ['0|1:1'])])
        b = ('b', ['S'], [(5, 'A', ['G'], ['0|1:1'])])
        result = stitch(data(a, b))
        self.assertEqual(len(result['blocks']), 1)
        self.assertEqual([m['flipped'] for m in result['blocks'][0]['members']], [False, False])

    def test_pairwise_unique_compatible_links_form_contradictory_cycle(self):
        a = ('a', ['S'], [(1, 'A', ['C'], ['1|0:9']), (3, 'G', ['T'], ['1|0:9'])])
        b = ('b', ['S'], [(1, 'A', ['C'], ['0|1:9']), (2, 'C', ['A'], ['0|1:9'])])
        c = ('c', ['S'], [(2, 'C', ['A'], ['0|1:9']), (3, 'G', ['T'], ['1|0:9'])])
        with self.assertRaisesRegex(ValueError, 'contradictory phase connections'):
            stitch(data(a, b, c))

    def test_consistently_parity_compatible_cycle_merges(self):
        a = ('a', ['S'], [(1, 'A', ['C'], ['1|0:9']), (3, 'G', ['T'], ['1|0:9'])])
        b = ('b', ['S'], [(1, 'A', ['C'], ['0|1:9']), (2, 'C', ['A'], ['1|0:9'])])
        c = ('c', ['S'], [(2, 'C', ['A'], ['0|1:9']), (3, 'G', ['T'], ['1|0:9'])])
        result = stitch(data(a, b, c))
        self.assertEqual(len(result['blocks']), 1)
        self.assertEqual([m['flipped'] for m in result['blocks'][0]['members']], [False, True, False])

    def test_shared_identity_conflict_rejected(self):
        a = ('a', ['S'], [(5, 'A', ['G'], ['1|0:1'])])
        b = ('b', ['S'], [(5, 'A', ['G'], ['0|0:1'])])
        with self.assertRaisesRegex(ValueError, 'incompatible phase blocks'):
            stitch(data(a, b))


class BlockIdentity(unittest.TestCase):
    def test_phase_set_groups_all_rows_of_one_fragment_sample_contig(self):
        # PS spans two rows; a second PS is a separate local block.
        a = ('a', ['S'], [(1, 'A', ['C'], ['1|0:3']), (5, 'A', ['G'], ['1|0:3']),
                          (9, 'A', ['T'], ['1|0:4'])])
        result = stitch(data(a))
        self.assertEqual(len(result['blocks']), 2)
        self.assertEqual([block['anchor']['phase_set'] for block in result['blocks']], [3, 4])
        self.assertEqual([len(block['alleles']) for block in result['blocks']], [2, 1])

    def test_same_phase_set_on_different_contigs_is_distinct(self):
        payload = {'fasta': '>c1\nAAAA\n>c2\nAAAA\n',
                   'fragments': [{'id': 'a', 'vcf': HEADER +
                                  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n'
                                  'c1\t1\t.\tA\tC\t.\tPASS\t.\tGT:PS\t1|0:5\n'
                                  'c2\t1\t.\tA\tG\t.\tPASS\t.\tGT:PS\t1|0:5\n'}]}
        result = stitch(payload)
        self.assertEqual(len(result['blocks']), 2)
        self.assertEqual([(b['contig'], b['anchor']['phase_set']) for b in result['blocks']],
                         [('c1', 5), ('c2', 5)])

    def test_anchor_is_smallest_fragment_id_then_numeric_phase_set(self):
        # phase_set 9 < 10 numerically while '10' < '9' lexicographically.
        a = ('b', ['S'], [(5, 'A', ['G'], ['1|0:9'])])
        b = ('a', ['S'], [(5, 'A', ['G'], ['1|0:10'])])
        result = stitch(data(a, b))
        self.assertEqual(result['blocks'][0]['anchor'], {'fragment_id': 'a', 'phase_set': 10})
        self.assertEqual([(m['fragment_id'], m['phase_set'], m['flipped'])
                          for m in result['blocks'][0]['members']],
                         [('a', 10, False), ('b', 9, False)])

    def test_conflict_detected_after_flip_without_order_dependence(self):
        a = ('a', ['S'], [(1, 'A', ['C'], ['1|0:1']), (3, 'G', ['T'], ['1|0:1'])])
        b = ('b', ['S'], [(1, 'A', ['C'], ['0|1:1']), (3, 'G', ['T'], ['1|0:1'])])
        with self.assertRaisesRegex(ValueError, 'incompatible phase blocks'):
            stitch(data(a, b))
        with self.assertRaisesRegex(ValueError, 'incompatible phase blocks'):
            stitch(data(b, a))


class SeparationAndMissingness(unittest.TestCase):
    def test_same_phase_set_across_fragments_without_overlap_stays_separate(self):
        a = ('a', ['S'], [(1, 'A', ['C'], ['1|0:7'])])
        b = ('b', ['S'], [(9, 'A', ['T'], ['0|1:7'])])
        result = stitch(data(a, b))
        self.assertEqual(len(result['blocks']), 2)
        self.assertEqual([block['anchor'] for block in result['blocks']],
                         [{'fragment_id': 'a', 'phase_set': 7}, {'fragment_id': 'b', 'phase_set': 7}])
        self.assertTrue(all(len(block['members']) == 1 and not block['members'][0]['flipped']
                            for block in result['blocks']))

    def test_partial_missing_complementary_but_both_directions_compatible(self):
        a = ('a', ['S'], [(5, 'A', ['G'], ['1|.:3'])])
        b = ('b', ['S'], [(5, 'A', ['G'], ['.|1:3'])])
        result = stitch(data(a, b))
        self.assertEqual(len(result['blocks']), 2)

    def test_all_missing_isolated_block_is_preserved(self):
        a = ('a', ['S'], [(5, 'A', ['G'], ['.|.:11'])])
        result = stitch(data(a))
        self.assertEqual(len(result['blocks']), 1)
        self.assertEqual(result['blocks'][0]['alleles'], [
            {'position': 5, 'ref': 'A', 'alt': 'G', 'copies': [None, None]}])
        self.assertEqual(result['blocks'][0]['members'],
                         [{'fragment_id': 'a', 'phase_set': 11, 'flipped': False}])

    def test_unknowns_do_not_override_known_values_in_merge(self):
        a = ('a', ['S'], [(1, 'A', ['C'], ['1|0:1']), (5, 'A', ['G'], ['1|.:1'])])
        b = ('b', ['S'], [(5, 'A', ['G'], ['1|0:1']), (9, 'A', ['T'], ['0|1:1'])])
        result = stitch(data(a, b))
        block = result['blocks'][0]
        self.assertEqual(len(block['members']), 2)
        copies = {(allele['position'], allele['alt']): allele['copies'] for allele in block['alleles']}
        self.assertEqual(copies[(5, 'G')], [1, 0])

    def test_different_samples_and_contigs_never_join(self):
        a = ('a', ['S1', 'S2'], [(5, 'A', ['G'], ['1|0:1', '0|1:1'])])
        b = ('b', ['S1', 'S3'], [(5, 'A', ['G'], ['1|0:1', '1|0:1'])])
        result = stitch(data(a, b))
        self.assertEqual(result['samples'], ['S1', 'S2', 'S3'])
        self.assertEqual(len(result['blocks']), 3)
        pairs = {(block['sample'], tuple(member['fragment_id'] for member in block['members']))
                 for block in result['blocks']}
        self.assertEqual(pairs, {('S1', ('a', 'b')), ('S2', ('a',)), ('S3', ('b',))})


class RepresentationAlignment(unittest.TestCase):
    def test_multi_alt_reordered_and_trimmable_representations_align(self):
        fragments = [
            ('a', ['S'], [(1, 'AC', ['AT', 'AG'], ['1|2:1']), (9, 'A', ['C'], ['1|0:1'])]),
            ('b', ['S'], [(2, 'C', ['G', 'T'], ['1|2:1'])]),
            ('c', ['S'], [(1, 'ACG', ['ATG'], ['1|0:1']), (9, 'A', ['C'], ['1|0:1'])]),
        ]
        result = stitch(data(*fragments))
        self.assertEqual(len(result['blocks']), 1)
        block = result['blocks'][0]
        copies = {(allele['position'], allele['ref'], allele['alt']): allele['copies']
                  for allele in block['alleles']}
        # Fragment a: T on copy0, G on copy1. Fragment b: G copy0, T copy1 -> flipped.
        # Fragment c: T on copy0 -> aligned with a.
        self.assertEqual(copies[(2, 'C', 'T')], [1, 0])
        self.assertEqual(copies[(2, 'C', 'G')], [0, 1])
        self.assertEqual(copies[(9, 'A', 'C')], [1, 0])
        self.assertEqual({m['fragment_id']: m['flipped'] for m in block['members']},
                         {'a': False, 'b': True, 'c': False})

    def test_multiallelic_other_alt_is_zero_not_unknown(self):
        # 1|2 for T,G: T -> (1,0), G -> (0,1); neither copy is null.
        result = stitch(data(('a', ['S'], [(2, 'C', ['T', 'G'], ['1|2:4'])])))
        copies = {allele['alt']: allele['copies'] for allele in result['blocks'][0]['alleles']}
        self.assertEqual(copies, {'T': [1, 0], 'G': [0, 1]})


def render(fragment, sample_order=None, alt_order=None, row_order=None):
    """Render a semantic fragment, remapping sample/ALT order and GT indexes together.

    fragment = {'id': str, 'samples': [str], 'ps': int,
                'rows': [(pos, ref, [alt names], {sample: (copy0, copy1)})]}
    Cells refer to allele names, so reordering ALT columns cannot desynchronise GT.
    """
    samples = list(sample_order if sample_order is not None else fragment['samples'])
    lines = [HEADER.rstrip('\n'),
             '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(samples)]
    rows = list(fragment['rows'])
    if row_order is not None:
        rows = [rows[i] for i in row_order]
    for pos, ref, alts, calls in rows:
        order = list(alt_order) if alt_order is not None else list(alts)
        if sorted(order) != sorted(alts):
            order = list(alts)
        index = {alt: i + 1 for i, alt in enumerate(order)}
        index[ref] = 0
        cells = []
        for sample in samples:
            copies = calls[sample]
            rendered = []
            for value in copies:
                rendered.append('.' if value is None else str(index[value]))
            cells.append('|'.join(rendered) + ':' + str(fragment['ps']))
        lines.append('\t'.join(['chr1', str(pos), '.', ref, ','.join(order),
                                '.', 'PASS', '.', 'GT:PS'] + cells))
    return {'id': fragment['id'], 'vcf': '\n'.join(lines) + '\n'}


def semantic_fragments():
    return [
        {'id': 'a', 'samples': ['S1', 'S2'], 'ps': 1, 'rows': [
            (1, 'AC', ['AT', 'AG'], {'S1': ('AT', 'AG'), 'S2': ('AT', 'AT')}),
            (9, 'A', ['C'], {'S1': ('C', 'A'), 'S2': (None, 'C')})]},
        {'id': 'b', 'samples': ['S2', 'S1'], 'ps': 1, 'rows': [
            (1, 'AC', ['AG', 'AT'], {'S1': ('AG', 'AT'), 'S2': ('AT', 'AT')}),
            (9, 'A', ['C'], {'S1': ('A', 'C'), 'S2': ('C', 'A')})]},
        {'id': 'c', 'samples': ['S1'], 'ps': 1, 'rows': [
            (1, 'ACG', ['ATG'], {'S1': ('ATG', 'ACG')})]},
    ]


class InputContract(unittest.TestCase):
    def test_reorders_do_not_change_result(self):
        plain = {'fasta': fasta(),
                 'fragments': [render(f) for f in semantic_fragments()]}
        permuted = {'fasta': fasta(),
                    'fragments': [render(f, sample_order=list(reversed(f['samples'])),
                                         alt_order=list(reversed(f['rows'][0][2])),
                                         row_order=list(reversed(range(len(f['rows'])))))
                                  for f in reversed(semantic_fragments())]}
        self.assertEqual(stitch(plain), stitch(permuted))

    def test_permuted_report_is_identical_except_source_hash(self):
        plain = {'fasta': fasta(),
                 'fragments': [render(f) for f in semantic_fragments()]}
        permuted = {'fasta': fasta(),
                    'fragments': [render(f, sample_order=list(reversed(f['samples'])),
                                         alt_order=list(reversed(f['rows'][0][2])),
                                         row_order=list(reversed(range(len(f['rows'])))))
                                  for f in reversed(semantic_fragments())]}
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            reports = []
            for name, payload in [('plain', plain), ('permuted', permuted)]:
                source = directory / f'{name}.json'
                source.write_text(json.dumps(payload))
                reports.append(run(source, True))
        self.assertNotEqual(reports[0]['source_sha256'], reports[1]['source_sha256'])
        for report in reports:
            report.pop('source_sha256')
        self.assertEqual(reports[0], reports[1])

    def test_duplicate_same_fragment_variant_row_rejected(self):
        payload = data(('a', ['S'], [(5, 'A', ['G'], ['1|0:1']), (5, 'A', ['G'], ['1|0:1'])]))
        with self.assertRaisesRegex(ValueError, 'duplicate variant'):
            stitch(payload)

    def test_extra_fragment_field_and_missing_fasta_rejected(self):
        with self.assertRaises(ValueError):
            stitch({'fragments': [{'id': 'a', 'vcf': HEADER}]})
        with self.assertRaises(ValueError):
            stitch({'fasta': fasta(), 'fragments': [{'id': 'a'}]})
        payload = data(('a', ['S'], [(5, 'A', ['G'], ['1|0:1'])]))
        payload['fragments'][0]['extra'] = 1
        with self.assertRaises(ValueError):
            stitch(payload)

    def test_cross_fragment_identity_is_evidence_not_collision(self):
        result = stitch(data(('a', ['S'], [(5, 'A', ['G'], ['1|0:1'])]),
                             ('b', ['S'], [(5, 'A', ['G'], ['1|0:1'])])))
        self.assertEqual(len(result['blocks']), 1)

    def test_within_fragment_normalized_collision_rejected(self):
        with self.assertRaisesRegex(ValueError, 'collision'):
            stitch(data(('a', ['S'], [(1, 'AC', ['AT'], ['1|0:1']),
                                      (2, 'C', ['T'], ['1|0:1'])])))

    def test_duplicate_fragment_id_rejected(self):
        with self.assertRaisesRegex(ValueError, 'duplicate fragment id'):
            stitch(data(('a', ['S'], [(5, 'A', ['G'], ['1|0:1'])]),
                        ('a', ['S'], [(9, 'A', ['G'], ['1|0:1'])])))

    def test_fragment_count_and_id_bounds(self):
        with self.assertRaisesRegex(ValueError, '1\\.\\.32'):
            stitch({'fasta': fasta(), 'fragments': []})
        with self.assertRaisesRegex(ValueError, '1\\.\\.32'):
            stitch({'fasta': fasta(), 'fragments': [
                {'id': f'f{i}', 'vcf': HEADER + '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n'}
                for i in range(33)]})
        with self.assertRaisesRegex(ValueError, 'invalid fragment id'):
            stitch(data(('bad/id', ['S'], [(5, 'A', ['G'], ['1|0:1'])])))

    def test_cell_format_rules(self):
        bad_cells = ['1/0:1', '1:1', '1|0', '1|0:1:2', '1.0:1', '1|0:0', '1|0:01', '1|0:-1']
        for cell in bad_cells:
            with self.subTest(cell=cell), self.assertRaises(ValueError):
                stitch(data(('a', ['S'], [(5, 'A', ['G'], [cell])])))

    def test_gt_out_of_range_rejected(self):
        with self.assertRaisesRegex(ValueError, 'out of range'):
            stitch(data(('a', ['S'], [(5, 'A', ['G'], ['2|0:1'])])))

    def test_format_must_be_gt_ps(self):
        payload = data(('a', ['S'], [(5, 'A', ['G'], ['1|0:1'])]))
        payload['fragments'][0]['vcf'] = payload['fragments'][0]['vcf'].replace('GT:PS', 'GT')
        with self.assertRaisesRegex(ValueError, 'unsupported'):
            stitch(payload)

    def test_header_only_vcf_defines_samples_without_blocks(self):
        vcf = HEADER + '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n'
        payload = {'fasta': fasta(), 'fragments': [{'id': 'a', 'vcf': vcf}]}
        result = stitch(payload)
        self.assertEqual(result['samples'], ['S1', 'S2'])
        self.assertEqual(result['source_records'], 0)
        self.assertEqual(result['blocks'], [])

    def test_sample_union_limit(self):
        first = [f'S{i}' for i in range(150)]
        second = [f'S{i}' for i in range(100, 250)]
        payload = data(('a', first, [(5, 'A', ['G'], ['1|0:1'] * 150)]),
                       ('b', second, [(9, 'A', ['G'], ['1|0:1'] * 150)]))
        with self.assertRaisesRegex(ValueError, 'sample union'):
            stitch(payload)

    def test_total_data_row_limit(self):
        rows = [(pos, 'A', ['C'], ['1|0:1']) for pos in range(1, 158)]
        payload = data(*[(f'f{i:02d}', ['S'], rows) for i in range(32)], seq='A' * 1000)
        with self.assertRaisesRegex(ValueError, 'total data rows'):
            stitch(payload)

    def test_outer_contract_fields(self):
        payload = data(('a', ['S'], [(5, 'A', ['G'], ['1|0:1'])]))
        payload['extra'] = True
        with self.assertRaises(ValueError):
            stitch(payload)
        payload = data(('a', ['S'], [(5, 'A', ['G'], ['1|0:1'])]))
        payload['fragments'][0]['vcf_extra'] = 'x'
        with self.assertRaises(ValueError):
            stitch(payload)

    def test_report_shape_and_sorting(self):
        result = stitch(data(('b', ['S2', 'S1'], [(9, 'A', ['G'], ['1|0:5', '0|1:2'])]),
                             ('a', ['S1'], [(5, 'A', ['G'], ['1|0:3'])])))
        self.assertEqual(set(result), {'samples', 'source_fragments', 'source_records', 'blocks'})
        self.assertEqual(result['samples'], ['S1', 'S2'])
        self.assertEqual([(b['sample'], b['contig'], b['anchor']['fragment_id'], b['anchor']['phase_set'])
                          for b in result['blocks']],
                         [('S1', 'chr1', 'a', 3), ('S1', 'chr1', 'b', 2), ('S2', 'chr1', 'b', 5)])
        for block in result['blocks']:
            self.assertEqual(set(block), {'sample', 'contig', 'anchor', 'members', 'alleles'})
            self.assertEqual(set(block['anchor']), {'fragment_id', 'phase_set'})
            self.assertEqual(block['members'],
                             sorted(block['members'], key=lambda m: (m['fragment_id'], m['phase_set'])))
            self.assertEqual(block['alleles'],
                             sorted(block['alleles'], key=lambda a: (a['position'], a['ref'], a['alt'])))
            for allele in block['alleles']:
                self.assertEqual(set(allele), {'position', 'ref', 'alt', 'copies'})
                self.assertEqual(len(allele['copies']), 2)
                self.assertTrue(all(copy in (0, 1, None) for copy in allele['copies']))


class MergeGuard(unittest.TestCase):
    def test_known_copy_conflict_during_merge_rejected(self):
        # Defensive guard: edge construction already rejects every pair that
        # shares an identity, so this path is only reachable directly.
        from variant_cohort.phase import _merge
        blocks = [
            {'sample': 'S', 'contig': 'chr1', 'fragment_id': 'a', 'phase_set': 1,
             'pairs': {('chr1', 1, 'A', 'C'): (1, 0)}},
            {'sample': 'S', 'contig': 'chr1', 'fragment_id': 'b', 'phase_set': 1,
             'pairs': {('chr1', 1, 'A', 'C'): (0, 0)}},
        ]
        with self.assertRaisesRegex(ValueError, 'conflicting phased allele copies'):
            _merge(blocks, [0, 1], {0: 0, 1: 0})

    def test_merge_fills_unknown_from_known(self):
        from variant_cohort.phase import _merge
        blocks = [
            {'sample': 'S', 'contig': 'chr1', 'fragment_id': 'a', 'phase_set': 1,
             'pairs': {('chr1', 1, 'A', 'C'): (None, 0)}},
            {'sample': 'S', 'contig': 'chr1', 'fragment_id': 'b', 'phase_set': 1,
             'pairs': {('chr1', 1, 'A', 'C'): (1, None)}},
        ]
        merged = _merge(blocks, [0, 1], {0: 0, 1: 0})
        self.assertEqual(merged['alleles'], [
            {'position': 1, 'ref': 'A', 'alt': 'C', 'copies': [1, 0]}])


class CliEndToEnd(unittest.TestCase):
    def _env(self):
        env = dict(os.environ, PYTHONPATH=str(SRC))
        return env

    def _cli(self, source, output, stitch_phase=True):
        command = [sys.executable, '-m', 'variant_cohort.cli',
                   '--source', str(source), '--output', str(output)]
        if stitch_phase:
            command.append('--stitch-phase')
        return subprocess.run(command, capture_output=True, env=self._env(), cwd=ROOT)

    def test_end_to_end_phase_stitch_report(self):
        payload = data(
            ('a', ['S'], [(1, 'A', ['C'], ['1|0:1']), (5, 'A', ['G'], ['1|0:1'])]),
            ('b', ['S'], [(5, 'A', ['G'], ['0|1:1']), (9, 'A', ['T'], ['1|0:1'])]),
        )
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source = directory / 'input.json'
            output = directory / 'report.json'
            source.write_text(json.dumps(payload))
            done = self._cli(source, output)
            self.assertEqual(done.returncode, 0, done.stderr)
            report = json.loads(output.read_text())
            self.assertEqual(set(report), {'schema_version', 'mode', 'source_sha256', 'result'})
            self.assertEqual(report['schema_version'], 2)
            self.assertEqual(report['mode'], 'phase_stitch')
            self.assertEqual(report['result']['source_fragments'], 2)
            self.assertEqual(report['result']['source_records'], 4)
            self.assertEqual(len(report['result']['blocks']), 1)
            self.assertEqual([m['flipped'] for m in report['result']['blocks'][0]['members']],
                             [False, True])

    def test_installed_phase_fixture_matches_library(self):
        report = run(ROOT / 'fixtures' / 'phase.json', True)
        self.assertEqual(report['schema_version'], 2)
        self.assertEqual(report['result'], stitch(json.loads((ROOT / 'fixtures' / 'phase.json').read_text())))

    def test_v2_rejects_v1_shaped_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            output = directory / 'report.json'
            output.write_text('keep')
            done = self._cli(ROOT / 'fixtures' / 'input.json', output)
            self.assertEqual(done.returncode, 2)
            self.assertEqual(output.read_text(), 'keep')

    def test_v1_entry_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            output = directory / 'report.json'
            done = self._cli(ROOT / 'fixtures' / 'input.json', output, stitch_phase=False)
            self.assertEqual(done.returncode, 0, done.stderr)
            report = json.loads(output.read_text())
            self.assertEqual(set(report), {'schema_version', 'source_sha256', 'result'})
            self.assertEqual(report['schema_version'], 1)
            self.assertEqual(set(report['result']), {'samples', 'source_records', 'alleles'})
            self.assertEqual(len(report['result']['alleles']), 2)

    def test_invalid_input_exits_2_and_preserves_output_bytes(self):
        payload = data(('a', ['S'], [(1, 'A', ['C'], ['1|0:9']), (3, 'G', ['T'], ['1|0:9'])]),
                       ('b', ['S'], [(1, 'A', ['C'], ['0|1:9']), (2, 'C', ['A'], ['0|1:9'])]),
                       ('c', ['S'], [(2, 'C', ['A'], ['0|1:9']), (3, 'G', ['T'], ['1|0:9'])]))
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source = directory / 'input.json'
            output = directory / 'report.json'
            source.write_text(json.dumps(payload))
            original = b'{"previous": true}\n'
            output.write_bytes(original)
            done = self._cli(source, output)
            self.assertEqual(done.returncode, 2)
            self.assertEqual(output.read_bytes(), original)

    def test_malformed_json_exits_2_and_preserves_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source = directory / 'input.json'
            output = directory / 'report.json'
            source.write_text('{"fasta": 1, "fasta": 2}')
            output.write_text('old')
            done = self._cli(source, output)
            self.assertEqual(done.returncode, 2)
            self.assertEqual(output.read_text(), 'old')

    def test_input_equals_output_and_symlink_rejected(self):
        payload = data(('a', ['S'], [(5, 'A', ['G'], ['1|0:1'])]))
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source = directory / 'input.json'
            source.write_text(json.dumps(payload))
            original = source.read_bytes()
            link = directory / 'link.json'
            link.symlink_to(source)
            for destination in (source, link):
                done = self._cli(source, destination)
                self.assertEqual(done.returncode, 2)
                self.assertEqual(source.read_bytes(), original)

    def test_success_replaces_atomically(self):
        payload = data(('a', ['S'], [(5, 'A', ['G'], ['1|0:1'])]))
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source = directory / 'input.json'
            output = directory / 'report.json'
            source.write_text(json.dumps(payload))
            output.write_text('stale')
            done = self._cli(source, output)
            self.assertEqual(done.returncode, 0, done.stderr)
            report = json.loads(output.read_text())
            self.assertEqual(report['schema_version'], 2)


if __name__ == '__main__':
    unittest.main()
