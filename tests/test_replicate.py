import copy
import json
import os
import random
import subprocess
import sys
import tempfile
import unittest
from itertools import combinations
from pathlib import Path

from variant_cohort.cli import run
from variant_cohort.replicate import audit_replicates

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'src'
HEADER = '##fileformat=VCFv4.2\n'
COLUMNS = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
FASTA = '>chr1\n' + 'ACGT' * 60 + '\n'
DEFAULT_RATIO = {'numerator': 1, 'denominator': 2}


def vcf_text(samples, rows, header=HEADER):
    """rows: (pos, ref, alts, cells) where cells align with samples."""
    lines = [header.rstrip('\n'), COLUMNS + '\t' + '\t'.join(samples)]
    for pos, ref, alts, cells in rows:
        lines.append('\t'.join(['chr1', str(pos), '.', ref,
                                ','.join(alts), '.', 'PASS', '.', 'GT']
                               + list(cells)))
    return '\n'.join(lines) + '\n'


def payload(samples, rows, minimum=1, ratio=None, fasta=FASTA):
    return {'fasta': fasta,
            'vcf': vcf_text(samples, rows),
            'minimum_informative_sites': minimum,
            'max_discordance': DEFAULT_RATIO if ratio is None else ratio}


def reorder(samples, rows, seed=0, shuffle_samples=False):
    """Permute record rows and ALT order (remapping GT).

    With shuffle_samples the sample columns are permuted too; callers must
    compare via named samples because left/right always follow lexicographic
    sample order.
    """
    sample_rng = random.Random(seed)
    sample_order = samples[:]
    if shuffle_samples:
        sample_rng.shuffle(sample_order)
    rows = copy.deepcopy(rows)
    row_rng = random.Random(seed + 101)
    row_order = list(range(len(rows)))
    row_rng.shuffle(row_order)
    out_rows = []
    for index in row_order:
        pos, ref, alts, cells = rows[index]
        alt_rng = random.Random(seed * 1009 + index)
        alt_perm = list(range(len(alts)))
        alt_rng.shuffle(alt_perm)
        new_alts = [alts[i] for i in alt_perm]
        # A genotype carrying old index i+1 is rewritten to the position of
        # that same allele sequence in the permuted ALT enumeration.
        remap = [0] + [1 + alt_perm.index(i) for i in range(len(alts))]
        remapped_cells = {}
        for name, cell in zip(samples, cells):
            sep = '|' if '|' in cell else '/'
            pieces = cell.replace('|', '/').split('/')
            remapped_cells[name] = sep.join(
                '.' if piece == '.' else str(remap[int(piece)])
                for piece in pieces)
        new_cells = [remapped_cells[name] for name in sample_order]
        out_rows.append((pos, ref, new_alts, new_cells))
    return sample_order, out_rows


def pair_map(result):
    return {(p['left'], p['right']): p for p in result['pairs']}


class ReplicateLibrary(unittest.TestCase):
    def test_report_shape_and_pair_coverage(self):
        samples = ['A', 'B', 'C', 'D']
        rows = [(1, 'A', ['C'], ['0/1'] * 4),
                (2, 'C', ['G'], ['1/1'] * 4)]
        result = audit_replicates(payload(samples, rows))
        self.assertEqual(set(result), {'samples', 'source_records',
                                       'parameters', 'pairs'})
        self.assertEqual(result['samples'], sorted(samples))
        self.assertEqual(result['source_records'], 2)
        self.assertEqual(result['parameters'],
                         {'minimum_informative_sites': 1,
                          'max_discordance': DEFAULT_RATIO})
        expected = list(combinations(sorted(samples), 2))
        self.assertEqual([(p['left'], p['right']) for p in result['pairs']],
                         expected)
        for pair in result['pairs']:
            self.assertEqual(set(pair),
                             {'left', 'right', 'complete_sites',
                              'informative_sites', 'concordant_sites',
                              'discordant_sites', 'status',
                              'mismatch_examples', 'omitted_mismatches'})
            self.assertEqual(pair['informative_sites'],
                             pair['concordant_sites'] + pair['discordant_sites'])

    def test_multi_alt_reorder_keeps_genotype_identity(self):
        samples = ['A', 'B']
        base = [(1, 'A', ['C', 'T', 'G'], ['0/1', '2/3']),
                (2, 'C', ['G', 'T'], ['1/2', '2/1'])]
        first = audit_replicates(payload(samples, base))
        for seed in range(6):
            reordered_samples, reordered_rows = reorder(samples, base, seed)
            again = audit_replicates(
                payload(reordered_samples, reordered_rows))
            self.assertEqual(again, first)
        pair = pair_map(first)[('A', 'B')]
        self.assertEqual(pair['informative_sites'], 2)
        self.assertEqual(pair['concordant_sites'], 1)
        self.assertEqual(pair['discordant_sites'], 1)
    def test_phase_separator_equivalence(self):
        samples = ['A', 'B']
        slash = audit_replicates(payload(
            samples, [(1, 'A', ['C', 'T'], ['0/1', '1/0']),
                      (2, 'C', ['G'], ['1|0', '0|1'])]))
        pipe = audit_replicates(payload(
            samples, [(1, 'A', ['C', 'T'], ['0|1', '1|0']),
                      (2, 'C', ['G'], ['1/0', '0/1'])]))
        self.assertEqual(pipe, slash)
        self.assertEqual(pair_map(slash)[('A', 'B')]['discordant_sites'], 0)

    def test_missing_calls_excluded_from_denominators(self):
        samples = ['A', 'B', 'C']
        rows = [
            (1, 'A', ['C'], ['0/1', './.', '0/1']),
            (2, 'C', ['G'], ['1/1', '0/1', '1/1']),
            (3, 'G', ['A'], ['0/1', '1/1', './.']),
            (4, 'T', ['N'], ['0/0', '0/0', '0/0']),
        ]
        result = pair_map(audit_replicates(payload(samples, rows)))
        ab = result[('A', 'B')]
        self.assertEqual(ab['complete_sites'], 3)
        self.assertEqual(ab['informative_sites'], 2)
        self.assertEqual(ab['discordant_sites'], 2)
        ac = result[('A', 'C')]
        self.assertEqual(ac['complete_sites'], 3)
        self.assertEqual(ac['informative_sites'], 2)
        bc = result[('B', 'C')]
        self.assertEqual(bc['complete_sites'], 2)
        self.assertEqual(bc['informative_sites'], 1)

    def test_reference_homozygous_sites_not_informative(self):
        samples = ['A', 'B']
        rows = [(1, 'A', ['C'], ['0/0', '0/0']),
                (2, 'C', ['G'], ['0/0', '0/0']),
                (3, 'G', ['A'], ['0/1', '0/1'])]
        result = audit_replicates(payload(samples, rows, minimum=2))
        pair = pair_map(result)[('A', 'B')]
        self.assertEqual(pair['complete_sites'], 3)
        self.assertEqual(pair['informative_sites'], 1)
        self.assertEqual(pair['concordant_sites'], 1)
        self.assertEqual(pair['status'], 'insufficient_evidence')

    def test_ratio_boundary_equality_is_candidate(self):
        samples = ['A', 'B']
        rows = []
        for pos in range(1, 5):
            ref = 'ACGT'[(pos - 1) % 4]
            alt = 'N'
            cells = (['0/1', '0/1'] if pos <= 2 else ['0/1', '1/1'])
            rows.append((pos, ref, [alt], cells))
        exact = audit_replicates(
            payload(samples, rows, ratio={'numerator': 1, 'denominator': 2}))
        self.assertEqual(pair_map(exact)[('A', 'B')]['status'], 'candidate')
        over = audit_replicates(
            payload(samples, rows, ratio={'numerator': 49, 'denominator': 100}))
        self.assertEqual(pair_map(over)[('A', 'B')]['status'], 'different')

    def test_one_extra_mismatch_flips_candidate_to_different(self):
        samples = ['A', 'B']
        # Two concordant informative sites, then discordant sites added.
        rows = [(pos, 'ACGT'[(pos - 1) % 4], ['N'],
                 ['0/1', '0/1']) for pos in (1, 2)]
        ratio = {'numerator': 1, 'denominator': 2}
        # 2 informative / 1 discordant: 1/2 is at the boundary -> candidate.
        rows_plus_one = rows + [(3, 'G', ['N'], ['0/1', '1/1'])]
        self.assertEqual(pair_map(audit_replicates(
            payload(samples, rows_plus_one, ratio=ratio)))[('A', 'B')]
            ['status'], 'candidate')
        # One more discordant site: 2/4 is at the boundary -> candidate.
        rows_plus_two = rows_plus_one + [(4, 'T', ['N'], ['0/1', '1/1'])]
        self.assertEqual(pair_map(audit_replicates(
            payload(samples, rows_plus_two, ratio=ratio)))[('A', 'B')]
            ['status'], 'candidate')
        # A fifth discordant site pushes 3/5 above 1/2 -> different.
        rows_plus_three = rows_plus_two + [(5, 'A', ['N'], ['0/1', '1/1'])]
        self.assertEqual(pair_map(audit_replicates(
            payload(samples, rows_plus_three, ratio=ratio)))[('A', 'B')]
            ['status'], 'different')

    def test_insufficient_evidence_status(self):
        samples = ['A', 'B']
        rows = [(1, 'A', ['C'], ['0/1', '1/1']),
                (2, 'C', ['G'], ['0/1', '0/1'])]
        result = audit_replicates(payload(samples, rows, minimum=3))
        self.assertEqual(pair_map(result)[('A', 'B')]['status'],
                         'insufficient_evidence')

    def test_zero_numerator_requires_perfect_match(self):
        samples = ['A', 'B']
        identical = payload(samples, [(1, 'A', ['C'], ['0/1', '1/0'])],
                            ratio={'numerator': 0, 'denominator': 1})
        self.assertEqual(pair_map(audit_replicates(identical))[('A', 'B')]
                         ['status'], 'candidate')
        different = payload(samples, [(1, 'A', ['C'], ['0/1', '1/1'])],
                            ratio={'numerator': 0, 'denominator': 1})
        self.assertEqual(pair_map(audit_replicates(different))[('A', 'B')]
                         ['status'], 'different')

    def test_empty_record_vcf_keeps_all_pairs_insufficient(self):
        samples = ['A', 'B', 'C']
        result = audit_replicates(payload(samples, []))
        self.assertEqual(result['source_records'], 0)
        self.assertEqual([(p['left'], p['right']) for p in result['pairs']],
                         list(combinations(sorted(samples), 2)))
        for pair in result['pairs']:
            self.assertEqual(pair['status'], 'insufficient_evidence')
            self.assertEqual(pair['complete_sites'], 0)
            self.assertEqual(pair['informative_sites'], 0)
            self.assertEqual(pair['mismatch_examples'], [])
            self.assertEqual(pair['omitted_mismatches'], 0)

    def test_mismatch_examples_canonical_order_and_truncation(self):
        samples = ['A', 'B']
        rows = []
        for pos in range(1, 12):
            ref = 'ACGT'[(pos - 1) % 4]
            alts = [a for a in ('N', 'G') if a != ref]
            rows.append((pos, ref, alts, ['0/1', '1/1']))
        result = pair_map(audit_replicates(payload(samples, rows)))['A', 'B']
        self.assertEqual(result['discordant_sites'], 11)
        self.assertEqual(len(result['mismatch_examples']), 8)
        self.assertEqual(result['omitted_mismatches'], 3)
        identities = [(e['contig'], e['position'], e['ref'], tuple(e['alts']))
                      for e in result['mismatch_examples']]
        self.assertEqual(identities, sorted(identities))
        self.assertEqual([e['position'] for e in result['mismatch_examples']],
                         list(range(1, 9)))
        first = result['mismatch_examples'][0]
        self.assertEqual(set(first), {'contig', 'position', 'ref', 'alts',
                                      'left_gt', 'right_gt'})
        self.assertEqual(first['left_gt'], sorted(first['left_gt']))
        self.assertEqual(first['right_gt'], sorted(first['right_gt']))

    def test_example_gt_uses_sequence_identities(self):
        samples = ['A', 'B']
        rows = [(1, 'A', ['T', 'C'], ['0/2', '1/1'])]
        result = pair_map(audit_replicates(payload(samples, rows)))['A', 'B']
        example = result['mismatch_examples'][0]
        self.assertEqual(example['alts'], ['C', 'T'])
        self.assertEqual(example['left_gt'], ['A', 'C'])
        self.assertEqual(example['right_gt'], ['T', 'T'])
        # Same biological genotype, different source ALT enumeration.
        swapped = payload(samples, [(1, 'A', ['C', 'T'], ['0/1', '2/2'])])
        swapped_example = pair_map(audit_replicates(swapped))[
            ('A', 'B')]['mismatch_examples'][0]
        self.assertEqual(swapped_example['left_gt'], ['A', 'C'])
        self.assertEqual(swapped_example['right_gt'], ['T', 'T'])

    def test_row_and_alt_permutation_invariance(self):
        samples = ['S1', 'S2', 'S3', 'S4']
        rows = [
            (1, 'A', ['C', 'T'], ['0/1', '1/1', '0/0', './1']),
            (2, 'C', ['G'], ['1|0', '0/1', '0/0', '1/1']),
            (3, 'G', ['A', 'N'], ['2/1', '0/2', './.', '1/1']),
            (4, 'T', ['N'], ['0/1', '0/0', '0/1', '0/0']),
        ]
        baseline = audit_replicates(payload(samples, rows))
        for seed in range(10):
            reordered_samples, reordered_rows = reorder(samples, rows, seed)
            self.assertEqual(reordered_samples, samples)
            variant = audit_replicates(
                payload(reordered_samples, reordered_rows))
            self.assertEqual(variant, baseline)

    def test_sample_column_permutation_invariance_by_named_samples(self):
        samples = ['S1', 'S2', 'S3', 'S4']
        rows = [
            (1, 'A', ['C', 'T'], ['0/1', '1/1', '0/0', './1']),
            (2, 'C', ['G'], ['1|0', '0/1', '0/0', '1/1']),
            (3, 'G', ['A', 'N'], ['2/1', '0/2', './.', '1/1']),
            (4, 'T', ['N'], ['0/1', '0/0', '0/1', '0/0']),
        ]
        baseline = pair_map(audit_replicates(payload(samples, rows)))
        for seed in range(10):
            reordered_samples, reordered_rows = reorder(
                samples, rows, seed, shuffle_samples=True)
            variant = pair_map(audit_replicates(
                payload(reordered_samples, reordered_rows)))
            self.assertEqual(set(variant), set(baseline))
            for key, pair in baseline.items():
                other = dict(variant[key])
                expected = dict(pair)
                # Example GTs stay attached to the named left/right samples.
                self.assertEqual(other, expected)

    def test_duplicate_canonical_identity_rejected_even_with_alt_permutation(self):
        samples = ['A', 'B']
        rows = [(1, 'A', ['C', 'T'], ['0/1', '1/1']),
                (1, 'A', ['T', 'C'], ['0/2', '2/2'])]
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            audit_replicates(payload(samples, rows))

    def test_haploid_extra_format_and_index_rules(self):
        samples = ['A', 'B']
        base = [(1, 'A', ['C'], ['0/1', '1/1'])]
        bad_cells = [
            ['1', '1/1'],
            ['.', '1/1'],
            ['0/1:30', '1/1'],
            ['0/2', '1/1'],
            ['0//1', '1/1'],
        ]
        for cells in bad_cells:
            with self.assertRaises(ValueError):
                audit_replicates(payload(samples, [(1, 'A', ['C'], cells)]))
        good = payload(samples, base)
        audit_replicates(good)

    def test_sample_and_scale_bounds(self):
        one = payload(['A'], [(1, 'A', ['C'], ['0/1'])])
        with self.assertRaisesRegex(ValueError, 'samples'):
            audit_replicates(one)
        many = payload(['S%d' % i for i in range(101)], [])
        with self.assertRaisesRegex(ValueError, 'samples'):
            audit_replicates(many)
        boundary = payload(['S%d' % i for i in range(100)], [])
        self.assertEqual(len(audit_replicates(boundary)['pairs']),
                         100 * 99 // 2)

    def test_parameter_validation(self):
        samples = ['A', 'B']
        rows = [(1, 'A', ['C'], ['0/1', '1/1'])]
        for minimum in (0, 5001, True, '2', 1.0):
            with self.assertRaises(ValueError):
                audit_replicates(payload(samples, rows, minimum=minimum))
        bad_ratios = [
            {'numerator': 2, 'denominator': 1},
            {'numerator': -1, 'denominator': 1},
            {'numerator': 0, 'denominator': 0},
            {'numerator': 0.5, 'denominator': 1},
            {'numerator': 1},
            {'numerator': 1, 'denominator': 2, 'extra': 3},
        ]
        for ratio in bad_ratios:
            with self.assertRaises(ValueError):
                audit_replicates(payload(samples, rows, ratio=ratio))

    def test_outer_contract_fields(self):
        data = payload(['A', 'B'], [(1, 'A', ['C'], ['0/1', '1/1'])])
        audit_replicates(data)
        for broken in (dict(data, extra=1), {k: v for k, v in data.items()
                                             if k != 'vcf'}):
            with self.assertRaises(ValueError):
                audit_replicates(broken)

    def test_reference_check_still_applies(self):
        data = payload(['A', 'B'], [(1, 'A', ['C'], ['0/1', '1/1'])],
                       fasta='>chr1\n' + 'C' * 80 + '\n')
        with self.assertRaisesRegex(ValueError, 'reference mismatch'):
            audit_replicates(data)


class CliEndToEnd(unittest.TestCase):
    def _env(self):
        return dict(os.environ, PYTHONPATH=str(SRC))

    def _cli(self, source, output, mode=True):
        command = [sys.executable, '-m', 'variant_cohort.cli',
                   '--source', str(source), '--output', str(output)]
        if mode:
            command.append('--audit-replicates')
        return subprocess.run(command, capture_output=True, env=self._env(),
                              cwd=ROOT)

    def test_end_to_end_fixture_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            output = directory / 'report.json'
            done = self._cli(ROOT / 'fixtures' / 'replicate.json', output)
            self.assertEqual(done.returncode, 0, done.stderr)
            report = json.loads(output.read_text())
            self.assertEqual(set(report), {'schema_version', 'mode',
                                           'source_sha256', 'result'})
            self.assertEqual(report['schema_version'], 5)
            self.assertEqual(report['mode'], 'replicate_audit')
            self.assertRegex(report['source_sha256'], r'[0-9a-f]{64}')
            result = report['result']
            statuses = {tuple(sorted((p['left'], p['right']))): p['status']
                        for p in result['pairs']}
            self.assertEqual(statuses[('S1', 'S2')], 'candidate')
            self.assertEqual(statuses[('S1', 'S3')], 'different')
            self.assertEqual(statuses[('S2', 'S3')], 'different')

    def test_library_and_cli_reports_match_except_hash_immutability(self):
        library = run(ROOT / 'fixtures' / 'replicate.json',
                      audit_replicates=True)
        raw = json.loads((ROOT / 'fixtures' / 'replicate.json').read_text())
        self.assertEqual(library['result'], audit_replicates(raw))
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / 'report.json'
            done = self._cli(ROOT / 'fixtures' / 'replicate.json', output)
            self.assertEqual(done.returncode, 0, done.stderr)
            self.assertEqual(json.loads(output.read_text()), library)

    def test_permuted_inputs_identical_report_except_hash(self):
        samples = ['A', 'B', 'C', 'D']
        rows = [
            (1, 'A', ['C', 'T'], ['0/1', '1/1', '0/0', './1']),
            (2, 'C', ['G'], ['1|0', '0/1', '0/0', '1/1']),
            (3, 'G', ['A', 'N'], ['2/1', '0/2', './.', '1/1']),
        ]
        baseline = audit_replicates(payload(samples, rows))
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            for seed in range(5):
                reordered_samples, reordered_rows = reorder(samples, rows, seed)
                source = directory / ('input%d.json' % seed)
                output = directory / ('report%d.json' % seed)
                source.write_text(json.dumps(
                    payload(reordered_samples, reordered_rows)))
                done = self._cli(source, output)
                self.assertEqual(done.returncode, 0, done.stderr)
                report = json.loads(output.read_text())
                self.assertEqual(report['result'], baseline)

    def test_invalid_input_exits_2_and_keeps_output_bytes(self):
        duplicate = payload(['A', 'B'],
                            [(1, 'A', ['C', 'T'], ['0/1', '0/2'])])
        duplicate['vcf'] += (
            'chr1\t1\t.\tA\tT,C\t.\tPASS\t.\tGT\t0/2\t0/1\n')
        bad_ref = payload(['A', 'B'], [(1, 'A', ['C'], ['0/1', '1/1'])])
        bad_ref['fasta'] = '>chr1\n' + 'C' * 80 + '\n'
        cases = [duplicate, bad_ref]
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            output = directory / 'report.json'
            for index, data in enumerate(cases):
                source = directory / ('bad%d.json' % index)
                source.write_text(json.dumps(data))
                original = b'{"previous": true}\n'
                output.write_bytes(original)
                done = self._cli(source, output)
                self.assertEqual(done.returncode, 2, done.stderr)
                self.assertEqual(output.read_bytes(), original)

    def test_input_equals_output_and_symlink_rejected(self):
        data = payload(['A', 'B'], [(1, 'A', ['C'], ['0/1', '1/1'])])
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source = directory / 'input.json'
            source.write_text(json.dumps(data))
            original = source.read_bytes()
            link = directory / 'link.json'
            link.symlink_to(source)
            for destination in (source, link):
                done = self._cli(source, destination)
                self.assertEqual(done.returncode, 2)
                self.assertEqual(source.read_bytes(), original)

    def test_atomic_replace_on_success(self):
        data = payload(['A', 'B'], [(1, 'A', ['C'], ['0/1', '1/1'])])
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source = directory / 'input.json'
            output = directory / 'report.json'
            source.write_text(json.dumps(data))
            output.write_text('stale')
            done = self._cli(source, output)
            self.assertEqual(done.returncode, 0, done.stderr)
            self.assertEqual(json.loads(output.read_text())['schema_version'], 5)

    def test_mode_flags_are_mutually_exclusive(self):
        data = payload(['A', 'B'], [(1, 'A', ['C'], ['0/1', '1/1'])])
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source = directory / 'input.json'
            output = directory / 'report.json'
            source.write_text(json.dumps(data))
            command = [sys.executable, '-m', 'variant_cohort.cli',
                       '--source', str(source), '--output', str(output),
                       '--audit-replicates', '--audit-pedigree']
            done = subprocess.run(command, capture_output=True,
                                  env=self._env(), cwd=ROOT)
            self.assertEqual(done.returncode, 2)

    def test_v5_rejects_other_mode_shaped_sources(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            output = directory / 'report.json'
            output.write_text('keep')
            done = self._cli(ROOT / 'fixtures' / 'input.json', output)
            self.assertEqual(done.returncode, 2)
            self.assertEqual(output.read_text(), 'keep')

    def test_legacy_modes_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            output = directory / 'report.json'
            done = self._cli(ROOT / 'fixtures' / 'input.json', output,
                             mode=False)
            self.assertEqual(done.returncode, 0, done.stderr)
            report = json.loads(output.read_text())
            self.assertEqual(report['schema_version'], 1)
            self.assertNotIn('mode', report)


if __name__ == '__main__':
    unittest.main()
