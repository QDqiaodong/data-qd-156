import json
import os
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from variant_cohort.cli import run
from variant_cohort.model import merge_cohorts

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'src'
HEADER = '##fileformat=VCFv4.2\n'
COLUMNS = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
FASTA = '>chr1\n' + 'ACGT' * 80 + '\n'


def vcf_text(samples, rows, chrom='chr1'):
    """rows: (pos, ref, alts, cells) where cells align with samples."""
    lines = [HEADER.rstrip('\n'), COLUMNS + '\t' + '\t'.join(samples)]
    for pos, ref, alts, cells in rows:
        lines.append('\t'.join([chrom, str(pos), '.', ref,
                                ','.join(alts), '.', 'PASS', '.', 'GT']
                               + list(cells)))
    return '\n'.join(lines) + '\n'


def empty_vcf(samples):
    return HEADER + COLUMNS + '\t' + '\t'.join(samples) + '\n'


def payload(cohorts, fasta=FASTA):
    """cohorts: list of (id, vcf_text)."""
    return {'fasta': fasta,
            'cohorts': [{'id': cid, 'vcf': text} for cid, text in cohorts]}


def remap_cell(cell, remap):
    """Rewrite one GT cell's numeric indexes through remap, keeping separator."""
    if '|' in cell:
        pieces, sep = cell.split('|'), '|'
    elif '/' in cell:
        pieces, sep = cell.split('/'), '/'
    else:
        pieces, sep = [cell], None
    out = ['.' if piece == '.' else str(remap[int(piece)]) for piece in pieces]
    return out[0] if sep is None else sep.join(out)


def permute(cohorts, seed):
    """cohorts: list of (id, samples, rows) -> permuted payload input.

    Permutes cohort order, sample columns, data rows and ALT enumerations,
    remapping GT indexes to preserve allele sequence identity.
    """
    rng = random.Random(seed)
    indexed = list(enumerate(cohorts))
    rng.shuffle(indexed)
    out = []
    for old_index, (cid, samples, rows) in indexed:
        sample_order = samples[:]
        rng.shuffle(sample_order)
        row_order = list(range(len(rows)))
        rng.shuffle(row_order)
        new_rows = []
        for row_index in row_order:
            pos, ref, alts, cells = rows[row_index]
            alt_perm = list(range(len(alts)))
            alt_rng = random.Random(seed * 100003 + old_index * 917 + row_index)
            alt_rng.shuffle(alt_perm)
            new_alts = [alts[i] for i in alt_perm]
            # Old index i+1 moves to the position of that allele sequence.
            remap = [0] + [1 + alt_perm.index(i) for i in range(len(alts))]
            cells_by_sample = dict(zip(samples, cells))
            new_cells = [remap_cell(cells_by_sample[name], remap)
                         for name in sample_order]
            new_rows.append((pos, ref, new_alts, new_cells))
        out.append((cid, vcf_text(sample_order, new_rows)))
    return payload(out)


def site_map(result):
    return {(s['contig'], s['position']): s for s in result['sites']}


class MergeLibrary(unittest.TestCase):
    def test_minimal_semantic_example_from_contract(self):
        batch_a = vcf_text(['S1'], [
            (2, 'C', ['T'], ['0/1']),
            (5, 'A', ['T'], ['0/1']),
        ])
        batch_b = vcf_text(['S2'], [
            (2, 'C', ['G'], ['1|0']),
        ])
        result = merge_cohorts(payload([('cA', batch_a), ('cB', batch_b)]))
        self.assertEqual(set(result), {'samples', 'source_cohorts',
                                       'source_records', 'sites'})
        self.assertEqual(result['samples'], ['S1', 'S2'])
        self.assertEqual(result['source_cohorts'], ['cA', 'cB'])
        self.assertEqual(result['source_records'], 3)
        sites = site_map(result)
        shared = sites[('chr1', 2)]
        self.assertEqual(shared['ref'], 'C')
        self.assertEqual(shared['alts'], ['G', 'T'])
        self.assertEqual(shared['source_cohorts'], ['cA', 'cB'])
        self.assertEqual(shared['genotypes'], [
            {'alleles': [0, 2], 'phased': False},
            {'alleles': [1, 0], 'phased': True},
        ])
        only_a = sites[('chr1', 5)]
        # S2's batch never observed this site: null, never an assumed 0/0.
        self.assertEqual(only_a['genotypes'], [
            {'alleles': [0, 1], 'phased': False}, None,
        ])
        self.assertEqual(only_a['source_cohorts'], ['cA'])

    def test_alt_union_remaps_by_sequence_identity(self):
        # Three batches enumerate different, partially overlapping ALT sets.
        a = vcf_text(['S1'], [(4, 'T', ['A', 'C'], ['0/1'])])
        b = vcf_text(['S2'], [(4, 'T', ['G', 'A'], ['2/1'])])
        c = vcf_text(['S3'], [(4, 'T', ['C'], ['1/1'])])
        result = merge_cohorts(payload([('b1', a), ('b2', b), ('b3', c)]))
        site = site_map(result)[('chr1', 4)]
        self.assertEqual(site['alts'], ['A', 'C', 'G'])
        # S1: REF + A -> [0,1]; S2 keeps source copy order 2/1 -> A then G
        # -> [1,3]; S3: C,C -> [2,2].
        self.assertEqual([g['alleles'] for g in site['genotypes']],
                         [[0, 1], [1, 3], [2, 2]])
        self.assertEqual(site['source_cohorts'], ['b1', 'b2', 'b3'])

    def test_phase_separator_copy_order_and_ploidy_preserved(self):
        rows = [(3, 'G', ['A', 'C', 'T'], ['2/1'])]
        a = vcf_text(['S1'], rows)
        b = vcf_text(['S2'], [(3, 'G', ['T', 'C', 'A'], ['1|3'])])
        c = vcf_text(['S3'], [(3, 'G', ['C'], ['1'])])          # haploid
        d = vcf_text(['S4'], [(3, 'G', ['A'], ['.'])])          # haploid missing
        result = merge_cohorts(payload([('k1', a), ('k2', b),
                                        ('k3', c), ('k4', d)]))
        genotypes = site_map(result)[('chr1', 3)]['genotypes']
        by_sample = dict(zip(result['samples'], genotypes))
        # Joint alts sorted: A=1, C=2, T=3.
        # S1 source 2/1 = C,A, kept in order -> [2,1], unphased.
        self.assertEqual(by_sample['S1'],
                         {'alleles': [2, 1], 'phased': False})
        # S2 source 1|3 = T,A, kept in order -> [3,1], phased.
        self.assertEqual(by_sample['S2'],
                         {'alleles': [3, 1], 'phased': True})
        self.assertEqual(by_sample['S3'],
                         {'alleles': [2], 'phased': False})
        self.assertEqual(by_sample['S4'],
                         {'alleles': [None], 'phased': False})

    def test_partial_missing_copies_keep_position(self):
        a = vcf_text(['S1'], [(6, 'C', ['G', 'T'], ['./2'])])
        b = vcf_text(['S2'], [(6, 'C', ['T'], ['1/.'])])
        result = merge_cohorts(payload([('x', a), ('y', b)]))
        genotypes = site_map(result)[('chr1', 6)]['genotypes']
        self.assertEqual(genotypes[0], {'alleles': [None, 2], 'phased': False})
        self.assertEqual(genotypes[1], {'alleles': [2, None], 'phased': False})

    def test_observed_hom_ref_differs_from_unobserved_null(self):
        a = vcf_text(['S1'], [(1, 'A', ['C'], ['0/0'])])
        b = vcf_text(['S2'], [(4, 'T', ['A'], ['0/1'])])
        result = merge_cohorts(payload([('a', a), ('b', b)]))
        sites = site_map(result)
        self.assertEqual(sites[('chr1', 1)]['genotypes'],
                         [{'alleles': [0, 0], 'phased': False}, None])
        self.assertEqual(sites[('chr1', 4)]['genotypes'],
                         [None, {'alleles': [0, 1], 'phased': False}])

    def test_compatible_same_ref_sites_merge(self):
        a = vcf_text(['S1'], [(2, 'C', ['T'], ['0/1'])])
        b = vcf_text(['S2'], [(2, 'C', ['T', 'G'], ['1/2'])])
        result = merge_cohorts(payload([('a', a), ('b', b)]))
        site = site_map(result)[('chr1', 2)]
        self.assertEqual(site['alts'], ['G', 'T'])
        self.assertEqual([g['alleles'] for g in site['genotypes']],
                         [[0, 2], [2, 1]])

    def test_sites_sorted_and_field_sets_exact(self):
        a = vcf_text(['S1'], [(8, 'T', ['A'], ['0/1']),
                              (2, 'C', ['G'], ['0/1'])])
        b = vcf_text(['S2'], [(2, 'C', ['A'], ['1/0'])])
        result = merge_cohorts(payload([('z', a), ('a', b)]))
        self.assertEqual([(s['contig'], s['position'], s['ref'])
                          for s in result['sites']],
                         [('chr1', 2, 'C'), ('chr1', 8, 'T')])
        for site in result['sites']:
            self.assertEqual(set(site), {'contig', 'position', 'ref', 'alts',
                                         'source_cohorts', 'genotypes'})
        self.assertEqual(result['source_cohorts'], ['a', 'z'])

    def test_ref_conflict_at_shared_site_rejected(self):
        # Both multi-base REF strings prefix-match the reference at the same
        # coordinate, so each batch passes the FASTA check alone, but the two
        # batches disagree on the REF string at the shared site.
        fasta = '>chr1\n' + 'ACGT' * 80 + '\n'
        a = vcf_text(['S1'], [(1, 'AC', ['AT'], ['0/1'])])
        b = vcf_text(['S2'], [(1, 'ACG', ['A'], ['0/1'])])
        with self.assertRaisesRegex(ValueError, 'REF conflict'):
            merge_cohorts(payload([('a', a), ('b', b)], fasta=fasta))

    def test_duplicate_site_within_cohort_rejected_even_with_different_alts(self):
        text = vcf_text(['S1'], [
            (2, 'C', ['T'], ['0/1']),
            (2, 'C', ['G'], ['0/1']),
        ])
        with self.assertRaisesRegex(ValueError, 'duplicate site'):
            merge_cohorts(payload([('a', text)]))

    def test_global_sample_name_overlap_rejected(self):
        a = vcf_text(['S1', 'S2'], [(1, 'A', ['C'], ['0/1', '0/0'])])
        b = vcf_text(['S2'], [(2, 'C', ['G'], ['0/1'])])
        with self.assertRaisesRegex(ValueError, 'shared across cohorts'):
            merge_cohorts(payload([('a', a), ('b', b)]))

    def test_empty_vcf_batches_keep_samples_and_add_no_sites(self):
        a = empty_vcf(['S1', 'S2'])
        b = vcf_text(['S3'], [(1, 'A', ['C'], ['0/1'])])
        c = empty_vcf(['S4'])
        result = merge_cohorts(payload([('a', a), ('b', b), ('c', c)]))
        self.assertEqual(result['samples'], ['S1', 'S2', 'S3', 'S4'])
        self.assertEqual(result['source_records'], 1)
        self.assertEqual(result['source_cohorts'], ['a', 'b', 'c'])
        self.assertEqual(len(result['sites']), 1)
        genotypes = result['sites'][0]['genotypes']
        self.assertEqual(genotypes, [None, None,
                                     {'alleles': [0, 1], 'phased': False},
                                     None])
        # All-header-only input is itself valid and yields an empty matrix.
        only_empty = merge_cohorts(payload([('a', a), ('c', c)]))
        self.assertEqual(only_empty['sites'], [])
        self.assertEqual(only_empty['source_records'], 0)
        self.assertEqual(only_empty['samples'], ['S1', 'S2', 'S4'])

    def test_invalid_gt_forms_rejected(self):
        for cell in ('0//1', '0/1/1', '0:1', '0/2', '-1/0', '0|1|0'):
            text = vcf_text(['S1'], [(1, 'A', ['C'], [cell])])
            with self.assertRaises(ValueError):
                merge_cohorts(payload([('a', text)]))
        # Haploid/diploid legal forms still parse.
        good = vcf_text(['S1'], [(1, 'A', ['C'], ['1'])])
        merge_cohorts(payload([('a', good)]))

    def test_scale_and_structure_limits(self):
        # 101 samples in one VCF exceeds the per-cohort cap.
        with self.assertRaisesRegex(ValueError, 'samples'):
            merge_cohorts(payload([('a', empty_vcf(
                ['S%03d' % i for i in range(101)]))]))
        # Global union over 200 across otherwise-valid batches.
        batches = [('c%d' % k, empty_vcf(
            ['s%03d' % (k * 100 + i) for i in range(100)])) for k in range(2)]
        batches.append(('c2', empty_vcf(['s200'])))
        with self.assertRaisesRegex(ValueError, '200'):
            merge_cohorts(payload(batches))
        # 33 cohorts exceeds the 32 cap.
        many = [('c%02d' % i, empty_vcf(['u%02d' % i])) for i in range(33)]
        with self.assertRaisesRegex(ValueError, 'cohorts'):
            merge_cohorts(payload(many))
        # 5001 total rows across two batches exceeds the joint cap while
        # neither single VCF exceeds its own 5000-row limit.
        fasta = '>chr1\n' + ('ACGT' * 1300) + '\n'
        rows_a = [(p, 'ACGT'[(p - 1) % 4], ['N'], ['0/1'])
                  for p in range(1, 2502)]
        rows_b = [(p, 'ACGT'[(p - 1) % 4], ['N'], ['0/1'])
                  for p in range(2502, 5002)]
        over = payload([('a', vcf_text(['S1'], rows_a)),
                        ('b', vcf_text(['S2'], rows_b))], fasta=fasta)
        with self.assertRaisesRegex(ValueError, '5000'):
            merge_cohorts(over)
        # Boundary: 5000 total rows is accepted.
        boundary = payload([('a', vcf_text(['S1'], rows_a)),
                            ('b', vcf_text(['S2'], rows_b[:-1]))],
                           fasta=fasta)
        self.assertEqual(merge_cohorts(boundary)['source_records'], 5000)

    def test_vcf_text_size_limit(self):
        padding = '##note=' + 'N' * 200 + '\n'
        text = (padding * 10000) + empty_vcf(['S1'])
        self.assertGreater(len(text), 2_000_000)
        with self.assertRaisesRegex(ValueError, 'VCF exceeds'):
            merge_cohorts(payload([('a', text)]))

    def test_cohort_id_and_outer_contract_validation(self):
        good = payload([('a', empty_vcf(['S1']))])
        merge_cohorts(good)
        for bad_id in ('', 'bad id', 'a/b', 'x' * 65):
            broken = payload([(bad_id, empty_vcf(['S1']))])
            with self.assertRaises(ValueError):
                merge_cohorts(broken)
        duplicate_id = payload([('a', empty_vcf(['S1'])),
                                ('a', empty_vcf(['S2']))])
        with self.assertRaisesRegex(ValueError, 'duplicate cohort'):
            merge_cohorts(duplicate_id)
        for broken in ({'fasta': FASTA},
                       {'fasta': FASTA, 'cohorts': []},
                       dict(good, extra=1),
                       {'fasta': FASTA, 'cohorts': 'nope'}):
            with self.assertRaises(ValueError):
                merge_cohorts(broken)
        bad_entry = {'fasta': FASTA,
                     'cohorts': [{'id': 'a', 'vcf': empty_vcf(['S1']),
                                  'extra': 1}]}
        with self.assertRaises(ValueError):
            merge_cohorts(bad_entry)

    def test_permutation_invariance(self):
        cohorts = [
            ('c1', ['S1', 'S2'], [
                (2, 'C', ['T', 'G'], ['0/1', '2|1']),
                (5, 'A', ['C'], ['./1', '0/0']),
            ]),
            ('c2', ['S3', 'S4'], [
                (2, 'C', ['G', 'A'], ['1/2', '1|0']),
                (8, 'T', ['A', 'C'], ['1/1', './.']),
            ]),
            ('c3', ['S5'], [
                (5, 'A', ['G', 'C'], ['2/1']),
            ]),
        ]
        baseline_text = payload(
            [(cid, vcf_text(samples, rows)) for cid, samples, rows in cohorts])
        baseline = merge_cohorts(baseline_text)
        for seed in range(12):
            variant = merge_cohorts(permute(cohorts, seed))
            self.assertEqual(variant, baseline)

    def test_duplicate_json_key_and_nonfinite_rejected_at_read(self):
        from variant_cohort.common import read
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / 'dup.json'
            source.write_text('{"fasta": ">x\\nA\\n", '
                              '"cohorts": [], "cohorts": []}')
            with self.assertRaisesRegex(ValueError, 'duplicate JSON key'):
                read(source)
            source.write_text('{"fasta": ">x\\nA\\n", "cohorts": [NaN]}')
            with self.assertRaisesRegex(ValueError, 'nonfinite'):
                read(source)


class CliEndToEnd(unittest.TestCase):
    def _env(self):
        return dict(os.environ, PYTHONPATH=str(SRC))

    def _cli(self, source, output, flag=True):
        command = [sys.executable, '-m', 'variant_cohort.cli',
                   '--source', str(source), '--output', str(output)]
        if flag:
            command.append('--merge-cohorts')
        return subprocess.run(command, capture_output=True, env=self._env(),
                              cwd=ROOT)

    def _cohorts_payload(self):
        return payload([
            ('cA', vcf_text(['S1', 'S2'], [
                (2, 'C', ['T'], ['0/1', '1|0']),
                (5, 'A', ['T'], ['0/0', '1']),
            ])),
            ('cB', vcf_text(['S3', 'S4'], [
                (2, 'C', ['G'], ['1|0', './1']),
                (8, 'T', ['C'], ['0/1', './.']),
            ])),
        ])

    def test_end_to_end_fixture_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / 'report.json'
            done = self._cli(ROOT / 'fixtures' / 'cohort_merge.json', output)
            self.assertEqual(done.returncode, 0, done.stderr)
            report = json.loads(output.read_text())
            self.assertEqual(set(report), {'schema_version', 'mode',
                                           'source_sha256', 'result'})
            self.assertEqual(report['schema_version'], 6)
            self.assertEqual(report['mode'], 'cohort_merge')
            self.assertRegex(report['source_sha256'], r'[0-9a-f]{64}')
            result = report['result']
            self.assertEqual(result['samples'], ['S1', 'S2', 'S3', 'S4'])
            self.assertEqual(result['source_cohorts'], ['batch-a', 'batch-b'])
            self.assertEqual(result['source_records'], 4)
            shared = site_map(result)[('chr1', 2)]
            self.assertEqual(shared['alts'], ['G', 'T'])
            self.assertEqual(shared['genotypes'][0],
                             {'alleles': [0, 2], 'phased': False})
            self.assertEqual(shared['genotypes'][2],
                             {'alleles': [1, 0], 'phased': True})

    def test_library_and_cli_reports_match(self):
        data = self._cohorts_payload()
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source = directory / 'input.json'
            output = directory / 'report.json'
            source.write_text(json.dumps(data))
            library = run(source, merge=True)
            self.assertEqual(library['result'], merge_cohorts(data))
            done = self._cli(source, output)
            self.assertEqual(done.returncode, 0, done.stderr)
            self.assertEqual(json.loads(output.read_text()), library)

    def test_permuted_inputs_identical_report_except_hash(self):
        cohorts = [
            ('c1', ['S1', 'S2'], [
                (2, 'C', ['T', 'G'], ['0/1', '2|1']),
                (5, 'A', ['C'], ['./1', '0/0']),
            ]),
            ('c2', ['S3', 'S4'], [
                (2, 'C', ['G', 'A'], ['1/2', '1|0']),
            ]),
        ]
        baseline = merge_cohorts(payload(
            [(cid, vcf_text(samples, rows)) for cid, samples, rows in cohorts]))
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            for seed in range(6):
                source = directory / ('in%d.json' % seed)
                output = directory / ('out%d.json' % seed)
                source.write_text(json.dumps(permute(cohorts, seed)))
                done = self._cli(source, output)
                self.assertEqual(done.returncode, 0, done.stderr)
                self.assertEqual(json.loads(output.read_text())['result'],
                                 baseline)

    def test_failures_exit_2_and_keep_output_bytes(self):
        conflict_fasta = '>chr1\n' + 'ACGT' * 80 + '\n'
        ref_conflict = payload([
            ('a', vcf_text(['S1'], [(1, 'AC', ['AT'], ['0/1'])])),
            ('b', vcf_text(['S2'], [(1, 'ACG', ['A'], ['0/1'])])),
        ], fasta=conflict_fasta)
        dup_site = payload([
            ('a', vcf_text(['S1'], [(2, 'C', ['T'], ['0/1']),
                                    (2, 'C', ['G'], ['0/1'])])),
        ])
        dup_sample = payload([
            ('a', empty_vcf(['S1'])),
            ('b', empty_vcf(['S1'])),
        ])
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            output = directory / 'report.json'
            for index, data in enumerate((ref_conflict, dup_site, dup_sample)):
                source = directory / ('bad%d.json' % index)
                source.write_text(json.dumps(data))
                original = b'{"previous": true}\n'
                output.write_bytes(original)
                done = self._cli(source, output)
                self.assertEqual(done.returncode, 2, done.stderr)
                self.assertEqual(output.read_bytes(), original)

    def test_input_equals_output_and_symlink_rejected(self):
        data = self._cohorts_payload()
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
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source = directory / 'input.json'
            output = directory / 'report.json'
            source.write_text(json.dumps(self._cohorts_payload()))
            output.write_text('stale')
            done = self._cli(source, output)
            self.assertEqual(done.returncode, 0, done.stderr)
            self.assertEqual(json.loads(output.read_text())['schema_version'], 6)

    def test_mode_flags_are_mutually_exclusive(self):
        data = self._cohorts_payload()
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source = directory / 'input.json'
            output = directory / 'report.json'
            source.write_text(json.dumps(data))
            for extra in ('--stitch-phase', '--project-assembly',
                          '--audit-pedigree', '--audit-replicates'):
                command = [sys.executable, '-m', 'variant_cohort.cli',
                           '--source', str(source), '--output', str(output),
                           '--merge-cohorts', extra]
                done = subprocess.run(command, capture_output=True,
                                      env=self._env(), cwd=ROOT)
                self.assertEqual(done.returncode, 2, extra)

    def test_v6_rejects_legacy_shaped_sources(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / 'report.json'
            output.write_text('keep')
            done = self._cli(ROOT / 'fixtures' / 'input.json', output)
            self.assertEqual(done.returncode, 2)
            self.assertEqual(output.read_text(), 'keep')

    def test_legacy_modes_v1_to_v5_unchanged(self):
        cases = [
            ('input.json', False, 1, None),
            ('phase.json', None, 2, 'phase_stitch'),
            ('projection.json', None, 3, 'assembly_projection'),
            ('pedigree.json', None, 4, 'pedigree_audit'),
            ('replicate.json', None, 5, 'replicate_audit'),
        ]
        flags = {
            2: '--stitch-phase',
            3: '--project-assembly',
            4: '--audit-pedigree',
            5: '--audit-replicates',
        }
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            for name, _, version, mode in cases:
                output = directory / ('report%d.json' % version)
                command = [sys.executable, '-m', 'variant_cohort.cli',
                           '--source', str(ROOT / 'fixtures' / name),
                           '--output', str(output)]
                if version >= 2:
                    command.append(flags[version])
                done = subprocess.run(command, capture_output=True,
                                      env=self._env(), cwd=ROOT)
                self.assertEqual(done.returncode, 0, done.stderr)
                report = json.loads(output.read_text())
                self.assertEqual(report['schema_version'], version)
                if mode is None:
                    self.assertNotIn('mode', report)
                else:
                    self.assertEqual(report['mode'], mode)


if __name__ == '__main__':
    unittest.main()
