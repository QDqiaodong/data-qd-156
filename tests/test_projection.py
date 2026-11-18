import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from variant_cohort.cli import run
from variant_cohort.project import project

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'src'
HEADER = '##fileformat=VCFv4.2\n'
COLUMNS = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'


def vcf(samples, rows, header=HEADER):
    """Build a GT-only VCF; rows are (contig, pos, ref, alts, cells)."""
    lines = [header.rstrip('\n'), COLUMNS + '\t' + '\t'.join(samples)]
    for contig, pos, ref, alts, cells in rows:
        lines.append('\t'.join([contig, str(pos), '.', ref, ','.join(alts),
                                '.', 'PASS', '.', 'GT'] + list(cells)))
    return '\n'.join(lines) + '\n'


def payload(source_fasta, target_fasta, vcf_text, chains):
    return {'source_fasta': source_fasta, 'target_fasta': target_fasta,
            'vcf': vcf_text, 'chains': chains}


def whole(chain_id, source_contig, target_contig, strand, length, source_start=1, target_start=1):
    return {'id': chain_id, 'source_contig': source_contig,
            'target_contig': target_contig, 'strand': strand,
            'blocks': [{'source_start': source_start,
                        'target_start': target_start, 'length': length}]}


def revcomp(text):
    return text.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


class NegativeStrandProjection(unittest.TestCase):
    def test_multiallelic_with_deletion_remaps_genotypes(self):
        # source revcomp is exactly the target, so projection is coordinate-free for ALT order.
        source = '>s1\nACGTTTGGCCAA\n'
        target = '>t1\nTTGGCCAAACGT\n'
        text = vcf(['S1', 'S2', 'S3', 'S4'], [
            ('s1', 2, 'CGTT', ['C', 'CGT'], ['1/2', '0/1', './1', '1'])])
        result = project(payload(source, target, text,
                                 [whole('m1', 's1', 't1', '-', 12)]))
        self.assertEqual(result['rejected'], [])
        entry = result['projected'][0]
        # revcomp(CGTT)=AACG, revcomp(C)=G, revcomp(CGT)=ACG; trimmed and sorted.
        self.assertEqual((entry['contig'], entry['position'], entry['ref'], entry['alts']),
                         ('t1', 8, 'AACG', ['ACG', 'G']))
        self.assertEqual(entry['genotypes'], [[1, 2], [0, 2], [2, None], [2]])

    def test_reverse_strand_N_stays_N(self):
        source = '>s\nAACNGT\n'
        target = '>t\nACNGTT\n'
        text = vcf(['S'], [('s', 3, 'CNGT', ['CAGT'], ['0/1'])])
        result = project(payload(source, target, text, [whole('m', 's', 't', '-', 6)]))
        entry = result['projected'][0]
        # revcomp(CNGT) = ACNG on the target; trimming leaves N as REF, never T.
        self.assertEqual((entry['position'], entry['ref'], entry['alts']), (3, 'N', ['T']))


class ReferenceSwap(unittest.TestCase):
    def test_target_reference_is_original_alt(self):
        source = '>s\nAAACCCGGGTTT\n'
        target = '>t\nAAATCCGGGTTT\n'
        text = vcf(['S1', 'S2'], [('s', 4, 'CCC', ['TCC', 'CCCC'], ['1/2', '0/1'])])
        result = project(payload(source, target, text,
                                 [whole('c', 's', 't', '+', 12)]))
        entry = result['projected'][0]
        self.assertEqual((entry['position'], entry['ref'], entry['alts']), (4, 'T', ['C', 'CC']))
        # original REF becomes an ALT; no allele may be dropped.
        self.assertEqual(len(entry['alts']), 2)
        self.assertEqual(entry['genotypes'], [[0, 2], [0, 1]])


class BlockCoverage(unittest.TestCase):
    def test_ref_spanning_adjacent_blocks_projects(self):
        source = '>s\nAACCGGTT\n'
        target = '>t\nAACCGGTT\n'
        blocks = [{'source_start': 1, 'target_start': 1, 'length': 4},
                  {'source_start': 5, 'target_start': 5, 'length': 4}]
        text = vcf(['S'], [('s', 3, 'CCGG', ['CCGA'], ['0/1'])])
        result = project(payload(source, target, text, [
            {'id': 'c', 'source_contig': 's', 'target_contig': 't',
             'strand': '+', 'blocks': blocks}]))
        entry = result['projected'][0]
        self.assertEqual((entry['position'], entry['ref'], entry['alts']), (6, 'G', ['A']))

    def test_target_gap_on_either_side_is_unmapped(self):
        source = '>s\nAACCGGTT\n'
        text = vcf(['S'], [('s', 3, 'CCGG', ['CCGA'], ['0/1'])])
        cases = [
            [{'source_start': 1, 'target_start': 1, 'length': 4},
             {'source_start': 5, 'target_start': 6, 'length': 3}],
            [{'source_start': 1, 'target_start': 2, 'length': 3},
             {'source_start': 5, 'target_start': 5, 'length': 4}],
        ]
        for blocks in cases:
            with self.subTest(blocks=blocks):
                target = '>t\n' + 'A' * 12 + '\n'
                result = project(payload(source, target, text, [
                    {'id': 'c', 'source_contig': 's', 'target_contig': 't',
                     'strand': '+', 'blocks': blocks}]))
                self.assertEqual(result['projected'], [])
                self.assertEqual(result['rejected'][0]['reason'], 'unmapped')

    def test_source_gap_between_blocks_is_unmapped(self):
        source = '>s\nAACCGGTT\n'
        target = '>t\nAACCGGTT\n'
        blocks = [{'source_start': 1, 'target_start': 1, 'length': 3},
                  {'source_start': 5, 'target_start': 5, 'length': 4}]
        text = vcf(['S'], [('s', 3, 'CCGG', ['CCGA'], ['0/1'])])
        result = project(payload(source, target, text, [
            {'id': 'c', 'source_contig': 's', 'target_contig': 't',
             'strand': '+', 'blocks': blocks}]))
        self.assertEqual(result['rejected'][0]['reason'], 'unmapped')

    def test_negative_strand_blocks_in_descending_target_order_project(self):
        source = '>s\nAACCGGTT\n'
        target = '>t\nAACCGGTT\n'
        blocks = [{'source_start': 1, 'target_start': 5, 'length': 4},
                  {'source_start': 5, 'target_start': 1, 'length': 4}]
        text = vcf(['S'], [('s', 3, 'CCGG', ['CCGA'], ['0/1'])])
        result = project(payload(source, target, text, [
            {'id': 'c', 'source_contig': 's', 'target_contig': 't',
             'strand': '-', 'blocks': blocks}]))
        # Source 3..6 maps to target 3..6 reversed; both REF and ALT are palindromic here.
        entry = result['projected'][0]
        self.assertEqual((entry['position'], entry['ref'], entry['alts']), (3, 'C', ['T']))

    def test_negative_strand_spanning_blocks_maps_with_reverse_coordinates(self):
        # Chain maps source 1..6 to target 7..12 and source 7..12 to target 1..6, i.e.
        # a full reverse-complement relation between the two assemblies.
        source = '>s\nACGTTTGGCCAA\n'
        target = '>t\nTTGGCCAAACGT\n'
        blocks = [{'source_start': 1, 'target_start': 7, 'length': 6},
                  {'source_start': 7, 'target_start': 1, 'length': 6}]
        text = vcf(['S'], [('s', 7, 'GGCC', ['GGCT'], ['0/1'])])
        result = project(payload(source, target, text, [
            {'id': 'c', 'source_contig': 's', 'target_contig': 't',
             'strand': '-', 'blocks': blocks}]))
        entry = result['projected'][0]
        # Source 7..10 lies wholly in the second block. On the - strand that block maps
        # source 7 to target 6 and source 10 to target 3, so the span is target 3..6.
        # revcomp(GGCC)=GGCC, revcomp(GGCT)=AGCC; trimming leaves G>A at position 3.
        self.assertEqual((entry['position'], entry['ref'], entry['alts']), (3, 'G', ['A']))
        self.assertEqual(entry['sources'][0]['position'], 7)

    def test_negative_strand_cross_block_gap_is_unmapped(self):
        source = '>s\nACGTTTGGCCAA\n'
        target = '>t\nTTGGCCAAACGTT\n'
        blocks = [{'source_start': 1, 'target_start': 8, 'length': 6},
                  {'source_start': 7, 'target_start': 1, 'length': 6}]
        # Source 5..8 straddles the 6/7 block boundary but the two target blocks are
        # separated (8..13 and 1..6), so the pair is not contiguous on the target.
        text = vcf(['S'], [('s', 5, 'TTGG', ['TTGT'], ['0/1'])])
        result = project(payload(source, target, text, [
            {'id': 'c', 'source_contig': 's', 'target_contig': 't',
             'strand': '-', 'blocks': blocks}]))
        self.assertEqual(result['projected'], [])
        self.assertEqual(result['rejected'][0]['reason'], 'unmapped')


class CandidateClassification(unittest.TestCase):
    def test_equivalent_chains_are_one_candidate(self):
        source = '>s\nAACCGG\n'
        target = '>t\nAACCGG\n'
        text = vcf(['S'], [('s', 3, 'C', ['T'], ['0/1'])])
        chains = [whole('z', 's', 't', '+', 6), whole('a', 's', 't', '+', 6)]
        result = project(payload(source, target, text, chains))
        self.assertEqual(result['rejected'], [])
        self.assertEqual(result['projected'][0]['chains'], ['a', 'z'])

    def test_distinct_usable_targets_are_ambiguous(self):
        source = '>s\nAACCGG\n'
        target = '>t1\nAACCGG\n>t2\nAATCGG\n'
        text = vcf(['S'], [('s', 3, 'C', ['T'], ['0/1'])])
        chains = [whole('c1', 's', 't1', '+', 6), whole('c2', 's', 't2', '+', 6)]
        result = project(payload(source, target, text, chains))
        self.assertEqual(result['projected'], [])
        self.assertEqual(result['rejected'][0]['reason'], 'ambiguous')

    def test_incompatible_target_is_not_ambiguous_against_a_usable_one(self):
        source = '>s\nAACCGG\n'
        target = '>good\nAACCGG\n>bad\nAAAAGG\n'
        text = vcf(['S'], [('s', 3, 'C', ['T'], ['0/1'])])
        chains = [whole('g', 's', 'good', '+', 6), whole('b', 's', 'bad', '+', 6)]
        result = project(payload(source, target, text, chains))
        self.assertEqual(result['projected'][0]['chains'], ['g'])

    def test_unmapped_and_reference_mismatch_are_distinct(self):
        source = '>s\nAACCGG\n>n\nAACCGG\n'
        target = '>t\nAAAAGG\n'
        rows = [('s', 3, 'C', ['T'], ['0/1']), ('n', 5, 'G', ['C'], ['0/1'])]
        text = vcf(['S'], rows)
        chains = [whole('c', 's', 't', '+', 6)]
        result = project(payload(source, target, text, chains))
        reasons = {(item['source']['contig'], item['source']['position']): item['reason']
                   for item in result['rejected']}
        self.assertEqual(reasons, {('s', 3): 'reference_mismatch', ('n', 5): 'unmapped'})

    def test_alt_of_different_length_cannot_become_reference(self):
        # target span is "C"; ALT "TT" is kept as ALT, never forced as the new REF.
        source = '>s\nAACCGG\n'
        target = '>t\nAATCGG\n'
        text = vcf(['S'], [('s', 3, 'C', ['TT', 'T'], ['1/2'])])
        result = project(payload(source, target, text, [whole('c', 's', 't', '+', 6)]))
        entry = result['projected'][0]
        self.assertEqual((entry['position'], entry['ref'], entry['alts']), (3, 'T', ['C', 'TT']))
        self.assertEqual(entry['genotypes'], [[0, 2]])


class ConvergenceAndMerging(unittest.TestCase):
    def test_two_sources_converging_with_identical_genotypes_merge(self):
        source = '>s\nAACCGG\n'
        target = '>t\nAACCGG\n'
        rows = [('s', 3, 'C', ['T'], ['0/1']), ('s', 2, 'AC', ['AT'], ['0/1'])]
        result = project(payload(source, target, vcf(['S'], rows),
                                 [whole('c', 's', 't', '+', 6)]))
        self.assertEqual(len(result['projected']), 1)
        entry = result['projected'][0]
        self.assertEqual((entry['contig'], entry['position'], entry['ref'], entry['alts']),
                         ('t', 3, 'C', ['T']))
        self.assertEqual([(s['contig'], s['position'], s['ref'], s['alts'])
                          for s in entry['sources']],
                         [('s', 2, 'AC', ['AT']), ('s', 3, 'C', ['T'])])

    def test_convergence_with_conflicting_genotypes_rejects_batch(self):
        source = '>s\nAACCGG\n'
        target = '>t\nAACCGG\n'
        rows = [('s', 3, 'C', ['T'], ['0/1']), ('s', 2, 'AC', ['AT'], ['1/1'])]
        with self.assertRaisesRegex(ValueError, 'conflicting genotypes'):
            project(payload(source, target, vcf(['S'], rows),
                            [whole('c', 's', 't', '+', 6)]))

    def test_duplicate_source_identity_rejected_even_with_reordered_alts(self):
        source = '>s\nAACCGG\n'
        target = '>t\nAACCGG\n'
        rows = [('s', 3, 'C', ['A', 'T'], ['1/2']), ('s', 3, 'C', ['T', 'A'], ['2/1'])]
        with self.assertRaisesRegex(ValueError, 'duplicate source identity'):
            project(payload(source, target, vcf(['S'], rows),
                            [whole('c', 's', 't', '+', 6)]))

    def test_same_identity_different_alt_enumeration_reports_sorted_alts(self):
        source = '>s\nAACCGG\n'
        target = '>t\nAATCGG\n'
        text = vcf(['S'], [('s', 3, 'C', ['T', 'A'], ['1/2'])])
        entry = project(payload(source, target, text,
                                [whole('c', 's', 't', '+', 6)]))['projected'][0]
        self.assertEqual(entry['sources'][0]['alts'], ['A', 'T'])


class GenotypeRules(unittest.TestCase):
    def test_phasing_separators_dropped_and_ploidy_kept(self):
        source = '>s\nAACCGG\n'
        target = '>t\nAATCGG\n'
        text = vcf(['S1', 'S2', 'S3'], [('s', 3, 'C', ['T'], ['1|0', '0/0', '1'])])
        entry = project(payload(source, target, text,
                                [whole('c', 's', 't', '+', 6)]))['projected'][0]
        self.assertEqual(entry['genotypes'], [[0, 1], [1, 1], [0]])

    def test_all_missing_diploid_genotype_stays_two_nulls(self):
        source = '>s\nAACCGG\n'
        target = '>t\nAATCGG\n'
        text = vcf(['S'], [('s', 3, 'C', ['T'], ['./.'])])
        entry = project(payload(source, target, text,
                                [whole('c', 's', 't', '+', 6)]))['projected'][0]
        self.assertEqual(entry['genotypes'], [[None, None]])

    def test_partial_missing_keeps_padding_and_sorts_known_first(self):
        source = '>s\nAACCGG\n'
        target = '>t\nAATCGG\n'
        text = vcf(['S'], [('s', 3, 'C', ['T'], ['./1'])])
        entry = project(payload(source, target, text,
                                [whole('c', 's', 't', '+', 6)]))['projected'][0]
        self.assertEqual(entry['genotypes'], [[0, None]])


class InputContract(unittest.TestCase):
    def test_exact_top_level_fields(self):
        source = '>s\nAACCGG\n'
        text = vcf(['S'], [('s', 3, 'C', ['T'], ['0/1'])])
        for mutation in ('extra', 'missing'):
            candidate = payload(source, source, text, [whole('c', 's', 's', '+', 6)])
            if mutation == 'extra':
                candidate['extra'] = True
            else:
                del candidate['chains']
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                project(candidate)

    def test_chain_field_rules(self):
        source = '>s\nAACCGG\n'
        text = vcf(['S'], [('s', 3, 'C', ['T'], ['0/1'])])
        base = whole('c', 's', 's', '+', 6)
        bad = []
        for key, value in [('strand', '*'), ('source_contig', 'missing'),
                           ('target_contig', 'missing'), ('id', 'bad/id'), ('id', '')]:
            item = json.loads(json.dumps(base))
            item[key] = value
            bad.append([item])
        for key in ('id', 'source_contig', 'target_contig', 'strand'):
            item = json.loads(json.dumps(base))
            del item[key]
            bad.append([item])
        for chains in bad:
            with self.subTest(chains=chains), self.assertRaises(ValueError):
                project(payload(source, source, text, chains))

    def test_duplicate_chain_id_rejected(self):
        source = '>s\nAACCGG\n'
        text = vcf(['S'], [('s', 3, 'C', ['T'], ['0/1'])])
        with self.assertRaisesRegex(ValueError, 'duplicate chain id'):
            project(payload(source, source, text,
                            [whole('c', 's', 's', '+', 6), whole('c', 's', 's', '+', 6)]))

    def test_chain_and_block_count_bounds(self):
        source = '>s\nAACCGG\n'
        text = vcf(['S'], [('s', 3, 'C', ['T'], ['0/1'])])
        with self.assertRaisesRegex(ValueError, '1\\.\\.128'):
            project(payload(source, source, text, []))
        with self.assertRaisesRegex(ValueError, '1\\.\\.128'):
            project(payload(source, source, text, [whole(f'c{i}', 's', 's', '+', 6)
                                                   for i in range(129)]))
        with self.assertRaisesRegex(ValueError, '1\\.\\.5000'):
            project(payload(source, source, text, [{
                'id': 'c', 'source_contig': 's', 'target_contig': 's',
                'strand': '+', 'blocks': []}]))

    def test_block_field_and_bound_rules(self):
        source = '>s\nAACCGG\n'
        text = vcf(['S'], [('s', 3, 'C', ['T'], ['0/1'])])
        base = {'source_start': 1, 'target_start': 1, 'length': 6}
        bad_blocks = [
            {'source_start': 0, 'target_start': 1, 'length': 1},
            {'source_start': True, 'target_start': 1, 'length': 1},
            {'source_start': 1, 'target_start': 1, 'length': 0},
            {'source_start': 1, 'target_start': 1, 'length': 7},
            {'source_start': 2, 'target_start': 2, 'length': 6},
            {'source_start': 1, 'target_start': 1},
        ]
        for block in bad_blocks:
            item = {'id': 'c', 'source_contig': 's', 'target_contig': 's',
                    'strand': '+', 'blocks': [block]}
            with self.subTest(block=block), self.assertRaises(ValueError):
                project(payload(source, source, text, [item]))
        with self.assertRaises(ValueError):
            project(payload(source, source, text, [{'id': 'c', 'source_contig': 's',
                                                    'target_contig': 's', 'strand': '+',
                                                    'blocks': [dict(base, extra=1)]}]))

    def test_block_ordering_rules(self):
        source = '>s\nAACCGGTT\n'
        text = vcf(['S'], [('s', 3, 'C', ['T'], ['0/1'])])
        bad = [
            # source start not strictly increasing / overlapping.
            {'strand': '+', 'blocks': [{'source_start': 5, 'target_start': 5, 'length': 2},
                                       {'source_start': 4, 'target_start': 7, 'length': 2}]},
            {'strand': '+', 'blocks': [{'source_start': 1, 'target_start': 1, 'length': 3},
                                       {'source_start': 3, 'target_start': 5, 'length': 2}]},
            # + strand target must strictly increase.
            {'strand': '+', 'blocks': [{'source_start': 1, 'target_start': 5, 'length': 2},
                                       {'source_start': 3, 'target_start': 4, 'length': 2}]},
            # - strand target must strictly decrease.
            {'strand': '-', 'blocks': [{'source_start': 1, 'target_start': 1, 'length': 2},
                                       {'source_start': 3, 'target_start': 5, 'length': 2}]},
            {'strand': '-', 'blocks': [{'source_start': 1, 'target_start': 5, 'length': 3},
                                       {'source_start': 4, 'target_start': 4, 'length': 2}]},
        ]
        for chain in bad:
            item = {'id': 'c', 'source_contig': 's', 'target_contig': 's'}
            item.update(chain)
            with self.subTest(chain=chain), self.assertRaises(ValueError):
                project(payload(source, source, text, [item]))

    def test_total_block_limit(self):
        source = '>s\n' + 'A' * 10000 + '\n'
        text = vcf(['S'], [('s', 1, 'A', ['C'], ['0/1'])])
        chain = {'id': 'c', 'source_contig': 's', 'target_contig': 's', 'strand': '+',
                 'blocks': [{'source_start': i * 2 + 1, 'target_start': i * 2 + 1, 'length': 1}
                            for i in range(5000)]}
        with self.assertRaisesRegex(ValueError, '20000'):
            project(payload(source, source, text, [dict(chain, id=f'c{i}') for i in range(5)]))

    def test_duplicate_target_interval_on_negative_strand_rejected(self):
        source = '>s\nAACCGGTT\n'
        text = vcf(['S'], [('s', 3, 'C', ['T'], ['0/1'])])
        item = {'id': 'c', 'source_contig': 's', 'target_contig': 's', 'strand': '-',
                'blocks': [{'source_start': 1, 'target_start': 5, 'length': 2},
                           {'source_start': 3, 'target_start': 5, 'length': 2}]}
        with self.assertRaises(ValueError):
            project(payload(source, source, text, [item]))

    def test_reference_must_be_independent_and_nonempty(self):
        source = '>s\nAACCGG\n'
        text = vcf(['S'], [('s', 3, 'C', ['T'], ['0/1'])])
        for target in ('', '>t\n'):
            with self.subTest(target=target), self.assertRaises(ValueError):
                project(payload(source, target, text, [whole('c', 's', 't', '+', 6)]))

    def test_source_fasta_vcf_mismatch_rejected(self):
        target = '>s\nAACCGG\n'
        text = vcf(['S'], [('s', 3, 'C', ['T'], ['0/1'])])
        with self.assertRaisesRegex(ValueError, 'reference mismatch'):
            project(payload('>s\nAATAGG\n', target, text, [whole('c', 's', 's', '+', 6)]))

    def test_header_only_vcf(self):
        text = HEADER + COLUMNS + '\tS1\tS2\n'
        result = project(payload('>s\nAACCGG\n', '>t\nAACCGG\n', text,
                                 [whole('c', 's', 't', '+', 6)]))
        self.assertEqual(result, {'samples': ['S1', 'S2'], 'source_records': 0,
                                  'projected': [], 'rejected': []})

    def test_chains_need_not_cover_every_contig(self):
        text = vcf(['S'], [('x', 1, 'A', ['C'], ['0/1'])])
        result = project(payload('>s\nAACCGG\n>x\nAACCGG\n', '>t\nAACCGG\n', text,
                                 [whole('c', 's', 't', '+', 6)]))
        self.assertEqual(result['rejected'][0]['reason'], 'unmapped')

    def test_every_header_sample_needs_a_gt_cell(self):
        text = HEADER + COLUMNS + '\tS1\tS2\n' + 's\t3\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\n'
        with self.assertRaises(ValueError):
            project(payload('>s\nAACCGG\n', '>t\nAACCGG\n', text,
                            [whole('c', 's', 't', '+', 6)]))

    def test_report_shape_and_sorting(self):
        text = vcf(['S2', 'S1'], [
            ('s', 7, 'T', ['A'], ['0/1', '1/1']),
            ('s', 3, 'C', ['T'], ['0/1', '0/0'])])
        result = project(payload('>s\nAACCGGTT\n', '>t\nAATCGGTT\n', text,
                                 [whole('c', 's', 't', '+', 8)]))
        self.assertEqual(set(result), {'samples', 'source_records', 'projected', 'rejected'})
        self.assertEqual(result['samples'], ['S1', 'S2'])
        self.assertEqual([(e['contig'], e['position']) for e in result['projected']],
                         [('t', 3), ('t', 7)])
        for entry in result['projected']:
            self.assertEqual(set(entry),
                             {'contig', 'position', 'ref', 'alts', 'genotypes', 'sources', 'chains'})
            for source in entry['sources']:
                self.assertEqual(set(source), {'contig', 'position', 'ref', 'alts'})
            for genotype in entry['genotypes']:
                self.assertEqual(len(genotype), 2)
                self.assertTrue(all(value is None or type(value) is int for value in genotype))


class PermutationInvariance(unittest.TestCase):
    def test_permuted_rows_samples_alts_chains_give_identical_result(self):
        source = '>s\nACGTTTGGCCAA\n'
        target = '>t\nTTGGCCAAACGT\n'
        plain = vcf(['S1', 'S2', 'S3'], [
            ('s', 2, 'CGTT', ['C', 'CGT'], ['1/2', '0/1', './1']),
            ('s', 8, 'GCC', ['G'], ['0/1', '1/1', '0/0'])])
        permuted = vcf(['S3', 'S1', 'S2'], [
            ('s', 8, 'GCC', ['G'], ['0/0', '0/1', '1/1']),
            ('s', 2, 'CGTT', ['CGT', 'C'], ['./2', '2/1', '0/2'])])
        chains = [whole('m', 's', 't', '-', 12)]
        self.assertEqual(
            project(payload(source, target, plain, chains)),
            project(payload(source, target, permuted, chains)))

    def test_chain_list_order_does_not_change_result(self):
        source = '>s\nAACCGG\n'
        text = vcf(['S'], [('s', 3, 'C', ['T'], ['0/1'])])
        chains = [whole('b', 's', 't', '+', 6), whole('a', 's', 't', '+', 6)]
        target = '>t\nAACCGG\n'
        first = project(payload(source, target, text, chains))
        second = project(payload(source, target, text, list(reversed(chains))))
        self.assertEqual(first, second)
        self.assertEqual(first['projected'][0]['chains'], ['a', 'b'])


class CliEndToEnd(unittest.TestCase):
    def _env(self):
        return dict(os.environ, PYTHONPATH=str(SRC))

    def _cli(self, source, output, *flags):
        command = [sys.executable, '-m', 'variant_cohort.cli',
                   '--source', str(source), '--output', str(output)]
        command.extend(flags)
        return subprocess.run(command, capture_output=True, env=self._env(), cwd=ROOT)

    def test_end_to_end_projection_report(self):
        data = payload('>s\nAACCGG\n', '>t\nAATCGG\n',
                       vcf(['S1', 'S2'], [('s', 3, 'C', ['T'], ['0/1', '1/1'])]),
                       [whole('c', 's', 't', '+', 6)])
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source = directory / 'input.json'
            output = directory / 'report.json'
            source.write_text(json.dumps(data))
            done = self._cli(source, output, '--project-assembly')
            self.assertEqual(done.returncode, 0, done.stderr)
            report = json.loads(output.read_text())
            self.assertEqual(set(report), {'schema_version', 'mode', 'source_sha256', 'result'})
            self.assertEqual(report['schema_version'], 3)
            self.assertEqual(report['mode'], 'assembly_projection')
            self.assertEqual(report['result']['samples'], ['S1', 'S2'])
            self.assertEqual(report['result']['source_records'], 1)
            entry = report['result']['projected'][0]
            self.assertEqual(entry['ref'], 'T')
            self.assertEqual(entry['alts'], ['C'])
            self.assertEqual(entry['genotypes'], [[0, 1], [0, 0]])

    def test_installed_projection_fixture_matches_library(self):
        path = ROOT / 'fixtures' / 'projection.json'
        report = run(path, project_assembly=True)
        self.assertEqual(report['schema_version'], 3)
        self.assertEqual(report['mode'], 'assembly_projection')
        self.assertEqual(report['result'],
                         project(json.loads(path.read_text())))

    def test_projection_and_stitch_are_mutually_exclusive(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            output = directory / 'report.json'
            done = self._cli(ROOT / 'fixtures' / 'projection.json', output,
                             '--project-assembly', '--stitch-phase')
            self.assertEqual(done.returncode, 2)
            self.assertFalse(output.exists())

    def test_projection_rejects_v1_and_v2_sources(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            output = directory / 'report.json'
            output.write_text('keep')
            for fixture in ('input.json', 'phase.json'):
                done = self._cli(ROOT / 'fixtures' / fixture, output, '--project-assembly')
                self.assertEqual(done.returncode, 2)
                self.assertEqual(output.read_text(), 'keep')

    def test_invalid_input_preserves_output_bytes(self):
        data = payload('>s\nAACCGG\n', '>t\nAACCGG\n',
                       vcf(['S'], [('s', 3, 'C', ['T'], ['0/1'])]),
                       [whole('c', 's', 't', '+', 6)])
        data['chains'][0]['blocks'][0]['length'] = 99
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source = directory / 'input.json'
            output = directory / 'report.json'
            source.write_text(json.dumps(data))
            original = b'{"previous": true}\n'
            output.write_bytes(original)
            done = self._cli(source, output, '--project-assembly')
            self.assertEqual(done.returncode, 2)
            self.assertEqual(output.read_bytes(), original)

    def test_input_equals_output_and_symlink_rejected(self):
        data = payload('>s\nAACCGG\n', '>t\nAACCGG\n',
                       vcf(['S'], [('s', 3, 'C', ['T'], ['0/1'])]),
                       [whole('c', 's', 't', '+', 6)])
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source = directory / 'input.json'
            source.write_text(json.dumps(data))
            original = source.read_bytes()
            link = directory / 'link.json'
            link.symlink_to(source)
            for destination in (source, link):
                done = self._cli(source, destination, '--project-assembly')
                self.assertEqual(done.returncode, 2)
                self.assertEqual(source.read_bytes(), original)

    def test_success_replaces_atomically(self):
        data = payload('>s\nAACCGG\n', '>t\nAACCGG\n',
                       vcf(['S'], [('s', 3, 'C', ['T'], ['0/1'])]),
                       [whole('c', 's', 't', '+', 6)])
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source = directory / 'input.json'
            output = directory / 'report.json'
            source.write_text(json.dumps(data))
            output.write_text('stale')
            done = self._cli(source, output, '--project-assembly')
            self.assertEqual(done.returncode, 0, done.stderr)
            self.assertEqual(json.loads(output.read_text())['schema_version'], 3)

    def test_permuted_source_report_identical_except_source_hash(self):
        source_fasta = '>s\nACGTTTGGCCAA\n'
        target_fasta = '>t\nTTGGCCAAACGT\n'
        chains = [whole('m', 's', 't', '-', 12)]
        plain = payload(source_fasta, target_fasta, vcf(['S1', 'S2'], [
            ('s', 2, 'CGTT', ['C', 'CGT'], ['1/2', '0/1']),
            ('s', 8, 'GCC', ['G'], ['0/1', '1/1'])]), chains)
        permuted = payload(source_fasta, target_fasta, vcf(['S2', 'S1'], [
            ('s', 8, 'GCC', ['G'], ['1/1', '0/1']),
            ('s', 2, 'CGTT', ['CGT', 'C'], ['0/2', '2/1'])]), list(reversed(chains)))
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            reports = []
            for name, data in [('plain', plain), ('permuted', permuted)]:
                source = directory / f'{name}.json'
                output = directory / f'{name}.report.json'
                source.write_text(json.dumps(data))
                done = self._cli(source, output, '--project-assembly')
                self.assertEqual(done.returncode, 0, done.stderr)
                reports.append(json.loads(output.read_text()))
        self.assertNotEqual(reports[0]['source_sha256'], reports[1]['source_sha256'])
        for report in reports:
            report.pop('source_sha256')
        self.assertEqual(reports[0], reports[1])

    def test_merge_conflict_exits_2_and_preserves_output(self):
        data = payload('>s\nAACCGG\n', '>t\nAACCGG\n',
                       vcf(['S'], [('s', 2, 'AC', ['AT'], ['0/1']),
                                   ('s', 3, 'C', ['T'], ['0/0'])]),
                       [whole('c', 's', 't', '+', 6)])
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source = directory / 'input.json'
            output = directory / 'report.json'
            source.write_text(json.dumps(data))
            original = b'{"previous": true}\n'
            output.write_bytes(original)
            done = self._cli(source, output, '--project-assembly')
            self.assertEqual(done.returncode, 2)
            self.assertIn(b'conflicting genotypes', done.stderr)
            self.assertEqual(output.read_bytes(), original)

    def test_legacy_modes_still_work(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            v1 = directory / 'v1.json'
            v2 = directory / 'v2.json'
            first = self._cli(ROOT / 'fixtures' / 'input.json', v1)
            second = self._cli(ROOT / 'fixtures' / 'phase.json', v2, '--stitch-phase')
            self.assertEqual(first.returncode, 0, first.stderr)
            self.assertEqual(second.returncode, 0, second.stderr)
            self.assertEqual(json.loads(v1.read_text())['schema_version'], 1)
            self.assertEqual(json.loads(v2.read_text())['schema_version'], 2)


if __name__ == '__main__':
    unittest.main()
