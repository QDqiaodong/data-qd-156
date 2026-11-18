import json
import os
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from variant_cohort.cli import run
from variant_cohort.common import read
from variant_cohort.model import qc_cohort

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'src'
HEADER = '##fileformat=VCFv4.2\n'
COLUMNS = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
FASTA_SEQ = 'ACGT' * 30
FASTA = '>chr1\n' + FASTA_SEQ + '\n'
POLICY = {
    'min_sample_call_rate': {'numerator': 1, 'denominator': 2},
    'min_site_call_rate': {'numerator': 1, 'denominator': 2},
    'max_sample_heterozygosity': {'numerator': 1, 'denominator': 2},
    'min_site_alt_count': 1,
}


def vcf_text(samples, rows, contigs=(('chr1', FASTA_SEQ),)):
    """rows: (contig, pos, ref, alts, cells) with cells aligned to samples."""
    lines = [HEADER.rstrip('\n'), COLUMNS + '\t' + '\t'.join(samples)]
    for contig, pos, ref, alts, cells in rows:
        lines.append('\t'.join([contig, str(pos), '.', ref,
                                ','.join(alts), '.', 'PASS', '.', 'GT']
                               + list(cells)))
    return '\n'.join(lines) + '\n'


def payload(vcf, policy=None, fasta=FASTA):
    return {'fasta': fasta,
            'vcf': vcf,
            'policy': dict(POLICY if policy is None else policy)}


def policy(**overrides):
    """Copy the base policy, replacing one threshold at a time."""
    value = json.loads(json.dumps(POLICY))
    value.update(overrides)
    return value


def full_fasta(contigs):
    return ''.join('>%s\n%s\n' % (name, seq) for name, seq in contigs)


def remap_cell(cell, remap):
    separator = '|' if '|' in cell else '/'
    pieces = cell.split(separator)
    return separator.join('.' if p == '.' else str(remap[int(p)])
                          for p in pieces)


def shift_representation(pos, ref, alts, sequence, direction):
    """Re-encode a record with a shared anchor base, preserving identity."""
    if direction == 'left' and pos > 1:
        extra = sequence[pos - 2]
        return pos - 1, extra + ref, [extra + alt for alt in alts]
    if direction == 'right':
        index = pos - 1 + len(ref)
        if index < len(sequence):
            extra = sequence[index]
            return pos, ref + extra, [alt + extra for alt in alts]
    return pos, ref, alts


def permute(samples, rows, sequences, seed):
    """Permute sample columns, row order, ALT order and record encoding.

    Genotypes are remapped so every allele keeps its sequence identity; the
    resulting VCF must therefore yield the same report modulo source_sha256.
    """
    rng = random.Random(seed)
    sample_order = samples[:]
    rng.shuffle(sample_order)
    row_order = list(range(len(rows)))
    rng.shuffle(row_order)
    out = []
    for index in row_order:
        contig, pos, ref, alts, cells = rows[index]
        local = random.Random(seed * 7919 + index)
        alt_perm = list(range(len(alts)))
        local.shuffle(alt_perm)
        new_alts = [alts[i] for i in alt_perm]
        remap = [0] + [1 + alt_perm.index(i) for i in range(len(alts))]
        new_pos, new_ref, new_alts = shift_representation(
            pos, ref, new_alts, sequences[contig],
            local.choice(['none', 'left', 'right', 'none']))
        by_sample = dict(zip(samples, cells))
        new_cells = [remap_cell(by_sample[name], remap)
                     for name in sample_order]
        out.append((contig, new_pos, new_ref, new_alts, new_cells))
    return vcf_text(sample_order, out, contigs=tuple(sequences.items()))


def site_index(result):
    return {(s['contig'], s['position'], s['ref'], s['alt']): s
            for s in result['sites']}


class QcLibrary(unittest.TestCase):
    def test_minimal_semantics_from_contract(self):
        rows = [('chr1', 2, 'C', ['T'], ['0/1', '0/0']),
                ('chr1', 5, 'A', ['G'], ['0/.', '1/1']),
                ('chr1', 9, 'A', ['C'], ['0/0', './.'])]
        result = qc_cohort(payload(vcf_text(['S1', 'S2'], rows),
                                   policy=policy(
                                       min_sample_call_rate={'numerator': 2,
                                                             'denominator': 3},
                                       max_sample_heterozygosity={
                                           'numerator': 1, 'denominator': 2},
                                       min_site_call_rate={'numerator': 0,
                                                           'denominator': 1},
                                       min_site_alt_count=0)))
        metrics = {m['sample']: m for m in result['sample_metrics']}
        self.assertEqual(result['samples'], ['S1', 'S2'])
        for name in ('S1', 'S2'):
            self.assertEqual(metrics[name]['total_sites'], 3)
            self.assertEqual(metrics[name]['called_sites'], 2)
            self.assertEqual(metrics[name]['call_rate'], '2/3')
            self.assertTrue(metrics[name]['retained'])
        # S1 has one heterozygous call; S2's second call is 1/1, which is a
        # homozygous ALT genotype and therefore not heterozygous.
        self.assertEqual(metrics['S1']['heterozygous_sites'], 1)
        self.assertEqual(metrics['S1']['heterozygosity'], '1/2')
        self.assertEqual(metrics['S2']['heterozygous_sites'], 0)
        self.assertEqual(metrics['S2']['heterozygosity'], '0/1')
        self.assertEqual(result['excluded_samples'], [])
        sites = site_index(result)
        # Third record is uncalled for S2, never an assumed 0/0.
        third = sites[('chr1', 9, 'A', 'C')]
        self.assertEqual(third['called_samples'], 1)
        self.assertEqual(third['call_rate'], '1/2')
        self.assertEqual(third['alt_count'], 0)
        self.assertEqual(third['alt_carriers'], 0)
        second = sites[('chr1', 5, 'A', 'G')]
        self.assertEqual((second['called_samples'], second['call_rate'],
                          second['alt_count'], second['alt_carriers']),
                         (1, '1/2', 2, 1))
        self.assertEqual(result['excluded_sites'], [])

    def test_sample_threshold_exclusion_rewrites_site_denominator(self):
        rows = [('chr1', 2, 'C', ['T'], ['0/1', './.']),
                ('chr1', 5, 'A', ['G'], ['0/.', '0/0']),
                ('chr1', 9, 'A', ['C'], ['0/0', './.'])]
        result = qc_cohort(payload(vcf_text(['S1', 'S2'], rows),
                                   policy=policy(
                                       min_sample_call_rate={'numerator': 2,
                                                             'denominator': 3},
                                       min_site_call_rate={'numerator': 0,
                                                           'denominator': 1},
                                       min_site_alt_count=0)))
        self.assertEqual(result['samples'], ['S1'])
        self.assertEqual([(e['sample'], e['reason'])
                          for e in result['excluded_samples']],
                         [('S2', 'low_call_rate')])
        self.assertEqual(result['excluded_samples'][0]['call_rate'], '1/3')
        self.assertEqual(result['excluded_samples'][0]['heterozygosity'], '0/1')
        sites = site_index(result)
        # Only S1 contributes; the two-input-sample denominator is wrong.
        second = sites[('chr1', 5, 'A', 'G')]
        self.assertEqual(second['called_samples'], 0)
        self.assertEqual(second['call_rate'], '0/1')
        self.assertEqual(second['alt_count'], 0)
        self.assertEqual(second['alt_carriers'], 0)
        self.assertEqual(second['genotypes'],
                         [{'alleles': [0, None], 'phased': False}])
        third = sites[('chr1', 9, 'A', 'C')]
        self.assertEqual((third['called_samples'], third['call_rate'],
                          third['alt_count']), (1, '1/1', 0))

    def test_zero_denominators_never_pass(self):
        text = HEADER + COLUMNS + '\tS1\tS2\n'
        result = qc_cohort(payload(text, policy=policy(
            min_sample_call_rate={'numerator': 0, 'denominator': 1},
            max_sample_heterozygosity={'numerator': 1, 'denominator': 1},
            min_site_call_rate={'numerator': 0, 'denominator': 1},
            min_site_alt_count=0)))
        self.assertEqual(result['source_records'], 0)
        self.assertEqual(result['samples'], [])
        self.assertEqual(result['sites'], [])
        self.assertEqual(result['excluded_sites'], [])
        for metric in result['sample_metrics']:
            self.assertEqual(metric['total_sites'], 0)
            self.assertEqual(metric['called_sites'], 0)
            self.assertIsNone(metric['call_rate'])
            self.assertEqual(metric['heterozygous_sites'], 0)
            self.assertIsNone(metric['heterozygosity'])
            self.assertFalse(metric['retained'])
        self.assertEqual(result['excluded_samples'], [
            {'sample': 'S1', 'reason': 'high_heterozygosity',
             'call_rate': None, 'heterozygosity': None},
            {'sample': 'S2', 'reason': 'high_heterozygosity',
             'call_rate': None, 'heterozygosity': None},
        ])

    def test_zero_called_sample_has_null_heterozygosity_and_fails(self):
        rows = [('chr1', 2, 'C', ['T'], ['./.', '0/1']),
                ('chr1', 5, 'A', ['G'], ['./.', '0/1'])]
        result = qc_cohort(payload(vcf_text(['S1', 'S2'], rows),
                                   policy=policy(
                                       min_sample_call_rate={'numerator': 0,
                                                             'denominator': 1},
                                       max_sample_heterozygosity={
                                           'numerator': 1, 'denominator': 1},
                                       min_site_call_rate={'numerator': 0,
                                                           'denominator': 1},
                                       min_site_alt_count=0)))
        metric = {m['sample']: m for m in result['sample_metrics']}['S1']
        # 0/2 is published in lowest terms.
        self.assertEqual(metric['call_rate'], '0/1')
        self.assertIsNone(metric['heterozygosity'])
        self.assertFalse(metric['retained'])
        self.assertEqual(result['samples'], ['S2'])
        self.assertEqual(result['excluded_samples'], [
            {'sample': 'S1', 'reason': 'high_heterozygosity',
             'call_rate': '0/1', 'heterozygosity': None}])

    def test_threshold_equality_is_retained(self):
        # Two called sites out of four: call rate 1/2 exactly at the minimum.
        # One heterozygous call out of two: heterozygosity 1/2 exactly at the
        # maximum.  The site matches exactly too.
        rows = [('chr1', 2, 'C', ['T'], ['0/1', '0/0']),
                ('chr1', 5, 'A', ['G'], ['0/0', '0/0']),
                ('chr1', 9, 'A', ['C'], ['./.', '0/0']),
                ('chr1', 13, 'A', ['C'], ['./.', '0/0'])]
        exact = qc_cohort(payload(vcf_text(['S1', 'S2'], rows),
                                  policy=policy(
                                      min_sample_call_rate={'numerator': 1,
                                                            'denominator': 2},
                                      max_sample_heterozygosity={
                                          'numerator': 1, 'denominator': 2},
                                      min_site_call_rate={'numerator': 1,
                                                          'denominator': 2},
                                      min_site_alt_count=0)))
        self.assertEqual(exact['samples'], ['S1', 'S2'])
        self.assertEqual(exact['excluded_samples'], [])
        self.assertEqual(exact['excluded_sites'], [])
        site = site_index(exact)[('chr1', 2, 'C', 'T')]
        self.assertEqual((site['call_rate'], site['alt_count']),
                         ('1/1', 1))
        # Site 9: called for S2 only -> call rate exactly 1/2, kept.
        exact_site = site_index(exact)[('chr1', 9, 'A', 'C')]
        self.assertEqual((exact_site['called_samples'],
                          exact_site['call_rate']), (1, '1/2'))

        stricter = qc_cohort(payload(vcf_text(['S1', 'S2'], rows),
                                     policy=policy(
                                         min_sample_call_rate={
                                             'numerator': 1,
                                             'denominator': 2},
                                         max_sample_heterozygosity={
                                             'numerator': 1,
                                             'denominator': 2},
                                         min_site_call_rate={
                                             'numerator': 2,
                                             'denominator': 3},
                                         min_site_alt_count=1)))
        # 1/2 < 2/3: site 9 now fails call rate; sites 5, 9 and 13 also have
        # alt_count 0, and the lexicographically smaller reason wins.
        self.assertEqual([(s['position'], s['reason'])
                          for s in stricter['excluded_sites']],
                         [(5, 'low_alt_count'), (9, 'low_alt_count'),
                          (13, 'low_alt_count')])

    def test_threshold_one_step_stricter_excludes(self):
        rows = [('chr1', 2, 'C', ['T'], ['0/1', '0/0']),
                ('chr1', 5, 'A', ['G'], ['0/0', '0/0']),
                ('chr1', 9, 'A', ['C'], ['./.', '0/0']),
                ('chr1', 13, 'A', ['C'], ['./.', '0/0'])]
        result = qc_cohort(payload(vcf_text(['S1', 'S2'], rows),
                                   policy=policy(
                                       min_sample_call_rate={'numerator': 3,
                                                             'denominator': 4},
                                       max_sample_heterozygosity={
                                           'numerator': 1, 'denominator': 2},
                                       min_site_call_rate={'numerator': 0,
                                                           'denominator': 1},
                                       min_site_alt_count=0)))
        self.assertEqual(result['samples'], ['S2'])
        self.assertEqual([(e['sample'], e['reason'])
                          for e in result['excluded_samples']],
                         [('S1', 'low_call_rate')])

    def test_dual_failure_reasons_are_fixed(self):
        rows = [('chr1', 2, 'C', ['T'], ['0/1', '0/1', '0/0', '0/0']),
                ('chr1', 5, 'A', ['G'], ['./.', '1/1', '0/0', '0/0']),
                ('chr1', 9, 'A', ['C'], ['./.', './.', '0/0', '0/0']),
                ('chr1', 13, 'A', ['C'], ['./.', './.', '0/0', '0/0'])]
        # S1: calls 1/4 (low) and heterozygosity 1/1 (high) -> both fail.
        # S2: calls 2/4 (passes exactly) and heterozygosity 1/2 (passes).
        # S3: all called, homozygous -> retained.  S4: all called -> retained.
        result = qc_cohort(payload(vcf_text(['S1', 'S2', 'S3', 'S4'], rows),
                                   policy=policy(
                                       min_sample_call_rate={'numerator': 1,
                                                             'denominator': 2},
                                       max_sample_heterozygosity={
                                           'numerator': 1, 'denominator': 2},
                                       min_site_call_rate={'numerator': 1,
                                                           'denominator': 2},
                                       min_site_alt_count=2)))
        self.assertEqual(result['samples'], ['S2', 'S3', 'S4'])
        self.assertEqual([(e['sample'], e['reason'])
                          for e in result['excluded_samples']],
                         [('S1', 'high_heterozygosity')])
        self.assertEqual(result['excluded_samples'][0]['call_rate'], '1/4')
        self.assertEqual(result['excluded_samples'][0]['heterozygosity'], '1/1')
        # Site 5 keeps two ALT copies from S2's 1/1 and passes.  Site 2 has
        # one ALT copy only; sites 9 and 13 have none.
        excluded = {(s['position'], s['reason']) for s in
                    result['excluded_sites']}
        self.assertEqual(excluded, {(2, 'low_alt_count'),
                                    (9, 'low_alt_count'),
                                    (13, 'low_alt_count')})
        self.assertEqual([(s['position'], s['alt']) for s in result['sites']],
                         [(5, 'G')])

    def test_dual_site_failure_uses_low_alt_count(self):
        rows = [('chr1', 2, 'C', ['T'], ['0/0', '0/0', '0/0']),
                ('chr1', 5, 'A', ['G'], ['0/1', './.', './.']),
                ('chr1', 9, 'A', ['C'], ['0/0', '0/0', '0/0'])]
        result = qc_cohort(payload(vcf_text(['S1', 'S2', 'S3'], rows),
                                   policy=policy(
                                       min_sample_call_rate={'numerator': 0,
                                                             'denominator': 1},
                                       max_sample_heterozygosity={
                                           'numerator': 1, 'denominator': 1},
                                       min_site_call_rate={'numerator': 2,
                                                           'denominator': 3},
                                       min_site_alt_count=2)))
        self.assertEqual(result['sites'], [])
        # Site 5 fails both tests (1 called of 3 retained, alt_count 1 < 2);
        # the lexicographically smaller reason wins.
        self.assertEqual([(s['position'], s['reason'])
                          for s in result['excluded_sites']],
                         [(2, 'low_alt_count'), (5, 'low_alt_count'),
                          (9, 'low_alt_count')])
        failed = [s for s in result['excluded_sites'] if s['position'] == 5][0]
        self.assertEqual((failed['called_samples'], failed['call_rate'],
                          failed['alt_count'], failed['alt_carriers']),
                         (1, '1/3', 1, 1))

    def test_phasing_separator_and_copy_order_preserved(self):
        rows = [('chr1', 2, 'C', ['T', 'G'], ['1|0', '0/1', '2|2']),
                ('chr1', 5, 'A', ['G'], ['./1', '1|.', '0/0'])]
        result = qc_cohort(payload(vcf_text(['S1', 'S2', 'S3'], rows),
                                   policy=policy(
                                       min_sample_call_rate={'numerator': 0,
                                                             'denominator': 1},
                                       max_sample_heterozygosity={
                                           'numerator': 1, 'denominator': 1},
                                       min_site_call_rate={'numerator': 0,
                                                           'denominator': 1},
                                       min_site_alt_count=0)))
        sites = site_index(result)
        # Row ALT order is T=1, G=2; each ALT becomes its own site.
        self.assertEqual(sites[('chr1', 2, 'C', 'T')]['genotypes'], [
            {'alleles': [1, 0], 'phased': True},
            {'alleles': [0, 1], 'phased': False},
            {'alleles': [0, 0], 'phased': True},
        ])
        self.assertEqual(sites[('chr1', 2, 'C', 'G')]['genotypes'], [
            {'alleles': [0, 0], 'phased': True},
            {'alleles': [0, 0], 'phased': False},
            {'alleles': [1, 1], 'phased': True},
        ])
        self.assertEqual(sites[('chr1', 5, 'A', 'G')]['genotypes'], [
            {'alleles': [None, 1], 'phased': False},
            {'alleles': [1, None], 'phased': True},
            {'alleles': [0, 0], 'phased': False},
        ])

    def test_alt_identity_remap_is_order_and_encoding_invariant(self):
        plain = vcf_text(['S1', 'S2'], [
            ('chr1', 2, 'C', ['T', 'G'], ['0/2', '1|2']),
            ('chr1', 5, 'A', ['G'], ['./1', '0/0']),
        ])
        # Same records: ALT enumeration swapped with remapped GT, and an
        # equivalent anchored encoding (left-extended REF/ALTs).
        swapped = vcf_text(['S1', 'S2'], [
            ('chr1', 2, 'C', ['G', 'T'], ['0/1', '2|1']),
            ('chr1', 4, 'TA', ['TG'], ['./1', '0/0']),
        ])
        relaxed = policy(min_sample_call_rate={'numerator': 0,
                                               'denominator': 1},
                         max_sample_heterozygosity={'numerator': 1,
                                                    'denominator': 1},
                         min_site_call_rate={'numerator': 0,
                                             'denominator': 1},
                         min_site_alt_count=0)
        baseline = qc_cohort(payload(plain, policy=relaxed))
        variant = qc_cohort(payload(swapped, policy=relaxed))
        self.assertEqual(variant, baseline)
        sites = site_index(baseline)
        self.assertEqual(sites[('chr1', 2, 'C', 'T')]['genotypes'], [
            {'alleles': [0, 0], 'phased': False},
            {'alleles': [1, 0], 'phased': True},
        ])
        self.assertEqual(sites[('chr1', 2, 'C', 'G')]['genotypes'], [
            {'alleles': [0, 1], 'phased': False},
            {'alleles': [0, 1], 'phased': True},
        ])
        self.assertEqual((sites[('chr1', 2, 'C', 'T')]['alt_count'],
                          sites[('chr1', 2, 'C', 'T')]['alt_carriers']),
                         (1, 1))
        self.assertEqual((sites[('chr1', 2, 'C', 'G')]['alt_count'],
                          sites[('chr1', 2, 'C', 'G')]['alt_carriers']),
                         (2, 2))
        self.assertEqual([(s['contig'], s['position'], s['ref'], s['alt'])
                          for s in baseline['sites']],
                         [('chr1', 2, 'C', 'G'), ('chr1', 2, 'C', 'T'),
                          ('chr1', 5, 'A', 'G')])

    def test_output_shape_and_ordering(self):
        rows = [('chr2', 3, 'G', ['A'], ['0/1', '0/0']),
                ('chr1', 8, 'T', ['C'], ['0/0', './.']),
                ('chr1', 2, 'C', ['T'], ['1/1', '0/1'])]
        contigs = (('chr1', FASTA_SEQ), ('chr2', FASTA_SEQ))
        result = qc_cohort(payload(vcf_text(['S2', 'S1'], rows,
                                            contigs=contigs),
                                   policy=policy(
                                       min_sample_call_rate={'numerator': 0,
                                                             'denominator': 1},
                                       max_sample_heterozygosity={
                                           'numerator': 1, 'denominator': 1},
                                       min_site_call_rate={'numerator': 0,
                                                           'denominator': 1},
                                       min_site_alt_count=0),
                                   fasta=full_fasta(contigs)))
        self.assertEqual(set(result), {'policy', 'source_records', 'samples',
                                       'sample_metrics', 'excluded_samples',
                                       'sites', 'excluded_sites'})
        self.assertEqual(result['source_records'], 3)
        self.assertEqual(result['samples'], ['S1', 'S2'])
        self.assertEqual([m['sample'] for m in result['sample_metrics']],
                         ['S1', 'S2'])
        for metric in result['sample_metrics']:
            self.assertEqual(set(metric), {'sample', 'total_sites',
                                           'called_sites', 'call_rate',
                                           'heterozygous_sites',
                                           'heterozygosity', 'retained'})
        for site in result['sites'] + result['excluded_sites']:
            self.assertEqual(set(site), {'contig', 'position', 'ref', 'alt',
                                         'called_samples', 'call_rate',
                                         'alt_count', 'alt_carriers',
                                         'genotypes', 'reason'}
                             if 'reason' in site else
                             {'contig', 'position', 'ref', 'alt',
                              'called_samples', 'call_rate', 'alt_count',
                              'alt_carriers', 'genotypes'})
            for genotype in site['genotypes']:
                self.assertEqual(set(genotype), {'alleles', 'phased'})
                self.assertEqual(len(genotype['alleles']), 2)
        self.assertEqual([(s['contig'], s['position'], s['ref'], s['alt'])
                          for s in result['sites']],
                         [('chr1', 2, 'C', 'T'), ('chr1', 8, 'T', 'C'),
                          ('chr2', 3, 'G', 'A')])
        self.assertEqual(result['excluded_sites'], [])
        self.assertEqual(result['policy'], policy(
            min_sample_call_rate={'numerator': 0, 'denominator': 1},
            max_sample_heterozygosity={'numerator': 1, 'denominator': 1},
            min_site_call_rate={'numerator': 0, 'denominator': 1},
            min_site_alt_count=0))

    def test_excluded_sites_sorted_by_reason_then_identity(self):
        rows = [('chr1', 13, 'A', ['C'], ['0/0', '0/0']),
                ('chr1', 5, 'A', ['G'], ['./.', '0/0']),
                ('chr1', 2, 'C', ['T'], ['./.', '0/0'])]
        result = qc_cohort(payload(vcf_text(['S1', 'S2'], rows),
                                   policy=policy(
                                       min_sample_call_rate={'numerator': 0,
                                                             'denominator': 1},
                                       max_sample_heterozygosity={
                                           'numerator': 1, 'denominator': 1},
                                       min_site_call_rate={'numerator': 1,
                                                           'denominator': 1},
                                       min_site_alt_count=1)))
        self.assertEqual([(s['position'], s['reason'])
                          for s in result['excluded_sites']],
                         [(2, 'low_alt_count'), (5, 'low_alt_count'),
                          (13, 'low_alt_count')])

    def test_empty_retained_set_has_deterministic_structure(self):
        rows = [('chr1', 2, 'C', ['T'], ['./.', '0/1']),
                ('chr1', 5, 'A', ['G'], ['./.', '0/1'])]
        result = qc_cohort(payload(vcf_text(['S1', 'S2'], rows),
                                   policy=policy(
                                       min_sample_call_rate={'numerator': 1,
                                                             'denominator': 2},
                                       max_sample_heterozygosity={
                                           'numerator': 2, 'denominator': 3},
                                       min_site_call_rate={'numerator': 0,
                                                           'denominator': 1},
                                       min_site_alt_count=0)))
        self.assertEqual(result['samples'], [])
        self.assertEqual([m['retained'] for m in result['sample_metrics']],
                         [False, False])
        self.assertEqual([(e['sample'], e['reason']) for e in
                          result['excluded_samples']],
                         [('S1', 'high_heterozygosity'),
                          ('S2', 'high_heterozygosity')])
        self.assertIsNone(result['excluded_samples'][0]['heterozygosity'])
        self.assertEqual(result['sites'], [])
        self.assertEqual(len(result['excluded_sites']), 2)
        for site in result['excluded_sites']:
            self.assertEqual(site['reason'], 'low_call_rate')
            self.assertEqual(site['called_samples'], 0)
            self.assertIsNone(site['call_rate'])
            self.assertEqual(site['alt_count'], 0)
            self.assertEqual(site['genotypes'], [])

    def test_empty_retained_set_reason_follows_dual_rule(self):
        # With a positive ALT floor both site conditions fail, so the fixed
        # lexicographic tie-break names low_alt_count, not low_call_rate.
        # Both samples have no complete call, so the retained set is empty.
        rows = [('chr1', 2, 'C', ['T'], ['./.', './1'])]
        result = qc_cohort(payload(vcf_text(['S1', 'S2'], rows),
                                   policy=policy(
                                       min_sample_call_rate={'numerator': 1,
                                                             'denominator': 1},
                                       max_sample_heterozygosity={
                                           'numerator': 1, 'denominator': 1},
                                       min_site_call_rate={'numerator': 0,
                                                           'denominator': 1},
                                       min_site_alt_count=1)))
        self.assertEqual(result['samples'], [])
        self.assertEqual([(s['position'], s['reason'])
                          for s in result['excluded_sites']],
                         [(2, 'low_alt_count')])

    def test_excluded_samples_sorted_by_reason_then_sample(self):
        # S1 and S4 fail the call-rate test; S3 fails only heterozygosity.
        # Per sample: S1 calls 4/4 with 3/4 heterozygous (high only),
        # S2 calls 4/4 with 1/4 heterozygous (retained), S3 calls 1/4 and is
        # heterozygous whenever called (both fail -> lexicographic winner),
        # S4 calls 1/4 with no heterozygous call (low call rate only).
        rows = [('chr1', 2, 'C', ['T'], ['0/1', '0/0', '0/1', './.']),
                ('chr1', 5, 'A', ['G'], ['0/1', '0/0', './.', './.']),
                ('chr1', 9, 'A', ['C'], ['0/1', '0/1', './.', '0/0']),
                ('chr1', 13, 'A', ['C'], ['0/0', '0/0', './.', './.'])]
        result = qc_cohort(payload(vcf_text(['S1', 'S2', 'S3', 'S4'], rows),
                                   policy=policy(
                                       min_sample_call_rate={'numerator': 1,
                                                             'denominator': 2},
                                       max_sample_heterozygosity={
                                           'numerator': 1, 'denominator': 2},
                                       min_site_call_rate={'numerator': 0,
                                                           'denominator': 1},
                                       min_site_alt_count=0)))
        self.assertEqual([(e['reason'], e['sample'])
                          for e in result['excluded_samples']],
                         [('high_heterozygosity', 'S1'),
                          ('high_heterozygosity', 'S3'),
                          ('low_call_rate', 'S4')])
        self.assertEqual(result['samples'], ['S2'])

    def test_permutation_invariance(self):
        contigs = (('chr1', FASTA_SEQ), ('chr2', 'TTGA' * 30))
        samples = ['S1', 'S2', 'S3', 'S4']
        rows = [
            ('chr1', 2, 'C', ['T', 'G'], ['0/1', '2|1', './.', '0/2']),
            ('chr1', 5, 'A', ['G'], ['./1', '0/0', '1/1', '1|0']),
            ('chr1', 9, 'A', ['C'], ['0/0', './.', '0/1', '1/1']),
            ('chr1', 13, 'A', ['C'], ['./.', '0/0', '0/0', '0/0']),
            ('chr2', 4, 'A', ['C', 'T'], ['1/2', '0/1', './.', '1|1']),
            ('chr2', 8, 'A', ['G'], ['1/1', './.', '0/0', './1']),
        ]
        sequences = dict(contigs)
        baseline_payload = payload(
            vcf_text(samples, rows, contigs=contigs),
            policy=policy(min_sample_call_rate={'numerator': 1,
                                                'denominator': 2},
                          max_sample_heterozygosity={'numerator': 2,
                                                     'denominator': 3},
                          min_site_call_rate={'numerator': 1,
                                              'denominator': 2},
                          min_site_alt_count=1),
            fasta=full_fasta(contigs))
        baseline = qc_cohort(baseline_payload)
        self.assertEqual(baseline['source_records'], 6)
        for seed in range(16):
            variant = qc_cohort(payload(
                permute(samples, rows, sequences, seed),
                policy=baseline_payload['policy'],
                fasta=full_fasta(contigs)))
            self.assertEqual(variant, baseline, 'seed %d' % seed)

    def test_normalized_identity_collision_rejected(self):
        text = vcf_text(['S1', 'S2'], [
            ('chr1', 2, 'C', ['T'], ['0/1', '0/0']),
            ('chr1', 1, 'AC', ['AT'], ['0/1', '0/0']),
        ])
        with self.assertRaisesRegex(ValueError, 'collision'):
            qc_cohort(payload(text))

    def test_duplicate_alt_enumeration_collision_rejected(self):
        # Two records whose ALT sets differ only in enumeration normalize to
        # the same identities: every ALT would be emitted twice.
        text = vcf_text(['S1', 'S2'], [
            ('chr1', 2, 'C', ['T', 'G'], ['0/1', '0/2']),
            ('chr1', 2, 'C', ['G', 'T'], ['0/1', '0/2']),
        ])
        with self.assertRaisesRegex(ValueError, 'collision'):
            qc_cohort(payload(text))

    def test_reference_mismatch_and_invalid_gt_rejected(self):
        bad_ref = vcf_text(['S1', 'S2'], [
            ('chr1', 2, 'A', ['T'], ['0/1', '0/0'])])
        with self.assertRaisesRegex(ValueError, 'reference mismatch'):
            qc_cohort(payload(bad_ref))
        for cell in ('0', '.', '0/1/1', '0//1', '0:1', '0/2', '-1/0', '0|1|1',
                     '0/1 0'):
            text = vcf_text(['S1', 'S2'], [
                ('chr1', 2, 'C', ['T'], [cell, '0/0'])])
            with self.assertRaises(ValueError, msg=cell):
                qc_cohort(payload(text))

    def test_outer_contract_and_policy_validation(self):
        good = payload(vcf_text(['S1', 'S2'], []))
        qc_cohort(good)
        for broken in ({'fasta': FASTA, 'vcf': good['vcf']},
                       {'fasta': FASTA, 'policy': dict(POLICY)},
                       dict(good, extra=1)):
            with self.assertRaises(ValueError):
                qc_cohort(broken)
        for bad_policy in (
                {'min_sample_call_rate': POLICY['min_sample_call_rate'],
                 'min_site_call_rate': POLICY['min_site_call_rate'],
                 'max_sample_heterozygosity':
                     POLICY['max_sample_heterozygosity']},
                dict(POLICY, min_site_alt_count=-1),
                dict(POLICY, min_site_alt_count=100001),
                dict(POLICY, min_site_alt_count=True),
                dict(POLICY, min_site_alt_count=1.0),
                dict(POLICY, min_sample_call_rate={'numerator': 0,
                                                   'denominator': 0}),
                dict(POLICY, min_sample_call_rate={'numerator': 2,
                                                   'denominator': 1}),
                dict(POLICY, min_site_call_rate={'numerator': -1,
                                                 'denominator': 2}),
                dict(POLICY, max_sample_heterozygosity={'numerator': 1,
                                                        'denominator': 2,
                                                        'extra': 0}),
                dict(POLICY, max_sample_heterozygosity={'numerator': True,
                                                        'denominator': 2}),
                dict(POLICY, extra_ratio={'numerator': 1, 'denominator': 2})):
            with self.assertRaises(ValueError, msg=str(bad_policy)):
                qc_cohort(payload(good['vcf'], policy=bad_policy))
        boundary = qc_cohort(payload(good['vcf'], policy=dict(
            POLICY, min_site_alt_count=100000)))
        self.assertEqual(boundary['policy']['min_site_alt_count'], 100000)

    def test_scale_limits(self):
        # 1 sample is below the 2..100 cohort range.
        one = HEADER + COLUMNS + '\tS1\n'
        with self.assertRaisesRegex(ValueError, 'samples'):
            qc_cohort(payload(one))
        # 101 samples exceeds the cap; 100 is accepted.
        with self.assertRaisesRegex(ValueError, 'samples'):
            qc_cohort(payload(HEADER + COLUMNS + '\t' + '\t'.join(
                'S%03d' % i for i in range(101)) + '\n'))
        qc_cohort(payload(HEADER + COLUMNS + '\t' + '\t'.join(
            'S%03d' % i for i in range(100)) + '\n'))
        # 17 ALTs exceed the per-row cap.
        many_alts = ['A' + 'C' * (index + 1) for index in range(17)]
        fasta = '>chr1\n' + 'C' * 40 + '\n'
        with self.assertRaises(ValueError):
            qc_cohort(payload(vcf_text(['S1', 'S2'], [
                ('chr1', 1, 'C', many_alts, ['0/1', '0/0'])]),
                fasta=fasta))
        # 5001 rows exceed the row cap; long reference sequence.
        long_fasta = '>chr1\n' + 'A' * 5001 + '\n'
        rows = [('chr1', pos, 'A', ['C'], ['0/1', '0/0'])
                for pos in range(1, 5001)]
        with self.assertRaisesRegex(ValueError, 'row'):
            qc_cohort(payload(vcf_text(['S1', 'S2'], rows + [
                ('chr1', 5001, 'A', ['C'], ['0/1', '0/0'])]),
                fasta=long_fasta))
        accepted = qc_cohort(payload(vcf_text(['S1', 'S2'], rows),
                                     fasta=long_fasta))
        self.assertEqual(accepted['source_records'], 5000)

    def test_vcf_size_limit(self):
        padding = '##note=' + 'N' * 200 + '\n'
        text = (padding * 10000) + HEADER + COLUMNS + '\tS1\tS2\n'
        self.assertGreater(len(text), 2_000_000)
        with self.assertRaisesRegex(ValueError, 'VCF exceeds'):
            qc_cohort(payload(text))

    def test_duplicate_json_key_and_nonfinite_rejected_at_read(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / 'dup.json'
            source.write_text('{"fasta": ">x\\nA\\n", "vcf": "", '
                              '"policy": {}, "policy": {}}')
            with self.assertRaisesRegex(ValueError, 'duplicate JSON key'):
                read(source)
            source.write_text('{"fasta": ">x\\nA\\n", "vcf": "", '
                              '"policy": NaN}')
            with self.assertRaisesRegex(ValueError, 'nonfinite'):
                read(source)


class QcCli(unittest.TestCase):
    def _env(self):
        return dict(os.environ, PYTHONPATH=str(SRC))

    def _cli(self, source, output, flag='--qc-cohort'):
        command = [sys.executable, '-m', 'variant_cohort.cli',
                   '--source', str(source), '--output', str(output)]
        if flag:
            command.append(flag)
        return subprocess.run(command, capture_output=True, env=self._env(),
                              cwd=ROOT)

    def _payload(self):
        rows = [('chr1', 2, 'C', ['T', 'G'], ['0/1', '2|1', '0/2', './.']),
                ('chr1', 5, 'A', ['G'], ['0/.', '1/1', '0/0', '1/1']),
                ('chr1', 9, 'A', ['C'], ['0/0', '.|.', '0/1', './.']),
                ('chr1', 13, 'A', ['C'], ['0/0', '0/0', '0/0', '0/1'])]
        return payload(vcf_text(['S1', 'S2', 'S3', 'S4'], rows),
                       policy=policy(
                           min_sample_call_rate={'numerator': 2,
                                                 'denominator': 3},
                           max_sample_heterozygosity={'numerator': 3,
                                                      'denominator': 4},
                           min_site_call_rate={'numerator': 1,
                                               'denominator': 2},
                           min_site_alt_count=1))

    def test_end_to_end_fixture_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / 'report.json'
            done = self._cli(ROOT / 'fixtures' / 'qc_cohort.json', output)
            self.assertEqual(done.returncode, 0, done.stderr)
            report = json.loads(output.read_text())
            self.assertEqual(set(report), {'schema_version', 'mode',
                                           'source_sha256', 'result'})
            self.assertEqual(report['schema_version'], 7)
            self.assertEqual(report['mode'], 'cohort_qc')
            self.assertRegex(report['source_sha256'], r'[0-9a-f]{64}')
            result = report['result']
            self.assertEqual(result['samples'], ['S1', 'S2', 'S3'])
            self.assertEqual(result['source_records'], 4)
            self.assertEqual([(e['sample'], e['reason'])
                              for e in result['excluded_samples']],
                             [('S4', 'low_call_rate')])
            self.assertEqual(
                [(s['position'], s['alt'], s['reason'])
                 for s in result['excluded_sites']],
                [(13, 'C', 'low_alt_count')])
            site = site_index(result)[('chr1', 2, 'C', 'T')]
            self.assertEqual((site['called_samples'], site['call_rate'],
                              site['alt_count'], site['alt_carriers']),
                             (3, '1/1', 2, 2))
            self.assertEqual(site['genotypes'], [
                {'alleles': [0, 1], 'phased': False},
                {'alleles': [0, 1], 'phased': True},
                {'alleles': [0, 0], 'phased': False},
            ])

    def test_library_and_cli_reports_match(self):
        data = self._payload()
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source = directory / 'input.json'
            output = directory / 'report.json'
            source.write_text(json.dumps(data))
            library = run(source, qc=True)
            self.assertEqual(library['result'], qc_cohort(data))
            done = self._cli(source, output)
            self.assertEqual(done.returncode, 0, done.stderr)
            self.assertEqual(json.loads(output.read_text()), library)

    def test_permuted_inputs_identical_report_except_hash(self):
        contigs = (('chr1', FASTA_SEQ),)
        samples = ['S1', 'S2', 'S3', 'S4']
        rows = [
            ('chr1', 2, 'C', ['T', 'G'], ['0/1', '2|1', './.', '0/2']),
            ('chr1', 5, 'A', ['G'], ['./1', '0/0', '1/1', '1|0']),
            ('chr1', 9, 'A', ['C'], ['0/0', './.', '0/1', '1/1']),
            ('chr1', 13, 'A', ['C'], ['./.', '0/0', '0/0', '0/0']),
        ]
        baseline = qc_cohort(payload(
            vcf_text(samples, rows, contigs=contigs),
            policy=policy(min_sample_call_rate={'numerator': 1,
                                                'denominator': 2},
                          max_sample_heterozygosity={'numerator': 2,
                                                     'denominator': 3},
                          min_site_call_rate={'numerator': 1,
                                              'denominator': 2},
                          min_site_alt_count=2)))
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            for seed in range(6):
                source = directory / ('in%d.json' % seed)
                output = directory / ('out%d.json' % seed)
                source.write_text(json.dumps(payload(
                    permute(samples, rows, dict(contigs), seed),
                    policy=policy(min_sample_call_rate={'numerator': 1,
                                                        'denominator': 2},
                                  max_sample_heterozygosity={'numerator': 2,
                                                             'denominator': 3},
                                  min_site_call_rate={'numerator': 1,
                                                      'denominator': 2},
                                  min_site_alt_count=2))))
                done = self._cli(source, output)
                self.assertEqual(done.returncode, 0, done.stderr)
                self.assertEqual(json.loads(output.read_text())['result'],
                                 baseline, 'seed %d' % seed)

    def test_failures_exit_2_and_keep_output_bytes(self):
        collision = vcf_text(['S1', 'S2'], [
            ('chr1', 2, 'C', ['T'], ['0/1', '0/0']),
            ('chr1', 1, 'AC', ['AT'], ['0/1', '0/0'])])
        bad_ref = vcf_text(['S1', 'S2'], [
            ('chr1', 2, 'A', ['T'], ['0/1', '0/0'])])
        bad_gt = vcf_text(['S1', 'S2'], [
            ('chr1', 2, 'C', ['T'], ['0/1/1', '0/0'])])
        one_sample = HEADER + COLUMNS + '\tS1\n'
        for index, text in enumerate((collision, bad_ref, bad_gt,
                                      one_sample)):
            with tempfile.TemporaryDirectory() as tmp:
                directory = Path(tmp)
                source = directory / ('bad%d.json' % index)
                output = directory / 'report.json'
                source.write_text(json.dumps(payload(text)))
                original = b'{"previous": true}\n'
                output.write_bytes(original)
                done = self._cli(source, output)
                self.assertEqual(done.returncode, 2, done.stderr)
                self.assertEqual(output.read_bytes(), original)

    def test_input_equals_output_and_symlink_rejected(self):
        data = self._payload()
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
            source.write_text(json.dumps(self._payload()))
            output.write_text('stale')
            done = self._cli(source, output)
            self.assertEqual(done.returncode, 0, done.stderr)
            self.assertEqual(json.loads(output.read_text())['schema_version'],
                             7)

    def test_mode_flags_are_mutually_exclusive(self):
        data = self._payload()
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source = directory / 'input.json'
            output = directory / 'report.json'
            source.write_text(json.dumps(data))
            for extra in ('--stitch-phase', '--project-assembly',
                          '--audit-pedigree', '--audit-replicates',
                          '--merge-cohorts'):
                command = [sys.executable, '-m', 'variant_cohort.cli',
                           '--source', str(source), '--output', str(output),
                           '--qc-cohort', extra]
                done = subprocess.run(command, capture_output=True,
                                      env=self._env(), cwd=ROOT)
                self.assertEqual(done.returncode, 2, extra)

    def test_v7_rejects_legacy_shaped_sources(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / 'report.json'
            output.write_text('keep')
            done = self._cli(ROOT / 'fixtures' / 'input.json', output)
            self.assertEqual(done.returncode, 2)
            self.assertEqual(output.read_text(), 'keep')

    def test_legacy_modes_v1_to_v6_unchanged(self):
        cases = [('input.json', None, 1, None),
                 ('phase.json', '--stitch-phase', 2, 'phase_stitch'),
                 ('projection.json', '--project-assembly', 3,
                  'assembly_projection'),
                 ('pedigree.json', '--audit-pedigree', 4, 'pedigree_audit'),
                 ('replicate.json', '--audit-replicates', 5,
                  'replicate_audit'),
                 ('cohort_merge.json', '--merge-cohorts', 6, 'cohort_merge')]
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            for name, flag, version, mode in cases:
                output = directory / ('report%d.json' % version)
                command = [sys.executable, '-m', 'variant_cohort.cli',
                           '--source', str(ROOT / 'fixtures' / name),
                           '--output', str(output)]
                if flag:
                    command.append(flag)
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
