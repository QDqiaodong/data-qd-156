"""Deterministic joint genotype matrix across disjoint VCF cohorts (schema v6).

Batches whose sample names never overlap are aligned on the exact site
identity ``(contig, position)``.  Records are never trimmed, left-aligned or
merged across different coordinates; a site only promises REF compatibility
for records sharing the exact same coordinate.  Every genotype copy is
remapped by allele sequence identity onto the lexicographic union of the
source ALT sequences, so source enumeration orders cannot leak into the
joint matrix.
"""
import re

from .common import require
from .vcf import cohort_records, reference

_ID = re.compile(r'[A-Za-z0-9_.-]{1,64}\Z')


def merge_cohorts(data):
    """Merge several disjoint-sample VCF batches into one genotype matrix."""
    require(data, 'fasta cohorts')
    ref = reference(data['fasta'])
    cohorts = data['cohorts']
    if not isinstance(cohorts, list) or not 1 <= len(cohorts) <= 32:
        raise ValueError('cohorts must contain 1..32 items')

    # First pass: validate every batch, enforce global sample and row budgets
    # and collect every source record keyed by its exact site identity.
    cohort_ids = set()
    all_samples = set()
    total_records = 0
    sites = {}
    for entry in cohorts:
        require(entry, 'id vcf')
        cohort_id = entry['id']
        if not isinstance(cohort_id, str) or not _ID.fullmatch(cohort_id):
            raise ValueError('invalid cohort id')
        if cohort_id in cohort_ids:
            raise ValueError('duplicate cohort id')
        cohort_ids.add(cohort_id)
        samples, rows = cohort_records(entry['vcf'], ref)
        overlap = all_samples.intersection(samples)
        if overlap:
            raise ValueError('sample name shared across cohorts: '
                             + min(overlap))
        all_samples.update(samples)
        if len(all_samples) > 200:
            raise ValueError('global sample union exceeds 200')
        total_records += len(rows)
        if total_records > 5000:
            raise ValueError('total data rows exceed 5000')
        local_sites = set()
        for chrom, pos, base, alts, calls in rows:
            identity = (chrom, pos)
            if identity in local_sites:
                raise ValueError('duplicate site within one cohort: '
                                 + '%s:%d' % identity)
            local_sites.add(identity)
            site = sites.get(identity)
            if site is None:
                site = {'ref': base, 'alt_set': set(), 'sources': []}
                sites[identity] = site
            elif site['ref'] != base:
                raise ValueError('REF conflict at shared site: '
                                 + '%s:%d' % identity)
            site['alt_set'].update(alts)
            site['sources'].append((cohort_id, samples, alts, calls))

    order = sorted(all_samples)
    sample_position = {name: i for i, name in enumerate(order)}
    out_sites = []
    for identity in sorted(sites):
        chrom, pos = identity
        site = sites[identity]
        joint_alts = sorted(site['alt_set'])
        genotypes = [None] * len(order)
        source_cohorts = []
        for cohort_id, samples, alts, calls in site['sources']:
            source_cohorts.append(cohort_id)
            # Sequence identity -> joint 1-based ALT index; REF stays 0 and
            # missing copies stay None.  Source numeric indexes are never
            # copied into the joint matrix.
            joint_index = {alt: i + 1 for i, alt in enumerate(joint_alts)}
            remap = [0] + [joint_index[alt] for alt in alts]
            for sample, (copies, phased) in zip(samples, calls):
                target = sample_position[sample]
                try:
                    alleles = [None if copy is None else remap[copy]
                               for copy in copies]
                except IndexError:
                    raise ValueError('internal GT remapping conflict')
                if genotypes[target] is not None:
                    raise ValueError('internal duplicate genotype slot')
                genotypes[target] = {'alleles': alleles, 'phased': phased}
        out_sites.append({
            'contig': chrom,
            'position': pos,
            'ref': site['ref'],
            'alts': joint_alts,
            'source_cohorts': sorted(source_cohorts),
            'genotypes': genotypes,
        })

    return {
        'samples': order,
        'source_cohorts': sorted(cohort_ids),
        'source_records': total_records,
        'sites': out_sites,
    }
