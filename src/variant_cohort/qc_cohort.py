"""Deterministic sample/site quality screening for one cohort VCF (schema v7).

The screen is strictly staged: samples are retained first using every data row
as the denominator, then site statistics are recomputed on the retained sample
subset only.  Nothing is back-filled after a site is dropped, and a normalized
allele identity is never merged, overwritten or emitted twice.

All threshold comparisons use exact integer cross-multiplication, so boundary
equality is retained.  Ratios are published as reduced ``numerator/denominator``
strings; a zero denominator is reported as ``null`` and can never satisfy a
threshold.
"""
from fractions import Fraction

from .common import require
from .vcf import qc_records, reference


def _ratio(value, field):
    require(value, 'numerator denominator')
    numerator = value['numerator']
    denominator = value['denominator']
    if type(numerator) is not int or type(denominator) is not int:
        raise ValueError(field + ' ratio must contain integers')
    if denominator < 1 or not 0 <= numerator <= denominator:
        raise ValueError('invalid ' + field)
    return numerator, denominator


def _policy(value):
    require(value, 'min_sample_call_rate min_site_call_rate '
                   'max_sample_heterozygosity min_site_alt_count')
    sample = _ratio(value['min_sample_call_rate'], 'min_sample_call_rate')
    site = _ratio(value['min_site_call_rate'], 'min_site_call_rate')
    hetero = _ratio(value['max_sample_heterozygosity'],
                    'max_sample_heterozygosity')
    alt_count = value['min_site_alt_count']
    if type(alt_count) is not int or not 0 <= alt_count <= 100000:
        raise ValueError('min_site_alt_count must be an integer 0..100000')
    return sample, site, hetero, alt_count


def _fraction(value, denominator):
    """Reduced integer fraction string, or None for a zero denominator."""
    if denominator == 0:
        return None
    reduced = Fraction(value, denominator)
    return '%d/%d' % (reduced.numerator, reduced.denominator)


def _identity_key(entry):
    return (entry['contig'], entry['position'], entry['ref'], entry['alt'])


def qc_cohort(data):
    """Return the schema-v7 cohort QC result for one diploid cohort VCF."""
    require(data, 'fasta vcf policy')
    policy = data['policy']
    sample_ratio, site_ratio, hetero_ratio, min_alt = _policy(policy)
    ref = reference(data['fasta'])
    header_samples, rows = qc_records(data['vcf'], ref)

    # Import lazily so that ``variant_cohort.model`` can re-export qc_cohort
    # without creating an import cycle.
    from .model import normalize

    # Canonicalize every source ALT once.  Identities are global to the batch:
    # two records that normalize onto the same (contig, position, ref, alt)
    # are a collision and reject the batch instead of being merged.
    seen = set()
    prepared = []
    for chrom, position, base, alts, calls in rows:
        identities = []
        for alt in alts:
            identity = (chrom,) + normalize(position, base, alt)
            if identity in seen:
                raise ValueError('normalized allele collision; '
                                 'reconcile source first')
            seen.add(identity)
            identities.append(identity)
        prepared.append((identities, calls))

    order = sorted(header_samples)
    header_slot = {name: index for index, name in enumerate(header_samples)}
    total_sites = len(prepared)

    # ---- sample stage: every data row is a denominator site -------------
    called = [0] * len(order)
    heterozygous = [0] * len(order)
    for identities, calls in prepared:
        for index, name in enumerate(order):
            copies, _ = calls[header_slot[name]]
            if copies[0] is not None and copies[1] is not None:
                called[index] += 1
                if copies[0] != copies[1]:
                    heterozygous[index] += 1

    call_num, call_den = sample_ratio
    hetero_num, hetero_den = hetero_ratio
    retained = []
    sample_metrics = []
    excluded_samples = []
    for index, name in enumerate(order):
        call_ok = (total_sites > 0 and
                   called[index] * call_den >= total_sites * call_num)
        hetero_ok = (called[index] > 0 and
                     heterozygous[index] * hetero_den
                     <= called[index] * hetero_num)
        entry = {
            'sample': name,
            'total_sites': total_sites,
            'called_sites': called[index],
            'call_rate': _fraction(called[index], total_sites),
            'heterozygous_sites': heterozygous[index],
            'heterozygosity': _fraction(heterozygous[index], called[index]),
            'retained': call_ok and hetero_ok,
        }
        sample_metrics.append(entry)
        if entry['retained']:
            retained.append(name)
        else:
            reasons = []
            if not call_ok:
                reasons.append('low_call_rate')
            if not hetero_ok:
                reasons.append('high_heterozygosity')
            excluded_samples.append({
                'sample': name,
                'reason': min(reasons),
                'call_rate': entry['call_rate'],
                'heterozygosity': entry['heterozygosity'],
            })
    excluded_samples.sort(key=lambda item: (item['reason'], item['sample']))

    # ---- site stage: recomputed on the retained samples only ------------
    retained_count = len(retained)
    site_num, site_den = site_ratio
    sites = []
    excluded_sites = []
    for identities, calls in prepared:
        for alt_index, identity in enumerate(identities, 1):
            called_samples = 0
            alt_count = 0
            alt_carriers = 0
            genotypes = []
            for name in retained:
                copies, phased = calls[header_slot[name]]
                complete = copies[0] is not None and copies[1] is not None
                if complete:
                    called_samples += 1
                carried = 0
                alleles = []
                for copy in copies:
                    if copy is None:
                        alleles.append(None)
                    elif copy == alt_index:
                        alleles.append(1)
                        carried += 1
                    else:
                        alleles.append(0)
                alt_count += carried
                if complete and carried:
                    alt_carriers += 1
                genotypes.append({'alleles': alleles, 'phased': phased})
            if alt_carriers > called_samples or alt_count < alt_carriers:
                raise ValueError('internal genotype count mismatch')

            contig, position, ref_base, alt_base = identity
            call_ok = (retained_count > 0 and
                       called_samples * site_den
                       >= retained_count * site_num)
            alt_ok = alt_count >= min_alt
            entry = {
                'contig': contig,
                'position': position,
                'ref': ref_base,
                'alt': alt_base,
                'called_samples': called_samples,
                'call_rate': _fraction(called_samples, retained_count),
                'alt_count': alt_count,
                'alt_carriers': alt_carriers,
                'genotypes': genotypes,
            }
            if call_ok and alt_ok:
                sites.append(entry)
            else:
                reasons = []
                if not call_ok:
                    reasons.append('low_call_rate')
                if not alt_ok:
                    reasons.append('low_alt_count')
                entry['reason'] = min(reasons)
                excluded_sites.append(entry)
    sites.sort(key=_identity_key)
    excluded_sites.sort(key=_identity_key)

    return {
        'policy': policy,
        'source_records': total_sites,
        'samples': retained,
        'sample_metrics': sample_metrics,
        'excluded_samples': excluded_samples,
        'sites': sites,
        'excluded_sites': excluded_sites,
    }
