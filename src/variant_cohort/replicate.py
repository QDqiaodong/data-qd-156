"""Offline pairwise replicate-candidate screening (schema v5).

The audit is a data-quality screen only.  It never infers relatedness,
identity or clinical meaning, produces no clusters and never takes a
transitive closure over pairwise candidate edges.

Each diploid GT is mapped onto the row's REF/ALT sequences and turned into
an unordered multiset of two sequence alleles, so phasing separators and the
source ALT enumeration order cannot change genotype identity.  Calls with a
missing copy are accepted but never participate in pairwise comparisons on
that record.
"""
from .common import require
from .vcf import reference, replicate_records

_MAX_EXAMPLES = 8


def _ratio(value):
    """Validate the max_discordance {numerator, denominator} object."""
    require(value, 'numerator denominator')
    numerator = value['numerator']
    denominator = value['denominator']
    if type(numerator) is not int or type(denominator) is not int:
        raise ValueError('discordance ratio must be integers')
    if denominator < 1 or not 0 <= numerator <= denominator:
        raise ValueError('invalid discordance ratio')
    return numerator, denominator


def audit_replicates(data):
    """Return the pairwise replicate-candidate screening result."""
    require(data, 'fasta vcf minimum_informative_sites max_discordance')
    minimum = data['minimum_informative_sites']
    if type(minimum) is not int or not 1 <= minimum <= 5000:
        raise ValueError('minimum_informative_sites must be 1..5000')
    numerator, denominator = _ratio(data['max_discordance'])
    ref = reference(data['fasta'])
    header_samples, rows = replicate_records(data['vcf'], ref)
    order = sorted(header_samples)
    n = len(order)
    # calls arrive aligned with the VCF header column order, which need not
    # be lexicographic; map header positions onto sorted pair positions.
    sorted_position = {name: i for i, name in enumerate(order)}
    header_to_sorted = [sorted_position[name] for name in header_samples]

    # Canonicalize every record once: identity uses lexicographically sorted
    # ALTs, so rows differing only in ALT enumeration share one identity.
    prepared = []
    seen = set()
    for chrom, pos, base, alts, calls in rows:
        sorted_alts = sorted(alts)
        identity = (chrom, pos, base, tuple(sorted_alts))
        if identity in seen:
            raise ValueError('duplicate variant record')
        seen.add(identity)
        sequences = [base] + alts
        genotypes = [None] * n
        complete_flags = [False] * n
        for header_pos, copies in enumerate(calls):
            target = header_to_sorted[header_pos]
            if any(copy is None for copy in copies):
                genotypes[target] = None
                complete_flags[target] = False
            else:
                genotypes[target] = tuple(
                    sorted(sequences[copy] for copy in copies))
                complete_flags[target] = True
        prepared.append((identity, sorted_alts, genotypes, complete_flags))

    # Iterate records in canonical identity order so retained mismatch
    # examples are already ordered by (contig, position, ref, alts).
    prepared.sort(key=lambda item: item[0])

    pairs = []
    for left_pos in range(n):
        for right_pos in range(left_pos + 1, n):
            left = order[left_pos]
            right = order[right_pos]
            complete_sites = 0
            informative_sites = 0
            concordant_sites = 0
            discordant_sites = 0
            examples = []
            for identity, sorted_alts, genotypes, complete_flags in prepared:
                if not (complete_flags[left_pos] and complete_flags[right_pos]):
                    continue
                complete_sites += 1
                left_gt = genotypes[left_pos]
                right_gt = genotypes[right_pos]
                if left_gt == (identity[2], identity[2]) and \
                        right_gt == (identity[2], identity[2]):
                    continue
                informative_sites += 1
                if left_gt == right_gt:
                    concordant_sites += 1
                else:
                    discordant_sites += 1
                    if len(examples) < _MAX_EXAMPLES:
                        chrom, pos, base, _ = identity
                        examples.append({
                            'contig': chrom,
                            'position': pos,
                            'ref': base,
                            'alts': list(sorted_alts),
                            'left_gt': list(left_gt),
                            'right_gt': list(right_gt),
                        })
            if informative_sites != concordant_sites + discordant_sites:
                raise ValueError('internal count mismatch')
            if informative_sites < minimum:
                status = 'insufficient_evidence'
            elif discordant_sites * denominator <= informative_sites * numerator:
                status = 'candidate'
            else:
                status = 'different'
            pairs.append({
                'left': left,
                'right': right,
                'complete_sites': complete_sites,
                'informative_sites': informative_sites,
                'concordant_sites': concordant_sites,
                'discordant_sites': discordant_sites,
                'status': status,
                'mismatch_examples': examples,
                'omitted_mismatches': discordant_sites - len(examples),
            })

    return {
        'samples': order,
        'source_records': len(rows),
        'parameters': {
            'minimum_informative_sites': minimum,
            'max_discordance': {'numerator': numerator,
                                'denominator': denominator},
        },
        'pairs': pairs,
    }
