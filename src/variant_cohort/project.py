"""Auditable projection of reference-checked VCF records onto another assembly (schema v3)."""
import bisect
import re

from .common import number, require
from .vcf import records, reference

_ID = re.compile(r'[A-Za-z0-9_.-]{1,64}\Z')
_COMPLEMENT = str.maketrans('ACGTN', 'TGCAN')


def _revcomp(sequence):
    """Reverse complement; N stays N."""
    return sequence.translate(_COMPLEMENT)[::-1]


def _parse_chains(value, source_ref, target_ref):
    """Validate chains and precompute their contiguous block maps."""
    if not isinstance(value, list) or not 1 <= len(value) <= 128:
        raise ValueError('chains must contain 1..128 items')
    chains = []
    ids = set()
    total_blocks = 0
    for item in value:
        require(item, 'id source_contig target_contig strand blocks')
        chain_id = item['id']
        if not isinstance(chain_id, str) or not _ID.fullmatch(chain_id):
            raise ValueError('invalid chain id')
        if chain_id in ids:
            raise ValueError('duplicate chain id')
        ids.add(chain_id)
        source_contig = item['source_contig']
        target_contig = item['target_contig']
        if not isinstance(source_contig, str) or source_contig not in source_ref:
            raise ValueError('unknown source contig')
        if not isinstance(target_contig, str) or target_contig not in target_ref:
            raise ValueError('unknown target contig')
        strand = item['strand']
        if strand not in ('+', '-'):
            raise ValueError('strand must be + or -')
        blocks = item['blocks']
        if not isinstance(blocks, list) or not 1 <= len(blocks) <= 5000:
            raise ValueError('blocks must contain 1..5000 items')
        total_blocks += len(blocks)
        if total_blocks > 20000:
            raise ValueError('total blocks exceed 20000')
        source_length = len(source_ref[source_contig])
        target_length = len(target_ref[target_contig])
        parsed = []
        previous = None
        for block in blocks:
            require(block, 'source_start target_start length')
            source_start = number(block['source_start'], 1)
            target_start = number(block['target_start'], 1)
            length = number(block['length'], 1)
            if source_start + length - 1 > source_length:
                raise ValueError('source block interval out of bounds')
            if target_start + length - 1 > target_length:
                raise ValueError('target block interval out of bounds')
            if previous is not None:
                previous_start, previous_target, previous_length = previous
                if source_start < previous_start + previous_length:
                    raise ValueError('source blocks must be strictly increasing and non-overlapping')
                if strand == '+' and target_start < previous_target + previous_length:
                    raise ValueError('target blocks must be strictly increasing and non-overlapping')
                if strand == '-' and target_start + length > previous_target:
                    raise ValueError('target blocks must be strictly decreasing and non-overlapping')
            previous = (source_start, target_start, length)
            parsed.append((source_start, target_start, length))
        chains.append({'id': chain_id, 'source_contig': source_contig,
                       'target_contig': target_contig, 'strand': strand,
                       'blocks': parsed, 'starts': [block[0] for block in parsed],
                       '_target_sequence': target_ref[target_contig]})
    return chains


def _mapped_span(chain, start, end):
    """Return the mapped target low coordinate, or None without full contiguous coverage.

    The whole source interval must be covered by consecutive blocks of one chain with
    no source gap and, when several blocks are used, no target gap either.
    """
    blocks = chain['blocks']
    starts = chain['starts']
    index = bisect.bisect_right(starts, start) - 1
    if index >= 0 and blocks[index][0] + blocks[index][2] - 1 >= start:
        position = index
    else:
        position = index + 1
    if position >= len(blocks) or blocks[position][0] > start:
        return None
    cursor = start
    selected = []
    previous = None
    index = position
    while cursor <= end:
        if index >= len(blocks):
            return None
        source_start, target_start, length = blocks[index]
        if source_start > cursor:
            return None
        if previous is not None:
            previous_source, previous_target, previous_length = previous
            if source_start != previous_source + previous_length:
                return None
            if chain['strand'] == '+':
                if target_start != previous_target + previous_length:
                    return None
            elif target_start + length != previous_target:
                return None
        selected.append(blocks[index])
        previous = blocks[index]
        cursor = source_start + length
        index += 1
    first = selected[0]
    last = selected[-1]
    span_length = end - start + 1
    if chain['strand'] == '+':
        low = first[1] + (start - first[0])
    else:
        low = last[1] + last[2] - 1 - (end - last[0])
    if low < 1 or low + span_length - 1 > len(chain['_target_sequence']):
        return None
    return low


def _trim_all(position, sequences):
    """Trim a common suffix, then a common prefix, keeping each sequence non-empty."""
    sequences = list(sequences)
    while (all(len(sequence) > 1 for sequence in sequences)
           and all(sequence[-1] == sequences[0][-1] for sequence in sequences)):
        sequences = [sequence[:-1] for sequence in sequences]
    while (all(len(sequence) > 1 for sequence in sequences)
           and all(sequence[0] == sequences[0][0] for sequence in sequences)):
        sequences = [sequence[1:] for sequence in sequences]
        position += 1
    return position, sequences


def _build_candidate(chain, low, target_ref, ref, alts, order, calls):
    """Build one usable candidate, or None when the target span is incompatible."""
    sequence = target_ref[chain['target_contig']]
    span = sequence[low - 1:low - 1 + len(ref)]
    if len(span) != len(ref):
        return None
    if chain['strand'] == '+':
        transformed_ref = ref
        transformed_alts = list(alts)
    else:
        transformed_ref = _revcomp(ref)
        transformed_alts = [_revcomp(alt) for alt in alts]

    if span == transformed_ref:
        pre = [transformed_ref] + transformed_alts
    else:
        if span not in transformed_alts:
            return None
        swap = transformed_alts.index(span)
        pre = [span] + [alt for index, alt in enumerate(transformed_alts) if index != swap] + [transformed_ref]
    if len(set(pre)) != len(pre):
        raise ValueError('inconsistent transformed allele set')

    old_sequences = [transformed_ref] + transformed_alts
    position, trimmed = _trim_all(low, pre)
    new_ref = trimmed[0]
    new_alts = sorted(trimmed[1:])
    if len(set(new_alts)) != len(new_alts) or new_ref in new_alts:
        raise ValueError('inconsistent normalized allele set')
    trimmed_by_sequence = {pre[index]: trimmed[index] for index in range(len(pre))}
    final_index = {allele: index + 1 for index, allele in enumerate(new_alts)}
    gt_map = []
    for allele in old_sequences:
        trimmed_allele = trimmed_by_sequence[allele]
        gt_map.append(0 if trimmed_allele == new_ref else final_index[trimmed_allele])

    genotypes = []
    for sample in order:
        copies = calls[sample]
        remapped = [None if copy is None else gt_map[copy] for copy in copies]
        remapped.sort(key=lambda value: (value is None, value))
        genotypes.append(tuple(remapped))
    key = (chain['target_contig'], position, new_ref, tuple(new_alts), tuple(genotypes))
    return {'contig': chain['target_contig'], 'position': position,
            'ref': new_ref, 'alts': new_alts, 'genotypes': genotypes, 'key': key}


def _source_key(source):
    return (source['contig'], source['position'], source['ref'], tuple(source['alts']))


def project(data):
    """Project VCF records onto target_fasta through validated chains."""
    require(data, 'source_fasta target_fasta vcf chains')
    source_ref = reference(data['source_fasta'])
    target_ref = reference(data['target_fasta'])
    samples, rows = records(data['vcf'], source_ref)
    chains = _parse_chains(data['chains'], source_ref, target_ref)
    order = sorted(samples)

    seen_sources = set()
    projected = {}
    rejected = []
    for chrom, pos, ref, alts, cell_values in rows:
        identity = (chrom, pos, ref, tuple(sorted(alts)))
        if identity in seen_sources:
            raise ValueError('duplicate source identity')
        seen_sources.add(identity)
        source = {'contig': chrom, 'position': pos, 'ref': ref,
                  'alts': sorted(alts)}
        calls = dict(zip(samples, cell_values))
        end = pos + len(ref) - 1
        classes = {}
        covered = False
        for chain in chains:
            if chain['source_contig'] != chrom:
                continue
            low = _mapped_span(chain, pos, end)
            if low is None:
                continue
            covered = True
            candidate = _build_candidate(chain, low, target_ref, ref, alts, order, calls)
            if candidate is None:
                continue
            entry = classes.get(candidate['key'])
            if entry is None:
                classes[candidate['key']] = (candidate, {chain['id']})
            else:
                entry[1].add(chain['id'])
        if not classes:
            rejected.append({'source': source,
                             'reason': 'reference_mismatch' if covered else 'unmapped'})
            continue
        if len(classes) > 1:
            rejected.append({'source': source, 'reason': 'ambiguous'})
            continue
        key, (candidate, chain_ids) = next(iter(classes.items()))
        allele_key = key[:4]
        genotypes = key[4]
        merged = projected.get(allele_key)
        if merged is None:
            projected[allele_key] = {'genotypes': genotypes,
                                     'chains': set(chain_ids),
                                     'sources': [source]}
        else:
            if merged['genotypes'] != genotypes:
                raise ValueError('conflicting genotypes for the same projected target')
            merged['chains'].update(chain_ids)
            merged['sources'].append(source)

    out = [{'contig': contig, 'position': position, 'ref': ref, 'alts': list(alts),
            'genotypes': [list(genotype) for genotype in entry['genotypes']],
            'sources': sorted(entry['sources'], key=_source_key),
            'chains': sorted(entry['chains'])}
           for (contig, position, ref, alts), entry in projected.items()]
    out.sort(key=lambda item: (item['contig'], item['position'], item['ref'], tuple(item['alts'])))
    rejected.sort(key=lambda item: _source_key(item['source']))
    return {'samples': order, 'source_records': len(rows),
            'projected': out, 'rejected': rejected}


analyze_projection = project
project_assembly = project
