"""Deterministic stitching of fragment-local VCF phase blocks (schema v2)."""
import re

from .common import require
from .model import normalize
from .vcf import phased_records, reference

_ID = re.compile(r'[A-Za-z0-9_.-]{1,64}\Z')


def _copy_value(copy, alt_index):
    """Map one GT copy onto the member indicator of a single ALT."""
    if copy is None:
        return None
    return 1 if copy == alt_index else 0


def _identities(chrom, position, base, alts, seen):
    """Normalized allele identities, rejecting collisions inside one fragment."""
    result = []
    for alt in alts:
        key = (chrom,) + normalize(position, base, alt)
        if key in seen:
            raise ValueError('normalized allele collision; reconcile source first')
        seen.add(key)
        result.append(key)
    return result


def _collect_blocks(parsed):
    blocks = []
    lookup = {}
    for fragment_id, samples, rows in parsed:
        seen = set()
        for chrom, position, base, alts, calls in rows:
            identities = _identities(chrom, position, base, alts, seen)
            for sample, (copies, phase_set) in zip(samples, calls):
                key = (sample, chrom, fragment_id, phase_set)
                block_id = lookup.get(key)
                if block_id is None:
                    block_id = len(blocks)
                    lookup[key] = block_id
                    blocks.append({'sample': sample, 'contig': chrom,
                                   'fragment_id': fragment_id, 'phase_set': phase_set,
                                   'pairs': {}})
                pairs = blocks[block_id]['pairs']
                for alt_index, identity in enumerate(identities, 1):
                    pairs[identity] = (_copy_value(copies[0], alt_index),
                                       _copy_value(copies[1], alt_index))
    return blocks


def _mark(bits, target, partners):
    """OR a partner bitset into every block named by ``bits``."""
    while bits:
        low = bits & -bits
        target[low.bit_length() - 1] |= partners
        bits ^= low


def _forced_edges(blocks):
    """Connect only pairs whose relative orientation is forced by shared alleles.

    Keep and swap compatibility are accumulated over every shared normalized
    allele identity: a direction survives only when no shared identity has a
    known disagreement.  Unknown copies never vote.  Exactly one surviving
    direction becomes a forced edge; two surviving directions leave the pair
    independent; none reject the batch.
    """
    groups = {}
    for block_id, block in enumerate(blocks):
        groups.setdefault((block['sample'], block['contig']), []).append(block_id)
    adjacency = {block_id: [] for block_id in range(len(blocks))}
    for group in groups.values():
        size = len(group)
        if size < 2:
            continue
        local = {block_id: index for index, block_id in enumerate(group)}
        index = {}
        for block_id in group:
            bit = 1 << local[block_id]
            for identity, (first, second) in blocks[block_id]['pairs'].items():
                index.setdefault(identity, []).append((bit, first, second))
        keep_bad = [0] * size
        swap_bad = [0] * size
        for members in index.values():
            zero0 = zero1 = one0 = one1 = 0
            for bit, first, second in members:
                if first == 0:
                    zero0 |= bit
                elif first == 1:
                    one0 |= bit
                if second == 0:
                    zero1 |= bit
                elif second == 1:
                    one1 |= bit
            if zero0 and one0:
                _mark(zero0, keep_bad, one0)
                _mark(one0, keep_bad, zero0)
            if zero1 and one1:
                _mark(zero1, keep_bad, one1)
                _mark(one1, keep_bad, zero1)
            if zero0 and one1:
                _mark(zero0, swap_bad, one1)
                _mark(one1, swap_bad, zero0)
            if zero1 and one0:
                _mark(zero1, swap_bad, one0)
                _mark(one0, swap_bad, zero1)

        for left in range(size):
            keep_bad[left] &= ~(1 << left)
            swap_bad[left] &= ~(1 << left)
        for left in range(size):
            left_keep = keep_bad[left]
            left_swap = swap_bad[left]
            if left_keep & left_swap:
                raise ValueError('incompatible phase blocks')
            forced_keep = left_swap & ~left_keep
            forced_swap = left_keep & ~left_swap
            if not forced_keep and not forced_swap:
                continue
            a = group[left]
            for bits, parity in ((forced_keep, 0), (forced_swap, 1)):
                while bits:
                    low = bits & -bits
                    right = low.bit_length() - 1
                    if right > left:
                        b = group[right]
                        adjacency[a].append((b, parity))
                        adjacency[b].append((a, parity))
                    bits ^= low
    return adjacency


def _orient(blocks, adjacency):
    orientation = {}
    components = []
    for root in range(len(blocks)):
        if root in orientation:
            continue
        orientation[root] = 0
        component = []
        stack = [root]
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbor, parity in adjacency[current]:
                expected = orientation[current] ^ parity
                if neighbor in orientation:
                    if orientation[neighbor] != expected:
                        raise ValueError('contradictory phase connections')
                else:
                    orientation[neighbor] = expected
                    stack.append(neighbor)
        components.append(component)
    return orientation, components


def _merge(blocks, component, orientation):
    anchor = min(component, key=lambda block_id: (
        blocks[block_id]['fragment_id'], blocks[block_id]['phase_set']))
    base = orientation[anchor]
    ordered = sorted(component, key=lambda block_id: (
        blocks[block_id]['fragment_id'], blocks[block_id]['phase_set']))
    merged = {}
    members = []
    for block_id in ordered:
        block = blocks[block_id]
        flipped = bool(orientation[block_id] ^ base)
        members.append({'fragment_id': block['fragment_id'],
                        'phase_set': block['phase_set'],
                        'flipped': flipped})
        for identity, pair in block['pairs'].items():
            if flipped:
                pair = (pair[1], pair[0])
            current = merged.get(identity)
            if current is None:
                merged[identity] = [pair[0], pair[1]]
                continue
            for slot in (0, 1):
                value = pair[slot]
                if value is None:
                    continue
                if current[slot] is None:
                    current[slot] = value
                elif current[slot] != value:
                    raise ValueError('conflicting phased allele copies')
    alleles = [{'position': identity[1], 'ref': identity[2], 'alt': identity[3],
                'copies': merged[identity]}
               for identity in sorted(merged, key=lambda item: (item[1], item[2], item[3]))]
    anchor_block = blocks[anchor]
    return {'sample': anchor_block['sample'],
            'contig': anchor_block['contig'],
            'anchor': {'fragment_id': anchor_block['fragment_id'],
                       'phase_set': anchor_block['phase_set']},
            'members': members,
            'alleles': alleles}


def stitch(data):
    """Build schema-v2 deterministic global blocks from phased fragments."""
    require(data, 'fasta fragments')
    ref = reference(data['fasta'])
    fragments = data['fragments']
    if not isinstance(fragments, list) or not 1 <= len(fragments) <= 32:
        raise ValueError('fragments must contain 1..32 items')

    parsed = []
    ids = set()
    sample_set = set()
    source_records = 0
    for item in fragments:
        require(item, 'id vcf')
        fragment_id = item['id']
        if not isinstance(fragment_id, str) or not _ID.fullmatch(fragment_id):
            raise ValueError('invalid fragment id')
        if fragment_id in ids:
            raise ValueError('duplicate fragment id')
        ids.add(fragment_id)
        samples, rows = phased_records(item['vcf'], ref)
        sample_set.update(samples)
        if len(sample_set) > 200:
            raise ValueError('sample union exceeds 200')
        source_records += len(rows)
        if source_records > 5000:
            raise ValueError('total data rows exceed 5000')
        parsed.append((fragment_id, samples, rows))

    blocks = _collect_blocks(parsed)
    adjacency = _forced_edges(blocks)
    orientation, components = _orient(blocks, adjacency)
    result_blocks = [_merge(blocks, component, orientation) for component in components]
    result_blocks.sort(key=lambda block: (
        block['sample'], block['contig'],
        block['anchor']['fragment_id'], block['anchor']['phase_set']))
    return {'samples': sorted(sample_set),
            'source_fragments': len(fragments),
            'source_records': source_records,
            'blocks': result_blocks}


analyze_phase = stitch
stitch_phase = stitch
