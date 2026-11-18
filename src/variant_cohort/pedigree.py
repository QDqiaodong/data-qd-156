"""Exact minimal-masking pedigree genotype consistency audit (schema v4).

Every VCF record is solved independently: samples are variable nodes, each
family is a factor node linked to its two parents and all its children, and
the resulting bipartite graph is required to be a forest.  A two-pass
(min, count) message propagation on each tree yields the exact global minimum
masking cost, the exact number of optimal complete-genotype assignments, and
the globally achievable genotype support of every sample.  Transmission
witnesses are not variables: only sample genotype tuples are counted, so
multiple inheritance paths behind one assignment never inflate the count.
"""
from .common import number, require
from .vcf import diploid_records, reference

_TABLES = {}


def _tables(n_alleles):
    """Unordered diploid genotype domain and Mendelian child sets."""
    cached = _TABLES.get(n_alleles)
    if cached is not None:
        return cached
    gts = [(a, b) for a in range(n_alleles) for b in range(a, n_alleles)]
    index = {gt: i for i, gt in enumerate(gts)}
    mend = []
    for a1, a2 in gts:
        row = []
        for b1, b2 in gts:
            children = {index[(x, y) if x <= y else (y, x)]
                        for x in (a1, a2) for y in (b1, b2)}
            row.append(sorted(children))
        mend.append(row)
    _TABLES[n_alleles] = (gts, mend)
    return _TABLES[n_alleles]


def _parse_families(value, sample_set):
    if not isinstance(value, list) or len(value) > 32:
        raise ValueError('families must contain 0..32 items')
    families = []
    seen_children = set()
    seen_pairs = set()
    for item in value:
        require(item, 'parents children')
        parents = item['parents']
        if not isinstance(parents, list) or len(parents) != 2:
            raise ValueError('parents must be exactly two sample names')
        first, second = parents
        for name in (first, second):
            if not isinstance(name, str) or name not in sample_set:
                raise ValueError('unknown parent sample')
        if first == second:
            raise ValueError('parents must be distinct samples')
        pair = frozenset((first, second))
        if pair in seen_pairs:
            raise ValueError('duplicate parent pair')
        seen_pairs.add(pair)
        children = item['children']
        if not isinstance(children, list) or not children:
            raise ValueError('children must be a non-empty array')
        if len(set(children)) != len(children):
            raise ValueError('duplicate child in family')
        for name in children:
            if not isinstance(name, str) or name not in sample_set:
                raise ValueError('unknown child sample')
            if name in pair:
                raise ValueError('child cannot be a parent of the same family')
            if name in seen_children:
                raise ValueError('sample is child of more than one family')
            seen_children.add(name)
        families.append((first, second, list(children)))
    return families


def _check_forest(families):
    """Reject any cycle in the sample/family bipartite graph."""
    parent = {}

    def find(node):
        if node not in parent:
            parent[node] = node
        root = node
        while parent[root] != root:
            root = parent[root]
        while parent[node] != root:
            parent[node], node = root, parent[node]
        return root

    for index, (p1, p2, children) in enumerate(families):
        family_node = ('family', index)
        for member in (p1, p2, *children):
            sample_node = ('sample', member)
            if find(family_node) == find(sample_node):
                raise ValueError('pedigree graph must be a forest')
            parent[find(family_node)] = find(sample_node)


def _structure(order, families):
    """Rooted forest over sample nodes 0..n-1 and family nodes n..n+F-1."""
    n = len(order)
    index = {name: i for i, name in enumerate(order)}
    fams = [(index[p1], index[p2], [index[c] for c in children])
            for p1, p2, children in families]
    adjacency = [[] for _ in range(n + len(fams))]
    for fi, (p1, p2, children) in enumerate(fams):
        node = n + fi
        for member in (p1, p2, *children):
            adjacency[node].append(member)
            adjacency[member].append(node)
    parent = [-1] * len(adjacency)
    children_of = [[] for _ in range(len(adjacency))]
    component = [-1] * len(adjacency)
    pre = []
    post = []
    components = 0
    for start in range(len(adjacency)):
        if component[start] != -1:
            continue
        component[start] = components
        queue = [start]
        pre.append(start)
        cursor = 0
        while cursor < len(queue):
            node = queue[cursor]
            cursor += 1
            for neighbor in adjacency[node]:
                if component[neighbor] != -1:
                    continue
                component[neighbor] = components
                parent[neighbor] = node
                children_of[node].append(neighbor)
                queue.append(neighbor)
                pre.append(neighbor)
        post.extend(reversed(queue))
        components += 1
    comp_samples = [[] for _ in range(components)]
    for sample in range(n):
        comp_samples[component[sample]].append(sample)
    return {'n': n, 'fams': fams, 'adjacency': adjacency, 'parent': parent,
            'children_of': children_of, 'pre': pre, 'post': post,
            'component': component, 'comp_samples': comp_samples}


def _solve(structure, unary, gts, mend):
    """Exact (min, count) two-pass propagation on the rooted forest.

    unary[v][g] is the masking cost of sample v taking genotype g.  Returns
    (minimum_cost, optimal_assignments, supports, required, possible) where
    supports[v] is the sorted list of genotype indices of sample v among
    globally optimal assignments and required/possible are sample indices.
    """
    n = structure['n']
    fams = structure['fams']
    adjacency = structure['adjacency']
    parent = structure['parent']
    children_of = structure['children_of']
    component = structure['component']
    comp_samples = structure['comp_samples']
    d = len(gts)
    pairs = d * d
    msg_cost = {}
    msg_count = {}

    def sample_message(v, target):
        cost = list(unary[v])
        count = [1] * d
        for neighbor in adjacency[v]:
            if neighbor == target:
                continue
            other_cost = msg_cost[(neighbor, v)]
            other_count = msg_count[(neighbor, v)]
            for g in range(d):
                cost[g] += other_cost[g]
                count[g] *= other_count[g]
        msg_cost[(v, target)] = cost
        msg_count[(v, target)] = count

    def family_messages(fam, targets):
        """Messages from one family factor to the given neighbor samples.

        Every message (neighbor -> fam) that exists at call time is an input;
        each target summarizes all inputs except the target's own subtree.
        Child-minimum tables are shared across targets of one call.
        """
        p1, p2, children = fams[fam - n]
        incoming = {}
        for neighbor in adjacency[fam]:
            key = (neighbor, fam)
            if key in msg_cost:
                incoming[neighbor] = (msg_cost[key], msg_count[key])
        child_min = {}
        child_ways = {}
        for child in children:
            if child not in incoming:
                continue
            child_cost, child_count = incoming[child]
            arr_cost = [0] * pairs
            arr_count = [1] * pairs
            for g1 in range(d):
                row = mend[g1]
                base = g1 * d
                for g2 in range(d):
                    best = None
                    ways = 0
                    for cg in row[g2]:
                        value = child_cost[cg]
                        if best is None or value < best:
                            best = value
                            ways = child_count[cg]
                        elif value == best:
                            ways += child_count[cg]
                    slot = base + g2
                    arr_cost[slot] = best
                    arr_count[slot] = ways
            child_min[child] = arr_cost
            child_ways[child] = arr_count
        total = [0] * pairs
        total_count = [1] * pairs
        for child in child_min:
            arr_cost = child_min[child]
            arr_count = child_ways[child]
            for slot in range(pairs):
                total[slot] += arr_cost[slot]
                total_count[slot] *= arr_count[slot]
        for target in targets:
            if target == p1 or target == p2:
                other = p2 if target == p1 else p1
                other_cost, other_count = incoming[other]
                cost = []
                count = []
                for g1 in range(d):
                    best = None
                    ways = 0
                    base = g1 * d
                    for g2 in range(d):
                        slot = base + g2
                        value = other_cost[g2] + total[slot]
                        candidate = other_count[g2] * total_count[slot]
                        if best is None or value < best:
                            best = value
                            ways = candidate
                        elif value == best:
                            ways += candidate
                    cost.append(best)
                    count.append(ways)
            else:
                first_cost, first_count = incoming[p1]
                second_cost, second_count = incoming[p2]
                if target in child_min:
                    excl = child_min[target]
                    excl_count = child_ways[target]
                else:
                    excl = [0] * pairs
                    excl_count = [1] * pairs
                cost = [None] * d
                count = [0] * d
                for g1 in range(d):
                    row = mend[g1]
                    base = g1 * d
                    head = first_cost[g1]
                    head_ways = first_count[g1]
                    for g2 in range(d):
                        slot = base + g2
                        value = head + second_cost[g2] + total[slot] - excl[slot]
                        candidate = (head_ways * second_count[g2]
                                     * (total_count[slot] // excl_count[slot]))
                        for cg in row[g2]:
                            if cost[cg] is None or value < cost[cg]:
                                cost[cg] = value
                                count[cg] = candidate
                            elif value == cost[cg]:
                                count[cg] += candidate
            msg_cost[(fam, target)] = cost
            msg_count[(fam, target)] = count

    for node in structure['post']:
        par = parent[node]
        if par == -1:
            continue
        if node < n:
            sample_message(node, par)
        else:
            family_messages(node, [par])
    for node in structure['pre']:
        if node < n:
            for child in children_of[node]:
                sample_message(node, child)
        else:
            family_messages(node, children_of[node])

    totals_cost = []
    totals_count = []
    for v in range(n):
        cost = list(unary[v])
        count = [1] * d
        for neighbor in adjacency[v]:
            other_cost = msg_cost[(neighbor, v)]
            other_count = msg_count[(neighbor, v)]
            for g in range(d):
                cost[g] += other_cost[g]
                count[g] *= other_count[g]
        totals_cost.append(cost)
        totals_count.append(count)
    comp_min = []
    comp_count = []
    for members in comp_samples:
        best = None
        ways = 0
        cost = totals_cost[members[0]]
        count = totals_count[members[0]]
        for g in range(d):
            if best is None or cost[g] < best:
                best = cost[g]
                ways = count[g]
            elif cost[g] == best:
                ways += count[g]
        comp_min.append(best)
        comp_count.append(ways)
    assignments = 1
    for ways in comp_count:
        assignments *= ways
    supports = []
    required = []
    possible = []
    for v in range(n):
        floor = comp_min[component[v]]
        support = [g for g in range(d) if totals_cost[v][g] == floor]
        supports.append(support)
        masked = [unary[v][g] > 0 for g in support]
        if all(masked):
            required.append(v)
        if any(masked):
            possible.append(v)
    return sum(comp_min), assignments, supports, required, possible


def audit(data):
    """Audit genotype consistency of a pedigree forest with minimal masking."""
    require(data, 'fasta vcf families weights')
    ref = reference(data['fasta'])
    samples, rows = diploid_records(data['vcf'], ref)
    sample_set = set(samples)
    weights = data['weights']
    if not isinstance(weights, dict) or set(weights) != sample_set:
        raise ValueError('weights must cover exactly the header samples')
    weight = {name: number(weights[name], 1, 1_000_000) for name in samples}
    families = _parse_families(data['families'], sample_set)
    _check_forest(families)
    order = sorted(samples)
    structure = _structure(order, families)
    records = []
    seen = set()
    for chrom, pos, base, alts, calls in rows:
        identity = (chrom, pos, base, tuple(sorted(alts)))
        if identity in seen:
            raise ValueError('duplicate variant record')
        seen.add(identity)
        sorted_alts = sorted(alts)
        rank = {alt: i + 1 for i, alt in enumerate(sorted_alts)}
        remap = [0] + [rank[alt] for alt in alts]
        gts, mend = _tables(len(alts) + 1)
        call_by_sample = dict(zip(samples, calls))
        unary = []
        for name in order:
            need = {}
            for copy in call_by_sample[name]:
                if copy is not None:
                    copy = remap[copy]
                    need[copy] = need.get(copy, 0) + 1
            w = weight[name]
            costs = []
            for a, b in gts:
                known = all((a == v) + (b == v) >= c for v, c in need.items())
                costs.append(0 if known else w)
            unary.append(costs)
        minimum, assignments, supports, required, possible = _solve(
            structure, unary, gts, mend)
        records.append({
            'contig': chrom,
            'position': pos,
            'ref': base,
            'alts': sorted_alts,
            'minimum_cost': minimum,
            'optimal_assignments': assignments,
            'required_masks': [order[v] for v in required],
            'possible_masks': [order[v] for v in possible],
            'genotypes': [[list(gts[g]) for g in support] for support in supports],
        })
    records.sort(key=lambda item: (item['contig'], item['position'],
                                   item['ref'], tuple(item['alts'])))
    return {'samples': order, 'source_records': len(rows), 'records': records}


audit_pedigree = audit
