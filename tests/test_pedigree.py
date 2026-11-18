import itertools
import json
import os
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from variant_cohort.cli import run
from variant_cohort.pedigree import audit

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'src'
HEADER = '##fileformat=VCFv4.2\n'
COLUMNS = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
FASTA = '>chr1\n' + 'ACGT' * 60 + '\n'


def vcf(samples, rows, header=HEADER):
    """Build a GT-only VCF; rows are (pos, ref, alts, cells) on chr1."""
    lines = [header.rstrip('\n'), COLUMNS + '\t' + '\t'.join(samples)]
    for pos, ref, alts, cells in rows:
        lines.append('\t'.join(['chr1', str(pos), '.', ref, ','.join(alts),
                                '.', 'PASS', '.', 'GT'] + list(cells)))
    return '\n'.join(lines) + '\n'


def payload(samples, rows, families, weights=None, fasta=FASTA):
    """rows are (pos, ref, alts, obs); obs maps sample -> known allele indexes."""
    if weights is None:
        weights = {s: 1 for s in samples}
    cells_rows = []
    for pos, ref, alts, obs in rows:
        cells = []
        for sample in samples:
            known = list(obs.get(sample, ()))
            cells.append('/'.join([str(a) for a in known] + ['.'] * (2 - len(known))))
        cells_rows.append((pos, ref, alts, cells))
    fams = [{'parents': [p1, p2], 'children': list(children)}
            for p1, p2, children in families]
    return {'fasta': fasta, 'vcf': vcf(samples, cells_rows),
            'families': fams, 'weights': weights}


def brute_force(samples, obs, families, weights, n_alleles):
    """Independent exhaustive enumeration of all genotype assignments."""
    gts = [(a, b) for a in range(n_alleles) for b in range(a, n_alleles)]
    order = sorted(samples)

    def matches(gt, known):
        return all(gt.count(v) >= list(known).count(v) for v in set(known))

    best = None
    count = 0
    optimal = []
    for assignment in itertools.product(gts, repeat=len(order)):
        value = dict(zip(order, assignment))
        ok = True
        for p1, p2, children in families:
            pool = {tuple(sorted((x, y)))
                    for x in set(value[p1]) for y in set(value[p2])}
            if any(value[c] not in pool for c in children):
                ok = False
                break
        if not ok:
            continue
        cost = sum(0 if matches(value[s], obs.get(s, ())) else weights[s]
                   for s in order)
        if best is None or cost < best:
            best, count, optimal = cost, 1, [value]
        elif cost == best:
            count += 1
            optimal.append(value)
    support = {s: set() for s in order}
    unmasked = {s: False for s in order}
    masked = {s: False for s in order}
    for value in optimal:
        for s in order:
            support[s].add(value[s])
            if matches(value[s], obs.get(s, ())):
                unmasked[s] = True
            else:
                masked[s] = True
    return {'minimum_cost': best, 'optimal_assignments': count,
            'support': {s: sorted(support[s]) for s in order},
            'required_masks': sorted(s for s in order if not unmasked[s]),
            'possible_masks': sorted(s for s in order if masked[s])}


def check_against_brute(case, samples, rows, families, weights):
    """Compare audit() with brute_force() row by row (alts given sorted)."""
    result = audit(payload(samples, rows, families, weights))
    case.assertEqual(len(result['records']), len(rows))
    by_position = {record['position']: record for record in result['records']}
    for pos, ref, alts, obs in rows:
        expected = brute_force(samples, obs, families, weights, len(alts) + 1)
        record = by_position[pos]
        case.assertEqual(record['minimum_cost'], expected['minimum_cost'], pos)
        case.assertEqual(record['optimal_assignments'],
                         expected['optimal_assignments'], pos)
        case.assertEqual(record['required_masks'], expected['required_masks'], pos)
        case.assertEqual(record['possible_masks'], expected['possible_masks'], pos)
        support = dict(zip(result['samples'], record['genotypes']))
        case.assertEqual(support,
                         {s: [list(gt) for gt in expected['support'][s]]
                          for s in result['samples']}, pos)
    return result


def random_forest(samples, rng, max_families):
    """Random valid family forest via incremental union-find."""
    parent = {}

    def find(node):
        if node not in parent:
            parent[node] = node
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    families = []
    used_children = set()
    used_pairs = set()
    target = rng.randint(0, min(max_families, len(samples)))
    attempts = 0
    while len(families) < target and attempts < 200:
        attempts += 1
        if len(samples) < 3:
            break
        p1, p2 = rng.sample(samples, 2)
        pair = frozenset((p1, p2))
        if pair in used_pairs:
            continue
        children = tuple(rng.sample(samples, rng.randint(1, min(3, len(samples)))))
        if any(c in used_children or c in (p1, p2) for c in children):
            continue
        saved = dict(parent)
        node = ('family', len(families))
        ok = True
        for member in [('sample', p1), ('sample', p2)] + [('sample', c) for c in children]:
            if find(node) == find(member):
                ok = False
                break
            parent[find(node)] = find(member)
        if not ok:
            parent.clear()
            parent.update(saved)
            continue
        families.append((p1, p2, children))
        used_children.update(children)
        used_pairs.add(pair)
    return families


class GivenExamples(unittest.TestCase):
    def test_spec_tie_example(self):
        samples = ['P', 'Q', 'C']
        result = audit(payload(samples, [(1, 'A', ['C'],
                                          {'P': (0, 0), 'Q': (0, 0), 'C': (0, 1)})],
                               [('P', 'Q', ('C',))]))
        self.assertEqual(result['samples'], ['C', 'P', 'Q'])
        self.assertEqual(result['source_records'], 1)
        record = result['records'][0]
        self.assertEqual((record['contig'], record['position'], record['ref'], record['alts']),
                         ('chr1', 1, 'A', ['C']))
        self.assertEqual(record['minimum_cost'], 1)
        self.assertEqual(record['optimal_assignments'], 5)
        self.assertEqual(record['required_masks'], [])
        self.assertEqual(record['possible_masks'], ['C', 'P', 'Q'])
        support = dict(zip(result['samples'], record['genotypes']))
        self.assertEqual(support['P'], [[0, 0], [0, 1], [1, 1]])
        self.assertEqual(support['Q'], [[0, 0], [0, 1], [1, 1]])
        self.assertEqual(support['C'], [[0, 0], [0, 1]])

    def test_transmission_witnesses_counted_once(self):
        samples = ['P', 'Q', 'C']
        result = audit(payload(samples, [(1, 'A', ['C'],
                                          {'P': (0, 1), 'Q': (0, 1), 'C': (0, 1)})],
                               [('P', 'Q', ('C',))]))
        record = result['records'][0]
        # Two inheritance paths behind C=0/1 must not double the assignment.
        self.assertEqual(record['minimum_cost'], 0)
        self.assertEqual(record['optimal_assignments'], 1)
        self.assertEqual(record['required_masks'], [])
        self.assertEqual(record['possible_masks'], [])

    def test_het_child_missing_has_three_assignments(self):
        samples = ['P', 'Q', 'C']
        result = audit(payload(samples, [(1, 'A', ['C'],
                                          {'P': (0, 1), 'Q': (0, 1), 'C': ()})],
                               [('P', 'Q', ('C',))]))
        record = result['records'][0]
        # C may be 0/0, 0/1 or 1/1; the 0/1 case has two witnesses, counted once.
        self.assertEqual(record['minimum_cost'], 0)
        self.assertEqual(record['optimal_assignments'], 3)
        support = dict(zip(result['samples'], record['genotypes']))
        self.assertEqual(support['C'], [[0, 0], [0, 1], [1, 1]])


class CrossGeneration(unittest.TestCase):
    def test_shared_missing_sample_forces_nonzero_joint_cost(self):
        samples = ['A', 'B', 'X', 'D', 'E']
        obs = {'A': (0, 0), 'B': (0, 0), 'X': (), 'D': (1, 1), 'E': (1, 1)}
        families = [('A', 'B', ('X',)), ('X', 'D', ('E',))]
        result = audit(payload(samples, [(1, 'A', ['C'], obs)], families))
        record = result['records'][0]
        # Each family alone admits cost zero; jointly impossible.
        self.assertEqual(record['minimum_cost'], 1)
        self.assertEqual(record['optimal_assignments'], 5)
        support = dict(zip(result['samples'], record['genotypes']))
        self.assertEqual(support['A'], [[0, 0], [0, 1], [1, 1]])
        self.assertEqual(support['B'], [[0, 0], [0, 1], [1, 1]])
        self.assertEqual(support['X'], [[0, 0], [0, 1]])
        self.assertEqual(support['D'], [[1, 1]])
        self.assertEqual(support['E'], [[0, 1], [1, 1]])
        self.assertEqual(record['required_masks'], [])
        self.assertEqual(record['possible_masks'], ['A', 'B', 'E'])
        check_against_brute(self, samples, [(1, 'A', ['C'], obs)], families,
                            {s: 1 for s in samples})

    def test_each_family_separately_has_zero_cost(self):
        obs = {'A': (0, 0), 'B': (0, 0), 'X': (), 'D': (1, 1), 'E': (1, 1)}
        for family, kept in [([('A', 'B', ('X',))], ['A', 'B', 'X']),
                             ([('X', 'D', ('E',))], ['X', 'D', 'E'])]:
            sub_obs = {s: obs[s] for s in kept}
            result = audit(payload(kept, [(1, 'A', ['C'], sub_obs)], family))
            self.assertEqual(result['records'][0]['minimum_cost'], 0)


class Siblings(unittest.TestCase):
    def test_two_children_jointly_constrain_parents(self):
        samples = ['P', 'Q', 'C1', 'C2']
        obs = {'P': (), 'Q': (), 'C1': (0, 0), 'C2': (1, 1)}
        families = [('P', 'Q', ('C1', 'C2'))]
        result = audit(payload(samples, [(1, 'A', ['C'], obs)], families))
        record = result['records'][0]
        # Only P=Q=0/1 can produce both a 0/0 and a 1/1 child.
        self.assertEqual(record['minimum_cost'], 0)
        self.assertEqual(record['optimal_assignments'], 1)
        support = dict(zip(result['samples'], record['genotypes']))
        self.assertEqual(support['P'], [[0, 1]])
        self.assertEqual(support['Q'], [[0, 1]])
        self.assertEqual(support['C1'], [[0, 0]])
        self.assertEqual(support['C2'], [[1, 1]])

    def test_sibling_conflict_spreads_mask_choice(self):
        samples = ['P', 'Q', 'C1', 'C2']
        obs = {'P': (0, 0), 'Q': (0, 0), 'C1': (0, 0), 'C2': (1, 1)}
        families = [('P', 'Q', ('C1', 'C2'))]
        check_against_brute(self, samples, [(1, 'A', ['C'], obs)], families,
                            {s: 1 for s in samples})
        result = audit(payload(samples, [(1, 'A', ['C'], obs)], families))
        record = result['records'][0]
        self.assertEqual(record['minimum_cost'], 1)
        # C2=1/1 is impossible under 0/0 x 0/0 unless a parent is masked.
        self.assertIn('C2', record['possible_masks'])


class TiedCostsAndCounting(unittest.TestCase):
    def test_isolated_missing_sample_counts_all_genotypes(self):
        result = audit(payload(['S'], [(1, 'A', ['C'], {'S': ()})], []))
        record = result['records'][0]
        self.assertEqual(record['minimum_cost'], 0)
        self.assertEqual(record['optimal_assignments'], 3)
        self.assertEqual(record['genotypes'], [[[0, 0], [0, 1], [1, 1]]])
        self.assertEqual(record['required_masks'], [])
        self.assertEqual(record['possible_masks'], [])

    def test_isolated_missing_samples_multiply(self):
        samples = ['S1', 'S2', 'S3']
        result = audit(payload(samples, [(1, 'A', ['C'], {})], []))
        self.assertEqual(result['records'][0]['optimal_assignments'], 27)

    def test_isolated_missing_multiallelic_combinatorics(self):
        result = audit(payload(['S1', 'S2'], [(1, 'A', ['C', 'G'], {})], []))
        # Six unordered diploid genotypes over three alleles, independently.
        self.assertEqual(result['records'][0]['optimal_assignments'], 36)
        self.assertEqual(result['records'][0]['genotypes'],
                         [[[0, 0], [0, 1], [0, 2], [1, 1], [1, 2], [2, 2]]] * 2)

    def test_partial_missing_support_and_count(self):
        result = audit(payload(['S'], [(1, 'A', ['C'], {'S': (0,)})], []))
        record = result['records'][0]
        self.assertEqual(record['optimal_assignments'], 2)
        self.assertEqual(record['genotypes'], [[[0, 0], [0, 1]]])

    def test_independent_families_counts_multiply(self):
        samples = ['P1', 'Q1', 'C1', 'P2', 'Q2', 'C2']
        obs = {'P1': (0, 0), 'Q1': (0, 0), 'C1': (0, 1),
               'P2': (0, 0), 'Q2': (0, 0), 'C2': (0, 1)}
        families = [('P1', 'Q1', ('C1',)), ('P2', 'Q2', ('C2',))]
        result = audit(payload(samples, [(1, 'A', ['C'], obs)], families))
        record = result['records'][0]
        self.assertEqual(record['minimum_cost'], 2)
        self.assertEqual(record['optimal_assignments'], 25)

    def test_weights_change_optimum_and_masks(self):
        samples = ['P', 'Q', 'C']
        obs = {'P': (0, 0), 'Q': (0, 0), 'C': (1, 1)}
        families = [('P', 'Q', ('C',))]
        # Masking the child costs 5, masking a parent 1: the optimum masks both
        # parents (one masked parent still cannot supply two allele-1 copies).
        weights = {'P': 1, 'Q': 1, 'C': 5}
        result = audit(payload(samples, [(1, 'A', ['C'], obs)], families, weights))
        record = result['records'][0]
        self.assertEqual(record['minimum_cost'], 2)
        self.assertEqual(record['optimal_assignments'], 4)
        self.assertEqual(record['required_masks'], ['P', 'Q'])
        self.assertEqual(record['possible_masks'], ['P', 'Q'])
        support = dict(zip(result['samples'], record['genotypes']))
        self.assertEqual(support['C'], [[1, 1]])
        self.assertEqual(support['P'], [[0, 1], [1, 1]])
        check_against_brute(self, samples, [(1, 'A', ['C'], obs)], families, weights)


class GenotypeMatching(unittest.TestCase):
    def test_homozygous_call_does_not_match_heterozygous_genotype(self):
        # 1/1 observation is incompatible with the 0/1 genotype.
        result = audit(payload(['S'], [(1, 'A', ['C'], {'S': (1, 1)})], []))
        record = result['records'][0]
        self.assertEqual(record['minimum_cost'], 0)
        self.assertEqual(record['genotypes'], [[[1, 1]]])
        self.assertEqual(record['optimal_assignments'], 1)

    def test_partial_missing_multiplicity(self):
        # 0/. matches 0/0 and 0/1 but not 0/2 or 1/1.
        result = audit(payload(['S'], [(1, 'A', ['C', 'G'], {'S': (0,)})], []))
        record = result['records'][0]
        self.assertEqual(record['genotypes'], [[[0, 0], [0, 1], [0, 2]]])
        self.assertEqual(record['optimal_assignments'], 3)

    def test_phased_diploid_cell_is_unphased_observation(self):
        samples = ['P', 'Q', 'C']
        data = payload(samples, [(1, 'A', ['C'], {})], [('P', 'Q', ('C',))])
        data['vcf'] = data['vcf'].replace('./.\t./.\t./.', '0|1\t0|1\t0|1')
        record = audit(data)['records'][0]
        self.assertEqual(record['minimum_cost'], 0)
        self.assertEqual(record['optimal_assignments'], 1)

    def test_multiallelic_family_against_brute_force(self):
        samples = ['P', 'Q', 'C']
        obs = {'P': (1, 2), 'Q': (0,), 'C': (2, 2)}
        families = [('P', 'Q', ('C',))]
        check_against_brute(self, samples, [(1, 'A', ['C', 'G', 'T'], obs)],
                            families, {'P': 2, 'Q': 3, 'C': 1})


class AltOrderAndRenumbering(unittest.TestCase):
    def test_alt_permutation_with_remapped_gt_gives_identical_record(self):
        samples = ['P', 'Q', 'C']
        families = [('P', 'Q', ('C',))]
        weights = {s: 1 for s in samples}
        plain = payload(samples, [(1, 'A', ['C', 'G'],
                                   {'P': (0, 1), 'Q': (1, 2), 'C': (2, 2)})],
                        families, weights)
        swapped = payload(samples, [(1, 'A', ['G', 'C'],
                                     {'P': (0, 2), 'Q': (2, 1), 'C': (1, 1)})],
                          families, weights)
        first = audit(plain)['records'][0]
        second = audit(swapped)['records'][0]
        self.assertEqual(first, second)
        self.assertEqual(first['alts'], ['C', 'G'])

    def test_duplicate_record_with_reordered_alts_rejected(self):
        samples = ['S']
        data = payload(samples, [(1, 'A', ['C', 'G'], {'S': (0, 1)}),
                                 (1, 'A', ['G', 'C'], {'S': (0, 2)})], [])
        with self.assertRaisesRegex(ValueError, 'duplicate variant record'):
            audit(data)

    def test_overlapping_distinct_records_analyzed_independently(self):
        samples = ['S']
        data = payload(samples, [(1, 'A', ['C'], {'S': (0, 1)}),
                                 (2, 'C', ['G'], {'S': (1, 1)})], [])
        result = audit(data)
        self.assertEqual(len(result['records']), 2)
        self.assertEqual([r['position'] for r in result['records']], [1, 2])


class BruteForceRandomized(unittest.TestCase):
    def test_random_small_cases_match_exhaustive_enumeration(self):
        rng = random.Random(20261118)
        trials = 0
        for _ in range(40):
            n_samples = rng.randint(1, 6)
            samples = [f'S{i}' for i in range(n_samples)]
            n_alts = rng.choice([1, 1, 2])
            n_alleles = n_alts + 1
            families = random_forest(samples, rng, 3)
            weights = {s: rng.randint(1, 4) for s in samples}
            obs = {}
            for s in samples:
                r = rng.random()
                if r < 0.35:
                    obs[s] = ()
                elif r < 0.6:
                    obs[s] = (rng.randrange(n_alleles),)
                else:
                    obs[s] = tuple(sorted(rng.choices(range(n_alleles), k=2)))
            alts = ['C', 'G'][:n_alts]
            check_against_brute(self, samples, [(1, 'A', alts, obs)],
                                families, weights)
            trials += 1
        self.assertEqual(trials, 40)

    def test_random_multiallelic_cases_match_brute_force(self):
        rng = random.Random(777)
        for _ in range(12):
            n_samples = rng.randint(3, 5)
            samples = [f'S{i}' for i in range(n_samples)]
            n_alleles = 4
            families = random_forest(samples, rng, 2)
            weights = {s: rng.randint(1, 3) for s in samples}
            obs = {}
            for s in samples:
                r = rng.random()
                if r < 0.4:
                    obs[s] = ()
                elif r < 0.7:
                    obs[s] = (rng.randrange(n_alleles),)
                else:
                    obs[s] = tuple(sorted(rng.choices(range(n_alleles), k=2)))
            check_against_brute(self, samples, [(1, 'A', ['C', 'G', 'T'], obs)],
                                families, weights)


class LargePedigrees(unittest.TestCase):
    def test_32_sample_chain_runs_and_counts_exactly(self):
        samples = [f's{i}' for i in range(32)]
        families = [(samples[2 * i - 2], samples[2 * i - 1], (samples[2 * i],))
                    for i in range(1, 16)]
        result = audit(payload(samples, [(1, 'A', ['C'], {})], families))
        record = result['records'][0]
        self.assertEqual(record['minimum_cost'], 0)
        # Independent chain count: genotype distribution after each family.
        gts = [(0, 0), (0, 1), (1, 1)]

        def mend(p, q):
            return {tuple(sorted((x, y))) for x in set(p) for y in set(q)}

        dist = {g: sum(1 for a in gts for b in gts if g in mend(a, b)) for g in gts}
        for _ in range(2, 16):
            dist = {g: sum(count * sum(1 for b in gts if g in mend(prev, b))
                           for prev, count in dist.items()) for g in gts}
        expected = sum(dist.values()) * 3  # s31 is an isolated missing sample
        self.assertEqual(record['optimal_assignments'], expected)
        self.assertTrue(all(g == [[0, 0], [0, 1], [1, 1]]
                            for g in record['genotypes']))
        self.assertEqual(record['required_masks'], [])
        self.assertEqual(record['possible_masks'], [])

    def test_wide_sibling_family_runs_and_counts_exactly(self):
        samples = ['P', 'Q'] + [f'c{i}' for i in range(30)]
        families = [('P', 'Q', tuple(samples[2:]))]
        result = audit(payload(samples, [(1, 'A', ['C'], {})], families))
        record = result['records'][0]
        gts = [(0, 0), (0, 1), (1, 1)]

        def mend(p, q):
            return {tuple(sorted((x, y))) for x in set(p) for y in set(q)}

        expected = sum(len(mend(p, q)) ** 30 for p in gts for q in gts)
        self.assertEqual(record['optimal_assignments'], expected)
        self.assertIsInstance(record['optimal_assignments'], int)
        self.assertEqual(record['minimum_cost'], 0)

    def test_32_sample_wide_family_with_calls(self):
        samples = ['P', 'Q'] + [f'c{i}' for i in range(30)]
        obs = {'P': (0, 1), 'Q': (0, 1)}
        obs.update({f'c{i}': (0, 0) for i in range(30)})
        families = [('P', 'Q', tuple(samples[2:]))]
        result = audit(payload(samples, [(1, 'A', ['C'], obs)], families))
        record = result['records'][0]
        # P=Q=0/1 can never yield a 0/0 child... it can: 0 from each parent.
        self.assertEqual(record['minimum_cost'], 0)
        self.assertEqual(record['optimal_assignments'], 1)


class InputContract(unittest.TestCase):
    def test_exact_top_level_fields(self):
        base = payload(['S'], [(1, 'A', ['C'], {'S': (0, 1)})], [])
        for mutation in ('extra', 'missing'):
            candidate = json.loads(json.dumps(base))
            if mutation == 'extra':
                candidate['extra'] = True
            else:
                del candidate['weights']
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                audit(candidate)

    def test_family_field_rules(self):
        samples = ['A', 'B', 'C', 'D']
        row = [(1, 'A', ['C'], {})]
        bad = [
            [{'parents': ['A'], 'children': ['C']}],                    # one parent
            [{'parents': ['A', 'B', 'D'], 'children': ['C']}],         # three parents
            [{'parents': ['A', 'A'], 'children': ['C']}],              # same parent twice
            [{'parents': ['A', 'B'], 'children': []}],                 # empty children
            [{'parents': ['A', 'B'], 'children': ['C', 'C']}],         # duplicate child
            [{'parents': ['A', 'B'], 'children': ['A']}],              # child is parent
            [{'parents': ['A', 'B'], 'children': ['ZZ']}],             # unknown child
            [{'parents': ['A', 'ZZ'], 'children': ['C']}],             # unknown parent
            [{'parents': ['A', 'B']}],                                 # missing children
            [{'parents': ['A', 'B'], 'children': ['C'], 'extra': 1}],  # extra field
        ]
        for families in bad:
            with self.subTest(families=families), self.assertRaises(ValueError):
                audit(payload(samples, row, families))

    def test_duplicate_parent_pair_rejected(self):
        samples = ['A', 'B', 'C', 'D']
        families = [('A', 'B', ('C',)), ('B', 'A', ('D',))]
        with self.assertRaisesRegex(ValueError, 'duplicate parent pair'):
            audit(payload(samples, [(1, 'A', ['C'], {})], families))

    def test_multiple_child_identity_rejected(self):
        samples = ['A', 'B', 'C', 'D', 'E']
        families = [('A', 'B', ('C',)), ('D', 'E', ('C',))]
        with self.assertRaisesRegex(ValueError, 'more than one family'):
            audit(payload(samples, [(1, 'A', ['C'], {})], families))

    def test_family_cycle_rejected(self):
        samples = ['A', 'B', 'C', 'D']
        families = [('A', 'B', ('C',)), ('B', 'C', ('D',))]
        with self.assertRaisesRegex(ValueError, 'forest'):
            audit(payload(samples, [(1, 'A', ['C'], {})], families))

    def test_indirect_family_cycle_rejected(self):
        samples = ['A', 'B', 'C', 'D', 'E', 'F']
        families = [('A', 'B', ('C',)), ('C', 'D', ('E',)), ('E', 'F', ('A',))]
        with self.assertRaisesRegex(ValueError, 'forest'):
            audit(payload(samples, [(1, 'A', ['C'], {})], families))

    def test_family_count_bound(self):
        # The 0..32 bound is a hard input limit and fires before item checks.
        families = [{'parents': ['A', 'B'], 'children': ['C']} for _ in range(33)]
        data = payload(['A', 'B', 'C'], [(1, 'A', ['C'], {})], [])
        data['families'] = families
        with self.assertRaisesRegex(ValueError, '0..32'):
            audit(data)
        # Zero families is allowed: every sample is an isolated variable.
        result = audit(payload(['A', 'B', 'C'], [(1, 'A', ['C'], {})], []))
        self.assertEqual(result['records'][0]['optimal_assignments'], 27)

    def test_weights_rules(self):
        samples = ['A', 'B']
        row = [(1, 'A', ['C'], {})]
        base = {s: 1 for s in samples}
        bad = [
            {'A': 1},                                  # missing sample
            {'A': 1, 'B': 1, 'ZZ': 1},                 # extra sample
            {'A': 0, 'B': 1},                          # zero
            {'A': -1, 'B': 1},                         # negative
            {'A': 1000001, 'B': 1},                    # above bound
            {'A': True, 'B': 1},                       # bool
            {'A': 1.5, 'B': 1},                        # non-integer
            {'A': '1', 'B': 1},                        # string
        ]
        for weights in bad:
            with self.subTest(weights=weights), self.assertRaises(ValueError):
                audit(payload(samples, row, [], weights=weights))
        ok = {'A': 1, 'B': 1000000}
        self.assertEqual(audit(payload(samples, row, [], weights=ok))['records'][0]
                         ['minimum_cost'], 0)

    def test_haploid_and_malformed_gt_rejected(self):
        for cell in ('0', '.', '1', '0/0/0', '0/', '/0', '0|', '0:.', '0/1/2'):
            with self.subTest(cell=cell), self.assertRaises(ValueError):
                data = payload(['S'], [(1, 'A', ['C'], {})], [])
                data['vcf'] = data['vcf'].replace('./.', cell)
                audit(data)

    def test_gt_allele_out_of_range_rejected(self):
        data = payload(['S'], [(1, 'A', ['C'], {})], [])
        data['vcf'] = data['vcf'].replace('./.', '2/2')
        with self.assertRaisesRegex(ValueError, 'out of range'):
            audit(data)

    def test_mode_bounds(self):
        # 33 samples
        samples = [f'S{i}' for i in range(33)]
        with self.assertRaisesRegex(ValueError, 'invalid samples'):
            audit(payload(samples, [(1, 'A', ['C'], {})], []))
        # 201 records
        rows = [(pos, 'A', ['C'], {}) for pos in range(1, 202)]
        with self.assertRaisesRegex(ValueError, 'row width/count'):
            audit(payload(['S'], rows, [], fasta='>chr1\n' + 'A' * 300 + '\n'))
        # 5 ALTs
        with self.assertRaisesRegex(ValueError, 'distinct sequence ALT'):
            audit(payload(['S'], [(1, 'A', ['C', 'G', 'T', 'CC', 'CA'], {})], []))

    def test_reference_mismatch_rejected(self):
        data = payload(['S'], [(1, 'A', ['C'], {})], [], fasta='>chr1\n' + 'T' * 100 + '\n')
        with self.assertRaisesRegex(ValueError, 'reference mismatch'):
            audit(data)

    def test_header_only_vcf_keeps_samples(self):
        text = HEADER + COLUMNS + '\tS2\tS1\n'
        data = {'fasta': FASTA, 'vcf': text, 'families': [],
                'weights': {'S1': 1, 'S2': 1}}
        result = audit(data)
        self.assertEqual(result, {'samples': ['S1', 'S2'], 'source_records': 0,
                                  'records': []})

    def test_report_shape_and_sorting(self):
        samples = ['B', 'A']
        rows = [(3, 'G', ['T'], {'A': (0, 1), 'B': ()}),
                (1, 'A', ['C'], {'A': (), 'B': (1, 1)})]
        result = audit(payload(samples, rows, []))
        self.assertEqual(set(result), {'samples', 'source_records', 'records'})
        self.assertEqual(result['samples'], ['A', 'B'])
        self.assertEqual(result['source_records'], 2)
        self.assertEqual([r['position'] for r in result['records']], [1, 3])
        for record in result['records']:
            self.assertEqual(set(record), {'contig', 'position', 'ref', 'alts',
                                           'minimum_cost', 'optimal_assignments',
                                           'required_masks', 'possible_masks',
                                           'genotypes'})
            self.assertIsInstance(record['position'], int)
            self.assertIsInstance(record['minimum_cost'], int)
            self.assertIsInstance(record['optimal_assignments'], int)
            self.assertGreater(record['optimal_assignments'], 0)
            self.assertEqual(len(record['genotypes']), 2)
            for support in record['genotypes']:
                self.assertEqual(support, sorted(support))
                for gt in support:
                    self.assertEqual(len(gt), 2)
                    self.assertLessEqual(gt[0], gt[1])


class PermutationInvariance(unittest.TestCase):
    def test_permuted_inputs_give_identical_result_except_hash(self):
        samples = ['P', 'Q', 'C', 'S']
        families = [('P', 'Q', ('C',))]
        weights = {'P': 1, 'Q': 2, 'C': 3, 'S': 4}
        rows = [(1, 'A', ['C'], {'P': (0, 0), 'Q': (0, 0), 'C': (0, 1), 'S': ()}),
                (2, 'C', ['G', 'T'], {'P': (1, 2), 'Q': (), 'C': (0, 2), 'S': (1, 1)})]
        plain = payload(samples, rows, families, weights)
        # Permute: sample columns, row order, ALT order (with GT remap),
        # family parent/children order.
        perm_samples = ['S', 'C', 'Q', 'P']
        perm_rows = [
            (2, 'C', ['T', 'G'], {'P': (2, 1), 'Q': (), 'C': (0, 1), 'S': (2, 2)}),
            (1, 'A', ['C'], {'P': (0, 0), 'Q': (0, 0), 'C': (0, 1), 'S': ()}),
        ]
        permuted = payload(perm_samples, perm_rows, [], weights)
        permuted['families'] = [{'parents': ['Q', 'P'], 'children': ['C']}]
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            reports = []
            for name, data in [('plain', plain), ('permuted', permuted)]:
                source = directory / f'{name}.json'
                source.write_text(json.dumps(data))
                reports.append(run(source, audit_pedigree=True))
        self.assertNotEqual(reports[0]['source_sha256'], reports[1]['source_sha256'])
        for report in reports:
            report.pop('source_sha256')
        self.assertEqual(reports[0], reports[1])

    def test_family_and_children_order_do_not_matter(self):
        samples = ['P', 'Q', 'C1', 'C2', 'R', 'S', 'T']
        obs = {'P': (0, 1), 'Q': (), 'C1': (0, 0), 'C2': (1, 1),
               'R': (0, 0), 'S': (0, 0), 'T': (0, 1)}
        first = payload(samples, [(1, 'A', ['C'], obs)],
                        [('P', 'Q', ('C1', 'C2')), ('R', 'S', ('T',))])
        second = payload(samples, [(1, 'A', ['C'], obs)],
                         [('S', 'R', ('T',)), ('Q', 'P', ('C2', 'C1'))])
        self.assertEqual(audit(first), audit(second))


class CliEndToEnd(unittest.TestCase):
    def _env(self):
        return dict(os.environ, PYTHONPATH=str(SRC))

    def _cli(self, source, output, *flags):
        command = [sys.executable, '-m', 'variant_cohort.cli',
                   '--source', str(source), '--output', str(output)]
        command.extend(flags)
        return subprocess.run(command, capture_output=True, env=self._env(), cwd=ROOT)

    def test_end_to_end_pedigree_audit_report(self):
        data = payload(['P', 'Q', 'C'],
                       [(1, 'A', ['C'], {'P': (0, 0), 'Q': (0, 0), 'C': (0, 1)})],
                       [('P', 'Q', ('C',))])
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source = directory / 'input.json'
            output = directory / 'report.json'
            source.write_text(json.dumps(data))
            done = self._cli(source, output, '--audit-pedigree')
            self.assertEqual(done.returncode, 0, done.stderr)
            report = json.loads(output.read_text())
            self.assertEqual(set(report), {'schema_version', 'mode',
                                           'source_sha256', 'result'})
            self.assertEqual(report['schema_version'], 4)
            self.assertEqual(report['mode'], 'pedigree_audit')
            self.assertEqual(set(report['result']),
                             {'samples', 'source_records', 'records'})
            record = report['result']['records'][0]
            self.assertEqual(record['minimum_cost'], 1)
            self.assertEqual(record['optimal_assignments'], 5)
            self.assertEqual(record['possible_masks'], ['C', 'P', 'Q'])

    def test_installed_pedigree_fixture_matches_library(self):
        path = ROOT / 'fixtures' / 'pedigree.json'
        report = run(path, audit_pedigree=True)
        self.assertEqual(report['schema_version'], 4)
        self.assertEqual(report['mode'], 'pedigree_audit')
        self.assertEqual(report['result'], audit(json.loads(path.read_text())))

    def test_audit_and_other_modes_are_mutually_exclusive(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            output = directory / 'report.json'
            for flags in [('--audit-pedigree', '--stitch-phase'),
                          ('--audit-pedigree', '--project-assembly')]:
                done = self._cli(ROOT / 'fixtures' / 'pedigree.json', output, *flags)
                self.assertEqual(done.returncode, 2)
                self.assertFalse(output.exists())

    def test_audit_rejects_other_mode_sources_and_preserves_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            output = directory / 'report.json'
            output.write_text('keep')
            for fixture in ('input.json', 'phase.json', 'projection.json'):
                done = self._cli(ROOT / 'fixtures' / fixture, output, '--audit-pedigree')
                self.assertEqual(done.returncode, 2)
                self.assertEqual(output.read_text(), 'keep')

    def test_invalid_input_preserves_output_bytes(self):
        data = payload(['P', 'Q', 'C'],
                       [(1, 'A', ['C'], {'P': (0, 0), 'Q': (0, 0), 'C': (0, 1)})],
                       [('P', 'Q', ('C',)), ('Q', 'C', ('P',))])  # cyclic
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source = directory / 'input.json'
            output = directory / 'report.json'
            source.write_text(json.dumps(data))
            original = b'{"previous": true}\n'
            output.write_bytes(original)
            done = self._cli(source, output, '--audit-pedigree')
            self.assertEqual(done.returncode, 2)
            self.assertEqual(output.read_bytes(), original)

    def test_input_equals_output_and_symlink_rejected(self):
        data = payload(['S'], [(1, 'A', ['C'], {'S': (0, 1)})], [])
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source = directory / 'input.json'
            source.write_text(json.dumps(data))
            original = source.read_bytes()
            link = directory / 'link.json'
            link.symlink_to(source)
            for destination in (source, link):
                done = self._cli(source, destination, '--audit-pedigree')
                self.assertEqual(done.returncode, 2)
                self.assertEqual(source.read_bytes(), original)

    def test_success_replaces_atomically(self):
        data = payload(['S'], [(1, 'A', ['C'], {'S': (0, 1)})], [])
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source = directory / 'input.json'
            output = directory / 'report.json'
            source.write_text(json.dumps(data))
            output.write_text('stale')
            done = self._cli(source, output, '--audit-pedigree')
            self.assertEqual(done.returncode, 0, done.stderr)
            self.assertEqual(json.loads(output.read_text())['schema_version'], 4)

    def test_legacy_modes_still_work(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            outputs = []
            for name, fixture, flags in [
                    ('v1', 'input.json', ()),
                    ('v2', 'phase.json', ('--stitch-phase',)),
                    ('v3', 'projection.json', ('--project-assembly',)),
                    ('v4', 'pedigree.json', ('--audit-pedigree',))]:
                output = directory / f'{name}.json'
                done = self._cli(ROOT / 'fixtures' / fixture, output, *flags)
                self.assertEqual(done.returncode, 0, (name, done.stderr))
                outputs.append(json.loads(output.read_text()))
            self.assertEqual([r['schema_version'] for r in outputs], [1, 2, 3, 4])


if __name__ == '__main__':
    unittest.main()
