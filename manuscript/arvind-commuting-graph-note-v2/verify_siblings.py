from collections import Counter, defaultdict
from itertools import product
import json

from verify import Ring


class PrimeRing(Ring):
    def __init__(self, prime, dual):
        super().__init__(prime, dual)
        self.q = prime

    def residue(self, value):
        return value % self.prime

    def nil(self, value):
        return self.prime * value

    def residue_mul(self, left, right):
        return left * right % self.prime

    def code(self, value):
        return value % self.prime, value // self.prime

    def from_code(self, residue, upper):
        return residue + self.prime * upper


class QuadraticRing:
    def __init__(self, dual):
        self.prime = 2
        self.q = 4
        self.size = 16
        self.dual = dual

    @staticmethod
    def coefficient_mul(left, right, modulus):
        left_low, left_high = left % modulus, left // modulus
        right_low, right_high = right % modulus, right // modulus
        low = (left_low * right_low - left_high * right_high) % modulus
        high = (left_low * right_high + left_high * right_low
                - left_high * right_high) % modulus
        return low + modulus * high

    def add(self, left, right):
        if self.dual:
            return left ^ right
        return (left % 4 + right % 4) % 4 + 4 * ((left // 4 + right // 4) % 4)

    def neg(self, value):
        if self.dual:
            return value
        return (-value % 4) + 4 * (-(value // 4) % 4)

    def mul(self, left, right):
        if not self.dual:
            return self.coefficient_mul(left, right, 4)
        low = self.coefficient_mul(left % 4, right % 4, 2)
        upper = (self.coefficient_mul(left % 4, right // 4, 2)
                 ^ self.coefficient_mul(left // 4, right % 4, 2))
        return low + 4 * upper

    def residue(self, value):
        if self.dual:
            return value % 4
        return value % 2 + 2 * ((value // 4) % 2)

    def nil(self, value):
        if self.dual:
            return 4 * value
        return 2 * (value % 2) + 8 * (value // 2)

    def residue_mul(self, left, right):
        return self.coefficient_mul(left, right, 2)

    def inverse(self, value):
        return next(candidate for candidate in range(self.size)
                    if self.mul(value, candidate) == 1)

    def code(self, value):
        if self.dual:
            return value % 4, value // 4
        return self.residue(value), (value % 4) // 2 + 2 * (value // 8)

    def from_code(self, residue, upper):
        if self.dual:
            return residue + 4 * upper
        return residue % 2 + 2 * (upper % 2) + 4 * (residue // 2) + 8 * (upper // 2)


def chart(ring, vector):
    first, second = vector
    if ring.residue(first):
        return 'I', first, ring.mul(ring.inverse(first), second)
    if ring.residue(second):
        coefficient = ring.mul(ring.inverse(second), first)
        return 'II', ring.code(coefficient)[1], second
    return 'III', ring.code(first)[1], ring.code(second)[1]


def decode(ring, coordinates):
    kind, first, second = coordinates
    if kind == 'I':
        return first, ring.mul(first, second)
    if kind == 'II':
        return ring.mul(second, ring.nil(first)), second
    return ring.nil(first), ring.nil(second)


def commutes(ring, left, right):
    return ring.mul(left[0], right[1]) == ring.mul(right[0], left[1])


def check_ring_laws(ring):
    elements = range(ring.size)
    for first in elements:
        assert ring.add(first, ring.neg(first)) == 0
        assert ring.mul(first, 1) == first
        assert ring.from_code(*ring.code(first)) == first
        for second in elements:
            assert ring.add(first, second) == ring.add(second, first)
            assert ring.mul(first, second) == ring.mul(second, first)
            for third in elements:
                assert ring.add(ring.add(first, second), third) == ring.add(first, ring.add(second, third))
                assert ring.mul(ring.mul(first, second), third) == ring.mul(first, ring.mul(second, third))
                assert ring.mul(first, ring.add(second, third)) == ring.add(ring.mul(first, second), ring.mul(first, third))
    for scalar in range(ring.q):
        for element in elements:
            assert ring.mul(ring.nil(scalar), element) == ring.nil(ring.residue_mul(scalar, ring.residue(element)))


def check_pair(square, dual):
    residue_size = square.q
    assert dual.q == residue_size
    vectors = list(product(range(square.size), repeat=2))
    mapping = {}
    blocks = defaultdict(set)
    nil_vectors = set()
    for vector in vectors:
        coordinates = chart(square, vector)
        assert decode(square, coordinates) == vector
        kind, first, second = coordinates
        if kind == 'I':
            target = kind, dual.from_code(*square.code(first)), dual.from_code(*square.code(second))
            blocks[kind, second].add(vector)
        elif kind == 'II':
            target = kind, first, dual.from_code(*square.code(second))
            blocks[kind, first].add(vector)
        else:
            target = coordinates
            nil_vectors.add(vector)
        mapping[vector] = decode(dual, target)
        assert chart(dual, mapping[vector]) == target
    assert len(set(mapping.values())) == len(vectors)
    assert len(blocks) == residue_size * (residue_size + 1)
    assert len(nil_vectors) == residue_size**2
    assert all(len(block) == residue_size**2 - residue_size for block in blocks.values())
    base_centralizers = Counter()
    zero_pairs = 0
    for left in vectors:
        centralizer_size = 0
        for right in vectors:
            source_commutes = commutes(square, left, right)
            assert source_commutes == commutes(dual, mapping[left], mapping[right])
            centralizer_size += source_commutes
        base_centralizers[centralizer_size] += 1
        zero_pairs += centralizer_size
    assert base_centralizers == Counter({residue_size**4: 1,
                                       residue_size**3: residue_size**2 - 1,
                                       residue_size**2: residue_size**4 - residue_size**2})
    for block in blocks.values():
        representative = next(iter(block))
        neighbours = {vector for vector in vectors if commutes(square, representative, vector)}
        assert neighbours - nil_vectors == block
        assert len(neighbours & nil_vectors) == residue_size
        assert all(commutes(square, left, right) for left in neighbours for right in neighbours)
    first_block, second_block = list(blocks.values())[:2]
    disjoint_edges = sorted(first_block)[:2] + sorted(second_block)[:2]
    assert len(set(disjoint_edges)) == 4
    for left_index, left in enumerate(disjoint_edges):
        for right_index, right in enumerate(disjoint_edges):
            if left_index != right_index:
                assert commutes(square, left, right) == (left_index // 2 == right_index // 2)
    for ring in (square, dual):
        nonzero_nil = ring.nil(1)
        path = [(1, 0), (nonzero_nil, 0), (0, nonzero_nil), (0, 1)]
        assert len(set(path)) == 4
        for left_index, left in enumerate(path):
            for right_index, right in enumerate(path):
                if left_index != right_index:
                    assert commutes(ring, left, right) == (abs(left_index - right_index) == 1)
        for element in range(ring.size):
            total = 0
            for unused in range(ring.prime):
                total = ring.add(total, element)
            if ring.dual:
                assert total == 0
            total = 0
            for unused in range(ring.prime**2):
                total = ring.add(total, element)
            assert total == 0
        total = 0
        for unused in range(ring.prime):
            total = ring.add(total, 1)
        assert (total == 0) == ring.dual
    expected = residue_size**6 + residue_size**5 - residue_size**3
    assert zero_pairs == expected
    return {'residue_field_size': residue_size,
            'group_order': residue_size**6,
            'ordered_base_pairs': len(vectors)**2,
            'zero_determinant_ordered_pairs': zero_pairs,
            'group_degree_multiplicities': {
                str(residue_size**6 - 1): residue_size**2,
                str(residue_size**5 - 1): residue_size**4 - residue_size**2,
                str(residue_size**4 - 1): residue_size**6 - residue_size**4},
            'clique_and_chromatic_number': residue_size**4,
            'conjugacy_class_number': residue_size**4 + residue_size**3 - residue_size,
            'primitive_blocks': len(blocks),
            'simplicial_primitive_blocks': True,
            'induced_P4': True,
            'induced_2K2': True,
            'derived_exponents': {'square': square.prime**2, 'dual': square.prime},
            'status': 'pass'}


def main():
    pairs = [(PrimeRing(prime, False), PrimeRing(prime, True)) for prime in [2, 3, 5]]
    pairs.insert(2, (QuadraticRing(False), QuadraticRing(True)))
    for pair in pairs:
        for ring in pair:
            check_ring_laws(ring)
    print(json.dumps({'arithmetic': 'exact finite-ring arithmetic',
                      'coverage': 'all base pairs, clique decomposition, degrees, conjugacy counts, induced P4 and 2K2, derived exponents for q=2,3,4,5',
                      'checks': [check_pair(*pair) for pair in pairs],
                      'status': 'pass'}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
