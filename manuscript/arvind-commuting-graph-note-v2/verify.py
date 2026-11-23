import argparse
from itertools import product
import json


class Ring:
    def __init__(self, prime, dual):
        self.prime = prime
        self.size = prime * prime
        self.dual = dual

    def add(self, left, right):
        if not self.dual:
            return (left + right) % self.size
        prime = self.prime
        return ((left % prime + right % prime) % prime
                + prime * ((left // prime + right // prime) % prime))

    def neg(self, value):
        if not self.dual:
            return -value % self.size
        prime = self.prime
        return (-value % prime) + prime * (-(value // prime) % prime)

    def mul(self, left, right):
        if not self.dual:
            return left * right % self.size
        prime = self.prime
        return ((left % prime) * (right % prime) % prime
                + prime * (((left % prime) * (right // prime)
                            + (left // prime) * (right % prime)) % prime))

    def inverse(self, value):
        assert value % self.prime != 0
        if not self.dual:
            return pow(value, -1, self.size)
        prime = self.prime
        inverse_residue = pow(value % prime, -1, prime)
        return (inverse_residue
                + prime * (-(value // prime) * inverse_residue**2 % prime))

    def determinant(self, left, right):
        return self.add(self.mul(left[0], right[1]),
                        self.neg(self.mul(right[0], left[1])))

    def triple_mul(self, left, right):
        return (self.add(left[0], right[0]), self.add(left[1], right[1]),
                self.add(self.add(left[2], right[2]), self.mul(left[0], right[1])))

    def chart(self, vector):
        first, second = vector
        prime = self.prime
        if first % prime:
            return ('I', first, self.mul(self.inverse(first), second))
        if second % prime:
            coefficient = self.mul(self.inverse(second), first)
            assert coefficient % prime == 0
            return ('II', coefficient // prime, second)
        return ('III', first // prime, second // prime)

    def decode(self, chart):
        kind, first, second = chart
        if kind == 'I':
            return first, self.mul(first, second)
        if kind == 'II':
            return self.mul(second, self.prime * first), second
        return self.prime * first, self.prime * second


def chart_incidence(left, right, prime):
    if left[0] == 'III' and right[0] != 'III':
        left, right = right, left
    if left[0] == 'II' and right[0] == 'I':
        left, right = right, left
    kinds = left[0], right[0]
    if kinds == ('I', 'I'):
        return left[2] == right[2]
    if kinds == ('I', 'II'):
        return False
    if kinds == ('I', 'III'):
        return right[2] == (left[2] % prime) * right[1] % prime
    if kinds == ('II', 'II'):
        return left[1] == right[1]
    if kinds == ('II', 'III'):
        return right[1] == 0
    assert kinds == ('III', 'III')
    return True


def check_prime(prime):
    assert prime >= 2 and all(prime % divisor for divisor in range(2, prime))
    square = Ring(prime, False)
    dual = Ring(prime, True)
    elements = range(prime**2)
    vectors = list(product(elements, repeat=2))
    charts = {vector: square.chart(vector) for vector in vectors}
    mapping = {vector: dual.decode(charts[vector]) for vector in vectors}
    assert len(set(mapping.values())) == len(vectors)
    assert mapping[(0, 0)] == (0, 0)
    for vector in vectors:
        assert square.decode(charts[vector]) == vector
        assert dual.chart(mapping[vector]) == charts[vector]
    commuting_pairs = 0
    for left in vectors:
        for right in vectors:
            square_zero = square.determinant(left, right) == 0
            dual_zero = dual.determinant(mapping[left], mapping[right]) == 0
            assert square_zero == dual_zero
            assert square_zero == chart_incidence(charts[left], charts[right], prime)
            commuting_pairs += square_zero
    for ring in (square, dual):
        for value in elements:
            residue = value % prime
            assert ring.mul(prime, value) == prime * residue
            for coefficient in range(prime):
                assert ring.mul(prime * coefficient, value) == prime * (coefficient * residue % prime)
            if residue:
                assert ring.mul(value, ring.inverse(value)) == 1
    assert all(dual.add(0, value) == value for value in elements)
    assert all(dual.add(value, dual.neg(value)) == 0 for value in elements)
    assert square.size == prime**2
    dual_axis_killed = []
    for value in elements:
        total = 0
        for unused in range(prime):
            total = dual.add(total, value)
        dual_axis_killed.append(total == 0)
    assert all(dual_axis_killed)
    assert prime % square.size != 0
    full_pairs = 0
    if prime <= 3:
        triples = list(product(elements, repeat=3))
        for left in triples:
            mapped_left = (*mapping[left[:2]], left[2])
            for right in triples:
                mapped_right = (*mapping[right[:2]], right[2])
                square_adjacent = (left != right and square.triple_mul(left, right)
                                   == square.triple_mul(right, left))
                dual_adjacent = (mapped_left != mapped_right
                                 and dual.triple_mul(mapped_left, mapped_right)
                                 == dual.triple_mul(mapped_right, mapped_left))
                assert square_adjacent == dual_adjacent
                full_pairs += 1
    return {'prime': prime, 'group_order': prime**6,
            'base_vectors': len(vectors), 'ordered_base_pairs': len(vectors)**2,
            'zero_determinant_ordered_pairs': commuting_pairs,
            'ordered_full_vertex_pairs': full_pairs,
            'derived_exponents': {'square': prime**2, 'dual': prime},
            'status': 'pass'}


def main():
    parser = argparse.ArgumentParser(description='Exact checks for the two Heisenberg commuting graphs.')
    parser.add_argument('--primes', nargs='+', type=int, default=[2, 3, 5])
    arguments = parser.parse_args()
    result = {'arithmetic': 'exact integer and finite-ring arithmetic',
              'coverage': 'chart bijections, determinant-zero equivalence, derived-axis exponents; full triple adjacency for primes at most 3',
              'checks': [check_prime(prime) for prime in arguments.primes],
              'status': 'pass'}
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
