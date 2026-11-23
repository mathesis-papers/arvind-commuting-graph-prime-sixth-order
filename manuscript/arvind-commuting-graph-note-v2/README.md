# Residue-field extensions of commuting-graph examples

The note extends the examples attributed to Eamonn O'Brien by constructing,
for every `q = p^f` with `p` prime and `f >= 1`, Heisenberg groups over
`R = (Z/p²Z)[X]/(h)` and `F_q[ε]/(ε²)` of order `q⁶`, with isomorphic full
commuting graphs but derived groups `(C_(p²))^f` and `(C_p)^(2f)`.
Here `h` is a monic lift of a degree-`f` irreducible polynomial over `F_p`.
The graphs are chordal and perfect, neither cographs, split nor threshold
graphs, and have explicit degree
distributions. Taking `f = 1` answers the second
sentence of Question 4 in Section 9, page 40 of Arvind, Cameron, Ma and
Maslova, *Aspects of the commuting graph*, arXiv:2305.07301v5.

## Files

- `note.tex`: complete proof, scientific references and verification scope.
- `note.pdf`: readable typeset note, supplied separately from the source archive.
- `metadata.json`: authors, license, covered problem and pinned scientific source.
- `coverage.json`: source locators and precise written, formal and finite-check scope.
- `verify.py`: dependency-free exact finite arithmetic checks.
- `checks.json`: expected output for the default prime range.
- `verify_siblings.py`: residue-field and graph-structure checks, including `F_4`.
- `sibling-checks.json`: expected output for `q = 2, 3, 4, 5`.

## Build

Use Tectonic **0.17.0** from the directory containing `note.tex`:

```sh
tectonic note.tex
```

Alternatively use a standard LaTeX installation:

```sh
pdflatex -interaction=nonstopmode -halt-on-error note.tex
pdflatex -interaction=nonstopmode -halt-on-error note.tex
```

The bibliography is embedded; no external bibliography or figure files are needed.

## Exact finite checks

Use Python 3.8 or later, without optimization flags:

```sh
python3 verify.py > actual-checks.json
diff -u checks.json actual-checks.json
python3 verify_siblings.py > actual-sibling-checks.json
diff -u sibling-checks.json actual-sibling-checks.json
```

The default checks use `p = 2, 3, 5`. They enumerate all ordered pairs of
base vectors, verify the chart description and graph bijection, and verify
the derived-axis exponent distinction. For `p = 2, 3`, they additionally
enumerate all ordered pairs of full group elements using the triple
multiplication law. Every operation is exact; no numerical tolerances are used.
The additional checker enumerates every ordered base-vector pair for
`q = 2, 3, 4, 5`, tests both rings' algebraic laws, the clique decomposition,
the degree and conjugacy-class formulas, an induced four-vertex path and
two induced disjoint edges.
For `q = 4` it uses `(Z/4Z)[X]/(X²+X+1)` and `F_4[ε]/(ε²)`.
The all-prime-power statement, product consequence and graph structure
have written proofs. The pinned Lean declarations linked in the note cover
the generic residue-coordinate determinant comparison and the concrete
prime-field (`f = 1`) existence and non-isoclinism results.

## Authors and license

Haobo Ma — ChronoAI Pte Ltd, Singapore.

Wenlin Zhang — National University of Singapore; ChronoAI Pte Ltd, Singapore.

CC BY 4.0 (`cc-by-4.0`).
