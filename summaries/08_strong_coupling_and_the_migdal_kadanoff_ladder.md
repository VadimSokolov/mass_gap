# Chapter 8. Strong Coupling and the Migdal Kadanoff Ladder

<!-- cite:begin -->
From *The Mass Gap: Thorin Measures, Brownian Motion and Yang Mills* by Nick Polson and Vadim Sokolov (2026), chapter 8. Cite as Polson and Sokolov (2026), *The Mass Gap: Thorin Measures, Brownian Motion and Yang Mills*, edition 1, https://doi.org/10.5281/zenodo.24546905. Book page: https://vsokolov.org/books/mass_gap/. Source: https://github.com/VadimSokolov/mass_gap.
<!-- cite:end -->

This chapter moves to four dimensional SU(2) lattice Yang Mills with the Wilson action and asks how far the gap can be controlled from strong coupling toward the continuum. It states the strong coupling character expansion and cluster expansion (the Osterwalder and Seiler theorem: a positive lattice mass and an area law), and explains why strong coupling control does not reach the continuum. It describes how glueball masses are measured by lattice Monte Carlo. It then reads the Migdal Kadanoff recursion, bond moving followed by decimation, as a Bayesian product followed by a convolution of laws, proves that each step preserves nonnegative character coefficients, so that the four dimensional lattice theory with the plaquette weight of any rung is reflection positive, and derives the coupling decrement through the mean of the Biane Pitman Yor variable, at leading order with the remainder not bounded (asymptotic freedom on the ladder). Finally it defines a heat kernel projected one dimensional map, proves its monotonicity from a monotone likelihood ratio property, and proves a first passage criterion for the correlation length. The gap criterion for the hierarchical lattice is conditional: (C) and (D) are open and (E) fails as stated, and even granted all three the remark bounds a correlation length, not the Thorin gap inequality.

## From the literature

- Character expansion, strong coupling cluster expansion, area law and positive lattice mass for small beta (Osterwalder and Seiler 1978), with the strong coupling mass series (Munster 1981).
- The Migdal recursion (Migdal 1975) and its explanation as bond moving followed by decimation, with bond moving giving one sided bounds on the free energy (Kadanoff 1975, 1976).
- Probability laws tied to theta and zeta functions, including the variable whose Laplace transform is sqrt(2q)/sinh sqrt(2q) (Biane, Pitman and Yor 2001).
- The duality between envelopes and mixtures (Polson and Scott 2016), and the Clebsch Gordan rule.
- Confinement on the Migdal Kadanoff hierarchical lattice (Ito 1985).
- Lattice Monte Carlo for glueballs (Creutz 1980; Kennedy and Pendleton 1985; Morningstar and Peardon 1999); asymptotic freedom (Gross and Wilczek 1973; Politzer 1973).

## New results

The new statements are the positivity lemma for the ladder, the monotonicity lemma for the projected map and the first passage criterion; none of them is a proof of a gap.

<!-- results:begin -->
- [Positivity of character coefficients is preserved by a Migdal Kadanoff step](../RESULTS.md#mk_step_preserves_positivity) (Lemma 8.3, page 128; proved): One step of the Migdal Kadanoff recursion keeps all character coefficients of the plaquette weight nonnegative, so the four dimensional lattice theory with the plaquette weight of any rung is reflection positive. New in this book.
- [Monotone likelihood ratio of the heat kernel makes the projected ladder map increasing](../RESULTS.md#mlr_implies_projected_map_increasing) (Lemma 8.7, page 134; proved): If the SU(2) heat kernel has monotone likelihood ratio in time, which the book proves elsewhere, then the one dimensional heat kernel projected Migdal Kadanoff map is increasing. New in this book.
- [First passage criterion for the projected ladder](../RESULTS.md#first_passage_criterion_projected_ladder) (Proposition 8.9, page 136; proved): If the projected ladder map lowers the coupling at every coupling and the decrement tends to a positive constant, every coupling reaches strong coupling in a number of steps linear in the coupling; this concerns the projected map and is not a gap theorem. New in this book.
<!-- results:end -->

## Key terms

strong coupling expansion, character expansion, cluster expansion, Wilson action, glueball, Migdal Kadanoff recursion, bond moving, decimation, hierarchical lattice, renormalization group ladder, convolution, reflection positivity
