# Chapter 5. Markov Chains, Mixing and Slice Variables

<!-- cite:begin -->
From *The Mass Gap: Thorin Measures, Brownian Motion and Yang Mills* by Nick Polson and Vadim Sokolov (2026), chapter 5. Cite as Polson and Sokolov (2026), *The Mass Gap: Thorin Measures, Brownian Motion and Yang Mills*, edition 1, https://doi.org/10.5281/zenodo.24546905. Book page: https://vsokolov.org/books/mass_gap/. Source: https://github.com/VadimSokolov/mass_gap.
<!-- cite:end -->

What limits the simulation of lattice gauge fields, and how do mixing, autocorrelation and slice variables relate to the mass gap? The chapter opens with Mark Kac: models as caricatures, Polya's recurrence theorem and marginal dimensions, Kac's lemma on mean return times, and the Feynman Kac formula linking Euclidean measures to the Hamiltonian gap. It then covers the Jerrum and Sinclair FPRAS for the ferromagnetic Ising partition function, canonical paths, the product estimator, and the relation between rapid mixing, log Sobolev constants and decay of correlations for spin systems. The Sokal section treats integrated autocorrelation time with Madras and Sokal windowing, Swendsen Wang and Edwards Sokal clusters, multigrid Monte Carlo, topological freezing in two dimensional U(1), and Lee Yang zeros. The last section covers slice variables, Higdon's partial decoupling, Besag's Markov random fields and a pseudolikelihood experiment for SU(2). It conjectures, with a sketch of proof, a log Sobolev bound for the sliced SU(2) link conditional.

## From the literature

- Feynman Kac formula (Kac 1949), models as caricatures (Kac 1969), Kac's lemma and the Ehrenfest urn, Polya's recurrence theorem (1921).
- FPRAS for the ferromagnetic Ising model (Jerrum and Sinclair 1993); spectral gap, log Sobolev constant and strong spatial mixing (Martinelli and Olivieri 1994).
- Log concave sampling in polynomial time (Applegate and Kannan 1991; Frieze, Kannan and Polson 1994, the first author's earlier work).
- Automatic windowing (Madras and Sokal 1988); Swendsen Wang algorithm (1987); Edwards Sokal representation (1988); multigrid Monte Carlo (Goodman and Sokal 1989).
- Lee Yang circle theorem (1952); zero free regions for Potts models (Sokal 2001).
- Slice variables and scale mixtures (Polson 1996, the first author's earlier work); auxiliary variables and partial decoupling (Higdon 1998).
- Hammersley Clifford theorem and the coding method (Besag 1974), pseudolikelihood (Besag 1975) and its consistency (Comets 1992).

## New results

<!-- results:begin -->
- [An FPRAS for the partition function of SU(2) lattice gauge theory](../RESULTS.md#op_fpras_su2_lattice_gauge) (Open problem 1, page 63; open): The book asks whether the partition function of SU(2) lattice gauge theory with the Wilson action can be approximated in time polynomial in the box size at every coupling. New in this book.
- [Slice variables give log concave conditionals](../RESULTS.md#slice_conditional_log_concave) (Proposition 5.2, page 69; proved): Slicing each quasi concave likelihood factor makes the conditional law of the parameter log concave, even when the posterior itself is multimodal. First stated in Polson (1996), Convergence of Markov chain Monte Carlo algorithms, in Bayesian Statistics 5, Oxford University Press.
- [Log Sobolev bound for the sliced SU(2) link conditional](../RESULTS.md#slice_cap_log_sobolev) (Conjecture 5.3, page 70; open): The book conjectures that, once the plaquette slices are cut at positive levels, the uniform law on the resulting intersection of spherical caps in SU(2) has a log Sobolev constant that does not depend on the coupling or the staples (a proof is only sketched). New in this book.
<!-- results:end -->

## Key terms

Markov chain Monte Carlo, mixing time, rapid mixing, spectral gap, log Sobolev inequality, critical slowing down, integrated autocorrelation time, Madras Sokal windowing, Swendsen Wang, Edwards Sokal, random cluster model, multigrid Monte Carlo, topological freezing, Kac lemma, Ehrenfest urn, Polya recurrence, Feynman Kac formula, FPRAS, canonical paths, product estimator, Lee Yang zeros, slice sampling, quasi concave, log concave, partial decoupling, Markov random field, Hammersley Clifford, pseudolikelihood, ICAR, SU(2) heat bath, geodesic convexity, Bakry Emery
