# Chapter 14. Inequalities Between the Ladder and the Lattice

<!-- cite:begin -->
From *The Mass Gap: Thorin Measures, Brownian Motion and Yang Mills* by Nick Polson and Vadim Sokolov (2026), chapter 14. Cite as Polson and Sokolov (2026), *The Mass Gap: Thorin Measures, Brownian Motion and Yang Mills*, edition 1, https://doi.org/10.5281/zenodo.24546905. Book page: https://vsokolov.org/books/mass_gap/. Source: https://github.com/VadimSokolov/mass_gap.
<!-- cite:end -->

This chapter compares the Migdal Kadanoff ladder with the true four dimensional SU(2) lattice gauge theory, the comparison that step 4 of the book's path to the mass gap needs. It tests the sandwich idea on one block spin step of scale factor two: the law of a blocked plaquette, the 2x2 Wilson loop holonomy, is measured by heat bath Monte Carlo on an 8^4 lattice at six couplings and compared with an independent plaquette caricature and with the Migdal Kadanoff bond moving caricature, in mean and in the likelihood ratio order. The measured law lies between them in mean at every coupling and in likelihood ratio order on most of the mass. The chapter then treats correlation inequalities: Griffiths and Ginibre inequalities for abelian models, a doubling inequality derived from reflection positivity and the Schwarz inequality that proves W22 >= W12^2 >= W11^4 for every compact gauge group, a rigorous heat bath (multihit, Rao Blackwell) upper bound for SU(2) loops, and Kadanoff's variational upper bound on the pressure that holds along the whole ladder. It reports a first empirical rung with pseudolikelihood fits of the blocked coupling, states a weak coupling sandwich conjecture, and lists what a proof of the Migdal Kadanoff end would need.

## From the literature

Griffiths' inequalities for Ising ferromagnets (Griffiths 1967); the general formulation of Griffiths inequalities (Ginibre 1970); the FKG inequality (Fortuin, Kasteleyn and Ginibre 1971); reflection positivity of the Wilson lattice measure and the loop inequalities it gives, among them W22 >= W12^2 >= W11^4 for every compact group (Osterwalder and Seiler 1978); convexity of the static quark potential from reflection positivity (Bachas 1986); the multihit estimator for Wilson loops and the perimeter bound it gives (Parisi, Petronzio and Rapuano 1983) and Rao Blackwellisation (Blackwell 1947); Migdal's and Kadanoff's recursions (Migdal 1975, Kadanoff 1975 and 1976); Kadanoff's variational bound on the pressure by bond moving; Monte Carlo renormalization group (Swendsen 1979).

## New results

Only the weak coupling sandwich conjecture is the authors' own; the other labelled statements of the chapter are known results presented in the book's form.

<!-- results:begin -->
- [Weak coupling sandwich for the 2x2 and 4x4 Wilson loops](../RESULTS.md#weak_coupling_sandwich_conjecture) (Conjecture 14.1, page 243; open): The book conjectures that at weak coupling the logarithms of the 2x2 and 4x4 Wilson loops lie strictly between the independent plaquette value and the Migdal Kadanoff value, a statement supported by four numerical cases and not proved. New in this book.
- [Likelihood ratio sandwich for the 2x2 Wilson loop](../RESULTS.md#op_likelihood_ratio_sandwich_2x2) (Open problem 23, page 245; open): Prove that the law of the 2x2 loop holonomy dominates the Migdal Kadanoff law built from the plaquette law in the likelihood ratio order, and that both ends of the sandwich hold in that order. New in this book.
- [Bare Migdal Kadanoff end of the Wilson loop sandwich](../RESULTS.md#op_bare_migdal_kadanoff_end) (Open problem 24, page 251; open): Prove that the expectation of the 2x2 Wilson loop in SU(2) lattice gauge theory is at most the fourth power of a Bessel ratio at the bond moving coupling, in every dimension and at every coupling. New in this book.
<!-- results:end -->

## Key terms

sandwich, Wilson loop, 2x2 loop, blocked plaquette, block spin, decimation, bond moving, Migdal Kadanoff, likelihood ratio order, monotone likelihood ratio, Griffiths inequality, Ginibre inequality, FKG, reflection positivity, correlation inequality, string tension, heat bath, multihit, Rao Blackwell, pseudolikelihood, Monte Carlo renormalization group, pressure bound, Kadanoff bound, lattice gauge theory, SU(2), comparison inequality
