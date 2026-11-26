# Chapter 9. The Ladder on Counts and Kingman's Coalescent

<!-- cite:begin -->
From *The Mass Gap: Thorin Measures, Brownian Motion and Yang Mills* by Nick Polson and Vadim Sokolov (2026), chapter 9. Cite as Polson and Sokolov (2026), *The Mass Gap: Thorin Measures, Brownian Motion and Yang Mills*, edition 1, https://doi.org/10.5281/zenodo.24546905. Book page: https://vsokolov.org/books/mass_gap/. Source: https://github.com/VadimSokolov/mass_gap.
<!-- cite:end -->

This chapter studies the exact Migdal Kadanoff ladder for SU(2) lattice gauge theory without the projection onto heat kernels. Heat kernels and Wilson laws are mixtures of the densities proportional to powers of one plus cosine of the angle, with a latent nonnegative integer count. Bond moving adds counts and decimation, the convolution of two such densities, thins them by an explicit positive kernel. The coalescence theorem shows this kernel is strictly positive, and the hypergeometric theorem identifies it as the law of the overlap of two random subsets, so that positivity is the Chu Vandermonde identity and the kernel is totally positive of order two. The harmonic law gives the mean count of a convolution, with the constant seven sixths, and hence the drift of minus thirteen sixteenths per rung from a point mass. The mean count is shown to be an observable, the expectation of the reciprocal of one minus cosine less two. The first rung theorem, stated on the class of laws with real zeros, gives the change of the mean count in terms of the variance to mean ratio, with explicit limits only for point masses and Wilson laws; the heat kernel value of minus one half is not proved. The second half shows that Brownian motion on SU(2) is a Wright Fisher diffusion and that heat flow is dual to the death process of Kingman's coalescent.

## From the literature

- the coalescent (Kingman 1982) and the random walks with spherical symmetry of Kingman (1963)
- the subadditive ergodic theorem (Kingman 1973)
- the expansion of the Wright Fisher transition density in Beta laws weighted by the ancestral process (Griffiths 1980, Tavare 1984)
- exact filters for the Wright Fisher diffusion through the dual process (Chaleyat Maurel and Genon Catalot 2009, Papaspiliopoulos and Ruggiero 2014)
- total positivity of birth and death transition functions (Karlin and McGregor 1959)
- the ancestral selection graph (Krone and Neuhauser 1997)
- the Ewens sampling formula (Ewens 1972)

## New results

The chapter's own results are the coalescence theorem, the hypergeometric kernel, the harmonic law, the first rung theorem, and the transfer of the Griffiths and Tavare duality to the SU(2) heat kernel; the open problems on log concavity and on the constant drift are also catalogued.

<!-- results:begin -->
- [Coalescence of Beta type densities under SU(2) convolution](../RESULTS.md#coalescence_of_counts) (Theorem 9.1, page 139; proved): Convolving two of the basic densities proportional to powers of one plus cosine on SU(2) gives a strictly positive mixture of such densities, with a known mean count. New in this book.
- [Hypergeometric law of coalescence](../RESULTS.md#coalescence_is_hypergeometric) (Theorem 9.4, page 141; proved): The coalescence kernel of the SU(2) convolution ladder is the hypergeometric law of the overlap of two random subsets, and its positivity is the Chu Vandermonde identity. New in this book.
- [Moments, total positivity and unnormalised form of the coalescence kernel](../RESULTS.md#hypergeometric_kernel_consequences) (Corollary 9.5, page 142; proved): The coalescence kernel has explicit mean and variance, is totally positive of order two, has an unnormalised closed form, and is the thinning of the entrance law of Kingman's death process. New in this book.
- [The exact Migdal Kadanoff ladder stays in the family of count mixtures](../RESULTS.md#mixture_family_closed_under_ladder) (Corollary 9.6, page 143; proved): Heat kernels and Wilson laws on SU(2) are mixtures of the basic densities, and the family of such mixtures is closed under products and convolutions, so every rung of the exact ladder stays in it. New in this book.
- [The mean count as an expectation of a reciprocal and a character sum](../RESULTS.md#mean_count_formula) (Corollary 9.7, page 145; proved): The mean count of a mixture equals the average of the reciprocal of one minus cosine, less two, and for the heat kernel it is twice the inverse coupling minus seven sixths to leading order. New in this book.
- [Harmonic law for the mean count of a convolution](../RESULTS.md#harmonic_law) (Theorem 9.8, page 147; proved): Under convolution the mean counts combine like resistances in parallel, with an explicit constant correction and an error of order one over the smallest index. New in this book.
- [Drift of the mean count from a point mass over one exact rung](../RESULTS.md#drift_from_point_mass) (Corollary 9.9, page 148; proved): One exact rung of the ladder started from a pure power density lowers the mean count by a constant, minus thirteen sixteenths for four copies and fourfold convolution. New in this book.
- [The latent count is an observable on the group](../RESULTS.md#count_is_an_observable) (Proposition 9.11, page 148; proved): The mean of the latent count of a link law equals the expectation of the reciprocal of the Wilson action density, less two. New in this book.
- [First rung theorem for laws with real zeros](../RESULTS.md#first_rung_theorem) (Theorem 9.13, page 150; proved): For link laws in the class with real zeros, one exact ladder rung changes the mean count by a constant depending only on the ratio of variance to mean, with explicit limits for point masses and Wilson laws; the heat kernel value is not proved. New in this book.
- [Finite variance estimator of the count on a lattice](../RESULTS.md#link_integrated_count_estimator) (Proposition 9.14, page 151; proved): Integrating out one link gives a bounded function of the configuration whose average is the count of a loop through that link. New in this book.
- [Laws with real zeros have a log concave count](../RESULTS.md#real_zero_laws_log_concave_count) (Proposition 9.15, page 154; proved): The latent count law of every link density with only real zeros beyond the edge of the group is log concave. New in this book.
- [Does the exact ladder preserve log concavity of the count](../RESULTS.md#ladder_preserves_log_concavity_problem) (Open problem 8, page 155; open): Open whether coalescence preserves log concavity of the latent count; true for laws of the real zero class, supported by exact computation otherwise. New in this book.
- [Explain the constant drift of the mean count](../RESULTS.md#constant_drift_problem) (Open problem 9, page 155; open): Open to prove the drift of minus one half per rung for heat kernels and at later rungs; the first rung from the real zero class is proved. New in this book.
- [Brownian motion on SU(2) is a Wright Fisher diffusion](../RESULTS.md#su2_heat_flow_is_wright_fisher) (Proposition 9.16, page 156; proved): In the variable one half of one plus cosine, the class angle of Brownian motion on SU(2) is the neutral two allele Wright Fisher diffusion with total mutation rate three, started from fixation. New in this book.
- [Kingman duality for the SU(2) heat kernel](../RESULTS.md#kingman_coalescent_duality) (Theorem 9.17, page 157; proved): Through the Wright Fisher form of Brownian motion on SU(2), the expansion of Griffiths and Tavare in Beta laws weighted by the ancestral process becomes an expansion of the SU(2) heat kernel in powers of one plus cosine, with Kingman's death process moving the weights. New in this book.
<!-- results:end -->

## Key terms

latent count, coalescence, hypergeometric distribution, Chu Vandermonde identity, total positivity, harmonic law, mean count, bond moving, decimation, Kingman coalescent, Wright Fisher diffusion, coalescent duality
