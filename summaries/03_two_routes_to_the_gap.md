# Chapter 3. Two Routes to the Gap

<!-- cite:begin -->
From *The Mass Gap: Thorin Measures, Brownian Motion and Yang Mills* by Nick Polson and Vadim Sokolov (2026), chapter 3. Cite as Polson and Sokolov (2026), *The Mass Gap: Thorin Measures, Brownian Motion and Yang Mills*, edition 1, https://doi.org/10.5281/zenodo.24546905. Book page: https://vsokolov.org/books/mass_gap/. Source: https://github.com/VadimSokolov/mass_gap.
<!-- cite:end -->

This chapter turns positivity into two translations of the mass gap into statements about generalized gamma convolutions (GGC), and tests them on solvable models. In Route A the offset spectral measure of a gauge invariant observable is itself a Thorin measure, so it defines a Yang Mills GGC, and the gap is the abscissa of exponential moments, the tail index of the GGC. In Route B the Polson 2 sin^2 lemma, from the authors' Riemann programme, transports the spectrum, regarded as a Levy measure, to a Thorin density built from a Levy exponent, namely the correlation at the offset minus the real part of the correlation continued to real time, divided by a multiple of real time. The gap is then the exact exponential decay rate of the Thorin densities in Euclidean time, with a monotone local rate. A uniform Thorin gap inequality is shown to be equivalent to a uniform spectral gap and to uniform exponential decay, with a free constant, and it passes to the continuum limit when the lattice offset measures converge weakly. A correlator with two isolated levels is not a GGC. The Ising chain with exact decimation, a massless field with power law correlations and the two dimensional Ising model (Onsager, Kaufman) illustrate a gap that survives, no gap, and a gap that closes at a critical point.

## From the literature

- Generalized gamma convolutions and Thorin measures (Thorin 1977); Bondesson's shot noise representation (Bondesson 1992); Kallen Lehmann spectral representation.
- Bernstein's theorem on completely monotone functions (Bernstein 1929).
- The one dimensional Ising model and its exact decimation (Ising chain, transfer matrix).
- Exact solution of the two dimensional Ising model, critical point and gap (Onsager 1944, Kaufman 1949).
- Glueball spectrum of lattice SU(3) (Morningstar and Peardon 1999); absence of a phase transition for SU(2) lattice gauge theory (Creutz 1980).

## New results

<!-- results:begin -->
- [The Yang Mills generalized gamma convolution of an observable](../RESULTS.md#yang_mills_ggc) (Proposition 3.1, page 30; proved): The offset spectral measure of a gauge invariant observable is a Thorin measure, so it defines a generalized gamma convolution whose Laplace transform is built from the zero momentum propagator. New in this book.
- [Mass gap equals the exponential moment threshold of the Yang Mills GGC](../RESULTS.md#gap_as_exponential_moment_threshold) (Theorem 3.2, page 31; proved): The mass gap of an observable is the exponential rate below which the associated generalized gamma convolution has finite exponential moments and above which it does not, so the gap is a tail index. New in this book.
- [The Polson 2 sin squared lemma](../RESULTS.md#polson_2sin2_lemma) (Lemma 3.5, page 33; proved): The compensated exponent of a single jump is a Thorin integral with an explicit nonnegative density built from a squared sine, the identity used to transport a Levy measure to a Thorin density. First stated in Polson and Sokolov (2026), Computing the Riemann Hypothesis, Volume I: Representations of Xi, Thorin Measures and the Reciprocal Xi Function, Zenodo.
- [The Thorin density is a Levy exponent](../RESULTS.md#thorin_density_is_levy_exponent) (Theorem 3.6, page 34; proved): The Thorin density obtained by transporting the spectrum is a Levy exponent divided by a multiple of real time: the exponent of a symmetric compound Poisson process, equal to the correlation at the offset minus the real part of the correlation continued to real time. New in this book.
- [Mass gap is the decay rate of the Thorin densities in Euclidean time](../RESULTS.md#gap_as_thorin_decay_rate) (Theorem 3.7, page 35; proved): The Thorin densities decay in Euclidean time at exactly the rate given by the gap, with a monotone local rate that decreases to the gap from above. New in this book.
- [Uniform Thorin gap inequality is equivalent to a uniform spectral gap, with constant one](../RESULTS.md#uniform_thorin_gap_equivalences) (Theorem 3.8, page 35; proved): A uniform inequality between Thorin densities at different Euclidean times is equivalent to a uniform spectral gap and to uniform exponential decay, and its constant can always be taken to be one. New in this book.
- [A uniform Thorin gap inequality gives a gap in the continuum limit](../RESULTS.md#continuum_gap_from_uniform_thorin_inequality) (Corollary 3.9, page 36; proved): If lattice offset measures converge weakly and satisfy the uniform Thorin gap inequality with one rate, the continuum theory has a gap at least that rate. New in this book.
- [A correlator with two isolated levels is not a GGC](../RESULTS.md#correlator_not_ggc_two_isolated_levels) (Proposition 3.10, page 37; proved): When the spectrum has two isolated levels, the normalised correlation in Euclidean time is not a generalized gamma convolution and not hyperbolically completely monotone. New in this book.
<!-- results:end -->

## Key terms

Route A, Route B, generalized gamma convolution, GGC, Thorin measure, Thorin density, Yang Mills GGC, Stieltjes function, exponential moment, abscissa of convergence, tail index, Levy exponent, compound Poisson process, 2 sin^2 lemma, Polson lemma, offset, Euclidean time, return amplitude, uniform Thorin gap inequality, mass gap, spectral gap, continuum limit, Ising chain, decimation, two dimensional Ising model, massless field, critical point, Riemann programme, hyperbolic completely monotone
