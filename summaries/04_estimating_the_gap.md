# Chapter 4. Estimating the Gap

<!-- cite:begin -->
From *The Mass Gap: Thorin Measures, Brownian Motion and Yang Mills* by Nick Polson and Vadim Sokolov (2026), chapter 4. Cite as Polson and Sokolov (2026), *The Mass Gap: Thorin Measures, Brownian Motion and Yang Mills*, edition 1, https://doi.org/10.5281/zenodo.24546905. Book page: https://vsokolov.org/books/mass_gap/. Source: https://github.com/VadimSokolov/mass_gap.
<!-- cite:end -->

This chapter treats the mass gap as a parameter to be estimated from a noisy Euclidean correlator. Recovering the spectral measure is the ill conditioned inverse Laplace transform problem (spectral reconstruction in lattice QCD), but the left edge of its support is well determined. The chapter builds a support constrained profile estimator: positive deconvolution by nonnegative least squares, a chi squared profile in the candidate gap, and a threshold rule, tested on a toy spectrum with three stable states and a two particle continuum, with a Bayesian version using a Laplace approximation and horseshoe priors. It then reviews Monte Carlo tools from Bayesian computation for lattice gauge theory: latent variable samplers, slice sampling heat bath, nested sampling and the density of states, the vertical likelihood and split sampling of Polson and Scott compared on a two dimensional SU(2) torus with exact partition function, and weighted Monte Carlo with the Yakowitz rate tested on thermodynamic integration. It can expose a wrong claim and measure a lattice gap, but cannot prove that a gap exists or reach the continuum limit.

## From the literature

- Backus Gilbert method (1968), maximum entropy, and smeared spectral densities for spectral reconstruction (Hansen, Lupo and Tantalo 2019).
- The EM algorithm for positive linear inverse problems (Vardi and Lee 1993), read as Monte Carlo inversion (Polson and Robert 1994); Bayesian nonparametric deconvolution (Padilla, Polson and Scott).
- Nonnegative least squares (Lawson and Hanson 1974); horseshoe priors (Carvalho, Polson and Scott 2010).
- Power posteriors (Friel and Pettitt 2008); the harmonic mean estimator (Newton and Raftery 1994); Troutt's vertical density sampler (1991); Wang and Landau sampling (2001); the 1/k ensemble (Hesselbo and Stinchcombe 1995); nested sampling (Skilling 2006).
- Vertical likelihood Monte Carlo and weighted slice sampling (Polson and Scott 2014, the first author's earlier work), which the chapter states and applies.
- Weighted Monte Carlo (Yakowitz, Krimmel and Szidarovszky 1978) and its mean square error rate of order n to the minus four (Masry and Cambanis 1990).

## New results

<!-- results:begin -->
No new results of this chapter are catalogued.
<!-- results:end -->

## Key terms

mass gap estimation, spectral reconstruction, inverse Laplace transform, positive deconvolution, nonnegative least squares, profile likelihood, support constrained estimator, Vardi Lee EM, Backus Gilbert, nested sampling, density of states, vertical likelihood, split sampling, weighted slice sampling, power posterior, thermodynamic integration, harmonic mean estimator, Wang Landau, partition function, Yakowitz rate, weighted Monte Carlo, trapezoid rule
