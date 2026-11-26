# Chapter 6. Hitting Times and Learning the Gap

<!-- cite:begin -->
From *The Mass Gap: Thorin Measures, Brownian Motion and Yang Mills* by Nick Polson and Vadim Sokolov (2026), chapter 6. Cite as Polson and Sokolov (2026), *The Mass Gap: Thorin Measures, Brownian Motion and Yang Mills*, edition 1, https://doi.org/10.5281/zenodo.24546905. Book page: https://vsokolov.org/books/mass_gap/. Source: https://github.com/VadimSokolov/mass_gap.
<!-- cite:end -->

A gap is a rate; this chapter reads it in two further ways and asks how it can be learned. First, the hazard rate of a hitting time (Roberts and Shortland): standard bounds place the spectral gap of a reversible diffusion within a factor two of Dirichlet eigenvalues of sets of measure at least one half; a reflection removes the factor; Brownian motion on SU(2) hitting the equator has an exact survival law with rate 3/2; and a hazard rate estimator of the gap uses positive data. Second, moderate deviations (Rubin and Sethuraman), the Cramer series, the one loop term of the Migdal Kadanoff ladder as curvature, and, under a stated hypothesis, de Finetti's theorem as a necessary test for a unique vacuum. Third, learning: sufficient statistics of lattice gauge fields, particle learning, and generative Bayesian computation (GBC) with neural and transformer generators. Fourth, quantum Bayesian computation: the Kogut Susskind Hamiltonian of one SU(2) plaquette as a tridiagonal matrix in characters, a proof that its vacuum has only real zeros, and the limits of quantum computation. None of this proves the gap in four dimensions.

## From the literature

- Hazard rate tangent approximation for hitting times (Roberts and Shortland 1995, 1997); Dirichlet eigenvalue bounds on the spectral gap, standard for Markov chains; Cheeger and Poincare inequalities; stability of Markov chains (Meyn and Tweedie 1993).
- Moderate deviations and Bayes risk efficiency (Rubin and Sethuraman 1965); Cramer series; large deviations of exchangeable vectors (Dinwoodie and Zabell 1992).
- Exponential families, Rao Blackwell, Pitman Koopman theorem, Cramer Rao bound.
- Particle learning (Carvalho, Johannes, Lopes and Polson 2010), generative Bayesian computation (Polson and Sokolov 2023) and quantum Bayesian computation (Polson, Sokolov and Xu 2023), the authors' earlier work.
- Hamiltonian lattice gauge theory (Kogut and Susskind 1975); phase estimation (Kitaev 1995); adiabatic quantum computation (Farhi, Goldstone, Gutmann and Sipser 2000); Peter and Weyl theorem; Laguerre and Polya theorem; Hurwitz's theorem.

## New results

<!-- results:begin -->
- [Gap of the odd sector as a Dirichlet eigenvalue](../RESULTS.md#reflection_dirichlet_odd_sector) (Proposition 6.2, page 76; proved): When a measure preserving reflection splits the state space into two halves, the lowest odd eigenvalue is the Dirichlet eigenvalue of the fixed set, which removes the factor of two in the hitting time bracket. New in this book.
- [Exact survival law of Brownian motion on SU(2) hitting the equator](../RESULTS.md#su2_equator_hitting_time) (Proposition 6.3, page 76; proved): The survival probability of Brownian motion on SU(2) until it reaches the equator is an explicit series whose hazard rate tends to the spectral gap 3/2. New in this book.
- [Gap of the one plaquette SU(2) Hamiltonian](../RESULTS.md#one_plaquette_hamiltonian_gap) (Proposition 6.6, page 91; proved): The single plaquette Kogut Susskind Hamiltonian is a tridiagonal matrix in the character basis, with gap 3 lambda/4 at strong coupling and tending to the square root of 2 at weak coupling. New in this book.
- [The vacuum of one plaquette has real zeros](../RESULTS.md#one_plaquette_vacuum_real_zeros) (Proposition 6.7, page 92; proved): The ground state of the one plaquette Hamiltonian, as a function of the cosine of the angle, has only real zeros at or below minus one, so its reciprocal is a Laplace transform of a generalized gamma convolution (Thorin form). New in this book.
- [Log concave count law of the one plaquette vacuum](../RESULTS.md#one_plaquette_vacuum_count_log_concave) (Corollary 6.9, page 94; proved): The vacuum density of one plaquette is a mixture of the count densities with a log concave mixing law whose mean is the count reported in the book's table. New in this book.
- [Real zeros of the vacuum of two plaquettes sharing a link](../RESULTS.md#two_plaquettes_real_zeros) (Open problem 2, page 96; open): Open whether the vacuum of two SU(2) plaquettes sharing a link, as a function of each plaquette variable, again has only real zeros below minus one. New in this book.
- [Register rule for the one plaquette gap](../RESULTS.md#register_rule_one_plaquette) (Open problem 3, page 96; open): Open whether the number of characters needed for the one plaquette gap grows like the square root of the vacuum count as the coupling goes to zero. New in this book.
- [Monotonicity of the Thorin atoms of the one plaquette vacuum](../RESULTS.md#vacuum_thorin_atoms_monotone) (Open problem 4, page 96; open): Open whether each Thorin atom of the one plaquette vacuum decreases monotonically toward one as the coupling goes to zero. New in this book.
<!-- results:end -->

## Key terms

hitting time, hazard rate, first passage time, Dirichlet eigenvalue, spectral gap, Poincare inequality, reversible diffusion, Brownian motion on SU(2), moderate deviations, Cramer series, exchangeability, de Finetti, sufficient statistic, exponential family, Fisher information, Cramer Rao bound, particle learning, generative Bayesian computation, GBC, transformer, quantum Bayesian computation, Kogut Susskind Hamiltonian, one plaquette, character truncation, real zeros, generalized gamma convolution, Thorin measure, phase estimation, adiabatic quantum computation
