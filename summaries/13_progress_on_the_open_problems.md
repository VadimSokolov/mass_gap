# Chapter 13. Progress on the Open Problems

<!-- cite:begin -->
From *The Mass Gap: Thorin Measures, Brownian Motion and Yang Mills* by Nick Polson and Vadim Sokolov (2026), chapter 13. Cite as Polson and Sokolov (2026), *The Mass Gap: Thorin Measures, Brownian Motion and Yang Mills*, edition 1, https://doi.org/10.5281/zenodo.24546905. Book page: https://vsokolov.org/books/mass_gap/. Source: https://github.com/VadimSokolov/mass_gap.
<!-- cite:end -->

This chapter reports the authors' progress on the open problems of the programme for the SU(2) Migdal Kadanoff renormalization ladder. The SU(2) heat kernel is a product of linear factors in the cosine of the angle (a Laguerre Polya class function) whose roots gamma_m(t) exceed one and move monotonically in time, and its reciprocal has a Thorin form with atoms at these roots. Hence the heat kernel has strict monotone likelihood ratio in time at every time and is strictly log concave in the cosine, and the projected map is monotone. Further results: an elementary large time bound; monotone likelihood ratio of the heat flow of the Wilson law exactly for kappa at most three halves; the weak coupling expansion of the ladder decrement for SU(2), and its first two orders for SU(N); a local heat kernel time; and total positivity of order two of convolution, proved for bond moved laws on the near hemisphere and failing, numerically, in a small cap at the antipode. The projected ladder contracts at every coupling: the time scaled score t d log K_t/dx decreases in t, by Hamilton's matrix Harnack inequality and an explicit small time analysis, so R(beta) is below beta for every beta and scale factor in four dimensions. Existence of the ladder string tension is proved under stated hypotheses and tested against lattice data; the exact ladder and the true lattice flow remain open.

## From the literature

- The Laguerre Polya class and multiplier sequences (Polya and Schur 1914), and the Polya class of the Riemann Xi function (Polya 1926).
- Monotone motion of zeros under the heat flow (de Bruijn 1950, Newman 1976).
- Total positivity of one dimensional diffusions (Karlin and McGregor 1959).
- The matrix Harnack inequality for the heat equation (Hamilton 1993).
- The Wilson action (Wilson 1974), the effective coupling (Parisi 1980), the Migdal Kadanoff recursion (Migdal 1975, Kadanoff 1976), lattice string tension data (Teper 1998, 1999) and the Biane, Pitman and Yor law (2001).

## New results

The headline results are the product formula and the contraction of the projected ladder.

<!-- results:begin -->
- [Product formula for the SU(2) heat kernel](../RESULTS.md#product_formula) (Theorem 13.1, page 213; proved): The SU(2) heat kernel is a product of linear factors in cos theta whose real roots all lie beyond the edge of the group and move monotonically in time. New in this book.
- [Monotone likelihood ratio in time for the SU(2) heat kernel](../RESULTS.md#likelihood_ratio_property) (Corollary 13.2, page 214; proved): The SU(2) heat kernel has strictly monotone likelihood ratio in time at every time, so the projected Migdal Kadanoff map is increasing at every coupling. New in this book.
- [Strict log concavity of the SU(2) heat kernel and its powers in cos theta](../RESULTS.md#heat_kernel_log_concavity) (Corollary 13.3, page 214; proved): The heat kernel, and every plaquette law made from it by bond moving, is strictly log concave as a function of the cosine of the angle. New in this book.
- [Thorin form of the SU(2) heat kernel](../RESULTS.md#heat_kernel_thorin_form) (Corollary 13.4, page 214; proved): The reciprocal of the SU(2) heat kernel, as a function of the cosine of the angle, is the Laplace transform of a generalized gamma convolution with unit atoms at the roots of the product formula. New in this book.
- [Elementary proof of monotone likelihood ratio at large time](../RESULTS.md#elementary_bound_large_time) (Proposition 13.5, page 215; proved): A direct bound on the character series shows the time monotonicity of the score for all times beyond 3.2140. New in this book.
- [Monotone likelihood ratio of the heat flow of the Wilson law](../RESULTS.md#wilson_law_likelihood_ratio_in_time) (Proposition 13.6, page 217; proved): The heat flow of the Wilson plaquette law has monotone likelihood ratio in time exactly when the coupling parameter is at most three halves. New in this book.
- [Weak coupling expansion of the SU(2) projected ladder](../RESULTS.md#weak_coupling_expansion_ladder) (Theorem 13.7, page 223; proved): The projected ladder map and the decrement of the coupling have an explicit weak coupling expansion, giving one quarter plus one over 48 beta at scale factor two. New in this book.
- [Weak coupling expansion of the projected ladder for SU(N)](../RESULTS.md#weak_coupling_expansion_sun) (Theorem 13.8, page 225; proved): For every N the first two orders of the ladder decrement are given in closed form, with two loop ratio N over 12 b squared. New in this book.
- [Bracketing a plaquette law between heat kernels in likelihood ratio order](../RESULTS.md#local_heat_kernel_time_bracket) (Proposition 13.9, page 227; proved): Any class function law with negative score lies between two heat kernel laws in likelihood ratio order, and the bracket survives bond moving and heat convolution. New in this book.
- [Total positivity of order two of convolution on SU(2)](../RESULTS.md#total_positivity_of_convolution) (Theorem 13.10, page 228; proved): A condition on the radial Laplacian ratio of a class function density makes convolution with it preserve likelihood ratio order. New in this book.
- [Heat kernels and Wilson laws satisfying the total positivity hypothesis](../RESULTS.md#tp2_heat_kernel_wilson) (Corollary 13.11, page 229; proved): Heat kernels always, and Wilson laws exactly up to kappa three halves, have totally positive convolution kernels, and the class is closed under convolution. New in this book.
- [Total positivity of bond moved heat kernel laws on the near hemisphere](../RESULTS.md#near_hemisphere_total_positivity) (Theorem 13.12, page 230; proved): Convolution with a bond moved heat kernel law has a totally positive kernel on angles whose sum is at most a right angle. New in this book.
- [Truncated monotonicity of a kernel that is totally positive on a box](../RESULTS.md#truncated_monotonicity) (Proposition 13.13, page 231; proved): Failure of total positivity outside a box costs an error bounded by the mass outside the box. New in this book.
- [Time scaled likelihood ratio of the SU(2) heat kernel](../RESULTS.md#time_scaled_likelihood_ratio) (Theorem 13.14, page 232; proved): The heat kernel score multiplied by time is strictly decreasing in time at every point, the key inequality behind the contraction of the ladder. New in this book.
- [The projected ladder contracts at every coupling](../RESULTS.md#the_ladder_contracts) (Corollary 13.17, page 235; proved): For the heat kernel projection of the Migdal Kadanoff map the coupling strictly decreases at every rung, for every coupling and scale factor in four dimensions. New in this book.
- [Existence of the string tension of the ladder](../RESULTS.md#string_tension_on_hierarchical_lattice_exists) (Proposition 13.18, page 236; proved): If the ladder reaches strong coupling with geometrically bounded higher coefficients, the rescaled string tension limit exists and is positive. New in this book.
<!-- results:end -->

## Key terms

product formula, SU(2) heat kernel, Thorin form, Laguerre Polya class, real zeros, monotone likelihood ratio, log concavity, Wilson plaquette law, weak coupling expansion, decrement, two loop coefficient, SU(N), total positivity of order two, bond moving, decimation, Migdal Kadanoff ladder, Hamilton matrix Harnack inequality, projected ladder, string tension, asymptotic scaling, effective coupling, hierarchical lattice, antipode cap
