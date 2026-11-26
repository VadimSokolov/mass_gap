# Chapter 10. Real Zeros

<!-- cite:begin -->
From *The Mass Gap: Thorin Measures, Brownian Motion and Yang Mills* by Nick Polson and Vadim Sokolov (2026), chapter 10. Cite as Polson and Sokolov (2026), *The Mass Gap: Thorin Measures, Brownian Motion and Yang Mills*, edition 1, https://doi.org/10.5281/zenodo.24546905. Book page: https://vsokolov.org/books/mass_gap/. Source: https://github.com/VadimSokolov/mass_gap.
<!-- cite:end -->

This chapter asks on what class of link laws the exact Migdal Kadanoff ladder for SU(2) acts well, and proposes a Lee Yang property: the density of a link is a product over zeros of one plus cosine of the angle divided by a number at least one, so all zeros in cosine are real and lie beyond the edge of the group. Convolution of class functions on SU(2) is written as a Hadamard product of Laurent coefficients through the circle, and Malo's theorem gives the preservation of real zeros by heat flow and by convolution when one factor has an antiderivative with real zeros and at least the same degree, so the coalescence kernel has real zeros and a Polya frequency sequence of reweighted weights. The Hartman Watson law is shown to be a generalized gamma convolution, with Thorin distribution function minus the phase of the Bessel function of imaginary order divided by pi, so the Wilson action is a Thorin mixture of heat kernels. Closure of the class under convolution, which would give a Lee Yang property for the whole ladder, is conjectured, with numerical evidence, and is left open.

## From the literature

- the Lee Yang property for zeros of partition functions (Lee and Yang 1952)
- Malo's theorem on Hadamard products of real rooted polynomials (Malo 1895)
- the Schur Szego composition theorem (Szego 1922)
- the characterisation of Polya frequency sequences (Aissen, Schoenberg, Whitney and Edrei 1952)
- generalized gamma convolutions and Thorin measures (Thorin 1977, Bondesson 1992)
- the Hartman Watson law (Hartman and Watson 1974) and its appearance in windings of planar Brownian motion (Yor 1980)
- complete monotonicity and infinite divisibility for Bessel functions (Kent 1978, Hartman 1976)
- reality of zeros under heat flow (de Bruijn 1950); the reality of the zeros of K_{iy} in y (Polya 1926)
- Neumann's product formula and Gegenbauer's addition theorem for Bessel functions (Watson 1944)
- operators preserving real zeros (Borcea and Branden 2009; Asano 1970); the Grace, Walsh and Szego theorem (Walsh 1922, Szego 1922)

## New results

The preservation of the class by heat flow proved here is also what Chapter 6 uses to show that the vacuum of one plaquette has real zeros.

<!-- results:begin -->
- [Convolution on SU(2) as a Hadamard product of Laurent coefficients](../RESULTS.md#su2_convolution_as_hadamard_product) (Lemma 10.2, page 163; proved): Convolution of class functions on SU(2) is a Hadamard product of the Laurent coefficients of one function multiplied by z minus its inverse and of an antiderivative of the other, which brings theorems on real zeros of Hadamard products to the group. New in this book.
- [Real zeros are preserved by convolution and heat flow on SU(2)](../RESULTS.md#real_zeros_under_convolution) (Theorem 10.4, page 164; proved): Convolution of class functions on SU(2) preserves real zeros beyond the edge of the group when one factor has such an antiderivative, and heat flow preserves the class of link densities with real zeros. New in this book.
- [The coalescence has real zeros](../RESULTS.md#coalescence_has_real_zeros) (Corollary 10.5, page 165; proved): The convolution of two basic densities is a polynomial with only real zeros beyond the edge, so its count is a sum of independent Bernoulli variables and the kernel weights form a Polya frequency sequence. New in this book.
- [Two Wilson laws convolve into the real zero class](../RESULTS.md#wilson_times_wilson) (Proposition 10.6, page 165; proved): By Gegenbauer's addition theorem the convolution of two Wilson laws on SU(2) is a normalised Bessel function of order one, so its zeros in cosine are explicit, real and beyond the edge of the group, and its count is a sum of independent Bernoulli variables. New in this book.
- [Thorin form of the Hartman Watson law](../RESULTS.md#thorin_form_hartman_watson) (Theorem 10.7, page 168; proved): The Hartman Watson law is a generalized gamma convolution with no drift, and its Thorin distribution function is minus the phase of the Bessel function of imaginary order divided by pi. New in this book.
- [The Wilson action is a Thorin mixture of heat kernels](../RESULTS.md#wilson_action_thorin_mixture_of_heat_kernels) (Corollary 10.8, page 169; proved): The Wilson law on SU(2) is a mixture of heat kernels over a random time that is a generalized gamma convolution, with Thorin measure the Hartman Watson one shifted by one half. New in this book.
- [Smoothed criterion for real zeros of a convolution](../RESULTS.md#smoothed_criterion_real_zeros) (Theorem 10.10, page 172; proved): The convolution of two polynomial laws with real zeros has real zeros whenever one of them convolved with the basic density of the other's degree has an antiderivative with real zeros. New in this book.
- [Smoothing conjecture for real zeros](../RESULTS.md#smoothing_conjecture) (Conjecture 10.11, page 173; open): Conjecturally, smoothing a polynomial law with real zeros by a basic density of no larger degree gives a law with a real zeroed antiderivative; it would imply closure of the class under convolution for polynomials. New in this book.
- [Radial flow identities for SU(2) convolution](../RESULTS.md#radial_flow_identities) (Theorem 10.12, page 174; proved): Convolution on SU(2) is a running average over the radius of convolutions on the circle of dilated laws, with explicit identities linking the two. New in this book.
- [Closure follows from the absence of a common zero along the flow](../RESULTS.md#closure_from_no_common_zero) (Theorem 10.14, page 175; proved): Convolution of two polynomial laws with real zeros has real zeros if the SU(2) and circle convolutions of their dilations never share a zero. New in this book.
- [A criterion by signs for real zeros of a convolution](../RESULTS.md#sign_criterion_real_zeros) (Theorem 10.15, page 176; proved): Real zeros of the convolution follow if all lower convolutions have alternating signs at the zeros of the circle convolution. New in this book.
- [Interlacing conjecture for circle and SU(2) convolution](../RESULTS.md#interlacing_conjecture) (Conjecture 10.16, page 176; open): Conjecturally the zeros of the circle convolution strictly interlace those of the SU(2) convolution; it would imply closure of the class under convolution. New in this book.
- [Conjecture that real zero laws are closed under convolution on SU(2)](../RESULTS.md#class_l_closed_under_convolution) (Conjecture 10.17, page 177; open): Conjecturally (supported by numerical evidence, not proved) the class of link laws with only real zeros beyond the edge of the group is closed under convolution, so the exact ladder has a Lee Yang property at every rung. New in this book.
- [A Malo type theorem for the SU(2) pairing](../RESULTS.md#malo_type_theorem_problem) (Open problem 10, page 178; open): Open problem to prove the smoothing or interlacing conjecture, or a Malo type theorem for the pairing of the circle lemma when neither factor has an antiderivative in the class. New in this book.
<!-- results:end -->

## Key terms

real zeros, Lee Yang property, Malo's theorem, Hadamard product, Polya frequency sequence, Schur Szego composition, generalized gamma convolution, Thorin measure, Hartman Watson law, heat flow, class L, closure under convolution
