# Chapter 11. Theta Kernels, Fermion Determinants and Howe Duality

<!-- cite:begin -->
From *The Mass Gap: Thorin Measures, Brownian Motion and Yang Mills* by Nick Polson and Vadim Sokolov (2026), chapter 11. Cite as Polson and Sokolov (2026), *The Mass Gap: Thorin Measures, Brownian Motion and Yang Mills*, edition 1, https://doi.org/10.5281/zenodo.24546905. Book page: https://vsokolov.org/books/mass_gap/. Source: https://github.com/VadimSokolov/mass_gap.
<!-- cite:end -->

This chapter asks which parts of a chain of theta lifts, used in an argument for the Riemann hypothesis announced in 2026 and read here through unpublished accounts, have a counterpart for lattice Yang Mills with gauge group SU(2); it gives no verdict on that argument. Plaquette laws with real zeros in the cosine of the angle (the Thorin class) are shown to be products of fermion determinants det(1+alpha U) with real hopping parameters, and their characters are theta lifts in the sense of Howe duality, given by Schur functions of the flavour group. Decimation, the convolution of two such laws, has an exact expansion over matchings of the complete graph on the flavours, which places the open closure conjecture for laws with real zeros next to the Heilmann and Lieb theorem. Three mechanisms of the expanded account that contain no number theory are translated: a one sided window criterion, which for a positive spectral measure becomes an elementary proposition with a monotone upper bound estimator of the gap built on negative binomial weights (a lattice Post and Widder kernel); a positive selector over channels; and an approximate tight frame lemma, given in a sharp form, which bounds averages such as string tensions and not gaps. The chapter ends with an inequality of type (M) for the renormalization ladder and open problems.

## From the literature

- Howe duality and the dual Cauchy identity for the pair of general linear groups (Howe 1989).
- The matching polynomial has only real zeros (Heilmann and Lieb 1972).
- The colour flavour transformation (Zirnbauer 1996) and its use for an induced lattice U(N) gauge theory (Budczies and Zirnbauer 2003).
- The Jacobi triple product, read as a free fermion.
- The inversion formula with gamma kernels (Post 1930, Widder 1941).
- The Harder and Narasimhan filtration (Harder and Narasimhan 1975), mentioned for the slope.
- The frame lemma, and its passage from a count to a mean, from an unpublished exposition written at the behest of M. O. Rubinstein (2026).

## New results

Only the statements about Yang Mills below are the authors' own.

<!-- results:begin -->
- [Laws with real zeros are products of SU(2) fermion determinants](../RESULTS.md#laws_with_real_zeros_are_determinants) (Proposition 11.1, page 180; proved): A plaquette law of the real zero class is a product of one flavour fermion determinants with real hopping parameters. New in this book.
- [Characters of a determinant law from the flavour group (Howe duality)](../RESULTS.md#characters_from_flavour_group) (Proposition 11.2, page 181; proved): The character coefficients of a product of fermion determinants are positive combinations of Schur functions of the hopping parameters. New in this book.
- [Convolution of class functions on SU(2) as an average over an ellipse](../RESULTS.md#convolution_ellipse_average) (Lemma 11.3, page 182; proved): The convolution of two class functions on SU(2) at cos theta equal to x is an area average over an ellipse in the plane of two cosines. New in this book.
- [Matching expansion for the convolution of determinant laws](../RESULTS.md#matching_expansion_of_decimation) (Theorem 11.4, page 182; proved): The convolution of two products of fermion determinants on SU(2) is an explicit sum over matchings of the flavours, a new exact formula for decimation. New in this book.
- [One sided windows and a monotone upper bound estimator for the gap](../RESULTS.md#windows_for_the_gap) (Proposition 11.5, page 185; proved): Adapting the window criterion of an unpublished account of the Riemann hypothesis argument, for a positive spectral measure a decay bound on a sparse set of times bounds the gap, and a negative binomial estimator gives monotone rigorous upper bounds. New in this book.
- [A sharp form of the frame lemma from an exposition of the Alpoge argument](../RESULTS.md#frame_lemma_sharp_form) (Proposition 11.7, page 187; proved): A sharpened constant in the frame lemma for sums of images of subspaces, attained for tight frames of isometries. New in this book.
- [An inequality of type (M) for the Migdal Kadanoff ladder](../RESULTS.md#op_type_m_inequality_for_ladder) (Open problem 11, page 191; open): Show that the exact blocked plaquette count stays below the ladder count plus a small error, until strong coupling is reached. New in this book.
- [Colour flavour transformation for SU(2) as a symplectic group](../RESULTS.md#op_symplectic_colour_flavour) (Open problem 13, page 193; open): Carry out the colour flavour transformation for SU(2) with its inner pairs and identify the mixing law. New in this book.
- [Bound the error per rung of the ladder below the drift](../RESULTS.md#op_slopes_error_per_rung) (Open problem 14, page 193; open): Bound the error per blocking step below the drift one half of the count, which with a bound on the first error would give the type (M) inequality. New in this book.
<!-- results:end -->

## Key terms

theta lift, theta kernel, Howe duality, dual Cauchy identity, Schur function, fermion determinant, hopping parameter, colour flavour transformation, matching expansion, matching polynomial, Heilmann and Lieb, real zeros, Thorin class, SU(2) decimation, symplectic group, window criterion, one sided window, negative binomial estimator, Post Widder, mass gap upper bound, approximate frame, frame lemma
