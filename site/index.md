# The Mass Gap <small>Thorin Measures, Brownian Motion and Yang Mills</small>

By Nick Polson and [Vadim Sokolov](https://orcid.org/0000-0002-6618-2965).

[![PDF](https://vsokolov.org/books/icon-pdf.svg)](https://github.com/VadimSokolov/mass_gap/releases/download/v1.0.0/mass_gap.pdf "PDF") [![Source on GitHub](https://vsokolov.org/books/icon-github.svg)](https://github.com/VadimSokolov/mass_gap/tree/v1.0.0 "Source on GitHub") [![DOI 10.5281/zenodo.24546905](https://zenodo.org/badge/DOI/10.5281/zenodo.24546905.svg)](https://doi.org/10.5281/zenodo.24546905)

## Description

The Mass Gap reads the Yang Mills mass gap problem through probability. Reflection positivity makes the spectral measure of a gauge field positive, and the book develops what that positivity buys: the gap as the exponential moment threshold of a generalized gamma convolution, and as the decay rate of its Thorin densities; Brownian motion on SU(2) and the exactly solvable theory in two dimensions; Monte Carlo estimation of the gap; and the Migdal Kadanoff renormalization ladder as a flow on laws of the plaquette. The new results are proved for the ladder, not for the four dimensional lattice. They include a product formula for the SU(2) heat kernel with monotone roots, whose reciprocal is the Laplace transform of a generalized gamma convolution; the contraction of the projected ladder at every coupling; a hypergeometric law for the coalescence of the latent count, with a harmonic law for its mean and a first rung theorem; the Thorin form of the Hartman Watson law; and the preservation of real zeros by convolution and the heat flow on SU(2). The book poses a conjecture, that link densities with real zeros beyond the edge of the group are closed under convolution, and reduces it to statements about one function at a time. Every computation is reported with its script, and theorems, leading order derivations, numerical results and open problems are kept apart. The mass gap in four dimensions is not proved here; the chapters on open problems say what is missing.

## How to cite

Polson and Sokolov (2026), *The Mass Gap: Thorin Measures, Brownian Motion and Yang Mills*, edition 1, https://doi.org/10.5281/zenodo.24546905.

```bibtex
@book{polsonsokolov2026mass_gap,
  author = {Polson, Nick and Sokolov, Vadim},
  title = {{The Mass Gap: Thorin Measures, Brownian Motion and Yang Mills}},
  year = {2026},
  publisher = {Independently published},
  doi = {10.5281/zenodo.24546905},
  url = {https://vsokolov.org/books/mass_gap/}
}
```

## Main results

- [Hypergeometric law of coalescence](https://vsokolov.org/books/mass_gap/results.html#coalescence_is_hypergeometric) (Theorem 9.4, page 141): The coalescence kernel of the SU(2) convolution ladder is the hypergeometric law of the overlap of two random subsets, and its positivity is the Chu Vandermonde identity.
- [Harmonic law for the mean count of a convolution](https://vsokolov.org/books/mass_gap/results.html#harmonic_law) (Theorem 9.8, page 147): Under convolution the mean counts combine like resistances in parallel, with an explicit constant correction and an error of order one over the smallest index.
- [First rung theorem for laws with real zeros](https://vsokolov.org/books/mass_gap/results.html#first_rung_theorem) (Theorem 9.13, page 150): For link laws in the class with real zeros, one exact ladder rung changes the mean count by a constant depending only on the ratio of variance to mean, with explicit limits for point masses and Wilson laws; the heat kernel value is not proved.
- [Real zeros are preserved by convolution and heat flow on SU(2)](https://vsokolov.org/books/mass_gap/results.html#real_zeros_under_convolution) (Theorem 10.4, page 164): Convolution of class functions on SU(2) preserves real zeros beyond the edge of the group when one factor has such an antiderivative, and heat flow preserves the class of link densities with real zeros.
- [Thorin form of the Hartman Watson law](https://vsokolov.org/books/mass_gap/results.html#thorin_form_hartman_watson) (Theorem 10.7, page 168): The Hartman Watson law is a generalized gamma convolution with no drift, and its Thorin distribution function is minus the phase of the Bessel function of imaginary order divided by pi.
- [Product formula for the SU(2) heat kernel](https://vsokolov.org/books/mass_gap/results.html#product_formula) (Theorem 13.1, page 213): The SU(2) heat kernel is a product of linear factors in cos theta whose real roots all lie beyond the edge of the group and move monotonically in time.
- [The projected ladder contracts at every coupling](https://vsokolov.org/books/mass_gap/results.html#the_ladder_contracts) (Corollary 13.17, page 235): For the heat kernel projection of the Migdal Kadanoff map the coupling strictly decreases at every rung, for every coupling and scale factor in four dimensions.
- [Thorin form of the SU(2) heat kernel](https://vsokolov.org/books/mass_gap/results.html#heat_kernel_thorin_form) (Corollary 13.4, page 214): The reciprocal of the SU(2) heat kernel, as a function of the cosine of the angle, is the Laplace transform of a generalized gamma convolution with unit atoms at the roots of the product formula.

All 97 results of the book (70 proved, 27 open) are on the page [Results by chapter](https://vsokolov.org/books/mass_gap/results.html), each with its statement, status and number.

***

<small>License: Text, figures and PDF: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Code: [MIT](https://opensource.org/license/mit).</small>

<small>This page is generated from the book's repository ([main](https://github.com/VadimSokolov/mass_gap/tree/main)). The numbers and pages are those of the PDF, [version 1.0.0](https://github.com/VadimSokolov/mass_gap/releases/tag/v1.0.0).</small>

<!-- Generated by bookrepo.py from book.yaml and results.yaml. Edit those files, not this one. -->
