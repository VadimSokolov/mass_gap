# Chapter 7. Two Dimensions and Other Exact Models

<!-- cite:begin -->
From *The Mass Gap: Thorin Measures, Brownian Motion and Yang Mills* by Nick Polson and Vadim Sokolov (2026), chapter 7. Cite as Polson and Sokolov (2026), *The Mass Gap: Thorin Measures, Brownian Motion and Yang Mills*, edition 1, https://doi.org/10.5281/zenodo.24546905. Book page: https://vsokolov.org/books/mass_gap/. Source: https://github.com/VadimSokolov/mass_gap.
<!-- cite:end -->

This chapter collects the models in which the mass gap can be computed exactly, as a calibration for the four dimensional problem. The main case is two dimensional Yang Mills theory with gauge group SU(2), on the lattice and in the continuum, where the gap is a statement about Brownian motion on the three sphere: the transfer matrix is diagonal in characters, with eigenvalues given by ratios of modified Bessel functions. The chapter proves a uniform Thorin gap inequality for the two dimensional lattice theory from the coarsest lattice down to the continuum, using exact decimation and a Bessel ratio inequality. It then gives a guided tour of the method in two dimensions, an identity between Riemann's kernel for the xi function and the SU(2) sphere partition function (the return density of Brownian motion on S^3), and short sections on dynamic versus static gaps, the renormalization group for statisticians, free fields and Kallen Lehmann spectral measures, and large N with free probability (Gross Witten Wadia transition, free unitary Brownian motion).

## From the literature

- Exact solution of two dimensional lattice Yang Mills by decimation, and the partition function on a closed surface as a sum over representations (Migdal 1975), explained geometrically by Witten (1991).
- Continuum two dimensional Yang Mills measure (Gross, King and Sengupta 1989; Sengupta 1997; Levy 2003), and the master field (Levy 2017).
- Bessel ratio inequalities (Amos 1974).
- Riemann's cosine transform for the xi function (Riemann 1859), and the work of Polya, de Bruijn and Newman on its deformations (Polya 1926; de Bruijn 1950; Newman 1976).
- Theta sum laws and the Mellin transform of xi (Biane, Pitman and Yor 2001).
- The large N limit ('t Hooft 1974), hierarchical models (Dyson 1969; Bleher and Sinai 1973; Ito 1985).
- The third order large N transition (Gross and Witten 1980; Douglas and Kazakov 1993), and free unitary Brownian motion (Biane 1997).
- Equivalences between strong spatial mixing, spectral gap and log Sobolev inequality (Dobrushin and Shlosman 1985, 1987; Stroock and Zegarlinski 1992; Martinelli and Olivieri 1994), the Bakry Emery criterion (Bakry and Emery 1985), and the strong coupling dynamics of Shen, Zhu and Zhu (2023).

## New results

The chapter proves the uniform two dimensional gap and the Riemann kernel identity below.

<!-- results:begin -->
- [Uniform gap for two dimensional SU(2) lattice Yang Mills from the coarsest lattice to the continuum](../RESULTS.md#two_dimensional_uniform_gap) (Theorem 7.2, page 102; proved): In two dimensional SU(2) lattice Yang Mills on a cylinder the gap in physical units stays above an explicit positive multiple of the continuum gap at every lattice spacing, which gives a uniform Thorin gap inequality with constant one. New in this book.
- [Riemann's kernel from the SU(2) sphere partition function](../RESULTS.md#riemann_kernel_su2_sphere_partition_function) (Proposition 7.3, page 107; proved): Riemann's kernel for the xi function is a combination of the return density of Brownian motion on the three sphere and its time derivative, that is, of the two dimensional SU(2) Yang Mills partition function on the sphere and its area derivative. New in this book.
- [Fisher zeros of the two dimensional sphere theory and the zeros of Xi](../RESULTS.md#op_fisher_zeros_sphere_theory) (Open problem 5, page 109; open): Study the complex zeros in the area of the SU(2) sphere partition function and relate them to the zeros of the xi function through the identity with Riemann's kernel. New in this book.
- [Route B Thorin densities across the Douglas Kazakov transition](../RESULTS.md#op_douglas_kazakov_thorin_densities) (Open problem 6, page 109; open): Follow the Thorin densities of the second route across the large N phase transition of two dimensional Yang Mills on the sphere, and ask for a number theoretic reading of the large N analogue of the identity with Riemann's kernel. New in this book.
<!-- results:end -->

## Key terms

two dimensional Yang Mills, SU(2), Brownian motion on the three sphere, heat kernel, character expansion, decimation, transfer matrix, modified Bessel function
