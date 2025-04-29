# BaCoN

Buffering between genes, where one gene can compensate for the loss of another gene, is fundamental for robust cellular functions. While experimentally testing all possible gene pairs is infeasible, gene buffering can be predicted genome-wide under the assumption that a gene’s buffering capacity depends on its expression level and its absence primes a severe fitness phenotype of the buffered gene. We developed BaCoN (Balanced Correlation Network), a post hoc unsupervised correction method that amplifies specific signals in expression-vs-fitness correlation networks. We quantified 147 million potential buffering relationships by associating CRISPR-Cas9-screening fitness effects with transcriptomic data across 1019 Cancer Dependency Map (DepMap) cell lines. BaCoN outperformed state-of-the-art methods, including multiple linear regression based on our compiled gene buffering prediction metrics. Combining BaCoN with batch correction or Cholesky data whitening further boosts predictive performance. We characterized 808 high-confidence buffering predictions and found that in contrast to buffering gene pairs overall, buffering paralogs were on different chromosomes. BaCoN performance increases with more screens and genes considered, making it a valuable tool for gene buffering predictions from the growing DepMap.

This is the official `R` package containing the functions required to compute a `BaCoN` matrix. 

# Quickstart

The quickest way to compute a `BaCoN` matrix is:

1. Install `devtools` from CRAN:

```{r}
if (!require("devtools", quietly = TRUE)) {
	install.packages("devtools")
	}
```

2. Install and load the `BaCoN` package into your workspace:

```{r}
devtools::install_github("billmannlab/BaCoN")

library(BaCoN)
```

3. Compute a BaCoN matrix:

```{r}
bacon_matrix <- BaCoN(correlation_matrix)
```

NOTE: 

In our manuscript, we suggest to combine `BaCoN` with [Cholesky whitening](https://doi.org/10.1080/00031305.2016.1277159) 
to achieve better results. 

Please refer to the [vignette](https://github.com/billmannlab/BaCoN/blob/main/vignettes/tutorial.pdf) for a tutorial. 

# Tutorial

A vignette describing how to predict buffering gene pairs using `BaCoN` can be found [here](https://github.com/billmannlab/BaCoN/blob/main/vignettes/tutorial.pdf). 

# Citation

To cite this package, please refer to:

- Rohde, T., Demirtas, Y., Süsser, S., Shaw, A., Kaulich, M., Billmann, M. `BaCoN` (Balanced Correlation Network) improves prediction of gene buffering. _Molecular Systems Biology_, (2025). DOI: [10.1038/s44320-025-00103-7](https://doi.org/10.1038/s44320-025-00103-7)

The `BaCoN` routine was used in the following manuscripts:

- Krieg, S., Rohde, T., Rausch, T. _et al._ Mitoferrin2 is a synthetic lethal target for chromosome 8p deleted cancers. _Genome Med_ **16**, 83 (2024). DOI: [10.1186/s13073-024-01357-w](https://doi.org/10.1186/s13073-024-01357-w)
