# BaCoN

Buffering between genes is fundamental for robust cellular functions. While experimentally testing all possible gene pairs is infeasible, gene buffering can be predicted genome-wide under the assumption that a gene’s buffering capacity depends on its expression level and the absence of this buffering capacity primes a severe fitness phenotype of the buffered gene. We developed BaCoN (Balanced Correlation Network), a post-hoc unsupervised correction method that amplifies specific signals in expression-vs-fitness effect correlation-based networks. We quantified 147 million potential buffering relationships by associating CRISPR-Cas9-screening fitness effects with transcriptomic data across 1019 Cancer Dependency Map (DepMap) cell lines. BaCoN outperformed state-of-the-art methods including multiple linear regression, based on our newly compiled metrics for gene buffering predictions. Combining BaCoN with batch correction or Cholesky data whitening further boosts predictive performance. We characterized a high-confidence list of 899 buffering predictions and found that while buffering genes overall are often syntenic, buffering paralogs are on different chromosomes. BaCoN performance increases with more screens and genes considered, making it a valuable tool for gene buffering predictions from the constantly growing DepMap.





This will be the R package containing the functions required to compute a BaCoN matrix. 

# Quickstart


The quickest way to compute a BaCoN matrix is:

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

# Citation

To cite this package, please refer to:
- Rohde, T., Demirtas, Y., Shaw, A., Billmann, M. BaCoN (Balanced Correlation Network) improves prediction of gene buffering. DOI: [10.1101/2024.07.01.601598](https://doi.org/10.1101/2024.07.01.601598)

The BaCoN routine was used in the following manuscripts:

- Krieg, S., Rohde, T., Rausch, T. _et al._ Mitoferrin2 is a synthetic lethal target for chromosome 8p deleted cancers. _Genome Med_ **16**, 83 (2024). DOI: [10.1186/s13073-024-01357-w](https://doi.org/10.1186/s13073-024-01357-w)
