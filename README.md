
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `{PhosView}`

<!-- badges: start -->

<!-- badges: end -->

## Disclaimer

PhosView has been developed and tested internally on small datasets
at the Proteomic Platform of the SFR Necker. As it is still in early development, updates will
be published as issues are identified and resolved.

Bug reports, feature requests, and feedback are welcome:
**[Open an issue](https://github.com/OMICShub-Imagine/PhosView/issues)**

*PhosView is a free R package and comes with ABSOLUTELY NO WARRANTY.*

## Quick Start

You can install the development version of `{PhosView}` like so:

### Installation

#### 1. Prerequisites

PhosView requires **R >= 4.1.0** (https://cran.r-project.org/bin/windows/base/) and RStudio (but is not mandatory)

#### 2. Install

```r
if (!require("remotes")) install.packages("remotes")

remotes::install_github("OMICShub-Imagine/PhosView")
```

### Run

You can launch the application by running:

``` r
library(PhosView)
PhosView::run_app()
```

## Tutorials

[How to use PhosView?](docs/index.html)


```r
library(PhosView)
```

## Citation

If you use PhosView in your research, please cite:

> Roger K, Ceccacci S, Guerrera C, Puig Lombardi ME (2026).
> *PhosView : An Integrated and Interactive Tool for Quantitative Proteomics and Phosphoproteomics Analysis.*
> R package version 1.0.0. Institut Imagine, Paris, France.

```r
citation("PhosView")
```

### Please also cite the underlying tools

PhosView wraps several published methods. 

| Tool | Citation |
|------|----------|
| **proteoDA** | Thurman TJ, Washam CL, Alkam D, Bird JT, Gies A, Dhusia K, et al. **proteoDA: a package for quantitative proteomics.** Journal of Open Source Software. 30 mai 2023;8(85):5184. [doi:10.21105/joss.05184](https://doi:10.21105/joss.05184) |
| **clusterProfiler** | Wu T, Hu E, Xu S, Chen M, Guo P, Dai Z, et al. **clusterProfiler 4.0: A universal enrichment tool for interpreting omics data.** The Innovation. août 2021;2(3):100141. [doi:10.1016/j.xinn.2021.100141](https://doi:10.1016/j.xinn.2021.100141) |
| **PhosR** | Kim HJ, Kim T, Xiao D, Yang P. **Protocol for the processing and downstream analysis of phosphoproteomic data with PhosR.** STAR Protocols. juin 2021;2(2):100585. [doi:10.1016/j.xpro.2021.100585](https://doi:10.1016/j.xpro.2021.100585) |
| **Pathview** | Luo W, Brouwer C. **Pathview: an R/Bioconductor package for pathway-based data integration and visualization.** Bioinformatics. 15 juill 2013;29(14):1830‑1. [doi:10.1093/bioinformatics/btt285](https://doi:10.1093/bioinformatics/btt285) |


## Authors

- **Kévin Roger** — Author, maintainer | Proteomics platform Necker, Université Paris Cité - Structure Fédérative de Recherche Necker, INSERM US24/CNRS UAR3633, Paris, France | Bioinformatics Core, SFR Necker – Institut Imagine, INSERM U1163, Paris, France 
- **Sara Ceccacci** — Contributor | Proteomics platform Necker, Université Paris Cité - Structure Fédérative de Recherche Necker, INSERM US24/CNRS UAR3633, Paris, France
- **Chiara Ida Guerrera** — Supervisor | Proteomics platform Necker, Université Paris Cité - Structure Fédérative de Recherche Necker, INSERM US24/CNRS UAR3633, Paris, France
- **Maria Emilia Puig Lombardi** — Supervisor | Bioinformatics Core, SFR Necker – Institut Imagine, INSERM U1163, Paris, France

## License

GPL (>= 3)

