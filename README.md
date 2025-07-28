# SubsRates

A pipeline to extract **terminal branch** dN, dS, and ω (omega) values from PAML (codeml) output, compute **root-to-tip distances** from phylogenetic trees, and summarize the results in a table and annotated tree plots.

It consists of:
- A **shell script** that parses codeml output.
- An **R script** that calculates root-to-tip distances and generates visualizations.


---

## Requirements

- **R** (≥ 4.0)
  - R packages:
    ```r
    install.packages(c("tidyverse", "ape"))
    remotes::install_github("ds8k/castor")
    ```

- **PAML** (v4.8 or compatible)

---

# Reference
This study is reported in the following paper:
Ceriotti, L. F., Gatica-Soria, L., & Sanchez-Puerta, M. V. (2022). Cytonuclear coevolution in a holoparasitic plant with highly disparate organellar genomes. Plant Molecular Biology, 109(6), 673-688. https://doi.org/10.1007/s11103-022-01266-9
