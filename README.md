# LinRace
LinRace is a method that integrates the lineage barcode and gene expression data using the asymmetric cell division model and infers cell lineage under a framework combining Neighbor Joining and maximum-likelihood heuristics.

![Figure1 (1)](https://user-images.githubusercontent.com/39555451/216690081-b1e437a3-ca60-4df4-9e43-ea6de0df3614.jpg)

## Installation
The R package for LinRace can be installed as:
`devtools::install_github("https://github.com/ZhangLabGT/LinRace")`

#### Please Cite

```
Xinhai Pan, Hechen Li, Pranav Putta, Xiuwei Zhang (2023). LinRace: single cell lineage reconstruction using paired lineage barcode and gene expression data. BioRxiv, 2023.04.12.536601. doi: https://doi.org/10.1101/2023.04.12.536601
```

## Usage
Example usages of LinRace can be refered in vignettes.
1. [Running LinRace on simulated datasets](vignettes/LinRace_test.Rmd)
2. [Comparing LinRace with other lineage reconstruction methods](vignettes/LinRace_compare.Rmd)
3. [Running LinRace on a mouse embryo dataset(embryo2 from M.Chan et al.)](vignettes/Mouse_test.Rmd)

## Results

The following figure shows an example application of LinRace on the ZF1-F3 sample (750 cells) from scGESTAULT datasets in [(1)](https://doi.org/10.1038/nbt.4103)]:
![LinRace_Figures (3)](https://github.com/ZhangLabGT/LinRace/assets/39555451/789a707d-7fe7-4e22-9c3f-d39d77234541)

a. LinRace reconstructed tree. The outer ring represents major cell type assignments and the inner colors on the edges show detailed intermediate cell type assignments. 
b. Properties of the reconstructed trees from different methods. The max depth means the maximum total edge length going from the root to any leaf cell. 
c. A detailed look of two GES (Gene Expression Subtrees) in LinRace reconstructed tree with inferred ancestral states.
