# Validation of TRM cell signature in hypersensitivity pneumonitis lung tissue
Almost as a bit of a chance discovery during my PhD, I found in some gene expression data (derived from [GSE47460](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47460) and available in [Genevestigator](https://genevestigator.com/)) that lung tissue from patients with hypersensitivity pneumonitis was enriched with genes that suggested that the lung contained far greater densities of T cells, specifically tissue-resident memory T cells, than lung tissue from patients with idiopathic pulmonary fibrosis or control healthy lung tissue. One reason that finding is significant is that it means that if we want to figure out more about the immunology of hypersensitivity pneumonitis we're better off directly sampling the lung than the blood, as tissue-resident T cells are less likely to recirculate in the peripheral blood.

I did need to prove it in a second dataset after having made the initial finding in the GSE47460 subset available in Genevestigator, so I used the gene count matrix available under [GSE150910](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150910). This repository contains the script for that validation analysis and the production of figure 1B in the letter that we submitted to AJRCCM ([Sendama et al., *Am J Respir Crit Care Med* 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10161757/)).
