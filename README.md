# A log-ratio biplot approach for exploring genetic relatedness based on identity by state

Relatedness:

- UN = unrelated
- 6TH = sixth degree relationships (second cousins once removed)
- 5TH = fifth degree relationships (second cousins)
- 4TH = fourth degree relationships (first cousins once removed)
- 3RD = third degree relationships (first cousins)
- 2ND = second degree relationships (half siblings, avuncular, grandparent-grandchild)
- 3/4S = three quarter siblings
- FS = full siblings
- PO = parent-offspring
- MZ = monozygotic twins

# Dependencies

- R (https://www.r-project.org/)
- PLINK 1.9 (https://www.cog-genomics.org/plink2/)

# logratio_biplot_pca() function

Input parameters:

| Parameter  | Description |
| ------------- | ------------- |
| data.in  | path of your PLINK files   |
| data.out  |  name of your dataset (no spaces allowed)  |
| chr  | chromosomes used for the analysis (default 1-22)  |
| geno  | missing call rate threshold (default 0) |
| hwe  | Hardy Weingberg Equilirbium mid p-value threshold (default 0.05)  |
| MAF  | minor allele frequency threshold (default 0.40)  |
| LD.window  | LD window size in kb (default 50) |
| LD.step  | number of variants to shift the window at each step (default 5)  |
| LD.correlation  | pairwise LD correlation (r<sup>2</sup>) threshold (default 0.2)  |
| rel.cutoff  | kinship coefficient threshold to detect potential unrelated individuals (default 0.025)  |
| mz.ibd2  | IBD2 threshold to detect potential MZ pairs (default 0.7)  |
| po.ibd1  | IBD1 threshold to detect potential PO pairs (default 0.7)  |
| relationships  | predicted relationships in the analysis (default UN, 6TH, 5TH, 4TH, 3RD, 2ND, 3/4S and FS)  |
| nsim.rel | number of simulations for each relationship used as training set (default 20)  |
| rel.colors  | colors for each relationship (default forestgreen, darkorange, deeppink2, cyan2, darkgoldenrod2, darkorchid2, black and dodgerblue2)  |
| peel.and.zoom  | peel and zoom log-ratio biplot approach (default TRUE) |

# Example of the CEU population from 1000Genomes project

```
setwd("github")
source("logratio_biplot_pca.R")

# logratio biplot with UN, 6TH, 5TH, 4TH, 3RD, 2ND, 3/4S and FS relationships and peel and zoom approach

my_logpca = logratio_biplot_pca(data.in="data/CEU",
                                data.out = "CEU_population",
                                nsim.rel = 20,
                                peel.and.zoom = T)

table(my_logpca$pairs.ids$prediction)
```


# References

- Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ (2015) Second-generation PLINK: rising to the challenge of larger and richer datasets. GigaScience, 4.
- Graffelman J, Moreno V (2013). The mid p-value in exact tests for Hardy-Weinberg equilibrium. Statistical applications in genetics and molecular biology, 12(4), 433-448.
- Graffelman J, Galván-Femenía I, de Cid R, Barceló-i-Vidal C. A log-ratio biplot approach for exploring genetic relatedness based on identity by state. 

# Peel and zoom approach

![alt text](https://github.com/ivangalvan/logratio_biplot_pca/blob/master/plots/CEU_population_logratio_biplot_pca_peel_and_zoom.png)

