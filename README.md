# A log-ratio biplot approach for exploring genetic relatedness based on identity by state

# Dependences:

- R
- PLINK 1.9 (https://www.cog-genomics.org/plink2/)

# logratio_biplot_pca() function

Input parameters

- data.in: the path of your PLINK files 
- data.out: the name of your dataset
- chr: chromosomes used for the analysis (default 1-22)
- geno:  missing call rate threshold (default 0)
- hwe: (default 0.05)
- MAF: minor allele frequency threshold (default 0.4)
- LD.window = 50
- LD.step = 5,
- LD.correlation = 0.2,
- rel.cutoff = 0.025,
- mz.ibd2 = 0.7,
- po.ibd1 = 0.7,
- nsim.rel = 100,
- relationships (default UN, 6TH, 5TH, 4TH, 3RD, 2ND, 3/4S and FS),
- rel.colors 
- peel.and.zoom 

# References


![alt text](https://github.com/ivangalvan/logratio_biplot_pca/blob/master/plots/CEU_population_logratio_biplot_pca_peel_and_zoom.png)

