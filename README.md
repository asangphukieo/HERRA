# HERRA
R script for Heritability Estimation using a Regularized Regression Approach (HERRA)
doi: 10.1371/journal.pone.0181269

- This is modified version of HERRA script from original R code (http://www.tau.ac.il/~gorfinem/).
  - add parallelization
  - add bootstrapping support for h^2
  - add random phenotype mode
  - run for each chromosome separately
  
# Requires:
  * Preprocess vcf file
    * vcftools
  * HERRA R script
    * glmnet: installed.packages("glmnet")
    * doMC: installed.packages("doMC")
    * boot: installed.packages("boot")
    * permute: installed.packages("permute")
  * Dataset
    * SNPs from Maela dataset of Streptococus pneumoniae genomes (snps.vcf.gz) (10.1371/journal.pgen.1004547)
    * phenotype data of drug resistance (resistances.phe)

# To use
  ### 1.  Preprocess vcf file
  - Filter out SNP containing missing value and minor allele frequency < 0.01 :

  ```
  vcftools --extract-FORMAT-info GT --out vcf_matrix_Maela --maf 0.01 --max-missing 1.0 --gzvcf $vcf
  ```
  
  - Replace genotype of '0/0' to 0, and '1/1' to 1 
  ```
  cut -f3- vcf_matrix_Maela.GT.FORMAT| sed s/"0\/0"/"0"/g | sed s/"1\/1"/"1"/g > vcf_matrix_Maela.txt
  ```
 
  ### 2. Run HERRA R script
  
  - 2.1 Run without bootstraping
```
R --slave -f HERRA.R --args vcf_matrix=vcf_matrix_Maela.txt phe_file=resistances.pheno use_cov=FALSE out=herra_allSNP_Maela boot=FALSE prevalence=0.1 random_phe=FALSE lambda_list=0.006, 0.008, 0.01, 0.0102, 0.0104
```

  - 2.2 Run with bootstraping
```
R --slave -f HERRA.R --args vcf_matrix=vcf_matrix_Maela.txt phe_file=resistances.pheno use_cov=FALSE core=30 out=herra_allSNP_Maela boot=TRUE n_boot=100 prevalence=0.1 random_phe=FALSE lambda_list=0.01
```

Argument|description
---|---
vcf_matrix| input matrix file (require)
phe_file| phenotype file (require)
use_cov| apply cavariate (default: TRUE)
cov_file| cavariate file (require if `use_cov=TRUE`)
out|output file name
random_phe| radom shuffling phenotype (default: FALSE)
core| cpu core for parallelization applicable for bootstrapping
boot| calculate bootstrapping support (default: TRUE)
n_boot| number of bootstrapping run (require if `boot=TRUE`, default: 100)
prevalence| disease prevalence for h^2 adjusted score (liability score) (default: 0.1)
lambda_list| shrinkage parameter value in variant screening step (default: "0.006, 0.008, 0.01, 0.0102, 0.0104"). If multiple shrinkage parameter values were used, the only maximum h^2 will be selected and reported.

  ### 3. output
  There are 2 files reported after running
1. log file containing h^2 score: herra_allSNP_Maela.out
```
grep -v "#" herra_allSNP_Maela.out
```

2. Histogram and quantile plot of bootstrapping result: herra_allSNP_Maela.pdf
![alt text](https://)

