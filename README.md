# CBioProfiler: a web and standalone pipeine for cancer biomarker and subtype characteration from molecular level to clinical settings

## Instroduction

**CBioProfiler** (**C**ancer **Bio**marker and subtype **Profiler**) was developed to facilitate researchers and clinicians to screen, characterize, annotate and translate cancer biomarkers and subtypes from molecular level to clinical settings more comfortably with graphical user interfaces (GUI), which will help implement targeted clinical diagnosis and treatment measures for different patients to achieve precision medicine.

CBioProfiler integrated a novel R package CuratedCancerPrognosisData that reviewed, curated and integrated the gene expression data and corresponding clinical data of 47,210 clinical samples from 268 gene expression studies of 43 common blood and solid tumors.

The whole pipeline of CBioProfiler includes two main pipelines: cancer biomarker pipeline and cancer subtype pipeline. The cancer biomarker pipeline includes 5 modules: (1) dimensionality reduction using three methods of weighted gene co-expression network analysis(WGCNA), univariate Cox proportional hazards regression model(CoxPH), differentially expressed gene (DEG) analysis, (2) benchmark experiment with 6 machine learning learners (Lasso, Ridge, Elastic net, Glmboost, Coxboost, Randomforest) using cross validation (CV) and nested cross validation (nCV) based on R package [mlr](https://cran.r-project.org/web/packages/mlr/index.html), (3) prediction model construction using Cox proportional hazards regression model and nomogram, (4) clinical annotation using a variety of clinical approaches (correlation with clinical features, Kaplan-Meier curve, CoxPH model, time-dependent ROC, most correlated genes, correlation with specific gene, gene expression in different groups, correlation with immune infiltration, correlation with stemness score, correlation with ESTIMATE score, correlation with immune checkpoint, correlation with IFN-gamma score, correlation with cytolytic activity, correlation with cancer pathway, correlation with metabolism pathway, correlation with hallmark signature, correlation with drug response), and (5) biological annotation using over-representation analysis (ORA) and gene set enrichment analysis (GSEA).
 
The subtype pipeline includes 3 modules: (1) data preprocessing (feature selection based on variance, median absolute deviation(MAD),CoxPH model, and principal component analysis (PCA), (2) subtype identification (integration of multiple unsupervised machine learning methods (K-means clustering (K-means), hierarchical clustering, partitioning around medoids (PAM) clustering, etc.) using two popular consensus clustering methods (ConsensusClusterPlus and M3C), (3) subtype evaluation and validation.

## Licence

Open source under GPLV3.0. Both 'CBioProfiler' software and curated cancer gene expression data 'CuratedCancerPrognosisData' are free for academic but ***non-commercial*** use.

## Important Notes

* Thanks for considering Web CBioProfiler for your study. Due to the limited computing power of the server, when multiple users use CBioProfiler at the same time, the response of the program may become slower, please be patient and only click the button once and wait until one step done.
* If you plan to use CBioProfiler for high-iterative nested cross validation calculations, we ***strongly recommend that you download the CBioProfiler source code and intall it to your R software for corresponding calculations***. This is caused by the limited computing power of the server. We apologize to you for this.
* The App will be disconnected from our server after a hour if there is no mouse action on the web browser.

## Citation

To be added


## Installation

CBioProfiler has reviewed, curated and integrated 45 common blood and solid tumor gene expression data and corresponding clinical data from GEO, TCGA, ICGC, TARGET, ArrayExpress and other public databases. These public data come from 47,210 clinical samples from 268 gene expression studies. If you use these public data when conducting your own research, we strongly recommend that you cite these public data. For a detailed introduction based on these public data, please refer to : https://liuxiaoping2020.github.io/CBioProfilerDatasource/

To promote the use of these public data, we developed an R package "CuratedCancerPrognosisData", which is a necessary condition for CBioProfiler to run public data.

To install CuratedCancerPrognosisData package you can download the associated source code to you local disk at https://zenodo.org/record/5728447#.Ya9vhsj1dk4 and install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) in accordance with your operaction system and R version and then:

**(1) In the R console:**

```{r setup, include=FALSE}
install.packages("path-to-CuratedCancerPrognosisData/CuratedCancerPrognosisData_1.0.tar.gz", repos = NULL, type="source")
```

**or (2) In the terminal line:**

```
R CMD INSTALL path-to-CuratedCancerPrognosisData/CuratedCancerPrognosisData_1.0.tar.gz
```

**After you have downloaded and installed "CuratedCancerPrognosisData", you can install " CBioProfiler" via:**

**（1）github:**

```
devtools::install_github("liuxiaoping2020/CBioProfiler")
```

**or (2) gitee:**

```
remotes::install_git("https://gitee.com/liuxiaoping2020/CBioProfiler")
```

## Running  CBioProfiler locally in R

```
library(CBioProfiler)
CBioProfiler()
```

## Docker image

CBioProfiler has been packaged as a Docker image, users can use CBioProfiler through Docker without any additional configuration. When they finished installing Docker on their personal computer. The Docker image for CBioProfiler can be found at https://hub.docker.com/r/liuxiaoping2020/CBioProfiler/tags

## Contact

If you have any questions regarding the installation and use of the CuratedCancerPrognosisData and CBioProfiler, please feel free to report them at https://github.com/liuxiaoping2020/CBioProfiler/issues, or you can contact the creator and maintainer of the package, Liu Xiaoping through liuxiaoping@whu.edu.cn.
