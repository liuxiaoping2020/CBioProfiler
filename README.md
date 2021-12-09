# CBioExplorer: a web and standalone application for screening, validation, and annotation of cancer survival related biomarkers from molecular level to clinical settings

## Instroduction

**CBioExplorer** (**C**ancer **Bio**marker **Explorer**) was developed to facilitate researchers and clinicians to screen,characterize, annotate and translate cancer biomarkers from molecular level to clinical settings more comfortably with graphical user interfaces (GUI). The whole pipeline of CBioExplorer includes data collection, data curation, dimensionality reduction using three methods of WGCNA, univariate Cox proportional hazards regression model, differentially expressed gene analysis, benchmark experiment with 6 machine learning learners (Lasso, Ridge, Elastic net, Glmboost, Coxboost, Randomforest) using cross validation and nested cross validation based on R package [mlr](https://cran.r-project.org/web/packages/mlr/index.html), prediction model construction using Cox proportional hazards regression model and nomogram, clinical annotation using a variety of clinical approaches, and biological annotation using over-representation analysis (ORA) and gene set enrichment analysis (GSEA).

## Licence

Open source under GPLV3.0. Both 'CBioExplorer' software and curated cancer gene expression data 'CuratedCancerPrognosisData' are free for academic but ***non-commercial*** use.

## Important Notes

* Thanks for considering Web CBioExplorer for your study. Due to the limited computing power of the server, when multiple users use CBioExplorer at the same time, the response of the program may become slower, please be patient and only click the button once and wait until one step done.
* If you plan to use CBioExplorer for high-iterative nested cross validation calculations, we ***strongly recommend that you download the CBioExplorer source code and intall it to your R software for corresponding calculations***. This is caused by the limited computing power of the server. We apologize to you for this.
* The App will be disconnected from our server after a hour if there is no mouse action on the web browser.

## Citation

To be added

## Installation

CBioExplorer has reviewed, curated and integrated 45 common blood and solid tumor gene expression data and corresponding clinical data from GEO, TCGA, ICGC, TARGET, ArrayExpress and other public databases. These public data come from 47,210 clinical samples from 268 gene expression studies. If you use these public data when conducting your own research, we strongly recommend that you cite these public data. For a detailed introduction based on these public data, please refer to : https://liuxiaoping2020.github.io/CBioExplorerDatasource/

To promote the use of these public data, we developed an R package "CuratedCancerPrognosisData", which is a necessary condition for CBioExplorer to run public data.

To install CuratedCancerPrognosisData package you can download the associated source code to you local disk at https://zenodo.org/record/5728447#.Ya9vhsj1dk4 and install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) in accordance with your operaction system and R version and then:

**(1) In the R console:**

```{r setup, include=FALSE}
install.packages("path-to-CuratedCancerPrognosisData/CuratedCancerPrognosisData_1.0.tar.gz", repos = NULL, type="source")
```

**or (2) In the terminal line:**

```
R CMD INSTALL path-to-CuratedCancerPrognosisData/CuratedCancerPrognosisData_1.0.tar.gz
```

**After you have downloaded and installed "CuratedCancerPrognosisData", you can install "CBioExplorer" via:**

**（1）github:**

```
devtools::install_github("liuxiaoping2020/CBioExplorer")
```

**or (2) gitee:**

```
remotes::install_git("https://gitee.com/liuxiaoping2020/CBioExplorer")
```

## Running CBioExplorer locally

```
library(CBioExplorer)
CBioExplorer()
```

## Docker image

CBioExplorer has been packaged as a Docker image, users can use CBioExplorer through Docker without any additional configuration. When they finished installing Docker on their personal computer, they can use the following code in the terminal to download and run CBioExplorer Docker image:

```
docker run -d -p 3838:3838 liuxiaoping2020/cbioexplorer
```

## Contact

If you have any questions regarding the installation and use of the CuratedCancerPrognosisData and CBioExplorer, please feel free to report them at https://github.com/liuxiaoping2020/CBioExplorer/issues, or you can contact the creator and maintainer of the package, Liu Xiaoping through liuxiaoping@whu.edu.cn.

