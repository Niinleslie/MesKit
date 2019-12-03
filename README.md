<img src="/vignettes/logo.png" height="80" width="240" />

- [MesKit]
  * [Introduction](#introduction)
  * [Installation](#installation)
  * [Usage](#usage)
  * [Shiny APP](#shiny-app)
  * [Use MesKit with Docker Container](#use-meskit-with-docker-container)
- [Designers](#designers)
- [Credit](#credit)
- [Maintainer](#maintainer)
- [Copyright](#copyright)
- [Citation](#citation)

## Introduction
Intra-tumor heterogeneity (ITH) is now thought to be a key factor contributing to the therapeutic failures and drug resistance, which have attracted increasing attention in the cancer research field. Here, we present an R package, MesKit, for characterizing cancer genomic ITH and inferring the tumor’s evolutionary history. MesKit provides a wide range of analysis including ITH evaluation, enrichment, signature, clone evolution analysis via implementation of well-established computational and statistical methods. 
The source code and documents are freely available through Github (https://github.com/Niinleslie/MesKit). We also developed a shiny application to provide easier analysis and visualization.


## Installation

#### Via Github(latest)

```R
install.packages("remotes")
remotes::install_github("Niinleslie/MesKit")
```

## Usage

For details, please visit https://github.com/Niinleslie/MesKit/blob/master/vignettes/MesKit.Rmd 

<div  align="center">    
<img src="/vignettes/overview.png" height="560" width="700" align = center/>
</div>
   


## Shiny APP

Our package can run with R shiny. To run it:

```R
pkg.suggested <- c('shiny', 'shinyjs','shinyBS','shinydashboard', 'shinyWidgets', 'shinycssloaders', 'DT','org.Hs.eg.db','BSgenome.Hsapiens.UCSC.hg19')
## if genomic reference version is hg38, change 'BSgenome.Hsapiens.UCSC.hg19' to 'BSgenome.Hsapiens.UCSC.hg38'

checkPackages <- function(pkg){
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package pkg needed for shiny app. Please install it.", call. = FALSE)
  }
}
lapply(pkg.suggested, checkPackages)
shiny::runApp(system.file("shiny", package = "MesKit"))
```

Also you can run shiny with:

```R
runMesKit()
```

## Use MesKit with Docker Container

To use MesKit with docker container, please visit https://github.com/Niinleslie/MesKit/blob/master/MesKit.docker.md


## Designers
* Jian Ren, renjian.sysu@gmail.com
* Qi Zhao, zhaoqi3@mail2.sysu.edu.cn
* Mengni Liu, liumn5@mail2.sysu.edu.cn
 
## Credit
This software were developed by:

* [Mengni Liu](liumn5@mail2.sysu.edu.cn), Sun Yat-sen university 
* [Jianyu Chen](chenjy327@mail2.sysu.edu.cn), Sun Yat-sen university 
* [Chengwei Wang](wangchw8@outlook.com), Sun Yat-sen university 

## Maintainer
Mengni Liu

## Copyright
Copyright © 2014-2018. RenLab from SYSUCC. All Rights Reserved
For more useful tools/applications, please go to renbal.org

## Citation
`Citation (from within R, enter citation("MesKit")):`

_MesKit: a tool kit for dissecting cancer evolution from multi-region derived tumor biopsies via somatic mutations_

