<img src="/vignettes/logo.png" height="80" width="240" /> 

# [M]()ulti-region [e]()xome [s]()equencing analysis tool [Kit]()

Intra-tumor heterogeneity (ITH) is now thought to be a key factor contributing to the therapeutic failures and drug resistance, which have attracted increasing attention in the cancer research field. Here, we present an R package, MesKit, for characterizing cancer genomic ITH and inferring the history of tumor evolution. MesKit provides a wide range of analysis including ITH evaluation, enrichment, signature, clone evolution analysis via implementation of well-established computational and statistical methods. 
The source code and documents are freely available through Github (https://github.com/Niinleslie/MesKit). We also developed a shiny application to provide easier analysis and visualization.


## Installation

```R
install.packages("remotes")
remotes::install_github("Niinleslie/MesKit")
```

## Usage
A complete documentation of MesKit can be found [here.](http://meskit.renlab.org/)

<div  align="left">    
<img src="/vignettes/MesKit_overview.png" height="500" width="600" />
</div>
   


## Shiny APP

For GUI-based analysis, users can use the following code to launch Shiny app build with the package.

```R
pkg.suggested <- c('shiny', 'shinyjs','shinyBS','shinydashboard', 'shinyWidgets', 'shinycssloaders', 'DT','org.Hs.eg.db','BSgenome.Hsapiens.UCSC.hg19')
## if genomic reference version is hg38, change 'BSgenome.Hsapiens.UCSC.hg19' to 'BSgenome.Hsapiens.UCSC.hg38'

# Install the required packages
checkPackages <- function(pkg){
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package pkg needed for shiny app. Please install it.", call. = FALSE)
  }
}
lapply(pkg.suggested, checkPackages)
# run shiny app from shiny package
shiny::runApp(system.file("shiny", package = "MesKit"))
```

Also you can run shiny with:

```R
runMesKit()
```

## Configure Shiny APP with Docker 

We provided a docker image for a quick configuration of shiny app bundle with shiny-server, please see the simple commands [here](https://github.com/Niinleslie/MesKit/blob/master/MesKit.docker.md).

## Authors
This software was mainly developed by:

* Mengni Liu, liumn5@mail2.sysu.edu.cn, Sun Yat-sen university 
* Chengwei Wang, wangchw8@outlook.com, Sun Yat-sen university 
* Jianyu Chen, chenjy327@mail2.sysu.edu.cn, Sun Yat-sen university 
* Xin Wang, wangx555@mail2.sysu.edu.cn, Sun Yat-sen university

## Supervised by 

* [Jian Ren](renjian@sysucc.org.cn) and [Qi Zhao](zhaoqi@sysucc.org.cn) from Bioinformatic Center of Sun Yat-sen University Cancer Center 

## Maintainer
[Mengni Liu](liumn5@mail2.sysu.edu.cn), Sun Yat-sen university  <br/>

## Copyright

Copyright © 2014-2019. RenLab from SYSUCC. All Rights Reserved<br/>
For more useful tools/applications, please go to [renlab.org](http://www.renlab.org)

## Citation
Citation (from within R, enter `citation("MesKit")):`
>
_MesKit: a tool kit for dissecting cancer evolution from multi-region derived tumor biopsies via somatic mutations_

