#' Prepare inputs for SCHISM
#' @description read tsv documents from PyClone and figure 
#' out the list of filtered clusters. According to those 
#' targeted clusters, filter targeted mutations and reorganize 
#' data as two outputs: clusterEstimates.tsv and 
#' mutation_to_cluster.tsv, which are the input of SCHISM.
#' 
#' @param dirClusterTsv specify the directory of cluster.tsv document
#' @param dirLociTsv specify the directory of loci.tsv document
#' @param dirOutput specify the directory of two output for SCHISM
#' @return clusterEstimates.tsv, mutation_to_cluster.tsv
#' 
#' @export prepareSchismInput
#' 
#' @examples
#' \dontrun{
#' prepareSchismInput(dir.cluster.tsv, dir.loci.tsv, dir.output)
#'}

## SCHISM method
#' Get outputs from SCHISM and draw the fishplot
#' @description read cluster.cellularity and GA.consensusTree documents 
#' from SCHISM. Construct fraction table and parents vector for drawing 
#' fishplot.
#' 
#' @param dir.cluster.cellularity specify the directory of cluster.cellularity 
#' document
#' @param dir.GA.consensusTree specify the directory of GA.consensusTree 
#' document
#' @return the fishplot showing the evolution relationship predicted by SCHISM
#' 
#' @export schism2Fishplot
#' @export schism2Timescape
#' @export inferByClonevol
#'
#' @examples
#' \dontrun{
#' schism2Fishplot(dir.cluster.cellularity, dir.GA.consensusTree)
#'}

## generate results for createFishPlotObjects
schism2Fishplot <- function(dirClusterCellularity, dirGAconsensusTree){
    ## get cellularity infomation
    clusterCellularity=read.table(dirClusterCellularity, sep="\t", 
                                   stringsAsFactors=FALSE, header=TRUE)
    
    ## get the list of samples and clusters in this function
    sampleLs=unique(clusterCellularity$sampleID)
    clusterLs=unique(clusterCellularity$clusterID)
    
    ## build the cancer celluarity/fraction table for fishplot
    frac.c=c()
    for (sampleName in sampleLs) {
        sample <- clusterCellularity[which(
            clusterCellularity$sampleID == sampleName), ]$cellularity
        frac.c <- c(frac.c,sample)
    }
    fracTable=matrix(frac.c*100, ncol=length(sampleLs))
    rownames(fracTable) <- clusterLs
    colnames(fracTable) <- sampleLs
    
    ## read the evolution relationship of different subclones
    GAconsensusTree=read.table(dirGAconsensusTree, sep="\t", 
                                stringsAsFactors=FALSE, header=TRUE)
    
    ## figure out clusters which will be first nodes of the evolution tree
    endsLs=unique(GAconsensusTree[which(
        !GAconsensusTree$parent %in% GAconsensusTree$child), 1])
    
    ## make sure the first node of the tree would be first in parents
    if (!rownames(fracTable)[1] %in% endsLs){
        fracEnd=fracTable[which(
            rownames(fracTable) == endsLs[1]),]
        fracTable <- fracTable[-which(
            rownames(fracTable) == endsLs[1]),]
        fracTable <- rbind2(fracEnd, fracTable)
        rownames(fracTable)[1] <- endsLs[1]
    }
    
    ## rearrange subclonal evolution relationship as a vector(parents)
    parents <- c(seq_along(clusterLs))
    for (cluster_name in clusterLs){
        parentCluster <- GAconsensusTree[which(
            GAconsensusTree$child == cluster_name), 1]
        if (cluster_name %in% endsLs){
            parents[which(rownames(fracTable) == cluster_name)] <- 0
        } else{
            parents[which(
                rownames(fracTable) == cluster_name)] <- which(
                    rownames(fracTable) == min(parentCluster))
        }
    }
    
    ## fishplot printing
    pdf('fish.pdf', width=8, height=5)
    fish=createFishObject(fracTable, parents, 
                          timepoints=c(seq_along(sampleLs)), 
                          fix.missing.clones=TRUE)
    if (length(clusterLs) > 10) {
        fish=setCol(fish, as.character(clusterLs))
    }
    fish=layoutClones(fish)
    fishPlot(fish, shape="spline", title.btm="PatientID", cex.title=0.5,
             vlines=seq(1, length(sampleLs)), vlab=sampleLs, pad.left=0.5)
    dev <- dev.off()
    message("Fishplot Done!")
}

## Timescape method
schism2Timescape <- function(dirClusterCellularity, 
                             dirGAconsensusTree, 
                             dirSampleInfo=NULL){
    ## get cellularity infomation
    clusterCellularity=read.table(dirClusterCellularity, 
                                   sep="\t", 
                                   stringsAsFactors=FALSE, 
                                   header=TRUE)
    names(clusterCellularity) <- c("timepoint", "clone_id", 
                                    "clonal_prev", "sd")
    clonalPrev <- clusterCellularity[c("timepoint", "clone_id", 
                                         "clonal_prev")]
    
    ## make sure the timepoint is specific for the real data.
    if (!is.null(dirSampleInfo)){
        ## read info file
        sampleInfoInput <- read.table(dirSampleInfo, quot="", 
                                        header=TRUE, fill=TRUE, 
                                        sep='', stringsAsFactors=FALSE)
        sampleInfoInput <- sampleInfoInput[order(sampleInfoInput$time),]
        timepoints <- sampleInfoInput$sample
        
        clonalPrevTemp <- data.frame()
        for (timepoint in timepoints){
            clonalPrevTimepoint <- clonalPrev[which(
                clonalPrev$timepoint == timepoint), ]
            clonalPrevTemp <- rbind(clonalPrevTemp, clonalPrevTimepoint)
        }
        clonalPrev <- clonalPrevTemp
    } 
    
    ## read the evolution relationship of different subclones
    GAconsensusTree=read.table(dirGAconsensusTree, sep="\t", 
                                stringsAsFactors=FALSE, header=TRUE)
    names(GAconsensusTree) <- c("source", "target")
    treeEdges <- GAconsensusTree
    
    timescape(clonalPrev, treeEdges, mutations="NA", clone_colours="NA",
              xaxis_title="Time Point", yaxis_title="Clonal Prevalence",
              phylogeny_title="Clonal Phylogeny", alpha=50,
              genotype_position="stack", perturbations="NA", sort=FALSE,
              show_warnings=TRUE, width=900, height=NULL)
    message("Timescape Plot Done!")
}

## clonevol method
inferByClonevol <- function(dir.cluster.tsv, plotOption="fishplot"){
    ## read mutations in cluster.tsv from PyClone and get targeted clusters
    cluster.tsv=read.table(dir.cluster.tsv, sep="\t",
                           stringsAsFactors=FALSE, header=TRUE)
    cluster_filter=cluster.tsv[!is.na(cluster.tsv$cluster_id) &
                                   cluster.tsv$size >= 5 &
                                   cluster.tsv$mean >= 0.1,]
    cluster_ls =unique(cluster_filter$cluster_id)
    
    ## read mutations within targeted clusters in loc.tsv file from PyClone
    loci.tsv=read.table(dir.loci.tsv, sep="\t",
                        stringsAsFactors=FALSE, header=TRUE)
    loci_filtered <- loci.tsv[which(loci.tsv$cluster_id %in% cluster_ls),]
    
    # loci_filtered <- na.omit(loci_filtered)
    sample.names <- unique(loci_filtered$sample_id)
    
    
    ## make sure the cluster id is continuous integer.
    for (cluster_newid in seq_along(ls.cluster_id)){
        loci_filtered$cluster_id[loci_filtered$cluster_id == ls.cluster_id[cluster_newid]] <- cluster_newid
    }
    
    ls.cluster_id <- unique(loci_filtered$cluster_id)
    
    ## separate each sample column
    dat.final <- data.frame()
    for (cluster in ls.cluster_id){
        dat.cluster <- loci_filtered[which(loci_filtered$cluster_id == cluster),]
        dat.samples <- data.frame(row.names=seq_along(dat.cluster[which(dat.cluster$sample_id == sample.names[1]),]$cellular_prevalence))
        for (sample in sample.names){
            col.sample <- data.frame(dat.cluster[which(dat.cluster$sample_id == sample),]$cellular_prevalence)
            names(col.sample) <- sample
            dat.samples <- cbind(dat.samples, col.sample)
        }
        dat.samples <- cbind(cluster, dat.samples)
        dat.final <- rbind(dat.final, dat.samples)
    }
    
    ##################### must be used ##################################
    consensusTree=infer.clonal.models(variants=dat.final,
                                      cluster.col.name='cluster',
                                      ccf.col.names=sample.names,
                                      cancer.initiation.model='monoclonal',
                                      subclonal.test='bootstrap',
                                      subclonal.test.model='non-parametric',
                                      founding.cluster=1,
                                      num.boots=1000,
                                      cluster.center='mean',
                                      min.cluster.vaf=0.01,
                                      # min probability of CCF is non-negative
                                      sum.p=0.05,
                                      # confidence interval estimate for CCF
                                      alpha=0.05)
    
    fishPlotInput <- generateFishplotInputs(results=consensusTree)
    
    if (plotOption == "fishplot"){
        fishes = createFishPlotObjects(fishPlotInput)
        for (i in seq_along(fishes)){
          pdf(paste('FISH', i, '.pdf', sep=""), width=8, height=5)
          fish = layoutClones(fishes[[i]])
          fish = setCol(fish,fishPlotInput$clonevol.clone.colors)
          if (length(sample.names) > 10) {
            fish=setCol(fish, as.character(sample.names))
          }
          fishPlot(fish,shape="spline", 
                   title.btm="Patient", 
                   cex.title=0.5,
                   vlines=seq(1, length(sample.names)), 
                   vlab=sample.names, pad.left=0.5)
          dev.off()
        }
        message("FishPlot Done!")
    } else if (plotOption == "timescape"){
        ## preparation for timscape input: treeEdges
        fishPlotTreeEdges <- fishPlotInput$parents[[1]]
        fishPlotTreeEdges <- 
            c(fishPlotTreeEdges[which(fishPlotTreeEdges == 0)], 
              fishPlotTreeEdges[which(fishPlotTreeEdges > 0)]-1)
        branchEdge <- data.frame(matrix(ncol = 2, nrow = 0))
        treeEdges <- data.frame(matrix(ncol = 2, nrow = 0))
        cloneList <- 0:(length(fishPlotTreeEdges)-1)
        for (counter in seq_along(fishPlotTreeEdges)){
            if (cloneList[counter] == 0){
                next()
            }
            branchEdge <- data.frame(fishPlotTreeEdges[counter], 
                                     cloneList[counter])
            treeEdges <- rbind(treeEdges, branchEdge)
        }
        colnames(treeEdges) <- c("source", "target")
        
        ## preparation for timscape input: clonalPrev
        fishPlotClonalPrev <- as.data.frame(fishPlotInput$cell.fractions[[1]])
        clonalPrev <- data.frame(matrix(ncol = 3, nrow = 0))
        clonalPrevTemp <- data.frame(matrix(ncol = 3, nrow = 0))
        for (timepoint in colnames(fishPlotClonalPrev)){
            for (cloneIdCounter in seq_along(cloneList)){
                clonalPrevTemp <- 
                    data.frame(timepoint, cloneList[cloneIdCounter], 
                               fishPlotClonalPrev[, which(colnames(fishPlotClonalPrev) == timepoint)][cloneIdCounter])
                clonalPrev <- rbind(clonalPrev, clonalPrevTemp)
            }
        }
        colnames(clonalPrev) <- c("timepoint", "clone_id", "clonal_prev")
        ## print 
        timescape(clonalPrev, treeEdges, mutations="NA", clone_colours="NA",
                  xaxis_title="Time Points", yaxis_title="Clonal Prevalence",
                  phylogeny_title="Clonal Phylogeny", alpha=50,
                  genotype_position="stack", perturbations="NA", sort=FALSE,
                  show_warnings=TRUE, width=900, height=NULL)
        message("Timescape Plot Done!")
    }
  }



