## 7 Visualization

### 7.1 MATH score
According to published researches, using MATH score to quantify tumor heterogeneity has certain clinical significance.
We can calculate MATH score through VAF of samples and present the results. Parameter `tsb` is set to select samples.

`MATH_score(maf_file, tsb = c("tsb1"))`


### 7.2 VAF plot
This function produces density distribution plot, which can show you the clusters of muatations with different Variant Allele Frequencies in all/selected samples. 
Sample(s), and whether to show MATH score can be determined by controlling `sample_option` and `show.MATH` arguments respectively.  

`VAF_plot(maf_file, sample_option = "OFA", theme_option = "aaas", file_format = "png", show.MATH = T)`


### 7.3 Shared/Private mutation (Upset)
Knowing how many matations are shared and owned privately by samples of one patient is helpful for understanding the tumor heterogeneity.
We can use `mut_shard_private` function to visualize the intersect mutations and their types in several samples of one patient, by producing a stack plot.
The parameter "show.num" can be set to determine whether to show the numbers of each muatations in the plot.

`mut.shared_private(maf_file, patientID = c("tsb1"), show.num = FALSE)`


### 7.4 Tumor_clone plot
If the cluster result of a patient and the information of locus are given, this function could visualize the outcomes by offering you a radar plot and a dotplot.
The plot could be changed according to the parameters `clone.min.mut` and `clone.min.aveCCF` you set.

`TumorClones_plot(patientID, ccf.dir = "", out.dir = "", clone.min.mut = 5, clone.min.aveCCF = 0.1)`


### 7.5 Fishplot
With information of cluster and locus, fishplot could provide customers with an intuitive and accurate representation of how an individual tumor is changing over time, which could make analysis easier.
In addition, fish plot may also find inches outside of cancer biology and could represent the changing landscapes of microbial populations.
A clusterEstimates.tsv and a mutation_to_cluster.tsv files would be produced after processing the data input,

(```)
	prepareSchismInput(dir.cluster.tsv, dir.loci.tsv, dir.output)
	write.table(clusterEstimates.tsv)
	write.table(mutation_to_cluster.tsv)
	schism2Fishplot(clusterEstimates.tsv, mutation_to_cluster.tsv)
(```)


### 7.6 GO anaysis
The GO database standardizes the gene products from functions, participating biological pathways and cell localization.
Through GO enrichment analysis, we can roughly understand where the differenal genes enrich, in what biological functions, pathways or cell localizations.
This function can offer a barplot and a dotplot to visualize the result of the GO enrichment. You can get all/seleted type of analysis by setting the `type` parameter.
Also, `pval` and `qval` can be controlled to meet different needs.

`GO_analysis(maf_file , patientID , type = "ALL", pval = 0.01, qval = 0.05)`


### 7.7 Pathway analysis





### 7.8 Mutational Signature
In order to get a phyligenetic tree, this function can define branches' set relationship by re-labeling their tumor sample barcode from the smallest set. 
And it will calcualte each branch's mutational signature weight according to cosmic reference and pick the maxium. 
Finally, a data frame includes each set/branch's mutational signature will be offered.

`Mutational_sigs_tree(maf_file, branch_file)`


### 7.9 NJtree plot
According to researches, constructing a cancer phyligenetic tree with mutation information could help people understand more about the heterogeneity of cancer in samples of certain petient. 
In a phyligenetic tree, the trunk represents shared mutations and the thinner branches represents private and low frequency muatations of different samples. Besides, the distances between branches show how close these samples are.
This function uses the method of Neighbor-joining to construct the phyligenetic tree.
This method minimizes the total distance of the phyligenetic tree by determining the nearest or adjacent paired classification units.
The maf_file would be sorted and transformed into an NJtree object before plotting NJtree. Apart from that, the trunk and branches are colored differently, in order to better show the signatures.
(```)
	maf <- read.Maf(patientID, dat.dir)
	njtree <- read.NJtree(maf_file, use.indel = T, use.ccf = F,mut.signature = TRUE, sig.min.mut.number = 50,)
	getMutSort(njtree)
	getPhyloTree(njtree)
	getNJtreeSignature(njtree)
(```)