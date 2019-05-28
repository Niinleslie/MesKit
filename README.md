# Meskit

## Introduction

## Documentation
1. Installation
2. Quick Start
3. Parameters
4. Output Interpretations
5. FAQ

### Quick Start 
#### Prepare input files
To analyse with our packages you need to provide a set of input files, including:
- Necessary details of samples;
- Mutation information (MAF files);
- Output files of Pyclone analysis (ccf files).

Here we give examples of details of samples and MAF files. And others you should prepare contain the output tables of Pyclone analysis (```clusters.tsv``` and ```loci.tsv``` for cluster and locus specific information).
 
##### Information of samples
> It should contain the sampleID, patientID, lesion and time.
 
 |  sample  |  patient |  lesion |  time  |
 ---- | ------ | ------ | ------
 | 311252-S | 311252 |  S      |   -   |
 | 311252-V |  311252  |  V  |     -   |
 |311252-TC1 | 311252 |  TC  |     -   |
 |311252-TC2 | 311252 |  TC  |     -   |
 
 
##### The MAF files
| Hugo_Symbol|	Chromosome | Start_Position |	End_Position |	Variant_Classification | Variant_Type |	Reference_Allele |	Tumor_Seq_Allele1 | Tumor_Seq_Allele2 |	Ref_allele_depth |	Alt_allele_depth |	VAF	| CDS_Change	| Protein_Change |	Tumor_Sample_Barcode |
|:-----| :------| :------ | :----- | :------ | :----- | :---- | :-----| :----- | :----- | :-------| :---- | :-----| :----- | :----- |
| LOC729737| 1 | 135207 |	135207	| RNA |	SNP |	C | C |	G |	40	| 4 | 0.0909 | NA |	NA | 311252-S |
|TTC34,ACTRT2| 1 | 2869474 | 2869474 |	IGR |INS | - | | CTCTCT |	43 | 8 | 0.1568 | NA |	NA | 311252-S |
|NBPF1|1 | 16908223 |	16908223 | Intron | SNP | T |	T |A|	142| 8 | 0.0533 | NA| NA | 311252-S|
|PRAMEF2 | 1 | 12921600 | 12921600 | Missense_Mutation | SNP | C |	C | T |73 |	3 |	0.0394 | c.C1391T |	p.P464L | 311252-S |
