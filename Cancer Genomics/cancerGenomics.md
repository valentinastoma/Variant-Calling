Maftools only requires somatic variants in Mutation Annotation Format
(MAF) and is independent of larger alignment files.

    library(devtools)

    ## Loading required package: usethis

    # install_github(repo = "PoisonAlien/maftools")
    library(maftools)

    # BSgenome.Hsapiens.UCSC.hg19 package
    # if (!requireNamespace("BiocManager", quietly = TRUE))
    #     install.packages("BiocManager")
    # BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
    # 
    # # NMF and pheatmap packages
    # install.packages(c("pheatmap", "NMF"))

Data
----

    #path to TCGA LAML MAF file
    laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') 
    #clinical information containing survival information and histology. This is optional
    laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') 

    laml = read.maf(maf = laml.maf, clinicalData = laml.clin)

    ## -Reading
    ## -Validating
    ## -Silent variants: 475 
    ## -Summarizing
    ## -Processing clinical data
    ## -Finished in 0.309s elapsed (0.820s cpu)

Accessing the information of the useful slots from MAF object (<span
class="citeproc-not-found" data-reference-id="Mayakonda">**???**</span>)

    laml

    ## An object of class  MAF 
    ##                    ID          summary  Mean Median
    ##  1:        NCBI_Build               37    NA     NA
    ##  2:            Center genome.wustl.edu    NA     NA
    ##  3:           Samples              193    NA     NA
    ##  4:            nGenes             1241    NA     NA
    ##  5:   Frame_Shift_Del               52 0.271      0
    ##  6:   Frame_Shift_Ins               91 0.474      0
    ##  7:      In_Frame_Del               10 0.052      0
    ##  8:      In_Frame_Ins               42 0.219      0
    ##  9: Missense_Mutation             1342 6.990      7
    ## 10: Nonsense_Mutation              103 0.536      0
    ## 11:       Splice_Site               92 0.479      0
    ## 12:             total             1732 9.021      9

    #Shows sample summry.
    head(getSampleSummary(laml))

    ##    Tumor_Sample_Barcode Frame_Shift_Del Frame_Shift_Ins In_Frame_Del
    ## 1:         TCGA-AB-3009               0               5            0
    ## 2:         TCGA-AB-2807               1               0            1
    ## 3:         TCGA-AB-2959               0               0            0
    ## 4:         TCGA-AB-3002               0               0            0
    ## 5:         TCGA-AB-2849               0               1            0
    ## 6:         TCGA-AB-2923               1               1            0
    ##    In_Frame_Ins Missense_Mutation Nonsense_Mutation Splice_Site total
    ## 1:            1                25                 2           1    34
    ## 2:            0                16                 3           4    25
    ## 3:            0                22                 0           1    23
    ## 4:            0                15                 1           5    21
    ## 5:            0                16                 1           2    20
    ## 6:            0                15                 3           0    20

    #Shows gene summary.
    head(getGeneSummary(laml))

    ##    Hugo_Symbol Frame_Shift_Del Frame_Shift_Ins In_Frame_Del In_Frame_Ins
    ## 1:        FLT3               0               0            1           33
    ## 2:      DNMT3A               4               0            0            0
    ## 3:        NPM1               0              33            0            0
    ## 4:        IDH2               0               0            0            0
    ## 5:        IDH1               0               0            0            0
    ## 6:        TET2              10               4            0            0
    ##    Missense_Mutation Nonsense_Mutation Splice_Site total MutatedSamples
    ## 1:                15                 0           3    52             52
    ## 2:                39                 5           6    54             48
    ## 3:                 1                 0           0    34             33
    ## 4:                20                 0           0    20             20
    ## 5:                18                 0           0    18             18
    ## 6:                 4                 8           1    27             17
    ##    AlteredSamples
    ## 1:             52
    ## 2:             48
    ## 3:             33
    ## 4:             20
    ## 5:             18
    ## 6:             17

Using ‘head’ for demonstation purposes in the kntited file and
conservation of space.

    #shows clinical data associated with samples
    head(getClinicalData(laml))

    ##    Tumor_Sample_Barcode FAB_classification days_to_last_followup
    ## 1:         TCGA-AB-2802                 M4                   365
    ## 2:         TCGA-AB-2803                 M3                   792
    ## 3:         TCGA-AB-2804                 M3                  2557
    ## 4:         TCGA-AB-2805                 M0                   577
    ## 5:         TCGA-AB-2806                 M1                   945
    ## 6:         TCGA-AB-2807                 M1                   181
    ##    Overall_Survival_Status
    ## 1:                       1
    ## 2:                       1
    ## 3:                       0
    ## 4:                       1
    ## 5:                       1
    ## 6:                       1

    #Shows all fields in MAF
    getFields(laml)

    ##  [1] "Hugo_Symbol"            "Entrez_Gene_Id"         "Center"                
    ##  [4] "NCBI_Build"             "Chromosome"             "Start_Position"        
    ##  [7] "End_Position"           "Strand"                 "Variant_Classification"
    ## [10] "Variant_Type"           "Reference_Allele"       "Tumor_Seq_Allele1"     
    ## [13] "Tumor_Seq_Allele2"      "Tumor_Sample_Barcode"   "Protein_Change"        
    ## [16] "i_TumorVAF_WU"          "i_transcript_name"

    #Writes maf summary to an output file with basename laml.
    write.mafSummary(maf = laml, basename = 'laml')

Visualization
-------------

This section covers multiple visualization methods that one can perfrom
with maftools. Specifically, plotmafSummary demonstrates number of
variants in each sample as a stacked barplot and variant types as a
boxplot summarized by Variant\_Classification.

    plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

![](cancerGenomics_files/figure-markdown_strict/unnamed-chunk-10-1.png)

### Drawing oncoplots

Oncoplots, or waterfall plots are also possible in maftools

    #oncoplot for top ten mutated genes.
    oncoplot(maf = laml, top = 10)

![](cancerGenomics_files/figure-markdown_strict/unnamed-chunk-11-1.png)

### Oncostrip

Oncostrip is used to visualize a set of genes, with the color strips
representing the mutations.

    oncostrip(maf = laml, genes = c('DNMT3A','NPM1', 'RUNX1'))

![](cancerGenomics_files/figure-markdown_strict/unnamed-chunk-12-1.png)

### Transition and transversion

This titv function returns a summary of trasnversion and transition
classified SNPs. Data can be presented as a bar plot as well.

    laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
    #plot titv summary
    plotTiTv(res = laml.titv)

![](cancerGenomics_files/figure-markdown_strict/unnamed-chunk-13-1.png)
\#\#\# Lollipop plot

These plots can be used to visualize amoni acid changes, whic stands for
the mutation spots on protein structure. Site preference for the
mutation on the spot of the protein can be mutational hot spots

    #lollipop plot for DNMT3A, which is one of the most frequent mutated gene in Leukemia.
    lollipopPlot(maf = laml, gene = 'DNMT3A', AACol = 'Protein_Change', showMutationRate = TRUE)

    ## 3 transcripts available. Use arguments refSeqID or proteinID to manually specify tx name.

    ##      HGNC refseq.ID protein.ID aa.length
    ## 1: DNMT3A NM_175629  NP_783328       912
    ## 2: DNMT3A NM_022552  NP_072046       912
    ## 3: DNMT3A NM_153759  NP_715640       723

    ## Using longer transcript NM_175629 for now.

    ## Removed 3 mutations for which AA position was not available

![](cancerGenomics_files/figure-markdown_strict/unnamed-chunk-14-1.png)

### Labelling points

We can label the plots iwth the amino acids - if everyhting is chosen,
all points will be annotated.

    lollipopPlot(maf = laml, gene = 'KIT', AACol = 'Protein_Change', labelPos = 816, refSeqID = 'NM_000222')

![](cancerGenomics_files/figure-markdown_strict/unnamed-chunk-15-1.png)

### Rainfall plots

Cancer genomes are charachterized by genomic loci with localized
hypermutations. rainfall plots are used to visualize inter variant
distance on a linera genomic scale

    brca <- system.file("extdata", "brca.maf.gz", package = "maftools")
    brca = read.maf(maf = brca, verbose = FALSE)

    rainfallPlot(maf = brca, detectChangePoints = TRUE, pointSize = 0.6)

    ## Processing TCGA-A8-A08B..

    ## Kataegis detected at:

    ##    Chromosome Start_Position End_Position nMuts Avg_intermutation_dist Size
    ## 1:          8       98129348     98133560     7               702.0000 4212
    ## 2:          8       98398549     98403536     9               623.3750 4987
    ## 3:          8       98453076     98456466     9               423.7500 3390
    ## 4:          8      124090377    124096810    22               306.3333 6433
    ## 5:         12       97436055     97439705     7               608.3333 3650
    ## 6:         17       29332072     29336153     8               583.0000 4081
    ##    Tumor_Sample_Barcode C>G C>T
    ## 1:         TCGA-A8-A08B   4   3
    ## 2:         TCGA-A8-A08B   1   8
    ## 3:         TCGA-A8-A08B   1   8
    ## 4:         TCGA-A8-A08B   1  21
    ## 5:         TCGA-A8-A08B   4   3
    ## 6:         TCGA-A8-A08B   4   4

![](cancerGenomics_files/figure-markdown_strict/unnamed-chunk-17-1.png)

### comparing mutation load against TCGA cohorts

TCGA contains over 30 different cancer cohorts and median mutation load
across them varies from as low as 7 per exome. tcgaCompare draws
distribution of variants compiled from over 10,000 WXS samples across 33
TCGA landmark cohorts (<span class="citeproc-not-found"
data-reference-id="Mayakonda">**???**</span>).

    laml.mutload = tcgaCompare(maf = laml, cohortName = 'Example-LAML')

    ## Performing pairwise t-test for differences in mutation burden..

![](cancerGenomics_files/figure-markdown_strict/unnamed-chunk-18-1.png)

### Plotting VAF

    plotVaf(maf = laml, vafCol = 'i_TumorVAF_WU')

![](cancerGenomics_files/figure-markdown_strict/unnamed-chunk-19-1.png)

### Genecloud

Visualization of the mutated genes - the size if proportional to the
number of samples it is mutated in.

    geneCloud(input = laml, minMut = 3)

![](cancerGenomics_files/figure-markdown_strict/unnamed-chunk-20-1.png)

9 - Analysis
------------

somaticInteractions function is used to estimate co ocurring or mutually
exclusive sets of genes through Fisher’s exeact test (<span
class="citeproc-not-found" data-reference-id="Mayakonda">**???**</span>)

    #exclusive/co-occurance event analysis on top 10 mutated genes. 
    somaticInteractions(maf = laml, top = 25, pvalue = c(0.05, 0.1))

![](cancerGenomics_files/figure-markdown_strict/unnamed-chunk-21-1.png)

    ##      gene1  gene2       pValue oddsRatio  00 11 01 10              Event
    ##   1: ASXL1  RUNX1 0.0001541586 55.215541 176  4 12  1       Co_Occurence
    ##   2:  IDH2  RUNX1 0.0002809928  9.590877 164  7  9 13       Co_Occurence
    ##   3:  IDH2  ASXL1 0.0004030636 41.077327 172  4  1 16       Co_Occurence
    ##   4:  FLT3   NPM1 0.0009929836  3.763161 125 17 16 35       Co_Occurence
    ##   5:  SMC3 DNMT3A 0.0010451985 20.177713 144  6 42  1       Co_Occurence
    ##  ---                                                                    
    ## 296: PLCE1  ASXL1 1.0000000000  0.000000 184  0  5  4 Mutually_Exclusive
    ## 297: RAD21  FAM5C 1.0000000000  0.000000 183  0  5  5 Mutually_Exclusive
    ## 298: PLCE1  FAM5C 1.0000000000  0.000000 184  0  5  4 Mutually_Exclusive
    ## 299: PLCE1  RAD21 1.0000000000  0.000000 184  0  5  4 Mutually_Exclusive
    ## 300:  EZH2  PLCE1 1.0000000000  0.000000 186  0  4  3 Mutually_Exclusive
    ##              pair event_ratio
    ##   1: ASXL1, RUNX1        4/13
    ##   2:  IDH2, RUNX1        7/22
    ##   3:  ASXL1, IDH2        4/17
    ##   4:   FLT3, NPM1       17/51
    ##   5: DNMT3A, SMC3        6/43
    ##  ---                         
    ## 296: ASXL1, PLCE1         0/9
    ## 297: FAM5C, RAD21        0/10
    ## 298: FAM5C, PLCE1         0/9
    ## 299: PLCE1, RAD21         0/9
    ## 300:  EZH2, PLCE1         0/7

### Detecting cancer driver genes based on positional clustering

oncodrive is based on algorithm oncodriveCLUST which was originally
implemented in Python. Concept is based on the fact that most of the
variants in cancer causing genes are enriched at few specific loci (aka
hot-spots). This method takes advantage of such positions to identify
cancer genes (Tamborero, Gonzalez-Perez, and Lopez-Bigas 2013).

    laml.sig = oncodrive(maf = laml, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')

    ## Warning in oncodrive(maf = laml, AACol = "Protein_Change", minMut = 5,
    ## pvalMethod = "zscore"): Oncodrive has been superseeded by OncodriveCLUSTL. See
    ## http://bg.upf.edu/group/projects/oncodrive-clust.php

    ## Estimating background scores from synonymous variants..

    ## Not enough genes to build background. Using predefined values. (Mean = 0.279; SD = 0.13)

    ## Estimating cluster scores from non-syn variants..

    ##   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   9%  |                                                                              |=========                                                             |  13%  |                                                                              |============                                                          |  17%  |                                                                              |===============                                                       |  22%  |                                                                              |==================                                                    |  26%  |                                                                              |=====================                                                 |  30%  |                                                                              |========================                                              |  35%  |                                                                              |===========================                                           |  39%  |                                                                              |==============================                                        |  43%  |                                                                              |=================================                                     |  48%  |                                                                              |=====================================                                 |  52%  |                                                                              |========================================                              |  57%  |                                                                              |===========================================                           |  61%  |                                                                              |==============================================                        |  65%  |                                                                              |=================================================                     |  70%  |                                                                              |====================================================                  |  74%  |                                                                              |=======================================================               |  78%  |                                                                              |==========================================================            |  83%  |                                                                              |=============================================================         |  87%  |                                                                              |================================================================      |  91%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%

    ## Comapring with background model and estimating p-values..

    ## Done !

    head(laml.sig)

    ##    Hugo_Symbol Frame_Shift_Del Frame_Shift_Ins In_Frame_Del In_Frame_Ins
    ## 1:        IDH1               0               0            0            0
    ## 2:        IDH2               0               0            0            0
    ## 3:        NPM1               0              33            0            0
    ## 4:        NRAS               0               0            0            0
    ## 5:       U2AF1               0               0            0            0
    ## 6:         KIT               1               1            0            1
    ##    Missense_Mutation Nonsense_Mutation Splice_Site total MutatedSamples
    ## 1:                18                 0           0    18             18
    ## 2:                20                 0           0    20             20
    ## 3:                 1                 0           0    34             33
    ## 4:                15                 0           0    15             15
    ## 5:                 8                 0           0     8              8
    ## 6:                 7                 0           0    10              8
    ##    AlteredSamples clusters muts_in_clusters clusterScores protLen   zscore
    ## 1:             18        1               18     1.0000000     414 5.546154
    ## 2:             20        2               20     1.0000000     452 5.546154
    ## 3:             33        2               32     0.9411765     294 5.093665
    ## 4:             15        2               15     0.9218951     189 4.945347
    ## 5:              8        1                7     0.8750000     240 4.584615
    ## 6:              8        2                9     0.8500000     976 4.392308
    ##            pval          fdr fract_muts_in_clusters
    ## 1: 1.460110e-08 1.022077e-07              1.0000000
    ## 2: 1.460110e-08 1.022077e-07              1.0000000
    ## 3: 1.756034e-07 8.194826e-07              0.9411765
    ## 4: 3.800413e-07 1.330144e-06              1.0000000
    ## 5: 2.274114e-06 6.367520e-06              0.8750000
    ## 6: 5.607691e-06 1.308461e-05              0.9000000

Size of the points proportional to the number of clusters found in the
gene

    plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = TRUE)

![](cancerGenomics_files/figure-markdown_strict/unnamed-chunk-23-1.png)

### Adding and summarizing pfam domains

Adds pfam domain information to the amino acid changes.

    laml.pfam = pfamDomains(maf = laml, AACol = 'Protein_Change', top = 10)

    ## Warning in pfamDomains(maf = laml, AACol = "Protein_Change", top = 10): Removed
    ## 50 mutations for which AA position was not available

![](cancerGenomics_files/figure-markdown_strict/unnamed-chunk-24-1.png)

    #Protein summary (Printing first 7 columns for display convenience)
    head(laml.pfam$proteinSummary[,1:7, with = FALSE], 10)

    ##       HGNC AAPos Variant_Classification  N total   fraction    DomainLabel
    ##  1: DNMT3A   882      Missense_Mutation 27    54 0.50000000  AdoMet_MTases
    ##  2:   IDH1   132      Missense_Mutation 18    18 1.00000000       PTZ00435
    ##  3:   IDH2   140      Missense_Mutation 17    20 0.85000000       PTZ00435
    ##  4:   FLT3   835      Missense_Mutation 14    52 0.26923077       PKc_like
    ##  5:   FLT3   599           In_Frame_Ins 10    52 0.19230769       PKc_like
    ##  6:  U2AF1    34      Missense_Mutation  7     8 0.87500000        zf-CCCH
    ##  7:   NRAS    61      Missense_Mutation  6    15 0.40000000 H_N_K_Ras_like
    ##  8:    KIT   816      Missense_Mutation  5    10 0.50000000       PTKc_Kit
    ##  9:   NRAS    13      Missense_Mutation  5    15 0.33333333 H_N_K_Ras_like
    ## 10:   FLT3   601           In_Frame_Ins  4    52 0.07692308       PKc_like

    #Domain summary (Printing first 3 columns for display convenience)
    head(laml.pfam$domainSummary[,1:3, with = FALSE], 15)

    ##         DomainLabel nMuts nGenes
    ##  1:        PKc_like    55      5
    ##  2:        PTZ00435    38      2
    ##  3:   AdoMet_MTases    33      1
    ##  4:           7tm_1    24     24
    ##  5:         COG5048    17     17
    ##  6: Cadherin_repeat    16     16
    ##  7:            Runt    16      1
    ##  8:             Dcm    15      1
    ##  9:  H_N_K_Ras_like    15      1
    ## 10:             P53    15      2
    ## 11:             FN3    12      9
    ## 12:            PTPc    10      5
    ## 13:         Tet_JBP    10      1
    ## 14:          bZIP_2    10      1
    ## 15:         COG1100     9      2

### Pan -Cancer comparison

    #MutsigCV results for TCGA-AML
    laml.mutsig <- system.file("extdata", "LAML_sig_genes.txt.gz", package = "maftools")
    pancanComparison(mutsigResults = laml.mutsig, qval = 0.1, cohortName = 'LAML', inputSampleSize = 200, label = 1)

    ## Significantly mutated genes in LAML (q < 0.1): 23

    ## Significantly mutated genes in PanCan cohort (q <0.1): 114

    ## Significantly mutated genes exclusive to LAML (q < 0.1):

    ##       gene pancan            q nMut log_q_pancan     log_q
    ##  1:  CEBPA  1.000 3.500301e-12   13   0.00000000 11.455895
    ##  2:   EZH2  1.000 7.463546e-05    3   0.00000000  4.127055
    ##  3: GIGYF2  1.000 6.378338e-03    2   0.00000000  2.195292
    ##  4:    KIT  0.509 1.137517e-05    8   0.29328222  4.944042
    ##  5:   PHF6  0.783 6.457555e-09    6   0.10623824  8.189932
    ##  6: PTPN11  0.286 7.664584e-03    9   0.54363397  2.115511
    ##  7:  RAD21  0.929 1.137517e-05    5   0.03198429  4.944042
    ##  8:  SMC1A  0.801 2.961696e-03    6   0.09636748  2.528460
    ##  9:   TET2  0.907 2.281625e-13   17   0.04239271 12.641756
    ## 10:    WT1  1.000 2.281625e-13   12   0.00000000 12.641756

![](cancerGenomics_files/figure-markdown_strict/unnamed-chunk-27-1.png)

    ##          gene   pancan            q nMut log_q_pancan    log_q
    ##   1:   ACVR1B 6.11e-02 1.000000e+00    0     1.213959  0.00000
    ##   2:     AKT1 2.68e-10 1.000000e+00    0     9.571865  0.00000
    ##   3:      APC 1.36e-13 1.000000e+00    0    12.866461  0.00000
    ##   4:    APOL2 7.96e-03 1.000000e+00    0     2.099087  0.00000
    ##   5: ARHGAP35 2.32e-12 1.000000e+00    1    11.634512  0.00000
    ##  ---                                                          
    ## 120:    U2AF1 4.07e-08 4.503311e-13    8     7.390406 12.34647
    ## 121:      VHL 2.32e-12 1.000000e+00    0    11.634512  0.00000
    ## 122:      WT1 1.00e+00 2.281625e-13   12     0.000000 12.64176
    ## 123:   ZNF180 8.60e-02 1.000000e+00    0     1.065502  0.00000
    ## 124:   ZNF483 2.37e-02 1.000000e+00    0     1.625252  0.00000

### Survival analysis

**Mutation in any given genes**

Survival analysis is an essential part of cohort based sequencing
projects.

    #Survival analysis based on grouping of DNMT3A mutation status
    mafSurvival(maf = laml, genes = 'DNMT3A', time = 'days_to_last_followup', Status = 'Overall_Survival_Status', isTCGA = TRUE)

    ## Looking for clinical data in annoatation slot of MAF..

    ## Number of mutated samples for given genes:

    ## DNMT3A 
    ##     48

    ## Removed 11 samples with NA's

    ## Median survival..

    ##     Group medianTime   N
    ## 1: Mutant        245  45
    ## 2:     WT        396 137

![](cancerGenomics_files/figure-markdown_strict/unnamed-chunk-28-1.png)

**Predict genesets associated with survival**

Identify set of genes which results in poor survival

    #Using top 20 mutated genes to identify a set of genes (of size 2) to predict poor prognostic groups
    prog_geneset = survGroup(maf = laml, top = 20, geneSetSize = 2, time = "days_to_last_followup", Status = "Overall_Survival_Status", verbose = FALSE)

    ## Removed 11 samples with NA's

    print(prog_geneset)

    ##     Gene_combination P_value    hr  WT Mutant
    ##  1:      FLT3_DNMT3A 0.00104 2.510 164     18
    ##  2:      DNMT3A_SMC3 0.04880 2.220 176      6
    ##  3:      DNMT3A_NPM1 0.07190 1.720 166     16
    ##  4:      DNMT3A_TET2 0.19600 1.780 176      6
    ##  5:        FLT3_TET2 0.20700 1.860 177      5
    ##  6:        NPM1_IDH1 0.21900 0.495 176      6
    ##  7:      DNMT3A_IDH1 0.29300 1.500 173      9
    ##  8:       IDH2_RUNX1 0.31800 1.580 176      6
    ##  9:        FLT3_NPM1 0.53600 1.210 165     17
    ## 10:      DNMT3A_IDH2 0.68000 0.747 178      4
    ## 11:      DNMT3A_NRAS 0.99200 0.986 178      4

Above results show a combination (N = 2) of genes which are associated
with poor survival (P &lt; 0.05).

    mafSurvGroup(maf = laml, geneSet = c("DNMT3A", "FLT3"), time = "days_to_last_followup", Status = "Overall_Survival_Status")

    ## Looking for clinical data in annoatation slot of MAF..

    ## Removed 11 samples with NA's

    ## Median survival..

    ##     Group medianTime   N
    ## 1: Mutant      242.5  18
    ## 2:     WT      379.5 164

![](cancerGenomics_files/figure-markdown_strict/unnamed-chunk-31-1.png)

### Comparing two cohorts

Detecting the mutation pattern in order to compare nultiple cohorts

    #Primary APL MAF
    primary.apl = system.file("extdata", "APL_primary.maf.gz", package = "maftools")
    primary.apl = read.maf(maf = primary.apl)

    ## -Reading
    ## -Validating
    ## --Non MAF specific values in Variant_Classification column:
    ##   ITD
    ## -Silent variants: 45 
    ## -Summarizing
    ## -Processing clinical data
    ## --Missing clinical data
    ## -Finished in 0.096s elapsed (0.346s cpu)

    #Relapse APL MAF
    relapse.apl = system.file("extdata", "APL_relapse.maf.gz", package = "maftools")
    relapse.apl = read.maf(maf = relapse.apl)

    ## -Reading
    ## -Validating
    ## --Non MAF specific values in Variant_Classification column:
    ##   ITD
    ## -Silent variants: 19 
    ## -Summarizing
    ## -Processing clinical data
    ## --Missing clinical data
    ## -Finished in 0.071s elapsed (0.277s cpu)

    #Considering only genes which are mutated in at-least in 5 samples in one of the cohort to avoid bias due to genes mutated in single sample.
    pt.vs.rt <- mafCompare(m1 = primary.apl, m2 = relapse.apl, m1Name = 'Primary', m2Name = 'Relapse', minMut = 5)
    print(pt.vs.rt)

    ## $results
    ##    Hugo_Symbol Primary Relapse         pval         or       ci.up      ci.low
    ## 1:         PML       1      11 1.529935e-05 0.03537381   0.2552937 0.000806034
    ## 2:        RARA       0       7 2.574810e-04 0.00000000   0.3006159 0.000000000
    ## 3:       RUNX1       1       5 1.310500e-02 0.08740567   0.8076265 0.001813280
    ## 4:        FLT3      26       4 1.812779e-02 3.56086275  14.7701728 1.149009169
    ## 5:      ARID1B       5       8 2.758396e-02 0.26480490   0.9698686 0.064804160
    ## 6:         WT1      20      14 2.229087e-01 0.60619329   1.4223101 0.263440988
    ## 7:        KRAS       6       1 4.334067e-01 2.88486293 135.5393108 0.337679367
    ## 8:        NRAS      15       4 4.353567e-01 1.85209500   8.0373994 0.553883512
    ## 9:      ARID1A       7       4 7.457274e-01 0.80869223   3.9297309 0.195710173
    ##         adjPval
    ## 1: 0.0001376942
    ## 2: 0.0011586643
    ## 3: 0.0393149868
    ## 4: 0.0407875250
    ## 5: 0.0496511201
    ## 6: 0.3343630535
    ## 7: 0.4897762916
    ## 8: 0.4897762916
    ## 9: 0.7457273717
    ## 
    ## $SampleSummary
    ##     Cohort SampleSize
    ## 1: Primary        124
    ## 2: Relapse         58

**Forest plots** Above results show two genes PML and RARA which are
highly mutated in Relapse APL compared to Primary APL. We can visualize
them.

    forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.1, color = c('royalblue', 'maroon'), geneFontSize = 0.8)

![](cancerGenomics_files/figure-markdown_strict/unnamed-chunk-34-1.png)

**Co-onco plots**

    genes = c("PML", "RARA", "RUNX1", "ARID1B", "FLT3")
    coOncoplot(m1 = primary.apl, m2 = relapse.apl, m1Name = 'PrimaryAPL', m2Name = 'RelapseAPL', genes = genes, removeNonMutated = TRUE)

![](cancerGenomics_files/figure-markdown_strict/unnamed-chunk-35-1.png)

**Lollipop plots-2**

Showing gene wise difference

    lollipopPlot2(m1 = primary.apl, m2 = relapse.apl, gene = "PML", AACol1 = "amino_acid_change", AACol2 = "amino_acid_change", m1_name = "Primary", m2_name = "Relapse")

    ## Gene: PML

    ## 9 transcripts available. Use arguments refSeqID or proteinID to manually specify tx name.

    ##    HGNC refseq.ID protein.ID aa.length
    ## 1:  PML NM_033238  NP_150241       882
    ## 2:  PML NM_002675  NP_002666       633
    ## 3:  PML NM_033249  NP_150252       585
    ## 4:  PML NM_033247  NP_150250       435
    ## 5:  PML NM_033239  NP_150242       829
    ## 6:  PML NM_033250  NP_150253       781
    ## 7:  PML NM_033240  NP_150243       611
    ## 8:  PML NM_033244  NP_150247       560
    ## 9:  PML NM_033246  NP_150249       423

    ## Using longer transcript NM_033238 for now.
    ## 9 transcripts available. Use arguments refSeqID or proteinID to manually specify tx name.

    ##    HGNC refseq.ID protein.ID aa.length
    ## 1:  PML NM_033238  NP_150241       882
    ## 2:  PML NM_002675  NP_002666       633
    ## 3:  PML NM_033249  NP_150252       585
    ## 4:  PML NM_033247  NP_150250       435
    ## 5:  PML NM_033239  NP_150242       829
    ## 6:  PML NM_033250  NP_150253       781
    ## 7:  PML NM_033240  NP_150243       611
    ## 8:  PML NM_033244  NP_150247       560
    ## 9:  PML NM_033246  NP_150249       423

    ## Using longer transcript NM_033238 for now.

![](cancerGenomics_files/figure-markdown_strict/unnamed-chunk-36-1.png)

Clinical enrichment analysis
----------------------------

takes in clinical features and performs enrichment analysis. Various
groupwise and pairwise comparison.

    fab.ce = clinicalEnrichment(maf = laml, clinicalFeature = 'FAB_classification')

    ## Sample size per factor in FAB_classification:

    ## 
    ## M0 M1 M2 M3 M4 M5 M6 M7 
    ## 19 44 44 21 39 19  3  3

    #Results are returned as a list. Significant associations p-value < 0.05
    fab.ce$groupwise_comparision[p_value < 0.05]

    ##    Hugo_Symbol Group1 Group2 n_mutated_group1 n_mutated_group2      p_value
    ## 1:        IDH1     M1   Rest         11 of 44         7 of 149 0.0002597371
    ## 2:        TP53     M7   Rest           3 of 3        12 of 190 0.0003857187
    ## 3:      DNMT3A     M5   Rest         10 of 19        38 of 174 0.0057610493
    ## 4:       CEBPA     M2   Rest          7 of 44         6 of 149 0.0117352110
    ## 5:       RUNX1     M0   Rest          5 of 19        11 of 174 0.0117436825
    ## 6:        NPM1     M5   Rest          7 of 19        26 of 174 0.0248582372
    ## 7:       CEBPA     M1   Rest          6 of 44         7 of 149 0.0478737468
    ##    OR_low   OR_high       fdr
    ## 1:      0 0.3926994 0.0308575
    ## 2:      0 0.1315271 0.0308575
    ## 3:      0 0.6406007 0.3072560
    ## 4:      0 0.6874270 0.3757978
    ## 5:      0 0.6466787 0.3757978
    ## 6:      0 0.8342897 0.6628863
    ## 7:      0 0.9869971 1.0000000

Above results shows IDH1 mutations are enriched in M1 subtype of
leukemia compared to rest of the cohort. Similarly DNMT3A is in M5,
RUNX1 is in M0

    plotEnrichmentResults(enrich_res = fab.ce, pVal = 0.05)

![](cancerGenomics_files/figure-markdown_strict/unnamed-chunk-39-1.png)

Drug -Gene interaction
----------------------

This plot shows potential druggable gene categories along with upto top
5 genes involved in them.

    dgi = drugInteractions(maf = laml, fontSize = 0.75)

![](cancerGenomics_files/figure-markdown_strict/unnamed-chunk-40-1.png)

    dnmt3a.dgi = drugInteractions(genes = "DNMT3A", drugs = TRUE)

    ## Number of claimed drugs for given genes:
    ##      Gene N
    ## 1: DNMT3A 7

    #Printing selected columns.
    dnmt3a.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]

    ##      Gene interaction_types    drug_name drug_claim_name
    ## 1: DNMT3A                                            N/A
    ## 2: DNMT3A                   DAUNORUBICIN    Daunorubicin
    ## 3: DNMT3A                     DECITABINE      Decitabine
    ## 4: DNMT3A                     IDARUBICIN      IDARUBICIN
    ## 5: DNMT3A                     DECITABINE      DECITABINE
    ## 6: DNMT3A         inhibitor   DECITABINE   CHEMBL1201129
    ## 7: DNMT3A         inhibitor  AZACITIDINE      CHEMBL1489

Oncogenic signaling pathway
---------------------------

Check for enrichment of known ocogenic singaling pathways in TCGA
cohorts/

    OncogenicPathways(maf = laml)

    ## Pathway alteration fractions

    ##        Pathway  N n_affected_genes fraction_affected
    ##  1:    RTK-RAS 85               18        0.21176471
    ##  2:      Hippo 38                7        0.18421053
    ##  3:      NOTCH 71                6        0.08450704
    ##  4:        MYC 13                3        0.23076923
    ##  5:        WNT 68                3        0.04411765
    ##  6:       TP53  6                2        0.33333333
    ##  7:       NRF2  3                1        0.33333333
    ##  8:       PI3K 29                1        0.03448276
    ##  9: Cell_Cycle 15                0        0.00000000
    ## 10:   TGF-Beta  7                0        0.00000000

![](cancerGenomics_files/figure-markdown_strict/unnamed-chunk-43-1.png)

Complete pathway visualization

    PlotOncogenicPathways(maf = laml, pathways = "RTK-RAS")

![](cancerGenomics_files/figure-markdown_strict/unnamed-chunk-44-1.png)

Mutational Signatures
---------------------

<!-- ```{r} -->
<!-- #Requires BSgenome object -->
<!-- library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE) -->
<!-- laml.tnm = trinucleotideMatrix(maf = laml, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19") -->
<!-- ``` -->
<!-- Analyze the differences in mutational patterns between APOBEC enriched and non-APOBEC enriched samples. -->
<!-- ```{r} -->
<!-- plotApobecDiff(tnm = laml.tnm, maf = laml, pVal = 0.2) -->
<!-- ``` -->
<!-- ### Signature analysis  -->
<!-- ```{r} -->
<!-- # encoutnered issues with parallel processing, neede to reinstall the following packages and and reload the libraries for NMF and foreach  -->
<!-- #install.packages(c( "foreach", "doParallel") ) -->
<!-- ``` -->
<!-- ```{r} -->
<!-- library(NMF) -->
<!-- library(foreach) -->
<!-- library(devtools) -->
<!-- # #install.packages('car') -->
<!-- # library(car) -->
<!-- ``` -->
<!-- ```{r} -->
<!-- laml.sign = estimateSignatures(mat = laml.tnm, nTry = 6, pConstant = 1) -->
<!-- laml.sig = extractSignatures(mat = laml.tnm, n = 3, pConstant = 1) -->
<!-- ``` -->
<!-- ```{r} -->
<!-- #Compare against original 30 signatures  -->
<!-- laml.og30.cosm = compareSignatures(nmfRes = laml.sig, sig_db = "legacy") -->
<!-- ``` -->
<!-- ```{r} -->
<!-- #Compate against updated version3 60 signatures -->
<!-- laml.v3.cosm = maftools::compareSignatures(nmfRes = laml.sig, sig_db = "SBS") -->
<!-- ``` -->
<!-- ```{r} -->
<!-- library(pheatmap) -->
<!-- pheatmap::pheatmap(mat = laml.og30.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures") -->
<!-- ``` -->
<!-- ```{r} -->
<!-- plotSignatures(nmfRes = laml.sig, title_size = 0.8) -->
<!-- ``` -->
<!-- ### Signature enrichment analysis  -->
<!-- Signatures can further be assigned to samples and enrichment analysis can be performd using signatureEnrichment funtion, which identifies mutations enriched in every signature identified. -->
<!-- ```{r} -->
<!-- laml.se = signatureEnrichment(maf = laml, sig_res = laml.sig) -->
<!-- ``` -->
<!-- ```{r} -->
<!-- plotEnrichmentResults(enrich_res = laml.se, pVal = 0.05) -->
<!-- ``` -->

References:
-----------

Tamborero, David, Abel Gonzalez-Perez, and Nuria Lopez-Bigas. 2013.
“OncodriveCLUST: Exploiting the Positional Clustering of Somatic
Mutations to Identify Cancer Genes.” *Bioinformatics* 29 (18): 2238–44.
<https://doi.org/10.1093/bioinformatics/btt395>.
