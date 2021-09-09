# Biv_ChIP

ChIP-Seq Data analysis 

**Step1 :** Generate the read tagdensity profiles for each promoter from bed formatted read alignments from H3K27me3 ChIP-Seq data in same manner as below

>perl genePromoterTagDensity.pl -i refSeqGenes-mm9 -b H3K27me3/H3K27me3_0hr_ES1546.bed -s /ddn/gs1/group/jothi/jothidata/UCSC/Data/mm9/mm9_chrom_lengths.txt -n 1 -F 200 -r 0 -w 2000 -o H3K27me3/H3K27me3_0hr_ES1546_TagDensity.txt &

    USAGE: genePromoterTagDensity.pl -i <file> -b <file> -s <file> -n <1/0> -d <1/0> -w <number> -o <file>
	-i      File containing list of genes information.      
	-b      *.bed file 1 with tags in BED format
        -s      *File containing chromsome lengths in UCSC format (https://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/mm9.chrom.sizes)
        -n      1 for normalization by number of tags, 0 otherwise
        -F      Average fragment length for tag shifting, e.g., 160, 200, 300
        -r      1 (-1/0) if same (opposite, both resp.) strand reads need to be considered; e.g., Gro-Seq data
        -w       number of base pairs to be considered on both sides of the TSS; e.g., 500
        -o      *name of the output file
File formats: 

	(1) File with gene information (-i) should be provided in the following format 
  	 name    chrom   strand  txStart txEnd   cdsStart        cdsEnd  exonCount       exonStarts      exonEnds        name2
  	 NM_028778       chr1    +       134212701       134230065       134212806       134228958       7       134212701,134221529,134224273,134224707,134226534,134227135,134227897,  134213049,134221650,134224425,134224773,134226654,134227268,134230065,    Nuak2
  	 .....
	(2)File with chromosome length information (-s) should be provided in the following format
   
  	 chr1    197195432  
  	 chr2    181748087
  	 .....
	
	(3)File with read alignment data (-b) should be provided in the BED format as follows
   
  	 chr1    3018108 3018158 U       0       -
  	 chr1    3001762 3001812 U       0       -
   	 .............
	   
The code outputs two files. For the given example-

	(1) H3K27me3_0hr_ES1546_TagDensity.txt -> Average tagdensity vector for each 100bp bin within the promoter window
	(2) H3K27me3_0hr_ES1546_TagDensity.txt_geneWise.txt-> For each gene promoter the tagdnesity vector for each of the 100bp non-overlapping bin 

**Step2:** Generate the read tagdensity profiles for each promoter from bed formatted read alignments from genomic-input data in same manner
>perl genePromoterTagDensity.pl -i refSeqGenes-mm9 -b INPUT/INPUT_0hr_ES1554.bed -s /ddn/gs1/group/jothi/jothidata/UCSC/Data/mm9/mm9_chrom_lengths.txt -n 1 -F 200 -r 0 -w 2000 -o INPUT/INPUT_0hr_ES1554_TagDensity_2kb.txt &

**Step3:** generate List of the promoters enriched for H3K27me3 as demonstrated in the example below 
>perl EnrichedGenePromoterReporting.pl -I ./INPUT/INPUT_0hr_ES1554_TagDensity_2kb.txt_geneWise.txt -C ./H3K27me3/H3K27me3_0hr_ES1546_TagDensity.txt_geneWise.txt -D 1 -F 3 -o H3K27me3_2kb_enriched_genes_0hr.txt

**Step4:** Repeat Steps 1-3 for H3K4me3 modification and generate List of the promoters enriched for H3K4me3 as demonstrated in the example below 

>perl genePromoterTagDensity.pl -i refSeqGenes-mm9 -b H3K4me3/H3K4me3_0hr_ES1530.bed -s /ddn/gs1/group/jothi/jothidata/UCSC/Data/mm9/mm9_chrom_lengths.txt -n 1 -F 200 -r 0 -w 500 -o H3K4me3/H3K4me3_0hr_ES1530_TagDensity.txt &

>perl genePromoterTagDensity.pl -i refSeqGenes-mm9 -b INPUT/INPUT_0hr_ES1554.bed -s /ddn/gs1/group/jothi/jothidata/UCSC/Data/mm9/mm9_chrom_lengths.txt -n 1 -F 200 -r 0 -w 500 -o INPUT/INPUT_0hr_ES1554_TagDensity_500bp.txt &

>perl EnrichedGenePromoterReporting.pl -I ./INPUT/INPUT_0hr_ES1554_TagDensity_500bp.txt_geneWise.txt -C ./H3K27me3/H3K4me3_0hr_ES1530_TagDensity.txt_geneWise.txt -D 1 -F 3 -o H3K4me3_500bp_enriched_genes_0hr.txt

**Step5:** Annotate four class chromatin state for a given timepoint as demostrated in the following example

>perl Annotate_Chromatin_state_github.pl -A H3K4me3_enriched_genes_0hr.txt -B H3K27me3_enriched_genes_0hr.txt -P refSeqGenes-mm9 -l mESC_0hr -o mESC_0hr_chromatin.txt
