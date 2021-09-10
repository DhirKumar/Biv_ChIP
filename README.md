# Biv_ChIP

ChIP-Seq Data analysis 

**Step1 :** Generate read/tag density profiles for each gene promoter from BED-formatted ChIP-Seq read alignments, as shown below. This script takes as input, among others, a file containing a set of genes (UCSC format), BED-formatted file containing ChIP-Seq reads (H3K4me3, H3K27me3, or Input) and a file containing chromosome lengths in UCSC format (e.g., https://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/mm9.chrom.sizes). Run this script two separate times to generate read density profiles for H3K4me3 and corresponding Input (using the same value for -w), and repeat this for H3K27me3 and corresponding Input (using the same value for -w) 

>USAGE: \>perl genePromoterTagDensity.pl -i \<input-file\> -b \<input-file\> -s \<input-file\> -n \<1/0\> -F \<int\> -r \<1/0\> -w \<number\> -o \<output-file\>

        -i      input file fontaining list of genes in UCSC format
        -b      input file containing ChIP-Seq reads/tags in BED format
        -s      input file containing chromsome (https://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/mm9.chrom.sizes)
        -n      1 for normalization by total number of reads in BED file, 0 otherwise
        -F      Average ChIP-Seq fragment length for read shifting, e.g., 160, 200, 300
        -r      1 (-1/0) if same (opposite, both resp.) strand reads need to be considered; e.g., Gro-Seq data
        -w      number of base pairs up- and down-stream of TSS, for promoter definition (in multiples of 100); e.g., 500
        -o      output file into which results will be stored

Input file formats: 

	(1) File containing list of genes (-i) should be provided in UCSC format (tab-separated)
	
  	 name    chrom   strand  txStart txEnd   cdsStart        cdsEnd  exonCount       exonStarts      exonEnds        name2
  	 NM_028778       chr1    +       134212701       134230065       134212806       134228958       7       134212701,134221529,134224273,134224707,134226534,134227135,134227897,  134213049,134221650,134224425,134224773,134226654,134227268,134230065,    Nuak2
  	 .....
	 .....


	(2) File containing ChIP-Seq read alignments (-b) should be provided in the BED format (tab-separated)
   
  	 chr1    3018108 3018158 U       0       -
  	 chr1    3001762 3001812 U       0       -
   	 .............


	 (3) File containing chromosome length information (-s) should be provided in the following format (tab-separated)
   
  	 chr1    197195432  
  	 chr2    181748087
  	 .....
	 .....	 


The code outputs two files. For example, if your output file name (-o) is H3K4me3_readDensity.txt, you get the following two output files:

	(1) H3K4me3_readDensity.txt containing average read density (per bp per read sequenced) within each of the 100 bp non-overlapping windows spanning +/-(-w) of TSS, in the following format:
	
	window_no	total_no_reads_within_that_window_for_all_genes	no_bp_spanning_that_window_for_all_genes	avg_read_density_per_bp_per_read_sequenced	sem_of_average
	
	
	(2) H3K4me3_readDensity.txt_geneWise.txt containing a matrix of read densitites for each of the 100bp windows for each gene promoter, in the following format:
	
	gene_ID	chr	txStart	txEnd	window_1	window_2	...	window_n
	

**Step2:** Generate a list of gene promoters enriched for H3K4me3, using the \*\_geneWise.txt" files for H3K4me3 and corresponding Input, as shown below:

>USAGE: \>perl EnrichedGenePromoterReporting.pl -I <input-file> -C <input-file> -D <int> -F <number> -o <out-file>

        -I      input file containing geneWise red-density matrix for genomic control input, as generated by genePromoterTagDensity.pl
        -C      input file containing geneWise red-density matrix for H3K4me3/H3K27me3, as generated by genePromoterTagDensity.pl
        -D      FDR threshold as a percentage between 0 to 100
        -F      Minimum FoldOverInput threshold (e.g. 3) for H3K4me3/H3K27me3 enrichment consideration
        -o      output file into which results will be stored


**Step3:** Repeat Steps 1-3 for H3K4me3 modification and generate List of the promoters enriched for H3K4me3 as demonstrated in the example below 

>perl genePromoterTagDensity.pl -i refSeqGenes-mm9 -b H3K4me3/H3K4me3_0hr_ES1530.bed -s mm9_chrom_lengths.txt -n 1 -F 200 -r 0 -w 500 -o H3K4me3/H3K4me3_0hr_ES1530_TagDensity.txt &

>perl genePromoterTagDensity.pl -i refSeqGenes-mm9 -b INPUT/INPUT_0hr_ES1554.bed -s mm9_chrom_lengths.txt -n 1 -F 200 -r 0 -w 500 -o INPUT/INPUT_0hr_ES1554_TagDensity_500bp.txt &

>perl EnrichedGenePromoterReporting.pl -I ./INPUT/INPUT_0hr_ES1554_TagDensity_500bp.txt_geneWise.txt -C ./H3K27me3/H3K4me3_0hr_ES1530_TagDensity.txt_geneWise.txt -D 1 -F 3 -o H3K4me3_500bp_enriched_genes_0hr.txt

**Step5:** Annotate four class chromatin state for a given timepoint as demostrated in the following example

>perl Annotate_Chromatin_state_github.pl -A H3K4me3_enriched_genes_0hr.txt -B H3K27me3_enriched_genes_0hr.txt -P refSeqGenes-mm9 -l mESC_0hr -o mESC_0hr_chromatin.txt

	USAGE: Annotate_Chromatin_state.pl -A <file> -B <file> -P <file> -o <file>
        -A      File containing gene/transcript-IDS enriched for H3K4me3
        -B      File containing gene/transcript-IDS enriched for H3K27me3
        -P      File containing promoter information with transcript IDs in the first column of the file
        -l      Column label for the cell type. e.g. mESC_0hr
        -o      *name of the output file
	
