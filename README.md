# Biv_ChIP

ChIP-Seq Data analysis 

Step1 : Generate the read tagdensity profiles for each promoter from bed formatted read alignments from H3K27me3 ChIP-Seq data in same manner as below
>perl genePromoterTagDensity.pl -i refSeqGenes-mm9 -b H3K27me3/H3K27me3_0hr_ES1546.bed -s /ddn/gs1/group/jothi/jothidata/UCSC/Data/mm9/mm9_chrom_lengths.txt -n 1 -F 200 -r 0 -w 2000 -o H3K27me3/H3K27me3_0hr_ES1546_TagDensity.txt &
The code outputs two files. For the given example-
(1) H3K27me3_0hr_ES1546_TagDensity.txt -> Average tagdensity vector for each 100bp bin within the promoter window
(2) H3K27me3_0hr_ES1546_TagDensity.txt_geneWise.txt-> For each gene promoter the tagdnesity vector for each of the 100bp non-overlapping bin 

Step2: Generate the read tagdensity profiles for each promoter from bed formatted read alignments from genomic-input data in same manner
>perl genePromoterTagDensity.pl -i refSeqGenes-mm9 -b INPUT/INPUT_0hr_ES1554.bed -s /ddn/gs1/group/jothi/jothidata/UCSC/Data/mm9/mm9_chrom_lengths.txt -n 1 -F 200 -r 0 -w 2000 -o INPUT/INPUT_0hr_ES1554_TagDensity_2kb.txt &

Step3: generate List of the promoters enriched for H3K27me3 as demonstrated in the example below 
>perl EnrichedGenePromoterFinding.pl -I ./INPUT/INPUT_0hr_ES1554_TagDensity_2kb.txt_geneWise.txt -C ./H3K27me3/H3K27me3_0hr_ES1546_TagDensity.txt_geneWise.txt -D 1 -F 3 -o H3K27me3_2kb_enriched_genes_0hr.txt

Step4: Repeat Steps 1-3 for H3K4me3 modification and generate List of the promoters enriched for H3K4me3 as demonstrated in the example below 
>perl genePromoterTagDensity.pl -i refSeqGenes-mm9 -b H3K4me3/H3K4me3_0hr_ES1530.bed -s /ddn/gs1/group/jothi/jothidata/UCSC/Data/mm9/mm9_chrom_lengths.txt -n 1 -F 200 -r 0 -w 500 -o H3K4me3/H3K4me3_0hr_ES1530_TagDensity.txt &
>perl PerGeneTagDensityDirectionalReads_K4.pl -i refSeqGenes-mm9 -b INPUT/INPUT_0hr_ES1554.bed -s /ddn/gs1/group/jothi/jothidata/UCSC/Data/mm9/mm9_chrom_lengths.txt -n 1 -F 200 -r 0 -w 500 -o INPUT/INPUT_0hr_ES1554_TagDensity_500bp.txt &
>perl EnrichedGenePromoterFinding.pl -I ./INPUT/INPUT_0hr_ES1554_TagDensity_500bp.txt_geneWise.txt -C ./H3K27me3/H3K4me3_0hr_ES1530_TagDensity.txt_geneWise.txt -D 1 -F 3 -o H3K4me3_500bp_enriched_genes_0hr.txt

Step5:Annotate four class chromatin state for a given timepoint as demostrated in the following example
>perl Annotate_Chromatin_state_github.pl -A H3K4me3_enriched_genes_0hr.txt -B H3K27me3_enriched_genes_0hr.txt -P refSeqGenes-mm9 -l mESC_0hr -o mESC_0hr_chromatin.txt
