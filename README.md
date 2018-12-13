# final.project
The sequences are from sorted skin cells that were taken from mouse wound margin and normal unwounded skin 


1.	Once logged in, I transferred to an interactive node on School of Medicine (SOM)
$ qrsh -q som -pe openmp 4-64
2.	Make a directory for the final project

$ cd  /data/users/gkashgar
$ mkdir final 
$ cd final

3.	Data where downloaded  from the link that was sent by the sequencing core  facility. 
4.	
$ wget http://hts.igb.uci.edu/gkashgar17091539

5.	Check downloading quality 

$ mdsum5* .txt.gz

6.	Unzip the .gz files 

$ gunzip *.txt.gz

7.	Obtain the indexed mouse genome (mm10). Since we will be using Bowtie2 to do the alignment, we will get the indexed genome from Bowtie2 website http://bowtie-bio.sourceforge.net/index.shtml. The files will need to be uncompressed and the result should be 4 files ending in *.bt2 and 2 files ending in *.rev.1/2.bt2. You can verify that they are there using the ls command.

$ wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/mm10.zip
$ unzip mm10.zip

8.	annotation for the mm10 genome. 
9.	
$ wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz
$ gunzip Mus_musculus_UCSC_mm10.tar.gz
$ tar –xvf Mus_musculus_UCSC_mm10.tar

10.	Alignment for sequencing to mouse genome 

$ module load bowtie2/2.2.3
$ module load tophat/2.0.12
$ module load samtools/0.1.19
$ tophat2 -T –G mm10.gtf -o c3_r1_thout mm10  4R042-L7-P1-ATTACTCG-TATAGCCT-Sequences.txt 
$ tophat2 -T –G mm10.gtf -o c2_r1_thout mm10 4R042-L7-P5-ATTACTCG-ATAGAGGC-Sequences.txt

11.	Genome assembly of expressed genes and transcripts using cufflinks 

$ module load cufflinks
$ cufflinks -o wt1_un_clout  wt1_un_thout /accepted_hits.bam
$ cufflinks -o wt1_w_clout  wt1_w_thout /accepted_hits.bam

12.	create a file assemblies.txt that lists the assembly file for each samples. 

$ touch assemblies.txt
$ nano assemblies.txt

./ wt1_un _clout/transcripts.gtf
./ wt1_w _clout/transcripts.gtf

13.	Run Cuffmerge to create a single merged transcriptome annotation.

$ cuffmerge -g genes.gtf -s ref.fa assemblies.txt

14.	Identify differentially expressed genes and transcripts using Cuffdiff. 

$ cuffdiff -o diff_out1 -b ref.fa -L c1,c2 –u merged_asm1/merged.gtf \ ./ wt1_un _thout/accepted_hits.bam,./ wt1_w _thout/accepted_hits.bam 


15.	Sort the bam files by name using the –n option.  
$ module load samtools/0.1.19
$ cd wt1_un _thout
$ samtools sort –n accepted_hits.bam accepted_hits.nsorted &
$ samtools sort accepted_hits.bam accepted_hits.sorted &
$ samtools index accepted_hits.sorted.bam &
$ cd wt1_w _thout
$ samtools sort –n accepted_hits.bam accepted_hits.nsorted &
$ samtools sort accepted_hits.bam accepted_hits.sorted &
$ samtools index accepted_hits.sorted.bam

16.	Generate a count matrix using HTseq-count.

$ module load enthought_python/7.3.2
$ htseq-count -f bam –s no wt1_un _thout/accepted_hits.sorted.bam mm10.gtf > htseq_out/ wt1_un.counts &
$ htseq-count -f bam –s no wt1_w _thout/accepted_hits.nsorted.bam mm10.gtf > htseq_out/ wt1_w.counts &

17.	Making a heatmap of the count matrix with DESeq2.

$ module load R/3.2.2
   $ R
   source(“http://bioconductor.org/biocLite.R”)
  it might ask if its okay to install so type Y
  biocLite(“DESeq2”)
  if it asks if you want to use a personal library type Y
  if updated is needed type n
   library(DESeq2)
   outputPrefix <- “d1_DESeq2”
   directory <- “/som/gkashgar/rna-seq/24h/d1”
   ATR_ d1 <- c(“wt1_un.counts”, “wt1_w.counts”)
  condition_d1 <- factor(substr(ATR_d1, 1, 2))
   table_d1 <- data.frame(sampleName=ATR_d1, fileName=ATR_d1, condition=condition_d1)
   table_d1
   allows you to view your table to verify that you set it up correctly.
   dds_d1 <- DESeqDataSetFromHTSeqCount(sampleTable=table_d1, directory=directory, design=~condition)
  colData(dds_d1)$condition <- factor(colData(dds_d1)$condition, levels = c(‘un’,‘w”))
 following step is the gut of DESeq2 analysis
  >ddsc1c2 <- DESeq(dds_d1)
   results_d1 <- results(ddsd1)
  order results by padj value (most significant to least)
   results_d1= subset(results_d1, padj<0.05)
   results_d1 <- results_d1[order(results_d1$padj), ]
 can view the results with the following command (optional)
  should see DataFrame of baseMean, log2Foldchange, stat, pval, padj
   head (results_d1)
 because this data was analyzed with count values and not with tuxedo(cuffdiff), we need to save these results and normalized reads to csv
  resultsdata_d1 <-merge(as.data.frame(results_d1), as.data.frame(counts(ddsd1,normalized =TRUE)), by=‘row.names’, sort = FALSE)
   names(resultsdata_d1)[1] <- ‘gene’
   head (resultsdata_d1)
   write.csv(resultsdata_d1, file = paste0(outputPrefix, “-d1_results_normalized.csv”))



