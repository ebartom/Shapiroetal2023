#!/bin/bash
#SBATCH -A b1042             ## account (unchanged)
#SBATCH -p genomics          ## "-p" instead of "-q"
#SBATCH -J ShapiroEtAl.codeTest         ## job name
#SBATCH --mail-type=FAIL,TIME_LIMIT_90
#SBATCH --mail-user=jason.shapiro@northwestern.edu
#SBATCH -o "%x.o%j"
#SBATCH -N 1                 ## number of nodes
#SBATCH -n 24                 ## number of cores
#SBATCH -t 24:00:00          ## walltime
#SBATCH --mem=50000


module purge all
module load python/anaconda3.6
module load deeptools
module load bedtools/2.29.1
export PATH=$PATH:/projects/b1025/tools/MACS-1.4.2/bin
export PYTHONPATH=/projects/b1025/tools/MACS-1.4.2/lib/python2.6/site-packages:$PYTHONPATH
module load R/3.3.3
# bedGraphToBigWig is a utility downloaded from the UCSC genome browser.
# Other used tools are either provided perl scripts or deeptools or bedtools utilities.

# The following code was used to generate Figure 6J-M, Extended data figure 5A-D, Extended data figure 6A, Extended data figure 8A-3.

# Primary analysis (alignment and peak calling) was done using Ceto, available at https://github.com/ebartom/NGSbartom, using hg19

#Generating ISPTZ files for POLII/H3K9me2 for  Control and DFO groups (for UCSC genome browser tracks).
# First generate bedgraphs from the bam files, subtracting input samples from the specific ChIP-seq files
# Bam files are not made available on Github, because they are too large.  Fastq files are available from GEO at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214019 and were aligned to hg19 using Bowtie.
echo "Generating ISPTZ tracks for PolII/H3K9me2"
bamCompare -b1 bam/293T_150DFO_H3K9me2_12hr.bam -b2 bam/293T_150DFO_Input_12hr.bam --operation subtract -o bedgraph/293T_150DFO_H3K9me2_12hr_InputSubtracted.bdg -of bedgraph &
bamCompare -b1 bam/293T_Control_H3K9me2_12hr.bam -b2 bam/293T_Control_Input_12hr.bam --operation subtract -o bedgraph/293T_Control_H3K9me2_12hr_InputSubtracted.bdg  -of bedgraph &
bamCompare -b1 bam/293T_150DFO_POLII_12hr.bam -b2 bam/293T_150DFO_Input_12hr.bam --operation subtract -o bedgraph/293T_150DFO_POLII_12hr_InputSubtracted.bdg  -of bedgraph &
bamCompare -b1 bam/293T_Control_POLII_12hr.bam -b2 bam/293T_Control_Input_12hr.bam --operation subtract -o bedgraph/293T_Control_POLII_12hr_InputSubtracted.bdg  -of bedgraph &

wait
# Next, take the input-subtracted bedgraphs above, and replace negative numbers with 0.  Then convert bedgraphs to bigwigs.
echo "Pushing negative values to zero, for input subtracted PolII/H3K9me2"
for sample in "293T_150DFO_H3K9me2_12hr_InputSubtracted" "293T_Control_H3K9me2_12hr_InputSubtracted" "293T_Control_POLII_12hr_InputSubtracted" "293T_150DFO_POLII_12hr_InputSubtracted"
 do
     awk '$4 < 0 {printf "%s\t%d\t%d\t%d\n",$1,$2,$3,0}' bedgraph/$sample.bdg > bedgraph/$sample.pushToZero.bdg
     awk '$4 >= 0' bedgraph/$sample.bdg >> bedgraph/$sample.pushToZero.bdg
     sort -k1,1 -k2,2n bedgraph/$sample.pushToZero.bdg > bedgraph/$sample.pushToZero.sorted.bdg
     /projects/p20742/tools/bin/bedGraphToBigWig bedgraph/$sample.pushToZero.sorted.bdg hg19.chrom.sizes bw/$sample.pushToZero.bw
done

# This ISPTZ analysis was repeated for publicly available data to generate ISPTZ files for shKDM3A3B and shC
# Publicly available data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127624 (KDM3A) and https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71885 (H3K9me2, KDM3B) were downloaded and aligned to generate bam files.
echo "Repeat ISPTZ for publicly available tracks."
bamCompare -b1 New_KDM3B/New_KDM3B/bam/H3K9me2_shC.bam -b2 New_KDM3B/New_KDM3B/bam/Input_shC.bam --pseudocount 1 --operation log2 -o bw/NewKDM3B_shC.H3K9me2_Log2FC_over_Input.bw &
bamCompare -b1 New_KDM3B/New_KDM3B/bam/H3K9me2_sh3A3B.bam -b2 New_KDM3B/New_KDM3B/bam/Input_sh3A3B.bam --pseudocount 1 --operation log2 -o bw/NewKDM3B_sh3A3B.H3K9me2_Log2FC_over_Input.bw &

bamCompare -b1 New_KDM3B/New_KDM3B/bam/H3K9me2_shC.bam -b2 New_KDM3B/New_KDM3B/bam/Input_shC.bam --operation subtract -o bedgraph/NewKDM3B_shC_H3K9me2.InputSubtracted.bdg -of bedgraph &
bamCompare -b1 New_KDM3B/New_KDM3B/bam/H3K9me2_sh3A3B.bam -b2 New_KDM3B/New_KDM3B/bam/Input_sh3A3B.bam --operation subtract -o bedgraph/NewKDM3B_sh3A3B_H3K9me2.InputSubtracted.bdg -of bedgraph &

wait

for sample in "NewKDM3B_shC_H3K9me2.InputSubtracted" "NewKDM3B_sh3A3B_H3K9me2.InputSubtracted"
do
    awk '$4 < 0 {printf "%s\t%d\t%d\t%d\n",$1,$2,$3,0}' bedgraph/$sample.bdg > bedgraph/$sample.pushToZero.bdg
    awk '$4 >= 0' bedgraph/$sample.bdg >> bedgraph/$sample.pushToZero.bdg
    sort -k1,1 -k2,2n bedgraph/$sample.pushToZero.bdg > bedgraph/$sample.pushToZero.sorted.bdg
    /projects/p20742/tools/bin/bedGraphToBigWig bedgraph/$sample.pushToZero.sorted.bdg hg19.chrom.sizes bw/$sample.pushToZero.bw
done

# Calculate logFC bigWigs, using a pseudocount of 1 so that there are no division by 0 errors.
echo "Calculate logFC bigwigs."
 bigwigCompare --bigwig1 bw/293T_150DFO_POLII_12hr_InputSubtracted.pushToZero.bw --bigwig2 bw/293T_Control_POLII_12hr_InputSubtracted.pushToZero.bw --pseudocount 1 --operation log2 -o bw/logFC.150DFOoverControl.POLII.ISPTZ.bw
 bigwigCompare --bigwig1 bw/293T_150DFO_H3K9me2_12hr_InputSubtracted.pushToZero.bw --bigwig2 bw/293T_Control_H3K9me2_12hr_InputSubtracted.pushToZero.bw --pseudocount 1 --operation log2 -o bw/logFC.150DFOoverControl.H3K9me2.ISPTZ.bw

 # Also for the publicly available data.
 bigwigCompare --bigwig1 bw/NewKDM3B_sh3A3B_H3K9me2.InputSubtracted.pushToZero.bw --bigwig2 bw/NewKDM3B_shC_H3K9me2.InputSubtracted.pushToZero.bw --pseudocount 1 --operation log2 -o bw/logFC.NewKDM3B.sh3A3BovershC.H3K9me2.ISPTZ.bw

# Generate subtracted bigwigs in addition to Log2FC
echo "Subtract bigwigs"
 bigwigCompare --bigwig1 bw/293T_150DFO_H3K9me2_12hr_InputSubtracted.pushToZero.bw --bigwig2 bw/293T_Control_H3K9me2_12hr_InputSubtracted.pushToZero.bw --pseudocount 1 --operation subtract -o bw/Subtract.150DFOminControl.H3K9me2.ISPTZ.bw
 bigwigCompare --bigwig1 bw/NewKDM3B_sh3A3B_H3K9me2.InputSubtracted.pushToZero.bw --bigwig2 bw/NewKDM3B_shC_H3K9me2.InputSubtracted.pushToZero.bw --pseudocount 1 --operation subtract -o bw/Subtract.NewKDM3B.sh3A3BminshC.H3K9me2.ISPTZ.bw
 
# Next, for the heatmaps, set up a gene list, and then plot the heatmaps.
 
#Generating genes with POLII Peaks
# In order to use "canonical" isoforms for each gene, we downloaded knownGene.gtf and knownCanonical.txt for hg19 from  http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/

# Retrieve the strand information from the hg19 known gene annotation and added it to the canonical genes.
awk -F "\t" '{print $7,$9}' bed/hg19.knownGene.gtf | awk -F "\"" '{print $1,$4}' | perl -pe "s/ gene_id /\t/g" | sort | uniq | awk '{print $2,$1}' | perl -pe "s/\ /\t/g" | sort  > bed/hg19.knownGene.justStrand.gtf
sort -k 5 bed/knownCanonical.txt  > bed/knownCanonical.geneSorted.txt
join -1 5 -2 1 bed/knownCanonical.geneSorted.txt bed/hg19.knownGene.justStrand.gtf | awk '{printf "%s\t%d\t%d\t%s\t%d\t%s\n",$2,$3,$4,$1,1,$7}' > bed/knownCanonical.withStrand.bed

# Identify all canonical transcripts that overlap a Pol II peak, and then remove any duplicates.
bedtools intersect -wa -a bed/knownCanonical.withStrand.bed -b bed/293T_*POLII*macsPeaks.bed | sort | uniq > bed/knownCanonical.genesWithPOLIIpeaks.bed

#Generating Heatmaps comparing  H3K9me2 KDM3A and KDM3B plusMin 5kb around POLII occupied genes --> Extended data figure 8A
# KDM3A data ( ENCFF387ROQ_KDM3A.bigWig ) downloaded from ENCODE
echo "Create heatmaps with H3K9me2 plus publicly available KDM3A and KDM3B"
computeMatrix reference-point --referencePoint TSS -S bw/logFC.150DFOoverControl.H3K9me2.ISPTZ.bw bw/ENCFF387ROQ_KDM3A.bigWig New_KDM3B/New_KDM3B/tracks/KDM3B_shC.bw -R bed/knownCanonical.genesWithPOLIIpeaks.bed -a 5000 -b 5000 -p 24 -out heatmaps/ISPTZ.allPOLII.boundgenes.Canonical.H3K9me2andKDM3Anew3Bsignal.5000.5000.TSSref.matrix.txt.gz

plotHeatmap --matrixFile heatmaps/ISPTZ.allPOLII.boundgenes.Canonical.H3K9me2andKDM3Anew3Bsignal.5000.5000.TSSref.matrix.txt.gz --outFileSortedRegions bed/ISPTZ.allPOLII.boundgenes.Canonical.H3K9me2andKDM3Anew3Bsignal.5000.5000.TSSref.heatmap.km3.bed --kmeans 3 -out heatmaps/ISPTZ.allPOLII.boundgenes.Canonical.H3K9me2andKDM3Anew3Bsignal.5000.5000.TSSref.heatmap.km3.pdf --samplesLabel "logFC_H3K9me2" "KDM3A" "KDM3B" --zMin -2 0 0 --zMax 2 20 1.5 --yMin -0 0 0 --yMax 0.4 8 0.4 --colorList 'blue,white,red' 'white,purple' 'white,green'


# Next generate correlation plot comparing the change in H3K9me2 after DFO treatment to knockdown of KDM3A and KDM3B --> Extended data figure 8C

# First make a bedfile with coordinates for 1 kb on either side of the TSS for each of the canonical genes that overlap PolII
echo "Pull out TSSs for canonical genes that overlap polII peaks"
awk '$6 == "+" {printf "%s\t%s\t%s\t%s\t1\t%s\n",$1,$2-1,$2,$4,$6}'  bed/knownCanonical.genesWithPOLIIpeaks.bed > bed/knownCanonical.genesWithPolIIpeaks.posStrand.1bp.upstream.bed
awk '$6 == "-" {printf "%s\t%s\t%s\t%s\t1\t%s\n",$1,$3,$3+1,$4,$6}'  bed/knownCanonical.genesWithPOLIIpeaks.bed > bed/knownCanonical.genesWithPolIIpeaks.negStrand.1bp.upstream.bed
cat bed/knownCanonical.genesWithPolIIpeaks.posStrand.1bp.upstream.bed bed/knownCanonical.genesWithPolIIpeaks.negStrand.1bp.upstream.bed |  sort -k 1,1 -k2,2n  > bed/knownCanonical.genesWithPOLIIpeaks.1bp.upstream.bed
bedtools slop -i bed/knownCanonical.genesWithPOLIIpeaks.1bp.upstream.bed -g hg19.chrom.sizes -l 999 -r 1000 > bed/knownCanonical.genesWithPOLIIpeaks.plusMin1kb.bed

# Then calculate the correlation within these regions for each sample.
echo "Make correlation table for subtracted files"
multiBigwigSummary BED-file \
		   --bwfiles \
		   bw/Subtract.150DFOminControl.H3K9me2.ISPTZ.bw \
		   bw/Subtract.NewKDM3B.sh3A3BminshC.H3K9me2.ISPTZ.bw \
       		   bw/Input_sh3A3B.bw \
 		   bw/Input_shC.bw \
		   bw/293T_Control_Input_12hr.bw \
		   bw/293T_150DFO_Input_12hr.bw \
		   --BED bed/knownCanonical.genesWithPOLIIpeaks.plusMin1kb.bed \
		   -p 24 -out heatmaps/InputandSubtract_Samples.plusMin1kb.bed.npz

echo "Make correlation table for log2FC files"
multiBigwigSummary BED-file \
		   --bwfiles \
		   bw/logFC.150DFOminControl.H3K9me2.ISPTZ.bw \
		   bw/logFC.NewKDM3B.sh3A3BminshC.H3K9me2.ISPTZ.bw \
       		   bw/Input_sh3A3B.bw \
 		   bw/Input_shC.bw \
		   bw/293T_Control_Input_12hr.bw \
		   bw/293T_150DFO_Input_12hr.bw \
		   --BED bed/knownCanonical.genesWithPOLIIpeaks.plusMin1kb.bed \
		   -p 24 -out heatmaps/InputandLog2FC_Samples.plusMin1kb.bed.npz

# and Plot
echo "Make plots from tables"
plotCorrelation --corData heatmaps/InputandLog2FC_Samples.plusMin1kb.bed.npz \
		--corMethod pearson \
		--whatToPlot heatmap \
		--removeOutliers --plotNumbers \
		-o heatmaps/InputandLog2FC_Samples.plusMin1kb.bed.heatmap.png

plotCorrelation --corData heatmaps/InputandSubtract_Samples.plusMin1kb.bed.npz \
		--corMethod pearson \
		--whatToPlot heatmap \
		--removeOutliers --plotNumbers \
		-o heatmaps/InputandSubtract_Samples.plusMin1kb.bed.heatmap.png


# Author: Patrick Ozark
#
# Original perl code for pause index written by Patrick Ozark and published in this paper:  https://pubmed.ncbi.nlm.nih.gov/28860207/
# Modified by Jason Shapiro and Elizabeth Bartom for iron chelation / MTOR project with Jason Shapiro, and included on Github 

module load bedtools
export PATH=$PATH:/projects/b1025/tools/MACS-1.4.2/bin
export PYTHONPATH=/projects/b1025/tools/MACS-1.4.2/lib/python2.6/site-packages:$PYTHONPATH
module load ngsplot
module load R/3.3.3
module load python/anaconda3.6


# Get full length transcripts for all known canonical genes with Pol II peaks, swapping uc knownCanonical IDs to ENST ids and removing
# transcripts shorter than 2000 bp.
filePrefix=final_transcripts_full_length
echo $filePrefix
# knownToEnsembl.txt and ensGene.txt files can be found here:  http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/
echo "Get the right genomic regions"
perl getFullLengthTranscripts.v4.pl bed/knownCanonical.genesWithPOLIIpeaks.bed knownToEnsembl.txt ensGene.txt >  bed/$filePrefix.all.bed

#Sort the bed file
sort -k 1,1 -k2,2n bed/$filePrefix.all.bed > bed/$filePrefix.sorted.bed

# Pull out the coordinates of the TSS's, first separating by strand.
awk -F "\t" '$6=="+" {printf "%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2-1,$2,$4,$5,$6}' bed/$filePrefix.sorted.bed > bed/$filePrefix.sorted.TSS.pos.bed 
awk -F "\t" '$6=="-" {printf "%s\t%s\t%s\t%s\t%s\t%s\n",$1,$3-1,$3,$4,$5,$6}' bed/$filePrefix.sorted.bed > bed/$filePrefix.sorted.TSS.neg.bed
cat bed/$filePrefix.sorted.TSS.pos.bed bed/$filePrefix.sorted.TSS.neg.bed > bed/$filePrefix.sorted.TSS.bed

# Remove overlapping genes and genes within 1 kb of each other (based on just TSS locations)
echo "Clean up the list; remove overlapping and duplicate transcripts"
perl removeOverlappingGenes.pl bed/$filePrefix.sorted.bed bed/$filePrefix.sorted.TSS.bed > bed/$filePrefix.nonOverlapping.bed

# Some EnsG ids have more than one "canonical transcript"; remove those
# (This appears to be the consequence of one:many and many:one relationships between the UC ids and the ENSG ids.)
awk -F "\t" '{print $4}' bed/$filePrefix.nonOverlapping.bed  | sort | uniq -c | sort -n | awk '$1 > 1 && $2 ne "" {print $2}' > bed/$filePrefix.multiTxGenes.txt
wc bed/$filePrefix.nonOverlapping.bed
wc bed/$filePrefix.multiTxGenes.txt
# These are 75 genes that have multiple transcripts for the same gene, despite all of the above screening.  Discarding them.
grep -v -f bed/$filePrefix.multiTxGenes.txt bed/$filePrefix.nonOverlapping.bed > bed/$filePrefix.nonOverlapping.uniq.bed

# Get final transcript coordinates in bed format for promoter (-100 to +300) and gene body (+300 to +2000) - output = prefix.promoter.bed and prefix.genebody.bed
echo "Get the coordinates for the promoter region vs gene body"
perl getPromoterGeneBody.v3.pl bed/$filePrefix.nonOverlapping.uniq.bed  

# Get final transcript coordinates in bed format for full gene body (+300 to TES/TTS) - output = prefix.fullgenebody.bed
perl getFullGeneBody.v3.pl bed/$filePrefix.nonOverlapping.uniq.bed

# Filter chrM from bedGraph files
grep -vwE "chrM" 293T_150DFO_PolII.ISPTZ.bdg > 293T_150DFO_PolII.ISPTZ.filtered.bdg
grep -vwE "chrM" 293T_Control_PolII.ISPTZ.bdg > 293T_Control_PolII.ISPTZ.filtered.bdg

# Sort the bed files to match the bedgraphs.
echo "Sorting bed files."
sort -k 1,1 -k2,2n $filePrefix.nonOverlapping.uniq.promoter.bed > $filePrefix.nonOverlapping.uniq.promoter.sorted.bed
sort -k 1,1 -k2,2n $filePrefix.nonOverlapping.uniq.genebody.bed > $filePrefix.nonOverlapping.uniq.genebody.sorted.bed
sort -k 1,1 -k2,2n $filePrefix.nonOverlapping.uniq.fullgenebody.bed > $filePrefix.nonOverlapping.uniq.fullgenebody.sorted.bed

# Get sum of normalized coverage for promoter, gene body, and full gene body for final transcripts
echo "Get normalized read counts for each class of gene region."
bedtools map -a $filePrefix.nonOverlapping.uniq.promoter.sorted.bed -b 293T_Control_PolII.ISPTZ.filtered.bdg -c 4 -o sum > norm.promoter.Control.counts.nonOverlapping.uniq.txt &
bedtools map -a $filePrefix.nonOverlapping.uniq.promoter.sorted.bed -b 293T_150DFO_PolII.ISPTZ.filtered.bdg -c 4 -o sum > norm.promoter.150DFO.counts.nonOverlapping.uniq.txt &

bedtools map -a $filePrefix.nonOverlapping.uniq.genebody.sorted.bed -b 293T_Control_PolII.ISPTZ.filtered.bdg -c 4 -o sum > norm.genebody.Control.counts.nonOverlapping.uniq.txt &
bedtools map -a $filePrefix.nonOverlapping.uniq.genebody.sorted.bed -b 293T_150DFO_PolII.ISPTZ.filtered.bdg -c 4 -o sum > norm.genebody.150DFO.counts.nonOverlapping.uniq.txt &

bedtools map -a $filePrefix.nonOverlapping.uniq.fullgenebody.sorted.bed -b 293T_Control_PolII.ISPTZ.filtered.bdg -c 4 -o sum > norm.fullgenebody.Control.counts.nonOverlapping.uniq.txt &
bedtools map -a $filePrefix.nonOverlapping.uniq.fullgenebody.sorted.bed -b 293T_150DFO_PolII.ISPTZ.filtered.bdg -c 4 -o sum > norm.fullgenebody.150DFO.counts.nonOverlapping.uniq.txt &

wait
wait
wait
wait
wait
wait

echo "Combine the data into one larger table"
for s in "Control" "150DFO"
do
    # Combine promoter and genebody in a single file
    join --nocheck-order -j 4 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.7 norm.promoter.$s.counts.nonOverlapping.uniq.txt  norm.genebody.$s.counts.nonOverlapping.uniq.txt | perl -pe "s/ /\t/g" > norm.promoter.genebody.$s.counts.nonOverlapping.uniq.txt
    join --nocheck-order -j 4 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.7 norm.promoter.genebody.$s.counts.nonOverlapping.uniq.txt  norm.fullgenebody.$s.counts.nonOverlapping.uniq.txt | perl -pe "s/ /\t/g" > norm.$s.counts.nonOverlapping.uniq.txt
done

#awk -F "\t" '{printf "%s\t%s\t%s\t%s.%s\t%f\t%f\t%f\t%f\n",$1,$2,$3,$4,$5,$6,$7,$8,($6+1)/($7+1)}'  norm.Control.counts.nonOverlapping.uniq.txt >  norm.Control.uniqcounts.plusRatio.txt

# These bed files can be loaded into UCSC genome browser to check the regions and scores
echo "Generate bed files for quality control"
for g in "promoter" "genebody" "fullgenebody"
do
    for s in "Control" "150DFO"
    do
	awk -F "\t" '{printf "%s\t%s\t%s\t%s.%s\t%d\t%s\n",$1,$2,$3,$4,$5,$7,$6}'  norm.$g.$s.counts.nonOverlapping.uniq.txt >  norm.$g.$s.counts.nonOverlapping.uniq.bed
    done
done

# Combine the counts for both conditions in one file.
join -j 5 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.7,2.8 norm.promoter.genebody.Control.counts.nonOverlapping.uniq.txt norm.promoter.genebody.150DFO.counts.nonOverlapping.uniq.txt | perl -pe "s/ /\t/g" > norm.promoter.genebody.Control.150DFO.counts.nonOverlapping.uniq.txt

# Calculate the log2FC ratio between DFO and Control for both promoters and genebody (with a pseudocount of 1)
echo "Calculating pausing statistics"
awk -F "\t" '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,log(($9+1)/($7+1))/log(2),log(($10+1)/($8+1))/log(2)}' norm.promoter.genebody.Control.150DFO.counts.nonOverlapping.uniq.txt > temp.txt
echo '#Chr Start Stop ENSGid ENSTid Strand PromoterControl GeneBodyControl Promoter150DFO GeneBody150DFO log2FC.Promoter log2FC.GeneBody' | perl -pe "s/ /\t/g" | cat  > PolR2A.Occupancy.GeneList.log2FC.txt
cat temp.txt >> PolR2A.Occupancy.GeneList.log2FC.txt

echo "These commands were originally done in Excel; now moved to shell"
# Sort the list by log2FC.GeneBody, and pull out the lines with values > 1.25*stdev(log2FC.GeneBody)+ mean(log2FC.GeneBody) = "0.776257", to be the "up" set.
sort -n -k 12 PolR2A.Occupancy.GeneList.log2FC.txt | grep -v ^# | awk -F "\t" '$12 > 0.776257 {print $4}' > PolR2A.Occupancy.GeneList.log2FC.up.txt

# Calculate Pausing index for both control and DFO.
awk -F "\t" '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,($7+1)/($8+1),($9+1)/($10+1)}' norm.promoter.genebody.Control.150DFO.counts.nonOverlapping.uniq.txt > temp.txt
echo '#Chr Start Stop ENSGid ENSTid Strand PromoterControl GeneBodyControl Promoter150DFO GeneBody150DFO PauseIndex.Control PauseIndex.150DFO' | perl -pe "s/ /\t/g" | cat > PolR2A.Occupancy.PausingIndex.txt
cat temp.txt >> PolR2A.Occupancy.PausingIndex.txt

# Calculate the Log2FC ratio of the pausing indices.
grep -v ^# PolR2A.Occupancy.PausingIndex.txt | awk -F "\t" '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$12/$11,log($12/$11)/log(2)}' | sort -nr -k 14 > temp.txt
echo '#Chr Start Stop ENSGid ENSTid Strand PromoterControl GeneBodyControl Promoter150DFO GeneBody150DFO PauseIndex.Control PauseIndex.150DFO PauseIndexRatio Log2FC.PauseIndex' | perl -pe "s/ /\t/g" | cat > PolR2A.Occupancy.PausingIndex.log2FC.txt
cat temp.txt >> PolR2A.Occupancy.PausingIndex.log2FC.txt

# Pull out the lines with values > 1.25*stdev(log2FC.PauseIndex)+mean(log2FC.PauseIndex) = "0.726543339", removing any genes that were in the up set, to give the paused set.
sort -n -k 14 PolR2A.Occupancy.PausingIndex.log2FC.txt | grep -v ^# | awk -F "\t" '$14 > 0.726543339 {print $4}' | grep -v -f PolR2A.Occupancy.GeneList.log2FC.up.txt > PolR2A.Occupancy.GeneList.log2FC.paused.txt

# Pull out the genes that are down in both promoter < mean(log2FC.Promoter) and then by the gene body < mean(log2FC.GeneBody) and then remove paused genes to get the "down" set.
sort -n -k 11 PolR2A.Occupancy.GeneList.log2FC.txt | grep -v ^# | awk -F "\t" '$11 < -0.4405576 && $12 < -0.24861 {print $4}' | grep -v -f PolR2A.Occupancy.GeneList.log2FC.paused.txt > PolR2A.Occupancy.GeneList.log2FC.down.txt

# Figure 6K was generated from the above data in Prism.

# Next make profile plot for Extended Data Figure 5A
# First make bed files for each group
echo "Generating bed files for different classes of affected genes"
grep -f PolR2A.Occupancy.GeneList.log2FC.up.txt bed/final_transcripts_full_length.sorted.TSS.bed | awk '{printf "%s\t%s\t%s\t%s.%s\t%d\t%s\n",$1,$2,$3,$4,$5,1000,$6}' > PolR2A.Occupancy.GeneList.log2FC.up.TSS.bed
grep -f PolR2A.Occupancy.GeneList.log2FC.paused.txt bed/final_transcripts_full_length.sorted.TSS.bed | awk '{printf "%s\t%s\t%s\t%s.%s\t%d\t%s\n",$1,$2,$3,$4,$5,1000,$6}' > PolR2A.Occupancy.GeneList.log2FC.paused.TSS.bed
grep -f PolR2A.Occupancy.GeneList.log2FC.down.txt bed/final_transcripts_full_length.sorted.TSS.bed | awk '{printf "%s\t%s\t%s\t%s.%s\t%d\t%s\n",$1,$2,$3,$4,$5,1000,$6}' > PolR2A.Occupancy.GeneList.log2FC.down.TSS.bed

# Then plot them.
echo "Plotting occupancy"
computeMatrix reference-point -S bw/logFC.150DFOoverControl.POLII.ISPTZ.bw -R PolR2A.Occupancy.GeneList.log2FC.up.TSS.bed PolR2A.Occupancy.GeneList.log2FC.paused.TSS.bed PolR2A.Occupancy.GeneList.log2FC.down.TSS.bed -a 2000 -b 100 -o PolR2a.Occupancy.allGeneLists.matrix.gz
plotProfile -m PolR2a.Occupancy.allGeneLists.matrix.gz --regionsLabel "Increased" "Paused" "Decreased" -o PolR2a.Occupancy.allGeneLists.plot.pdf

# To Run GSEA we need common gene names not ENST values
echo "Extracting mapping between uc ID and cgn and ENSG and ENST to convert to a list of ranked CGNs"
# Download the conversion from UCSC using this command:
mysql   --user=genome   -N   --host=genome-mysql.cse.ucsc.edu   -A   -D hg19   -e "select ensGene.name, name2, value from ensGene, ensemblToGeneName where ensGene.name = ensemblToGeneName.name" >   queryOutput.txt
awk '{print $2,$3}' queryOutput.txt | perl -pe "s/ /\t/g"  | sort | uniq > EnsgToCGN.txt

# Sort the data table by ENST value.
sort -k 4  PolR2A.Occupancy.GeneList.log2FC.txt > PolR2A.Occupancy.GeneList.log2FC.ensgSorted.txt

# # Combine the data table with the gene names and sort by absolute value of log2FC in the gene body to create a preranked file for GSEA.
join -2 1 -1 4 PolR2A.Occupancy.GeneList.log2FC.ensgSorted.txt EnsgToCGN.txt | sort | uniq | awk '{printf "%s\t%s\n",$13,sqrt($12*$12)}' | sort -k 2 -nr > logFC.DFOvsCtrl.rankedGeneBodyCounts.rnk

# # Download GSEA to your laptop and load the .rnk file.
# # Select run GSEA preranked and select the .rnk file as your input.
# # Set "Collapse/Remap" to "No_Collapse"
# # Run the analysis!
