# pwd
cd /home/rad/users/gaurav/projects/workflows/nfatacseq

# Parameters for the script
species="mouse"
user="roman"
projName="organoidsFign"
outdir="/media/rad/HDD1/atacseq"
jobdir="/home/rad/users/gaurav/projects/workflows/nfatacseq"
projDir="${outdir}/${user}/${projName}"

#########################################################################################
# 2) PERFORM CUSTOM PEAKCALLING
#########################################################################################
# 2.1) Perform cutoff analysis to get the threshold for the peaks
# Output paramerts for downstream analysis
projDir="${outdir}/${user}/${projName}"
bamDir=${projDir}/results/bwa/mergedLibrary
analysisDir="${projDir}/analysis"; mkdir -p ${analysisDir}
qvalueDir="0_1" # qvalue = 0.1
qvalue="${qvalueDir//_/.}" 
customPeaksDir="${analysisDir}/customPeaks/mergedLibrary"; mkdir -p ${customPeaksDir}
peaksOutDir=${customPeaksDir}/q${qvalueDir}

# Custom peakcalling with calculated threshold using the --cutoff-analysis parameter
# The output columns are: 
  # 1) p-score -log10(pvalue) cutoff, 
  # 2) q-score -log10(qvalue) cutoff, 
  # 3) number of peaks called, 
  # 4) total length in basepairs of peaks called, 
  # 5) average length of peak in basepair
# --nolambda in the absence of a control is only valid for broad peaks, so diffuse histone modification ChIPs. If you apply it on sharp peak ChIPs like transcription factors, you deactivate the local background estimation that aims to decide if a peak is true or false based on the signal landscape surrounding the putative peak.

for b in ${bamDir}/*.bam;
do
  bname=$(basename ${b} .mLb.clN.sorted.bam)
  /home/rad/miniconda/bin/macs2 callpeak -f BAM --seed=39751 -g mm  --keep-dup all  --cutoff-analysis --nolambda --broad --broad-cutoff  ${qvalue} --outdir ${peaksOutDir} -t ${b} -n ${bname}
done

# From the cutoff analysis, we identified a cutoff of qvalue=0.01
qvalueDir="0_01" ; # qvalue=0.01 or qscore ~ 1.96
qvalue="${qvalueDir//_/.}" 
customPeaksDir="${analysisDir}/customPeaks/mergedLibrary"; mkdir -p ${customPeaksDir}
peaksOutDir=${customPeaksDir}/q${qvalueDir}
for b in ${bamDir}/*.bam;
do
  bname=$(basename ${b} .mLb.clN.sorted.bam)
  macs2 callpeak -f BAM --seed=39751 -g mm  --keep-dup all  --nolambda --broad -q ${qvalue} --broad-cutoff  ${qvalue} --outdir ${peaksOutDir} -t ${b} -n ${bname}
done

# 2.2) Generate consensus peaks region file and raw count matrix
consensusDir="${peaksOutDir}/consensus"
echo -e "- ${projDir}\n- ${peaksOutDir}\n- ${consensusDir}"
bash /home/rad/users/gaurav/projects/workflows/nfchipseq/scripts/generate_rawCount_mergedPeak_files.sh ${projDir} broadPeak ${peaksOutDir} ${consensusDir}

#########################################################################################
# 3) DOWNSTREAM ANALYSIS
#########################################################################################
# 3.1) Output paramerts for downstream analysis
projDir="${outdir}/${user}/${projName}"
bamDir=${projDir}/results/bwa/mergedLibrary
analysisDir="${projDir}/analysis"; mkdir -p ${analysisDir}
qvalueDir="0_01" ; # qvalue=0.01 or qscore ~ 1.96
qvalue="${qvalueDir//_/.}" 
customPeaksDir="${analysisDir}/customPeaks/mergedLibrary"; mkdir -p ${customPeaksDir}
peaksOutDir=${customPeaksDir}/q${qvalueDir}
consensusDir="${peaksOutDir}/consensus"

bname="${projName}_consensus_peaks"
consensusPeaksBed="${consensusDir}/${bname}.bed"
rawCountsTxtFile="${consensusDir}/${bname}_rawCounts.txt"
peaksAnnTxtFile="${consensusDir}/${bname}_annotation.txt"
origAnnFile="${consensusDir}/interimFiles/${projName}_consensus_peaks.mLb.clN.boolean.txt"

# 3.2) Parse concensus raw matrix and boolean matrix to get annoation files
echo "bash scripts/parse_nfatacseq_consensus_peaks_annotation.sh ${species} ${user} ${projName} ${consensusDir} ${bamDir} ${origAnnFile} ${jobdir}"
bash scripts/parse_nfatacseq_consensus_peaks_annotation.sh ${species} ${user} ${projName} ${consensusDir} ${bamDir} ${origAnnFile} ${jobdir}

# Output files are:
echo "- Consensus bed file: ${consensusPeaksBed}"
echo "- Raw peaks count   : ${rawCountsTxtFile}"
echo "- Peaks annotation  : ${peaksAnnTxtFile}"

# 3.2) Identify differential accessibility regions (DARs) between different biological conditions
# 3.2.1) Get the counts matrix file
countsMatrixFile="${consensusDir}/${bname}_rawCounts.mat"
cut -f4- ${rawCountsTxtFile} > ${countsMatrixFile}

# 3.2.2) Get the count files
darDir="/media/rad/HDD1/atacseq/roman/organoidsFign/analysis/customPeaks/mergedLibrary/q0_01/consensus/DARs"
infile="/media/rad/HDD1/atacseq/roman/organoidsFign/analysis/customPeaks/mergedLibrary/q0_01/consensus/organoidsFign_consensus_peaks_rawCounts.mat"; dos2unix ${infile}
countsDIR="${darDir}/counts"; mkdir -p ${countsDIR}
ncols=`expr $(awk '{print NF; exit}' ${infile})`;
for ((i=2;i<=ncols;i++)); do
 bname=`head -1 ${infile}|cut -f${i}`
 ofname="${countsDIR}/${bname}_Counts.txt"
 echo "processing $bname..."
 cut -f1,${i} ${infile} --output-delimiter=$'\t' | grep -v ${bname} > ${ofname}
done

# 3.2.3) Get the control and treatment files
DEseq2IN="${darDir}/deseq2IN"; mkdir -p ${DEseq2IN}
# Get the DEseq2IN files
cd ${countsDIR}
ls *SW*       > ${DEseq2IN}/FignKrasExpressed.txt
ls *{V5,Pdk}* > ${DEseq2IN}/KrasExpressed.txt
cd - 

# 3.2.4) Get the DARs
Rscript scripts/R_RUVseq_DESeq2_DAR_analysis.R ${countsDIR} ${DEseq2IN}/KrasExpressed.txt ${DEseq2IN}/FignKrasExpressed.txt 'KrasExpressed,FignKrasExpressed' ${darDir}/FignKrasExpressed_over_KrasExpressed FALSE gene TRUE TRUE FALSE

# 3.2.5) Get heatmaps of the DARs
l2fc=0.25
bm=5
padj=0.05
python scripts/draw_features_heatmap.py -if=${darDir}/results_DEseq2/FignKrasExpressed_over_KrasExpressed_DE_RESULTS.txt -of=${darDir}/heatmap/ignKrasExpressed_over_KrasExpressed_DARs.txt -lf=${l2fc} -bm=${bm} -pj=${padj} -xc -yc -sm




############################################
# 3.2.6) Get the DARs
sampleAnnotation="${consensusDir}/${bname}_sample_annotation.txt"

R
# Load libraries
cat("- Loading libraries...\n")
suppressPackageStartupMessages(library("DESeq2", warn.conflicts=FALSE, quietly=TRUE))
suppressPackageStartupMessages(library("edgeR" , warn.conflicts=FALSE, quietly=TRUE))
suppressPackageStartupMessages(library("RUVSeq", warn.conflicts=FALSE, quietly=TRUE))
suppressPackageStartupMessages(library("vsn", warn.conflicts=FALSE, quietly=TRUE))
suppressPackageStartupMessages(library("RColorBrewer", warn.conflicts=FALSE, quietly=TRUE))
suppressPackageStartupMessages(library("BiocParallel", warn.conflicts=FALSE, quietly=TRUE))
register(MulticoreParam(4))

# Input and output file
countsMatrixFile  <- "/media/rad/HDD1/atacseq/roman/organoidsFign/analysis/customPeaks/mergedLibrary/q0_01/consensus/organoidsFign_consensus_peaks_rawCounts.mat"
outputBaseFileName<- "organoidsFign_consensus_peaks"
sampleAnnFile     <- "/media/rad/HDD1/atacseq/roman/organoidsFign/analysis/customPeaks/mergedLibrary/q0_01/consensus/organoidsFign_consensus_peaks_sample_annotation.txt"
outputdir         <- "/media/rad/HDD1/atacseq/roman/organoidsFign/analysis/customPeaks/mergedLibrary/q0_01/consensus/DARs"; system(paste0("mkdir -p ", outputdir))
plotsdir          <- paste(outputdir, "/plots", sep=''); system(paste("mkdir -p", plotsdir, sep=' '));

# Get basenames
basefilename      <- tools::file_path_sans_ext(basename(countsMatrixFile))
extfilename       <- tools::file_ext(countsMatrixFile)

# Get input data 
cat("\n- Getting input data\n")
sampleCountsDF    <- read.table(countsMatrixFile, row.names=1, sep="\t", header=TRUE, quote = "")
coldata           <- read.table(sampleAnnFile, row.names=1, sep="\t", header=FALSE, quote = "")
colnames(coldata) <- "condition"
coldata$condition <- factor(coldata$condition)

# Create RLE plots for matrix
cat("\n- Getting matrix of counts data from counts matrix file...\n")
dds    <- DESeqDataSetFromMatrix(countData = sampleCountsDF, colData = coldata, design = ~ condition)

# Prefiltering
keep <- rowSums(counts(dds)) >= 50
dds <- dds[keep,]

# DAR analysis
dds <- DESeq(dds)
res <- results(dds)
res

# Plot counts
plotCountsFile <- paste0(plotsdir, "/", outputBaseFileName, "_count_plots.png", sep='') 
png(filename=plotCountsFile, height=600, width=600, bg="white", res=100)
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
dev.off()

# Shrinkage estimators
resultsNames(dds)
# [1] "Intercept"                                   
# [2] "condition_KrasExpressed_vs_FignKrasExpressed"

# Because we are interested in treated vs untreated, we set 'coef=2'
resLFC     <- lfcShrink(dds, coef="condition_KrasExpressed_vs_FignKrasExpressed", type="apeglm")
resNorm    <- lfcShrink(dds, coef=2, type="normal")
resAsh     <- lfcShrink(dds, coef=2, type="ashr")
MAplotFile <- paste0(plotsdir, "/", outputBaseFileName, "_MA_plot.png", sep='') 
png(filename=MAplotFile, height=600, width=2000, bg="white", res=100)
par(mfrow=c(1,3), mar=c(4,4,2,1))
plotMA(resLFC,  main="apeglm")
plotMA(resNorm, main="normal")
plotMA(resAsh,  main="ashr")
dev.off()

# Extracting transformed values
ntd <- normTransform(dds)
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

meanSDplotFile <- paste0(plotsdir, "/", outputBaseFileName, "_sizefactor_meanSD_plot.png", sep='') 
png(filename=meanSDplotFile, height=600, width=600, bg="white", res=100)
par(mfrow=c(1,3), mar=c(4,4,2,1))
meanSdPlot(assay(ntd), main="SizeFactor")
dev.off()

meanSDplotFile <- paste0(plotsdir, "/", outputBaseFileName, "_VSD_meanSD_plot.png", sep='') 
png(filename=meanSDplotFile, height=600, width=600, bg="white", res=100)
par(mfrow=c(1,3), mar=c(4,4,2,1))
meanSdPlot(assay(vsd))
dev.off()

meanSDplotFile <- paste0(plotsdir, "/", outputBaseFileName, "_RLD_meanSD_plot.png", sep='') 
png(filename=meanSDplotFile, height=600, width=600, bg="white", res=100)
par(mfrow=c(1,3), mar=c(4,4,2,1))
meanSdPlot(assay(rld))
dev.off()

# Heatmap of the sample-to-sample distances
samCorplotFile <- paste0(plotsdir, "/", outputBaseFileName, "_sample_correlation_plot.png", sep='') 
png(filename=samCorplotFile, height=1200, width=1200, bg="white", res=100)
sampleDists                <- dist(t(assay(rld)))
sampleDistMatrix           <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, rownames(sampleDistMatrix), sep="-")
colnames(sampleDistMatrix) <- paste(vsd$condition, rownames(sampleDistMatrix), sep="-")
colors                     <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

PCAplotFile <- paste0(plotsdir, "/", outputBaseFileName, "_PCA_plot.png", sep='') 
png(filename=PCAplotFile, height=1200, width=1400, bg="white", res=100)
plotPCA(vsd, intgroup=c("condition"))
dev.off()