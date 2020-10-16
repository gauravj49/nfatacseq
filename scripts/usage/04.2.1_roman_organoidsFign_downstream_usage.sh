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

# 3.2.6.1) Get the raw counts file for the DARs
python -c "import sys; import pandas as pd; import datatable as dt; input_file=sys.argv[1]; input_bed=sys.argv[2]; outtxt_file=sys.argv[3]; infileDT = dt.fread(input_file, sep='\t', header=True, nthreads=16); bedDT = dt.fread(input_bed, sep='\t', header=True, nthreads=16); infileDT = dt.fread(input_file, sep='\t', header=True, nthreads=16); bedDT = dt.fread(input_bed, sep='\t', header=True, nthreads=16); infileDF = infileDT.to_pandas(); infileDF.set_index('PeakID', inplace=True); bedDF = bedDT.to_pandas(); bedDF.set_index('feature', inplace=True); infileDF = infileDF[infileDF.index.isin(bedDF.index)]; infileDF.reset_index(level=0, inplace=True); infileDFColumns = infileDF.columns.tolist(); infileDFColumns.insert(3,infileDFColumns.pop(0)); infileDF = infileDF[infileDFColumns]; infileDF.to_csv(outtxt_file, header=True, index=False, sep='\t', float_format='%.0f');" /media/rad/HDD1/atacseq/roman/organoidsFign/analysis/customPeaks/mergedLibrary/q0_01/consensus/organoidsFign_consensus_peaks_rawCounts.txt /media/rad/HDD1/atacseq/roman/organoidsFign/analysis/customPeaks/mergedLibrary/q0_01/consensus/DARs/heatmap/FignKrasExpressed_over_KrasExpressed_significant_filtered_DARs.txt /media/rad/HDD1/atacseq/roman/organoidsFign/analysis/customPeaks/mergedLibrary/q0_01/consensus/DARs/heatmap/FignKrasExpressed_over_KrasExpressed_significant_filtered_rawCounts.txt

# 3.2.6.2) Get the annotation file for the DARs
python -c "import sys; import pandas as pd; import datatable as dt; input_file=sys.argv[1]; input_bed=sys.argv[2]; outtxt_file=sys.argv[3]; infileDT = dt.fread(input_file, sep='\t', header=True, nthreads=16); bedDT = dt.fread(input_bed, sep='\t', header=True, nthreads=16); infileDT = dt.fread(input_file, sep='\t', header=True, nthreads=16); bedDT = dt.fread(input_bed, sep='\t', header=True, nthreads=16); infileDF = infileDT.to_pandas(); infileDF.set_index('PeakID', inplace=True); bedDF = bedDT.to_pandas(); bedDF.set_index('feature', inplace=True); infileDF = infileDF[infileDF.index.isin(bedDF.index)]; infileDF.reset_index(level=0, inplace=True); infileDFColumns = infileDF.columns.tolist(); infileDFColumns.insert(3,infileDFColumns.pop(0)); infileDF = infileDF[infileDFColumns]; infileDF.to_csv(outtxt_file, header=True, index=False, sep='\t', float_format='%.0f');" /media/rad/HDD1/atacseq/roman/organoidsFign/analysis/customPeaks/mergedLibrary/q0_01/consensus/organoidsFign_consensus_peaks_annotation.txt /media/rad/HDD1/atacseq/roman/organoidsFign/analysis/customPeaks/mergedLibrary/q0_01/consensus/DARs/heatmap/FignKrasExpressed_over_KrasExpressed_significant_filtered_DARs.txt /media/rad/HDD1/atacseq/roman/organoidsFign/analysis/customPeaks/mergedLibrary/q0_01/consensus/DARs/heatmap/FignKrasExpressed_over_KrasExpressed_significant_filtered_annotation.txt

# --------------------------------------------------------

# 4.1) Integration of Chipseq and ATACseq samples

# 4.1.1) Copy the relevant peaks files in the subfolders
customPeaksDir="${consensusDir}/customConsensus"; mkdir -p ${customPeaksDir}
mkdir -p ${customPeaksDir}/{PK,PKF}/peaks
cp -rv ${peaksOutDir}/*{V5,Pdk}*.broadPeak ${customPeaksDir}/PK/peaks
cp -rv ${peaksOutDir}/*SW*.broadPeak ${customPeaksDir}/PKF/peaks

for pd in PK PKF;
do 
  customSamplesPeaksDir="${customPeaksDir}/${pd}"
  peaksDir=${customSamplesPeaksDir}/peaks
  echo -e "- ${projDir}\n- ${peaksDir}\n- ${customSamplesPeaksDir}"
  bash /home/rad/users/gaurav/projects/workflows/nfchipseq/scripts/generate_rawCount_mergedPeak_files.sh ${projDir} broadPeak ${peaksDir} ${customSamplesPeaksDir}

  origAnnFile="${customSamplesPeaksDir}/interimFiles/${projName}_consensus_peaks.mLb.clN.boolean.txt"
  scriptDir="/home/rad/users/gaurav/projects/workflows/nfatacseq"
  bash /home/rad/users/gaurav/projects/workflows/nfatacseq/scripts/parse_nfatacseq_consensus_peaks_annotation.sh ${species} ${user} ${projName} ${customSamplesPeaksDir} ${bamDir} ${origAnnFile} ${scriptDir}
  echo ""
done

# Intersect with the chip data
pkconsensusbed="/media/rad/HDD1/atacseq/roman/organoidsFign/analysis/customPeaks/mergedLibrary/q0_01/consensus/customConsensus/PK/organoidsFign_consensus_peaks.bed"
pkfconsensusbed="/media/rad/HDD1/atacseq/roman/organoidsFign/analysis/customPeaks/mergedLibrary/q0_01/consensus/customConsensus/PKF/organoidsFign_consensus_peaks.bed"

for chipbed in $(ls ${customPeaksDir}/*_ChipSeq.txt);
do
  cbname=$(basename ${chipbed} _ChipSeq.txt)
  intersectBed -a ${chipbed} -b ${pkconsensusbed}  -wo -f 0.5 > ${customPeaksDir}/50pc_covered_${cbname}_PK_organoidsFign_consensus_regions.txt
  intersectBed -a ${chipbed} -b ${pkfconsensusbed} -wo -f 0.5 > ${customPeaksDir}/50pc_covered_${cbname}_PKF_organoidsFign_consensus_regions.txt
  
  intersectBed -a ${chipbed} -b ${pkconsensusbed}  -wo        > ${customPeaksDir}/atleast_1bp_overlap_${cbname}_PK_organoidsFign_consensus_regions.txt
  intersectBed -a ${chipbed} -b ${pkfconsensusbed} -wo        > ${customPeaksDir}/atleast_1bp_overlap_${cbname}_PKF_organoidsFign_consensus_regions.txt
done

