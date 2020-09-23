# pwd
cd /home/rad/users/gaurav/projects/workflows/nfatacseq

# Parameters for the script
species="mouse"
user="miguel"
projName="TMColCancer"
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

# From the cutoff analysis, we identified a cutoff of pvalue=0.002
pvalueDir="0_002" ; # pvalue=0.002 or pscore ~ 2.70
pvalue="${pvalueDir//_/.}" 
customPeaksDir="${analysisDir}/customPeaks/mergedLibrary"; mkdir -p ${customPeaksDir}
peaksOutDir=${customPeaksDir}/p${pvalueDir}
for b in ${bamDir}/*.bam;
do
  bname=$(basename ${b} .mLb.clN.sorted.bam)
  macs2 callpeak -f BAM --seed=39751 -g mm  --keep-dup all  --nolambda --broad -p ${pvalue} --broad-cutoff  ${pvalue} --outdir ${peaksOutDir} -t ${b} -n ${bname}
done

# 2.2) Generate consensus peaks region file and raw count matrix
consensusDir="${peaksOutDir}/consensus"
echo -e "- ${projDir}\n- ${peaksOutDir}\n- ${consensusDir}"
bash /home/rad/users/gaurav/projects/workflows/nfchipseq/scripts/generate_rawCount_mergedPeak_files.sh ${projDir} broadPeak ${peaksOutDir} ${consensusDir}

# 2.3) 
#########################################################################################
# 3) DOWNSTREAM ANALYSIS
#########################################################################################
# 3.1) Output paramerts for downstream analysis
projDir="${outdir}/${user}/${projName}"
bamDir=${projDir}/results/bwa/mergedLibrary
analysisDir="${projDir}/analysis"; mkdir -p ${analysisDir}
pvalueDir="0_002" ; # pvalue=0.002 or pscore ~ 2.70
pvalue="${pvalueDir//_/.}" 
customPeaksDir="${analysisDir}/customPeaks/mergedLibrary"; mkdir -p ${customPeaksDir}
peaksOutDir=${customPeaksDir}/p${pvalueDir}
consensusDir="${peaksOutDir}/consensus"

bname="${projName}_consensus_peaks"
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
