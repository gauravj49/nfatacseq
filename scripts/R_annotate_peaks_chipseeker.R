#! /usr/bin/env Rscript

# ****************************************************************
# GOAL  : Annotate ATACseq peaks
# USAGE : Rscript R_annotate_peaks_chipseeker.R <input_tab_file> <output_annotation_file>
# ****************************************************************

################ load libraries ###############
suppressPackageStartupMessages(library("argparse"))

## Main logic
main <- function(){    

	##### Parse command line arguments ############
	args        <- check_options()
	inputfile   <- args$inputfile
	outputfile  <- args$outputfile
	species     <- args$species
	macs2Peaks  <- args$macs2Peaks
	addChr      <- args$addChr

  cat("- Input command line arguments ...\n")
  cat(paste0("\t- inputfile   = ",inputfile,"\n"))
  cat(paste0("\t- outputfile  = ",outputfile,"\n"))
  cat(paste0("\t- species     = ",species,"\n"))
  cat(paste0("\t- macs2Peaks  = ",macs2Peaks,"\n"))
  cat(paste0("\t- addChr      = ",addChr,"\n"))
	################ Load libraries ###############
	load_libraries()

  ################ Species AnnDBs ###############
  # Special mention: https://support.bioconductor.org/p/125609/
  # hub <- AnnotationHub()

  if(species == 'mouse'){
    # query(hub, c("Mus Musculus","EnsDb"))
    # txdb <- hub[['AH78811']]
    txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene
    andb <- "org.Mm.eg.db"
  }else if (species == 'human') {
    # query(hub, c("Homo sapiens","EnsDb"))
    # txdb <- hub[['AH78783']]
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    andb <- "org.Hs.eg.db"
  }
  
  ################ Main logic ###################
  # Annotate each sample peaks with ChIPseeker
	## Get the input data
	cat("- Reading input file ...\n")
  if(macs2Peaks){
    origPeakDT <- fread(inputfile, header=FALSE, sep="\t")
    setnames(origPeakDT, c("PeakChrom", "PeakStart", "PeakEnd", "PeakID", "PeakScore","PeakStrand","PeakSignalValue","PeakPvalue","PeakQvalue","PeakPointSource"))
    if (addChr){
      origPeakDT[,PeakChrom := paste0('chr',PeakChrom)]
    }
    peaksBedDT <- origPeakDT[,1:4]

  }else{
    origPeakDT <- fread(inputfile, header=TRUE, sep="\t")
    peaksBedDT <- origPeakDT[,c("PeakChrom", "PeakStart", "PeakEnd", "PeakID")]
  }
  print(str(as.data.frame(peaksBedDT)))
  # peaksBedDT[,PeakChrom := paste0('chr',PeakChrom)]
  peaksBed   <- makeGRangesFromDataFrame(as.data.frame(peaksBedDT))

  # Annotate regions
  cat("\n\t- Annotate regions\n")
  peakAnno <- annotatePeak(peaksBed, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb=andb)

  # Convert to dataframe 
  peakAnnoDT <- as.data.table(peakAnno)
  # print(peakAnnoDT)
  # Rename the header
  # Original header: > colnames(as.data.table(peakAnno))
  # [1]  "seqnames"      "start"         "end"           "width"        
  # [5]  "strand"        "annotation"    "geneChr"       "geneStart"    
  # [9]  "geneEnd"       "geneLength"    "geneStrand"    "geneId"       
  # [13] "transcriptId"  "distanceToTSS" "ENTREZID"      "SYMBOL"       
  # [17] "GENENAME" 
  setnames(peakAnnoDT, c("PeakChrom", "PeakStart", "PeakEnd", "PeakWidth", "PeakStrand", "DetailedGenomicAnnotation", "GeneChrom", "GeneStart", "GeneEnd", "GeneLength", "GeneStrand", "GeneID", "TranscriptID", "DistanceToTSS", "EntrezID", "GeneName", "GeneDesc"))

  # Copy the DetailedGenomicAnnotation column as GenomicAnnotation column
  peakAnnoDT[,GenomicAnnotation:=DetailedGenomicAnnotation]

  # Replace the detailed annotation to the abstract annotation
  peakAnnoDT[DetailedGenomicAnnotation %like%   'exon 1 '   , GenomicAnnotation:='ExonFirst']
  peakAnnoDT[!(DetailedGenomicAnnotation %like% 'exon 1 ')  , GenomicAnnotation:='ExonOther']
  peakAnnoDT[DetailedGenomicAnnotation %like%   'intron 1 ' , GenomicAnnotation:='IntronFirst']
  peakAnnoDT[!(DetailedGenomicAnnotation %like% 'intron 1 '), GenomicAnnotation:='IntronOther']
  peakAnnoDT[DetailedGenomicAnnotation=='Distal Intergenic' , GenomicAnnotation:='IntergenicDistal']
  peakAnnoDT[DetailedGenomicAnnotation=="3' UTR"            , GenomicAnnotation:='ThreeUTR']
  peakAnnoDT[DetailedGenomicAnnotation=="5' UTR"            , GenomicAnnotation:='FiveUTR' ]
  peakAnnoDT[DetailedGenomicAnnotation=='Downstream (1-2kb)', GenomicAnnotation:='DownstreamProximal']
  peakAnnoDT[DetailedGenomicAnnotation=='Downstream (<1kb)' , GenomicAnnotation:='DownstreamBasal']
  peakAnnoDT[DetailedGenomicAnnotation=='Downstream (2-3kb)', GenomicAnnotation:='DownstreamDistal']
  peakAnnoDT[DetailedGenomicAnnotation=='Promoter (1-2kb)'  , GenomicAnnotation:='PromoterProximal']
  peakAnnoDT[DetailedGenomicAnnotation=='Promoter (<=1kb)'  , GenomicAnnotation:='PromoterBasal']
  peakAnnoDT[DetailedGenomicAnnotation=='Promoter (2-3kb)'  , GenomicAnnotation:='PromoterDistal']

  # Reorder the columns
  setcolorder(peakAnnoDT, c("PeakChrom", "PeakStart", "PeakEnd", "PeakWidth", "PeakStrand", "GenomicAnnotation", "DetailedGenomicAnnotation", "GeneChrom", "GeneStart", "GeneEnd", "GeneLength", "GeneStrand", "GeneID", "TranscriptID", "DistanceToTSS", "EntrezID", "GeneName", "GeneDesc"))

  # peakAnnoDT[,unique(GenomicAnnotation)]
  #  [1] "IntronFirst"        "IntronOther"        "PromoterBasal"     
  #  [4] "ThreeUTR"           "IntergenicDistal"   "PromoterDistal"    
  #  [7] "PromoterProximal"   "DownstreamProximal" "DownstreamBasal"   
  # [10] "DownstreamDistal"   "FiveUTR"

  # Merge the orginal data table with annotation data table
  # There are few entires in the annotation data table that ...
  # ... were not present but the output datatable should be of same size as input
  # Create Temporary ID for merging
  origPeakDT[,mergeID:=paste0(PeakChrom,PeakStart,PeakEnd)] 
  peakAnnoDT[,mergeID:=paste0(PeakChrom,PeakStart,PeakEnd)] 

  if(macs2Peaks){
    origPeakDT <- origPeakDT[,c("PeakID", "PeakScore","PeakSignalValue","PeakPvalue","PeakQvalue","PeakPointSource","mergeID")]
    # Merge on mergeID
    mergedPeaksDT <- merge(peakAnnoDT, origPeakDT,all.y=T, by="mergeID")
  }else{
    # Merge on common ids and keep all the entires from the first data table
    mergedPeaksDT <- merge(peakAnnoDT, origPeakDT,all.y=T)
  }

  # Remove additional columns
  mergedPeaksDT[, c('mergeID') :=NULL]

  # Move PeakID column to the 4th position
  setcolorder(mergedPeaksDT, c("PeakChrom", "PeakStart", "PeakEnd", "PeakID"))

  # Save results in the output file
  fwrite(mergedPeaksDT, outputfile, sep = "\t", quote=F)

  # Sort alphanumerically with PeakID
  system(paste0("sort -k1,1V -k2,2g -k3,3g ", outputfile, " -o ", outputfile))

  cat(paste0("\t- ",outputfile,"\n"))

  # Print annotation database log information
  cat("\n\t- Annotation database log information\n")
  print(txdb)
}

##################### USER DEFINIED FUNCTIONS #########################
# Load Libraries
load_libraries <- function(){
	# Load libraries at start
	# if (!requireNamespace("BiocManager", quietly=TRUE))
  #     install.packages("BiocManager")
  # BiocManager::install("annotatr")
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(ChIPseeker))
  suppressPackageStartupMessages(library(genomation))
  suppressPackageStartupMessages(library(UpSetR))
  suppressPackageStartupMessages(library(AnnotationHub))
  suppressPackageStartupMessages(library(org.Mm.eg.db))
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  suppressPackageStartupMessages(library(TxDb.Mmusculus.UCSC.mm10.ensGene))
  suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
}

check_options <- function(){
	# Description of the script
	desc <- sprintf("
	----------------- SAMPLE USAGE ------------------
		- Rscript scripts/R_annotate_peaks_chipseeker.R -if=/media/rad/HDD1/atacseq/anja/rep_tALLcellLineMm/analysis/rep_tALLcellLineMm_analysis_consensus_peaks_annotation.tab -of=/media/rad/HDD1/atacseq/anja/rep_tALLcellLineMm/analysis/rep_tALLcellLineMm_analysis_consensus_peaks_annotation.txt -sp=mouse
	-------------------------------------------------
	CONTACT: 
		Gaurav Jain
		gaurav.jain@tum.de
	-------------------------------------------------\n
	")
	# create parser object
	parser <- ArgumentParser(description=cat(desc))

	# Add arguments 
	parser$add_argument("-if", "--inputfile"  , dest="inputfile"  , help="*Input mac2 peaks tab file"   , type="character", required=TRUE)
	parser$add_argument("-of", "--outputfile" , dest="outputfile" , help="*Output annotation file"      , type="character", required=TRUE)
	parser$add_argument("-sp", "--species"    , dest="species"    , help="*Species for annotation"      , type="character", required=TRUE, choices=(c('human', 'mouse')))
	parser$add_argument("-mf", "--macs2peaks" , dest="macs2Peaks", help="If set, the input file is macs2 peaks", action='store_true')
	parser$add_argument("-ac", "--addchr"     , dest="addChr"    , help="If set, add 'chr' to the first column", action='store_true')

	# Print the help message only if no arguments are supplied
	if(length(commandArgs(TRUE))==0){
		cat(desc)
		parser$print_help()
		quit()
	}

	# get command line options, if help option encountered print help and exit,
	# otherwise if options not found on command line then set defaults,
	args <- parser$parse_args()
	return(args)
}

## Call the main function in the end
main()

