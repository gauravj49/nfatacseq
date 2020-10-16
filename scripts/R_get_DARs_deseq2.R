# GOAL  : Identification of differential accessibility regions (DARs) between different biological conditions
# USAGE : 


# Command Arguments
warnings         <- warnings();
args             <- commandArgs(TRUE);
countsMatrixFile <- args[1]                           # 
conditionVals    <- unlist(strsplit(args[2], ','));   # 
outputRlePlot    <- args[3]

cat("- Loading libraries...\n")
suppressPackageStartupMessages(library("DESeq2", warn.conflicts=FALSE, quietly=TRUE))
suppressPackageStartupMessages(library("edgeR" , warn.conflicts=FALSE, quietly=TRUE))
suppressPackageStartupMessages(library("RUVSeq", warn.conflicts=FALSE, quietly=TRUE))
suppressPackageStartupMessages(library("BiocParallel", warn.conflicts=FALSE, quietly=TRUE))
register(MulticoreParam(4))

main <- function() {
	cat("\n- Getting input data\n")

	# Get input data 
	sampleCountsDF <- read.table(countsMatrixFile, row.names=1, sep="\t", header=T, quote = "")


	cat("\n- Getting input data...\n")
	# Apply filter
	filter         <- apply(sampleCountsDF, 1, function(x) length(x[x>1])>=2)
	filtered       <- sampleCountsDF[filter,]

	cat("\n- Getting matrix of counts data from counts matrix file...\n")
	# Create RLE plots for matrix
	matset         <- newSeqExpressionSet(as.matrix(filtered))

	# Get basenames
	outputDir       <- dirname(outputRlePlot)
	basefilename    <- tools::file_path_sans_ext(basename(outputRlePlot))
	extfilename     <- tools::file_ext(outputRlePlot)








	# Get RLE plots for same scale
	rlePlotsDir     <- paste(outputDir, "/rle_plots", sep=''); system(paste("mkdir -p", rlePlotsDir, sep=' '));
	ylimvals        <- c(-5,5)
	cat("\n- Getting RLE plots for same scale (",ylimvals,")...\n")
	RLEplotFile     <- paste(rlePlotsDir, "/", basefilename, ".png", sep='') 
	png(filename=RLEplotFile, height=900, width=1000, bg="white", res=100)
	par(mar=c(20,4,1,1));
	plotRLE(matset, outline=FALSE,col='royalblue',las=2, cex.axis=0.5, ylim=ylimvals, main = tools::file_path_sans_ext(basename(RLEplotFile)), ylab = "RLE")
	dev.off()

	cat("\n- Getting zoomedin RLE plots...\n")
	# Get zoomedin RLE plots
	zrlePlotsDir    <- paste(outputDir, "/rle_plots/zoomedIn", sep=''); system(paste("mkdir -p", zrlePlotsDir, sep=' '));
	zRLEplotFile    <- paste(zrlePlotsDir, "/", basefilename, "_zoomedIN.png", sep='') 
	png(filename=zRLEplotFile, height=900, width=1000, bg="white", res=100)
	par(mar=c(20,4,1,1));
	plotRLE(matset, outline=FALSE,col='royalblue',las=2, cex.axis=0.5,main = tools::file_path_sans_ext(basename(zRLEplotFile)), ylab = "RLE")
	dev.off()
}

# Run the main function
main()
