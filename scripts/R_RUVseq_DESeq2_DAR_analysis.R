# GOAL  : Find Differentially Expressed Genes/Transcripts/SmallRNAs
# USAGE : Rscript scripts/buildBoundaryModel.R input/H1ES/hESC_combined.domain.bed input/chrom_coordinates_HG19.bed 20000 1
# SOURCE: 
#		- https://harshinamdar.wordpress.com/2014/11/11/quick-tutorial-on-deseq2/
#	 	- https://sites.google.com/site/princetonhtseq/tutorials/rna-seq
# Contact: Gaurav Jain (gaurav.jain@dzne.de)

# Command Arguments
warnings         <- warnings();
args             <- commandArgs(TRUE);
inputDir         <- args[1]
ControlFiles     <- as.vector(read.table(args[2], header=FALSE, sep="\t", na.strings="NA", quote = "")$V1);
TreatmentFiles   <- as.vector(read.table(args[3], header=FALSE, sep="\t", na.strings="NA", quote = "")$V1);
conditionVals    <- unlist(strsplit(args[4], ','));
outputFile       <- args[5];
erccSpikeINs     <- args[6];
rnaClass         <- args[7];  # mirna or allncrna or genes
normalizedReads  <- args[8];  # If set, also get normalized reads
filterGenes      <- args[9];  # If set FALSE, do not filter reads during ruvseq (Default: TRUE)
saveRdata        <- args[10]; # If set, save all data in the logs folder. This can restore all variables for debugging purpose later


main <- function() {
	Sys.setenv("DISPLAY"=":0")
	options(warn=-1)
	cat("1) Input Control Files  : \n")
	print(ControlFiles)

	cat("\n2) Input Treatment Files: \n")
	print(TreatmentFiles)

	cat("\n3) Condition Values     : \n")
	print(conditionVals)

	# Create the output directories and files
	cat("\n4) Create the output directories and files ...\n")
	outputDir       <- dirname(outputFile);
	resultsDir      <- paste(outputDir, "/results_DEseq2", sep='')
	resultsOutFile  <- paste(resultsDir, "/", basename(outputFile), "_DE_RESULTS.txt", sep='')
	pcaDir          <- paste(outputDir, "/pca", sep='')
	digPlotDirPng   <- paste(outputDir, "/diagnosis_plots", sep='')
	digPlotDirPdf   <- paste(outputDir, "/diagnosis_plots/pdf", sep='')
	system(paste("mkdir -p", pcaDir, resultsDir, digPlotDirPng, digPlotDirPdf, sep=' '))

	# load libraries
	cat("\t4.1) Load DESeq2 library ...\n")
	suppressPackageStartupMessages(library("DESeq2", warn.conflicts=FALSE, quietly=TRUE))

	# Setup the sample table before filtering
	cat("\t4.2) Setup the sample table before filtering ...\n")
	sampleFiles            <- c(ControlFiles, TreatmentFiles)
	sampleCondition        <- c(rep(conditionVals[1], length(ControlFiles)), rep(conditionVals[2], length(TreatmentFiles)))
	sampleTable            <- data.frame(sampleName = sampleFiles, fileName = sampleFiles, condition = sampleCondition)
	dds                    <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = inputDir, design= ~ condition)
	colData(dds)$condition <- factor(colData(dds)$condition, levels=conditionVals)
	sampleCountsDF         <- as.data.frame(counts(dds))

	# Remove rnaClass (mirna or allncrna) from column names
	oldColNames            <- names(sampleCountsDF)
	newColNames			   <- gsub(paste("_",rnaClass,"Counts.txt", sep=''), "" , oldColNames)
	names(sampleCountsDF)  <- newColNames

	cat("\t4.3) Load RColorBrewer library and create colors for plots ...\n")
	suppressPackageStartupMessages(library(RColorBrewer, warn.conflicts=FALSE, quietly=TRUE))
	(mycols <- brewer.pal(8, "Dark2")[1:length(unique(factor(sampleCondition)))])

	# RUN RUVseq
	cat("\n5) Remove unwanted variation using RUVseq ...\n")
	# Create the pdf filename
	pdffile<-paste(pcaDir, "/", basename(outputFile), "_ruvseq_pca.pdf",sep='')
	# Start PDF device driver to save output to figure.pdf
	pdf(pdffile)
	# Run RUVseq and get the set object
	if (erccSpikeINs){
		set3<- run_ruvseq(conditionVals, sampleCountsDF, TRUE, filterGenes)
	}else{
		set3<- run_ruvseq(conditionVals, sampleCountsDF, FALSE, filterGenes)
	}

	cat("\n6) Run differential expression analysis with DESeq2 ...\n")
	cat("\t6.1) Create a coldata frame with unwanted variations factors and instantiate the DESeqDataSet ...\n")
	fdds <- DESeqDataSetFromMatrix(countData = counts(set3), colData = pData(set3), design = ~ W_1 + x)
	# fdds <- DESeq(fdds, fitType="local")
	# fdds <- DESeq(fdds, fitType="mean")
	fdds <- DESeq(fdds)
	
	cat("\t6.2) Regularized log transformation for clustering/heatmaps, etc ...\n")
	# Transform data to log space and visualize samples
	rld <- rlogTransformation(fdds, blind = TRUE)

	cat("\n7) Get results and write to output files ...\n")
	cat("\t7.1) Get results and order it according to padjusted values ...\n")
	# Get results and write to output files
	res <- results(fdds)
	res <- res[order(res$padj), ]

	# Save data results and normalized reads to csv!
	cat("\t7.2) Save data results and normalized reads to csv ...\n")
	resdata <- merge(as.data.frame(res), as.data.frame(counts(fdds,normalized=T)), by='row.names',sort=F)
	names(resdata)[1] <- 'feature'
	write.table(resdata, file = resultsOutFile, row.names = F, sep = '\t', quote = F)

	cat("\n8) Plotting diagnosis plots\n")
	cat("\t8.1) Plotting PCA plots...\n")
	plotpca(rld, pcaDir, outputFile, mycols)

	cat("\t8.2) Plotting volcano plots...\n")
	tryCatch({
	plot_volcanoPlot(resdata, digPlotDirPng, digPlotDirPdf, outputFile)
	},
	error=function(e){
		cat("\n\t- Error in plotting volcano plots")
		print(e)
	},
	warning=function(w){
		print(w)
	})
	
	cat("\t8.3) Plotting MA plots...\n")
	tryCatch({
	plot_maPlot(resdata, digPlotDirPng, digPlotDirPdf, outputFile)
	},
	error=function(e){
		cat("\n\t- Error in plotting MA plots")
		print(e)
	},
	warning=function(w){
		print(w)
	})

	cat("\t8.4) Plotting pvalues plots...\n")
	tryCatch({
	plot_pvalPlot(res, digPlotDirPng, digPlotDirPdf, outputFile)
	},
	error=function(e){
		cat("\n\t- Error in plotting pvalue plots:")
		print(e)
	},
	warning=function(w){
		print(w)
	})

	cat("\t8.5) Plotting independent filtering plots...\n")
	tryCatch({
		plot_indfilPlot(res, digPlotDirPng, digPlotDirPdf, outputFile)
	},
	error=function(e){
		cat("\n\t- Error in plotting independent filtering plots")
		print(e)
	},
	warning=function(w){
		print(w)
	})

	tryCatch({
		cat("\t8.6) Plotting sample correlation heatmaps...\n")
	plot_sampleDistplot(fdds, rld, digPlotDirPng, digPlotDirPdf, outputFile)
	},
	error=function(e){
		cat("\n\t- Error in plotting sample correlation heatmaps")
		print(e)
	},
	warning=function(w){
		print(w)
	})


	cat("\n ***********************************\n")
	cat("- Output Results Files          : ", resultsOutFile, "\n")
	options(warn=0)

	# Get the normalized data in separate folders with RLE plots comparing SFN, VST and RAW reads
	if (normalizedReads){
		cat("\n9) Get normalized reads ...\n")
		get_normalized_reads(fdds, sampleCountsDF, outputFile)
	}

	# Print session info to the session's log file
	logdir <- paste(outputDir,"/logs", sep='')
	system(paste("mkdir -p ", logdir, sep=''))
	session_logfile <- paste(logdir,"/session_info_",basename(outputFile),".log" ,sep='');
	print_session_info(session_logfile)

	## Save all data to a Rdata file for debugging purpose
	# # Source: http://stat.ethz.ch/R-manual/R-devel/library/base/html/load.html
	# # To load later in R interactive environment i.e. restore the saved values to the current environment
	# local({
	# 	attach(session_rdata)
	# 	ls()
	# })
	if (saveRdata){
		session_rdata <- paste(logdir,"/session_data_",basename(outputFile),".rda" ,sep='');
		save(list = ls(all = TRUE), file= session_rdata)
	}
}

########## USER DEFINED FUNCTIONS #################
# Run RUVseq
run_ruvseq <- function(conditionVals, sampleCountsDF, erccSpikeINs=FALSE, filterGenes=TRUE){
	suppressPackageStartupMessages(library("edgeR" , warn.conflicts=FALSE, quietly=TRUE))
	suppressPackageStartupMessages(library("RUVSeq", warn.conflicts=FALSE, quietly=TRUE))

	################ First pass RUVseq ############
	# par: sets the bottom, left, top and right margins

	# Filter out non-expressed genes, by requiring more than 1 reads in at least two samples for each gene.
	cat("\t5.1) Filter out non-expressed genes ...\n")
	if (filterGenes){
		filter        <- apply(sampleCountsDF, 1, function(x) length(x[x>1])>=2)
		filtered      <- sampleCountsDF[filter,]
	}else{
		filtered      <- sampleCountsDF
	}

	filteredGenes <- rownames(filtered)
	#	genes <- filteredGenes[grep("^ENS", filteredGenes)]
	if (erccSpikeINs){
		spikes <- rownames(filtered)[grep("ERCC", rownames(filtered))]
	}
	genes <- filteredGenes
	cat("\t\t- Total filtered genes:", length(genes), "\n")

	# Get RUVseq Expression set
	cat("\t5.2) Get RUVseq Expression set ...\n")
	#***********************************************************************
	# IMPORTANT: Need to create the sampleCondition again ... sometimes the control and treatment is switched which gives inverted log2fc
	# Here very important part is to recreate the factor as they might switch the levels .... 
	
	# # WRONG LEVELS:
	#  [1] HypCL HypCL HypCL HypCL HypCL HypCL Hyp4h Hyp4h Hyp4h Hyp4h Hyp4h Hyp4h
	#  	 Levels: Hyp4h HypCL

	# # CORRECT LEVELS:
	#  [1] HypCL HypCL HypCL HypCL HypCL HypCL Hyp4h Hyp4h Hyp4h Hyp4h Hyp4h Hyp4h
	#  	 Levels: HypCL Hyp4h
	#***********************************************************************
	
	sampleCondition        <- c(rep(conditionVals[1], length(ControlFiles)), rep(conditionVals[2], length(TreatmentFiles)))
	x   <- factor(sampleCondition, levels=conditionVals)
	set <- newSeqExpressionSet(as.matrix(filtered), phenoData = data.frame(x, row.names=colnames(filtered)))

	# The boxplots of relative log expression (RLE = log-ratio of read count to median read count across sample) and plots of principal components (PC) in Figure 1 reveal a clear need for betwen-sample normalization.
	colors <- brewer.pal(3, "Set2")
	par(mar=c(20,4,1,1));       # bottom, left, top and right margins
	plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x], las=2)
	par(mar=c(5.1,4.1,4.1,2.1)) # reset to defaults (not necessarilly current values)
	plotPCA(set, col=colors[x], cex=0.9)

	if (erccSpikeINs){
		# 2.2 RUVg: Estimating the factors of unwanted variation using control genes
		# To estimate the factors of unwanted variation, we need a set of negative control genes, i.e., genes that can be assumed not to be influenced by the covariates of interest (in the case of the zebrafish dataset, the Gallein treatment). In many cases, such a set can be identified, e.g., housekeeping genes or spike-in controls. If a good set of negative controls is not readily available, one can define a set of “in-silico empirical” controls as in Section 2.4.
		# Here, we use the ERCC spike-ins as controls and we consider k = 1 factors of unwanted variation.
		cat("\t5.3) Set1 : RUVg normalization based on spike-in controls ...\n")
		set1 <- RUVg(set, spikes, k=1)
		par(mar=c(20,4,1,1));       # bottom, left, top and right margins
		plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=colors[x], las=2)
		par(mar=c(5.1,4.1,4.1,2.1)) # reset to defaults (not necessarilly current values)
		plotPCA(set1, col=colors[x], cex=0.9)
	}else{
		cat("\t5.3) Set1: Not used here because of no spike-in controls ...\n")
	}

	tryCatch({
			# 2.4 Empirical control genes
			# If no genes are known a priori not to be influenced by the covariates of interest, one can obtain a set of “in-silico empirical” negative controls, e.g., least significantly DE genes based on a first-pass DE analysis performed prior to RUVg normalization
			design      <- model.matrix(~x, data=pData(set))
			y           <- DGEList(counts=counts(set), group=x)
			y           <- calcNormFactors(y, method="upperquartile")
			y           <- estimateGLMCommonDisp(y, design)
			y           <- estimateGLMTagwiseDisp(y, design)
			fit         <- glmFit(y, design)
			lrt         <- glmLRT(fit, coef=2)
			top         <- topTags(lrt, n=nrow(set))$table
			toplen      <- length(rownames(top))
			maxempGenes <- round(toplen/3)
			empirical   <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:maxempGenes]))]
	
			# Here, we consider all but the top 5,000 genes as ranked by edgeR p-values.
			cat("\t5.4) Set2: RUVg - Empirical control genes ...\n")
			cat("\t\t - Here, we consider all but the top ",maxempGenes," genes as ranked by edgeR p-values\n")
			set2 <- RUVg(set, empirical, k=1)
			par(mar=c(20,4,1,1));       # bottom, left, top and right margins
			plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[x], las=2)
			par(mar=c(5.1,4.1,4.1,2.1)) # reset to defaults (not necessarilly current values)
			plotPCA(set2, col=colors[x], cex=0.9)
		},
		error=function(e){
			cat("\n\t5.4) Set2: Error in performing RUVseq from empirical control genes:")
			print(e)
		},
		warning=function(w){
			cat("\n\t5.4) Set2: Warning in performing RUVseq from empirical control genes:")
			print(w)
		}
	)

	cat("\t5.5) Set3: RUVs - Estimating the factors of unwanted variation using replicate samples ...\n")
	# 3 RUVs: Estimating the factors of unwanted variation using replicate samples
	# As an alternative approach, one can use the RUVs method to estimate the factors of unwanted variation using replicate/negative control samples for which the covariates of interest are constant
	# First, we need to construct a matrix specifying the replicates. In the case of the zebrafish dataset, we can consider the three treated and the three control samples as replicate groups. This information is passed to RUVs in the following way.
	differences <- makeGroups(sampleCondition, conditionVals)
	set3 <- RUVs(set, genes, k=1, differences)
	par(mar=c(20,4,1,1));       # bottom, left, top and right margins
	plotRLE(set3, outline=FALSE, ylim=c(-4, 4), col=colors[x], las=2)
	par(mar=c(5.1,4.1,4.1,2.1)) # reset to defaults (not necessarilly current values)
	plotPCA(set3, col=colors[x], cex=0.9)

	par(new=TRUE)
	# Turn off device driver (to flush output to PNG/PDF file)
	dev.off()
	if (erccSpikeINs){
		return(set1)
	} else {
		return(set3)
	}
}

# Make a matrix suitable for use with RUVSeq methods such as RUVs().
makeGroups <- function(xs, conditionVals) {
	#***********************************************************************
	# IMPORTANT: Need to create the sampleCondition again ... sometimes the control and treatment is switched which gives inverted log2fc
	# Here very important part is to recreate the factor as they might switch the levels .... 
	
	# # WRONG LEVELS:
	#  [1] HypCL HypCL HypCL HypCL HypCL HypCL Hyp4h Hyp4h Hyp4h Hyp4h Hyp4h Hyp4h
	#  	 Levels: Hyp4h HypCL

	# # CORRECT LEVELS:
	#  [1] HypCL HypCL HypCL HypCL HypCL HypCL Hyp4h Hyp4h Hyp4h Hyp4h Hyp4h Hyp4h
	#  	 Levels: HypCL Hyp4h
	#***********************************************************************
	

	# Each row in the returned matrix corresponds to a set of replicate samples.
	# The number of columns is the size of the largest set of replicates; rows for
	# smaller sets are padded with -1 values.
	# 
	# @param xs A vector indicating membership in a group.
	# @seealso RUVSeq::RUVs
	# @example 
	#  makeGroups(c("A","B","B","C","C","D","D","D","A"))
	#       [,1] [,2] [,3]
	#  [1,]    1    9   -1
	#  [2,]    2    3   -1
	#  [3,]    4    5   -1
	#  [4,]    6    7    8
	xs <- factor(xs, levels=conditionVals)
	groups <- matrix(-1, nrow = length(levels(xs)), ncol = max(table(xs)))
	for (i in 1:length(levels(xs))) {
		idxs <- which(xs == levels(xs)[i])
		groups[i,1:length(idxs)] <- idxs
	}
	groups
}

changeLevels <- function(dtable, column, newlevs, recycling=FALSE) {
    # Source : https://github.com/rsaporta/pubR/blob/gitbranch/dt_changeLevels.R
    # Source2: https://stackoverflow.com/questions/14634964/how-does-one-change-the-levels-of-a-factor-column-in-a-data-table#14635190

    # get name of table
    dtName <- as.character(match.call()[[2]])

    # We are getting the table from the calling environment,  get(dtName, envir=parent.frame(1))
    #   and making changes to [[column]]

    # TODO:  Implement generic function
    # First, check to make sure it's a dt. 
    dtClass <- class(get(dtName, envir=parent.frame(1)))
    if (!"data.table" %in% dtClass) {
      stop ("dtable is not a data.table.  Please use `levels<-` for data.frame")
    }

    # Check for appropriate length on new levels
    oldLevsLength <- length(levels(get(dtName, envir=parent.frame(1))[[column]]))
    newLevsLength <- length(newlevs)

    # If levels is wrong size, we will give feedback as to the difference in size. 
    diffs <- newLevsLength - oldLevsLength
    
    # newlevs is too big 
    if (diffs > 0) {
      stop ("New levels has ", diffs, " level", ifelse(diffs==1, "", "s"), " too many." )
    
    # newlevs is too small
    } else if (diffs < 0) {
      stop ("New levels has ", -diffs, " level", ifelse(diffs==-1, "", "s"), " too few." )
    
    # this clause should never hit, but putting it in just to be safe. 
    } else if (diffs != 0) { 
      stop ("Something went wrong. Unsure what.")
    }

    # change levels using setattr();  passing through the return value.  
    return(setattr( get(dtName, envir=parent.frame(1))[[column]], 
              "levels",newlevs))

 
	##-----------------------------  EXAMPLE  -----------------------------###
	#   library(data.table)
	#   mydt <- data.table(id=1:6, value=as.factor(c("A", "A", "B", "B", "B", "C")), test=c(2, 2, 3, 4, 5, 6), key="id") 
	#
	#   newLevs.good    <- c("X", "Y", "Z")
	#   newLevs.tooFew  <- c("P", "Q")
	#   newLevs.tooMany <- c("R", "S", "T", "U", "W")
	#   originalLevs    <- c("A", "B", "C") 
	# 
	#
	#   changeLevels(mydt, "value", newLevs.good);  mydt
	#   changeLevels(mydt, "value", originalLevs);  mydt
	#   changeLevels(mydt, "value", newLevs.tooFew);  mydt
	#   changeLevels(mydt, "value", newLevs.tooMany);  mydt 
	##-----------
  }


# Draw PCA plots
plotpca <- function(rld, pcaDir, outputFile, mycols){
	# Create the png filename
	pngfile  <- paste(pcaDir, "/", basename(outputFile), "_pca.png",sep='');
	# Start PNG device driver to save output to figure.png
	png(filename=pngfile, height=8,width=8, res=300, units="in");
	# Plot PCA without labels
	rld_pca(rld, colors=mycols, intgroup="x", main=paste("PCA (", conditionVals[2],", ", conditionVals[1],")", sep=""), pcaDir=pcaDir)

	# Create the pdf filename
	pdffile<-paste(pcaDir, "/", basename(outputFile), "_pca.pdf",sep='')
	# Start PDF device driver to save output to figure.pdf
	pdf(pdffile)
	# Plot PCA with labels
	rld_pca(rld, colors=mycols, intgroup="x", main=paste("PCA (", conditionVals[2],", ", conditionVals[1],")", sep=""), plotlabels=TRUE, pcaDir=pcaDir)
}

rld_pca <- function (rld, intgroup = "condition", ntop = 500, plotlabels=FALSE, colors=NULL, legendpos="topright", main="PCA Biplot", textcx=1, pcaDir, ...) {
	suppressPackageStartupMessages(require(genefilter, warn.conflicts=FALSE, quietly=TRUE))
	suppressPackageStartupMessages(require(calibrate, warn.conflicts=FALSE, quietly=TRUE))
	suppressPackageStartupMessages(require(RColorBrewer, warn.conflicts=FALSE, quietly=TRUE))
	rv = rowVars(assay(rld))
	select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
	pca = prcomp(t(assay(rld)[select, ]))
	fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
	if (is.null(colors)) {
		if (nlevels(fac) >= 3) {
			colors = brewer.pal(nlevels(fac), "Paired")
		} else {
			colors = c("black", "red")
		}
	}
	
	cexsize = 2.0
	if (plotlabels){
		# Get the names
		rnames          <- rownames(as.data.frame(pca$x))
		nnum            <- seq(1:length(rnames))
		rarray          <- setNames(nnum, rnames)
		pcIDs           <- as.data.frame(as.matrix(rarray))

		# Rename first column as snow
		names(pcIDs)    <- c("sno")
		
		# Put rownames as second column
 		pcIDs <- cbind(data.frame(pcIDs, row.names=NULL), rownames(pcIDs))

 		# Rename the second column
		colnames(pcIDs)[2] <- 'sampleNames'

		# Save the dataframe in the output file
		txtFiledir      <- paste(pcaDir,"/fileIds", sep=''); system(paste("mkdir -p ", txtFiledir, sep=''));
		sampleIDfile    <- paste(txtFiledir, "/", basename(outputFile), "_fileIds.txt",sep='');
		write.table(pcIDs, file = sampleIDfile, row.names = F, sep = '\t', quote = F)
		
	
		# Put legend on bottom 1/8th of the chart
		layout(matrix(c(1,2,2)), heights=c(length(rnames)*6, length(rnames)*2))  

		# circle size
		cexsize = 1.0 
	}

	# Plot the pca
	pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
	pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
	pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
	pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
	
	plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], col=NULL, pch=21, cex=cexsize, xlab=pc1lab, ylab=pc2lab, main=main, ...)

	if (plotlabels){
		# Add text to the points
		with(as.data.frame(pca$x), textxy(PC1, PC2, labs=nnum, cex=0.7))

		# setup for no margins on the legend
		# c(bottom, left, top, right)
		par(mar=c(0, 0, 0, 0))
		plot.new()

		legend('center','groups', paste(nnum, rnames, sep=" = "), xpd=TRUE, ncol=2, cex=0.5, text.col=colors[fac], bty = "n")
		legend("top", legend=levels(fac), col=colors, pch=20, xpd=TRUE, ncol=2, cex=0.5, bty = "n")

		# Restore default clipping rect
		par(mar=c(5, 4, 4, 2) + 0.1)
	} else{
		# Expand right side of clipping rect to make room for the legend
		legend("top",                   # Location of legend 
			xpd = TRUE,                      # Allow drawing outside plot area
			xjust = 0,                       # Left justify legend box on x
			yjust = .5,                      # Center legend box on y
			legend = levels(fac),            # legend element labels
			col = colors,                    # Legend Element colors
			pch = 20,                        # Legend Element Styles
			cex = 0.5,
			ncol=2
		)

	}

	par(new=TRUE)
	# Turn off device driver (to flush output to PNG file)
	dev.off()
}

## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="top", labelsig=TRUE, textcx=1, ...) {
	with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
	with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
	with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
	with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
	if (labelsig) {
		require(calibrate)
		with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
	}
	legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR < ",sigthresh,sep=""), paste("abs(LogFC) > ",lfcthresh,sep=""), "both "), pch=20, col=c("red","orange","green"), ncol=3)

	par(new=TRUE)
	# Turn off device driver (to flush output to PNG/PDF file)
	dev.off()
}

plot_volcanoPlot <- function(resdata, digPlotDirPng, digPlotDirPdf, outputFile){
	mtitle <- paste("Volcano Plot (", conditionVals[2],", ", conditionVals[1],")", sep="")
	
	# Plot scaled volcano plots 
	# Create the png filename and start PNG device driver to save output to figure.png
	pngfile  <- paste(digPlotDirPng, "/", basename(outputFile), "_volcanoplot.png",sep='');
	png(filename=pngfile, height=8,width=8, res=300, units="in");
	volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2), ylim=c(0, 15), labelsig=FALSE, main=mtitle)
	
	# Create the pdf filename and start PDF device driver to save output to figure.pdf
	pdffile<-paste(digPlotDirPdf, "/", basename(outputFile), "_volcanoplot.pdf",sep='')
	pdf(pdffile)
	volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2), ylim=c(0, 15), labelsig=FALSE, main=mtitle)

	# Plot auto scaled volcano plots 
	# Create the png filename and start PNG device driver to save output to figure.png
	pngfile  <- paste(digPlotDirPng, "/", basename(outputFile), "_autoscaled_volcanoplot.png",sep='');
	png(filename=pngfile, height=8,width=8, res=300, units="in");
	volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, labelsig=FALSE, main=mtitle)
	
	# Create the pdf filename and start PDF device driver to save output to figure.pdf
	pdffile<-paste(digPlotDirPdf, "/", basename(outputFile), "_autoscaled_volcanoplot.pdf",sep='')
	pdf(pdffile)
	volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, labelsig=FALSE, main=mtitle)
}

## MA plot
## Could do with built-in DESeq2 function:
## DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
maplot <- function (res, main="MA Plot", thresh=0.05, labelsig=TRUE, textcx=1, ...) {
	with(res, plot(baseMean + 0.1, log2FoldChange, pch=20, cex=.5, log="x", main=main, ...))
	with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
	if (labelsig) {
		require(calibrate)
		with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
	}
	par(new=TRUE)
	# Turn off device driver (to flush output to PNG/PDF file)
	dev.off()
}

plot_maPlot <- function(resdata, digPlotDirPng, digPlotDirPdf, outputFile){
	mtitle <- paste("MA Plot (", conditionVals[2],", ", conditionVals[1],")", sep="")

	# Create the png filename
	pngfile  <- paste(digPlotDirPng, "/", basename(outputFile), "_maplot.png",sep='');
	# Start PNG device driver to save output to figure.png
	png(filename=pngfile, height=8,width=8, res=300, units="in");
	# Plot PCA without labels
	maplot(resdata, main=mtitle, labelsig=FALSE)

	# Create the pdf filename
	pdffile<-paste(digPlotDirPdf, "/", basename(outputFile), "_maplot.pdf",sep='')
	# Start PDF device driver to save output to figure.pdf
	pdf(pdffile)
	# Plot PCA with labels
	maplot(resdata, main=mtitle, labelsig=FALSE)
}

## Examine plot of p-values
pvalplot <- function (res, main="p-values", ...) {
	#hist(res$pvalue, breaks=50, col="grey", main=main)
	hist(res$pvalue[res$baseMean > 1], breaks=50, col="grey", main=main)
	par(new=TRUE)
	# Turn off device driver (to flush output to PNG/PDF file)
	dev.off()
}

plot_pvalPlot <- function(res, digPlotDirPng, digPlotDirPdf, outputFile){
	mtitle <- paste("Pvalues Plot (", conditionVals[2],", ", conditionVals[1],")", sep="")
	
	# Create the png filename
	pngfile  <- paste(digPlotDirPng, "/", basename(outputFile), "_pvalplot.png",sep='');
	# Start PNG device driver to save output to figure.png
	png(filename=pngfile, height=8,width=8, res=300, units="in");
	# Plot PCA without labels
	pvalplot(res, main=mtitle)

	# Create the pdf filename
	pdffile<-paste(digPlotDirPdf, "/", basename(outputFile), "_pvalplot.pdf",sep='')
	# Start PDF device driver to save output to figure.pdf
	pdf(pdffile)
	# Plot PCA with labels
	pvalplot(res, main=mtitle)
}

## Examine independent filtering
indfilplot <- function (res, main="Independent filtering",...) {
	plot(metadata(res)$filterNumRej, type="b", xlab="quantiles of 'baseMean'", ylab="number of rejections", main=main)
	par(new=TRUE)
	# Turn off device driver (to flush output to PNG/PDF file)
	dev.off()
}

plot_indfilPlot <- function(res, digPlotDirPng, digPlotDirPdf, outputFile){
	mtitle <- paste("Independent Filtering Plot (", conditionVals[2],", ", conditionVals[1],")", sep="")

	# Create the png filename
	pngfile  <- paste(digPlotDirPng, "/", basename(outputFile), "_independent_filtering_plot.png",sep='');
	# Start PNG device driver to save output to figure.png
	png(filename=pngfile, height=8,width=8, res=300, units="in");
	# Plot PCA without labels
	indfilplot(res, main=mtitle)

	# Create the pdf filename
	pdffile<-paste(digPlotDirPdf, "/", basename(outputFile), "_independent_filtering_plot.pdf",sep='')
	# Start PDF device driver to save output to figure.pdf
	pdf(pdffile)
	# Plot PCA with labels
	indfilplot(res, main=mtitle)
}

## Heatmap of sample-to-sample distances using the Poisson Distance
sampleDistplot <- function (fdds, rld, main="Independent filtering",...) {
	# Another option for calculating sample distances is to use the Poisson Distance (Witten 2011), implemented in the PoiClaClu package. This measure of 
	# dissimilarity between counts also takes the inherent variance structure of counts into consideration when calculating the distances between samples.
	# The PoissonDistance function takes the original count matrix (not normalized) with samples as rows instead of columns, so we need to transpose the 
	# counts in dds.
	suppressMessages(library("PoiClaClu"))
	suppressMessages(library("pheatmap"))
	poisd <- PoissonDistance(t(counts(fdds)))

	# In order to plot the sample distance matrix with the rows/columns arranged by the distances in our distance matrix, we manually provide sampleDists 
	# to the clustering_distance argument of the pheatmap function. Otherwise the pheatmap function would assume that the matrix contains the data values 
	# themselves, and would calculate distances between the rows/columns of the distance matrix, which is not desired.
	samplePoisDistMatrix <- as.matrix( poisd$dd)

	# Change the row names of the distance matrix to contain Treatment and Control instead of sample ID, so that we have all this information in view when looking at the heatmap.
	rownames(samplePoisDistMatrix) <- paste( attr(rld,"colData")$condition, rownames(attr(rld,"colData")),sep="-")
	colnames(samplePoisDistMatrix) <- NULL
	# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
	colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
	# colors <- colorRampPalette(rev(c('gold','darkorange','darkred')))(256)
	# colors <- colorRampPalette(c('darkblue','blue','lightblue','white','orange','red','darkred'))(1024)
	pheatmap(samplePoisDistMatrix,
         clustering_distance_rows=poisd$dd,
         clustering_distance_cols=poisd$dd,
         col=colors,
         fontsize=8,
         main=main)

	#par(new=TRUE)
	# Turn off device driver (to flush output to PNG/PDF file)
	dev.off()
}

plot_sampleDistplot <- function(fdds, rld, digPlotDirPng, digPlotDirPdf, outputFile){
	mtitle <- paste("Sample Correlation (Poisson Distance) Plot (", conditionVals[2],", ", conditionVals[1],")", sep="")

	# Create the png filename
	pngfile  <- paste(digPlotDirPng, "/", basename(outputFile), "_sampleCorrelation_poissonDistance_plot.png",sep='');
	# Start PNG device driver to save output to figure.png
	png(filename=pngfile, height=8,width=16, res=300, units="in");
	# Plot PCA without labels
	sampleDistplot(fdds, rld, main=mtitle)

	# Create the pdf filename
	pdffile<-paste(digPlotDirPdf, "/", basename(outputFile), "_sampleCorrelation_poissonDistance_plot.pdf",sep='')
	# Start PDF device driver to save output to figure.pdf
	pdf(pdffile)
	# Plot PCA with labels
	sampleDistplot(fdds, rld, main=mtitle)
}

# Print session info as log file formatted in tabular format
# Source: https://stackoverflow.com/questions/21967254/how-to-write-a-reader-friendly-sessioninfo-to-text-file
print_session_info <- function(session_logfile){
	suppressPackageStartupMessages(library("devtools"))
	suppressPackageStartupMessages(library("knitr"))

	# Get all the session info to the variable
	my_session_info <- devtools::session_info()

	# Print it in the tabular format using knitr
	writeLines(text = {
	    paste(sep = "\n", collapse = "",
	          paste0(rep("-", 80), collapse = ""),
	          paste(paste0(rep("-", 32), collapse = ""),
	                "R environment",
	                paste0(rep("-", 33), collapse = "")),
	          paste0(rep("-", 80), collapse = ""),
	          paste(knitr::kable(data.frame(setting = names(my_session_info$platform),
	                                  value = as.character(my_session_info$platform))), collapse = "\n"),
	          paste0(rep("-", 80), collapse = ""),
	          paste(paste0(rep("-", 35), collapse = ""),
	                "packages",
	                paste0(rep("-", 35), collapse = "")),
	          paste0(rep("-", 80), collapse = ""),
	          paste(knitr::kable(my_session_info$packages), collapse = "\n")
	          )
	}, con = session_logfile)
}

# Get normalized reads 
get_normalized_reads <- function(fdds, sampleCountsDF, outputFile){
	# 1) RLE plots for raw data
	# 2) Get size factor normalized reads (SFN)
	# 3) Get variance stabilizing transformation (VST) normalized reads

	# Get basenames
	outputDir       <- paste(dirname(outputFile),"/snormalized_counts", sep='')
	basefilename    <- tools::file_path_sans_ext(basename(outputFile))
	extfilename     <- tools::file_ext(outputFile)
	rlePlotsDir     <- paste(outputDir, "/rle_plots", sep=''); system(paste("mkdir -p", rlePlotsDir, sep=' '));
	ylimvals        <- c(-5,5)

	# Get annotation table
	sampleFiles     <- c(ControlFiles, TreatmentFiles)
	sampleCondition <- c(rep(conditionVals[1], length(ControlFiles)), rep(conditionVals[2], length(TreatmentFiles)))
	sampleDiagnosis <- c(rep(0, length(ControlFiles)), rep(1, length(TreatmentFiles)))
	newSampleTable  <- data.frame(sampleName = names(sampleCountsDF), condition = sampleCondition, diagnosis = sampleDiagnosis, row.names=names(sampleCountsDF))

	# 1) # Create RLE plots for RAW counts
	cat("\t9.1) Create RLE plots for RAW counts ...\n")
	filter          <- apply(sampleCountsDF, 1, function(x) length(x[x>1])>=2)
	filtered        <- sampleCountsDF[filter,]
	set             <- newSeqExpressionSet(as.matrix(filtered))
	rawRLEplotFile  <- paste(rlePlotsDir, "/", basefilename, "_raw.png", sep='') 
	png(filename=rawRLEplotFile, height=900, width=1000, bg="white", res=100)
	par(mar=c(20,4,1,1));
	plotRLE(set, outline=FALSE,col='navy',las=2, cex.axis=0.5, ylim=ylimvals)
	dev.off()

	# 2) Get quantile normalized reads (QNR)
	cat("\t9.2) Get quantile normalized reads (QNR) ...\n")
	suppressPackageStartupMessages(library("preprocessCore", warn.conflicts=FALSE, quietly=TRUE))
	qnrCountsDF       <- normalize.quantiles(as.matrix(sampleCountsDF))
	qnrfilter         <- apply(qnrCountsDF, 1, function(x) length(x[x>1])>=2)
	qnrfiltered       <- qnrCountsDF[qnrfilter,]
	qnrset            <- newSeqExpressionSet(as.matrix(qnrfiltered))

	qnrrawRLEplotFile <- paste(rlePlotsDir, "/", basefilename, "_qnr.png", sep='') 
	png(filename=qnrrawRLEplotFile, height=900, width=1000, bg="white", res=100)
	par(mar=c(20,4,1,1));
	plotRLE(qnrset, outline=FALSE,col='deepskyblue',las=2, cex.axis=0.5, ylim=ylimvals)
	dev.off()

	# Get QNR related directories and filenames
	qnrnormalizedDir      <- paste(outputDir,'/qnr/',conditionVals[2],'_',conditionVals[1], sep='')
	qnrtnormalizedDir     <- paste(qnrnormalizedDir, "/transposed", sep='')
	indqnrnormalizedDir   <- paste(qnrnormalizedDir, "/qnr_normalized_counts", sep='')
	qnrnormalizedOutFile  <- paste(qnrnormalizedDir , "/", basefilename, "_qnr_normalized.txt", sep='')
	qnrtnormalizedOutFile <- paste(qnrtnormalizedDir, "/", basefilename, "_featureInCols_qnr_normalized.txt", sep='')
	system(paste("mkdir -p", qnrnormalizedDir, qnrtnormalizedDir, indqnrnormalizedDir, sep=' '))

	# Convert rownames as feature column to save it in the csv file
	qnrrownames <- rownames(sampleCountsDF[qnrfilter,])
	qnrcolnames <- names(sampleCountsDF[qnrfilter,])
	colnames(qnrfiltered) <- qnrcolnames
	rownames(qnrfiltered) <- qnrrownames
	qnrfiltered <- cbind(feature = rownames(qnrfiltered), qnrfiltered)
	rownames(qnrfiltered) <- 1:nrow(qnrfiltered)
	
	# Save the counts file to output csv file
	write.table(qnrfiltered   , file = qnrnormalizedOutFile , row.names = F, sep = '\t', quote = F)
	write.table(t(qnrfiltered), file = qnrtnormalizedOutFile, col.names = F, sep = '\t', quote = F)

	# Get the normalized individual count files
	indvcmd <- paste("bash scripts/get_individual_countFiles.sh", qnrnormalizedOutFile, indqnrnormalizedDir, rnaClass, 1, 2, sep=' ')
	system(indvcmd)

	# Get the annotation dataframe by merging the sampleTable and qnrfiltered normalized counts
	annDF           <- merge(newSampleTable,t(qnrfiltered), by = "row.names")
	rownames(annDF) <- annDF[,1]
	annDF[,1]       <- NULL
	# Moeve diagnosis column to the end
	# https://stackoverflow.com/questions/18339370/reordering-columns-in-a-large-dataframe
	annDF <- annDF[c(setdiff(names(annDF), c("diagnosis")), c("diagnosis"))]
	qnrAnnOutFile   <- paste(qnrnormalizedDir, "/", basefilename, "_annotation_qnr_normalized.txt", sep='')
	write.table(annDF, file = qnrAnnOutFile, row.names = F, sep = '\t', quote = F)

	# 3) Get size factor normalized reads (SFN)
	cat("\t9.3) Get size factor normalized reads (SFN) ...\n")
	sfnCountsDF    <- as.data.frame(counts(fdds,normalized=T))
	# Create RLE plots for SFN counts
	sfnfilter      <- apply(sfnCountsDF, 1, function(x) length(x[x>1])>=2)
	sfnfiltered    <- sfnCountsDF[sfnfilter,]
	sfnset         <- newSeqExpressionSet(as.matrix(sfnfiltered))
	sfnRLEplotFile <- paste(rlePlotsDir, "/", basefilename, "_sfn.png", sep='') 
	png(filename=sfnRLEplotFile, height=900, width=1000, bg="white", res=100)
	par(mar=c(20,4,1,1));
	plotRLE(sfnset, outline=FALSE,col='lightskyblue',las=2, cex.axis=0.5, ylim=ylimvals)
	dev.off()

	# 4) Get variance stabilizing transformation (VST) normalized reads
	cat("\t9.4) Get variance stabilizing transformation (VST) normalized reads ...\n")
	# fitType="mean", a VST is applied for Negative Binomial distributed counts, 'k', with a fixed dispersion, 'a': ( 2 asinh(sqrt(a k)) - log(a) - log(4) )/log(2).
	vst <- varianceStabilizingTransformation(fdds, blind=FALSE, fitType="mean")

	# vsd is now the normalized log2-transformed data
	vsd <- assay(vst)  

	# Convert rownames as feature column to save it in the csv file
	vsd <- cbind(feature = rownames(vsd), vsd)
	rownames(vsd) <- 1:nrow(vsd)

	# Get VST related directories and filenames
	vstnormalizedDir      <- paste(outputDir,'/vst/',conditionVals[2],'_',conditionVals[1], sep='')
	vsttnormalizedDir     <- paste(vstnormalizedDir, "/transposed", sep='')
	indvstnormalizedDir   <- paste(vstnormalizedDir, "/vst_normalized_counts", sep='')
	vstnormalizedOutFile  <- paste(vstnormalizedDir , "/", basefilename, "_vst_normalized.txt", sep='')
	vsttnormalizedOutFile <- paste(vsttnormalizedDir, "/", basefilename, "_featureInCols_vst_normalized.txt", sep='')
	system(paste("mkdir -p", vstnormalizedDir, vsttnormalizedDir, indvstnormalizedDir, sep=' '))

	# Save the counts file to output csv file 
	# This step is just so that I can get the vsd object to dataframe
	write.table(vsd   , file = vstnormalizedOutFile , row.names = F, sep = '\t', quote = F)
	
	# Create RLE plots for VST counts
	vsd            <- read.table(vstnormalizedOutFile, row.names=1, sep="\t", header=T, quote = "", check.names=FALSE)

	# Remove the minimum to get 0s and 0s again
	vsdMin         <- min(vsd)
	vsd            <- vsd - vsdMin
	vstfilter      <- apply(vsd, 1, function(x) length(x[x>1])>=2) 
	vstfiltered    <- vsd[vstfilter,]
	vstset         <- newSeqExpressionSet(as.matrix(vstfiltered))
	rlePlotsDir    <- paste(outputDir, "/rle_plots", sep=''); system(paste("mkdir -p", rlePlotsDir, sep=' '));
	vstRLEplotFile <- paste(rlePlotsDir, "/", basefilename, "_vst.png", sep='') 
	png(filename=vstRLEplotFile, height=900, width=1000, bg="white", res=100)
	par(mar=c(20,4,1,1));
	plotRLE(vstset, outline=FALSE,col='dodgerblue',las=2, cex.axis=0.5, ylim=ylimvals)
	dev.off()

	vstZoomRLEplotFile <- paste(rlePlotsDir, "/", basefilename, "_vst_zoomedIN.png", sep='') 
	png(filename=vstZoomRLEplotFile, height=900, width=1000, bg="white", res=100)
	par(mar=c(20,4,1,1));
	plotRLE(vstset, outline=FALSE,col='dodgerblue',las=2, cex.axis=0.5)
	dev.off()

	# Before processing further
	# Get the annotation dataframe by merging the sampleTable and vsd normalized counts
	annDF           <- merge(newSampleTable,t(vsd), by = "row.names")
	rownames(annDF) <- annDF[,1]
	annDF[,1]       <- NULL
	# Moeve diagnosis column to the end
	# https://stackoverflow.com/questions/18339370/reordering-columns-in-a-large-dataframe
	annDF <- annDF[c(setdiff(names(annDF), c("diagnosis")), c("diagnosis"))]
	vstAnnOutFile   <- paste(vstnormalizedDir, "/", basefilename, "_annotation_vst_normalized.txt", sep='')
	write.table(annDF, file = vstAnnOutFile, row.names = F, sep = '\t', quote = F)

	# Convert rownames as feature column to save it in the csv file
	vsd <- cbind(feature = rownames(vsd), vsd)
	rownames(vsd) <- 1:nrow(vsd)

	# Save the counts file to output csv file
	write.table(vsd   , file = vstnormalizedOutFile , row.names = F, sep = '\t', quote = F)
	write.table(t(vsd), file = vsttnormalizedOutFile, col.names = F, sep = '\t', quote = F)

	# Get the normalized individual count files
	indvcmd <- paste("bash scripts/get_individual_countFiles.sh", vstnormalizedOutFile, indvstnormalizedDir, rnaClass, 1, 2, sep=' ')
	system(indvcmd)
}

# Run the main function
main()
