# mutantAligner.R -- turn TRIP Screen Mutant Strain NGS data into relative abundance data to model growth defects

# most recent updates:  Sept 2019 - package up as stand alone scripts to post to ISB's TB resource
#			May 2022  - version for Magda's RIF data
#			Aug 2022  - merge in some improvements from Shuyi's version at SCRI
#						changing notation from old "TFOE" acronym to "TRIP" like the publication
#						https://www.nature.com/articles/s41564-020-00810-x
#						changing plot formats to PDF
#			Nov 2022  - removing support for Bowtie1, cleaning & trimming for GitHub upload
#			Dec 2022  - adjust paths to let these scripts be called from any user path
#			Dec 2022  - making the 'sharedDayZero' argument more flexible
#			Mar 2023  - adding normalization by uninduced read abundance


# if not done already by a caller, define the path to the TRIP_Screen installion 
if ( ! exists( "TRIP.PATH")) TRIP.PATH <- "~/TRIP_Screen"	# set this as needed 

library( DuffyNGS)
require( Biostrings)
multicore.setup(1)	# for debugging, only use 1 core.  Once stable, up it as wanted
setCurrentSpecies( "MT_H37")


# the Bowtie alignment executable.  Using only Bowtie2 as of 2022
BOWTIE2_PROGRAM <- Sys.which( "bowtie2")
if ( ! nchar( BOWTIE2_PROGRAM)) stop( "Unable to find 'bowtie2' on the search path")

# locations of the genones and target indexes, and the set of mutants
TRIP.BowtieIndex <- file.path( TRIP.PATH, "BowtieIndexes/MTb.TRIP.MutantGenome")
TRIP.GenomeFile <- file.path( TRIP.PATH, "CustomMutantGenome/MTb.TRIP.MutantGenome.fasta")
TRIP.CompleteGeneFile <- file.path( TRIP.PATH, "CustomMutantGenome/TRIP.Complete207.Primers.txt")

# folder layout where to read from and write to
ResultsFolder <- "./Results"

# the required fields in the sample key file
SampleKeyColumns <- c("SampleID", "Experiment", "Day", "Replicate", "ATc_Induced", "Class",
				"SubClass1", "SubClass2", "MutantPool", "File1", "File2", "FastqPath")
MetaDataColumns <- SampleKeyColumns[ 1:9]
ModelMetaDataColumns <- c( MetaDataColumns, "Run", "ReadCount")

# fix what day numbers will be used as the baseline timepoint for modeling growth defects over time.
DAY_ZERO <- c( 0, 1)

# allow setting a maximum number of raw reads to align per sample
MAX_RAW_READS <- 2000000

# crop low abunance values for math stability.  Samples below this threshold excluded from modelling.
# RPMHEG = Reads Per Million per Hundred Expected Genes  (i.e number of mutants in the pool)
MIN_LOG2_RPMHEG <- 2
MIN_READ_PAIRS <- 20000

checkX11( bg="white", width=10, height=7)
par( mai=c(0.6,0.4,0.6, 0.4))

# allow a TRIM3 field in the Sample Key file to override the command line argument,
# so we can process many different runs at once.  The mate pair nature of this data,
# where we use local alignment to find the anchor in one mate and the TF gene in the 
# other, is based on using relatively short read lengths, like 50bp to 75bp.  Very
# long reads tend to see both anchor and TF in both mate reads, which breaks the logic
# used for finding pairs.  Therefore long reads must be trimmed from the 3' end to yield
# the expected shorter 50-75bp lenghts.


# top level function for the standard TRIP Screen:  (Transcriptional Regulator Induced Phenotype) 
# Trys to run all 5 main steps of the workflow:
#   1. Bowtie alignment of each mate pair file separately
#   2. Count up and categorize pairs, generate Pair Pie plots
#   3. Calculate RPMHEG abundance
#   4. Write out count data for all samples, generate alignment summary stats and overview plot
#   5. If all above sets succeed, then Run the final growth decay modeling step

# By default, the pipeline will use previously generate data results from previous steps, and only generate 
# those data if they are missing.  If you want to force recalculation, use the 'doXXXX' flags to override these defaults.

# Most problems are due to mistakes or omissions in the 'SampleKey' file.  It can take some effort to organize 
# the file to match what the pipeline expects.

# Sharing Day Zero samples:  typically there will not be baseline samples for every subset
# of samples to be modeled.  The 'sharedDayZero' argument is provided to let you specify a vector of
# character strings that are names from the 'Class/SubClass1/SubClass2' columns, to help the tool decide which day zero
# samples can be used by all other classes of experiments/drugs, etc.   Not very elegant, can be
# adjusted as needed for future experiments.


do.all <- function( sampleKeyFile, doBowtie=FALSE, dropZeroGenes=TRUE, 
				verbose=FALSE, experiment=NULL, doAlignStats=doBowtie, doPairs=FALSE, 
				makePies=FALSE, makeROC=FALSE, doMODEL=TRUE, sharedDayZero=NULL,
				trim5=0, trim3=0, max.reads=MAX_RAW_READS) {
	
	# set the needed Genome & Target index for TRIP
	BowtieIndex <<- TRIP.BowtieIndex
	GenomeFile <<- TRIP.GenomeFile
	CompleteGeneFile <<- TRIP.CompleteGeneFile

	# get all the Rv Gene names from the complete genome
	genome <- loadFasta( GenomeFile, verbose=F)
	genes <- grep( "^Rv|Empty", genome$desc, value=T)
	nGenes <- length(genes)
	cat( "\nN_Genes in Mutant Genome file:    ", nGenes)

	# load that sample key
	allSamples <- loadSampleKeyFile( sampleKeyFile)
	nSamples <- nrow( allSamples)

	# allow specifying just one experiment at a time, for testing, etc.
	if ( ! is.null( experiment)) {
		allSamples <- subset( allSamples, Experiment == experiment)
		nSamples <- nrow( allSamples)
	}

	# the 'Run Name' is extracted from the name of the Samples Key filename for the results folder name
	runName <- extractRunName( sampleKeyFile)
	results.path <- createResultsFolder( sampleKeyFile)

	# some storage is for the entire run... even when more than one experiment is in the Sample Key file
	outStats <<- matrix( nrow=nSamples, ncol=8)
	rownames(outStats) <<- allSamples$SampleID
	colnames(outStats) <<- c( "StartingReads", "ValidMate1", "Pct_ValidMate1", "ValidMate2", 
				"Pct_ValidMate2", "ValidPairs", "Pct_ValidPairs", "N_Genes")

	# allow more than one 'Experiment' per Sample Key file, and process each one separately
	expFactor <- factor( allSamples$Experiment)
	tapply( 1:nrow(allSamples), expFactor, function(x) {

		# do the subset of samples for this experiment
		samples <- allSamples[x, ]
		nSamples <- nrow( samples)

		# get the details about this experiment
		myExperiment <- samples$Experiment[1]
		cat( "\n\nProcessing Experiment:   ", myExperiment)

		# make a big matrix to hold all the data
		out <- matrix( 0, nrow=nGenes, ncol=nSamples)
		colnames(out) <- samples$SampleID
		rownames(out) <- genes
		poolNames <- nExpectedGenes <- vector( length=nSamples)
		#poolLists <- vector( mode="list")

		# do the alignments on multiple cores at once, if set to allow
		# see:  multicore.setup()
		doOneMutantAlign <- function( i) {
			sid <- samples$SampleID[i]
			# allow a TRIM3 field in the Sample Key file to override the command line argument,
			if ( "TRIM3" %in% colnames(samples)) {
				myTrim3 <- as.integer( samples$TRIM3[i])
				if ( ! is.na( myTrim3)) trim3 <- myTrim3
			}

			mutantAlign( sid, samples$File1[i], samples$File2[i], fastqPath=samples$FastqPath[i],
					doBowtie=doBowtie, results.path=results.path, 
					trim5=trim5, trim3=trim3, max.reads=max.reads, verbose=verbose)
		}
		multicore.lapply( 1:nSamples, FUN=doOneMutantAlign)
		
		# make folders of "AlignStats" images that summarize Bowtie success rates
		if (doBowtie && doAlignStats) {
			doMutantAlignStats( samples$SampleID, results.path=results.path)
		}

		# next, let's tabulate all the read pairs, for each sample on multiple cores at once...
		multicore.lapply( samples$SampleID, FUN=tabulateMutantAlignments.TRIP, doBowtie=doBowtie, 
				doPairs=doPairs, results.path=results.path, max.reads=max.reads)

		# now we can do the final summary and pie plots.  This step must be done live,
		# can not be taken from previous stored results...
		for ( i in 1:nSamples) {
			sid <- samples$SampleID[i]
			# stick the alignment results into the overall run table, not just for this one experiment
			rowptr <- x[i]
			ans <- NULL
			ans <- summarizeMutantAlignments.TRIP( sid, rowptr=rowptr, results.path=results.path, makePie=makePies)
			if ( is.null(ans)) next

			# put these counts where they go for this experiment
			where <- match( names(ans), genes, nomatch=0)
			out[ where, i] <- ans[ where > 0]

			# get the details about this pool
			myPool <- samples$MutantPool[i]
			poolGenes <- getMutantPoolGenes( myPool)
			nGenesThisPool <- length( poolGenes)
			poolNames[i] <- myPool
			nExpectedGenes[i] <- length(poolGenes)
		}

		# geneate ROC curves about the distribution of read counts
		if (makeROC) {
			rocAns <- makeROCimage( out, poolGenes=poolGenes, text="Raw Counts", 
							experiment=myExperiment, results.path=results.path)
			ROC.rawCount.cutoff <- rocAns$cutpoint
		} else {
			ROC.rawCount.cutoff <- 10
		}

		sumPerGene <- apply( out, MARGIN=1, FUN=sum)
		nDetect <- sum( sumPerGene > ROC.rawCount.cutoff)
		cat( "\n\nN_Detected Genes over", nSamples, "samples:  ", nDetect, "\tROC cut: ", ROC.rawCount.cutoff)
		if ( nDetect < nGenes && dropZeroGenes) {
			toDrop <- which( sumPerGene == 0)
			out <- out[ -toDrop, ]
			cat( "\nDropped genes with zero reads:  ", length( toDrop))
		}

		# normalize to RPMHEG units and write the various result files
		writeResultTables( out, samples, results.path=results.path, nExpectedGenes=nExpectedGenes, 
					poolGenes=poolGenes, experiment=myExperiment, makeROC=makeROC, runName=runName) 

		cat( "\nFinished Experiment:  ", myExperiment, "\n")

	})  ## done with tapply() of all experiments in the run

	# write out the final alignment stats overview table too..
	os <- data.frame( "SampleID"=rownames(outStats), outStats, stringsAsFactors=F)
	outfile <- file.path( results.path, "AlignmentStatistics.txt")
	write.table( os, outfile, sep="\t", quote=F, row.names=F)

	plotAlignmentOverview( outStats, allSamples, results.path=results.path)

	# when we are all done, try to run the growth defect code as a last step
	if ( doMODEL) {
		model.All.TRIP.Samples( sampleKeyFile, sharedDayZero=sharedDayZero, experiment=experiment, normalized=FALSE)
		model.All.TRIP.Samples( sampleKeyFile, sharedDayZero=sharedDayZero, experiment=experiment, normalized=TRUE)
	}
	
	cat( "\nDone.\n")
}


# function to do Bowtie alignment for both mate files of one sample
mutantAlign <- function( sid, file1, file2, fastqPath=".", doBowtie=TRUE, 
				results.path=".", trim5=0, trim3=0, max.reads=MAX_RAW_READS, verbose=TRUE) {

	# verify FASTQ files are found
	file1 <- file.path( fastqPath, file1)
	file2 <- file.path( fastqPath, file2)
	if ( ! file.exists( file1)) {
		cat("\n\nError:  Mate 1 Fastq file not found: ", sid, "\tFile: ", file1)
		stop()
	}
	if ( ! file.exists( file2)) {
		cat("\n\nError:  Mate 2 Fastq file not found: ", sid, "\tFile: ", file2)
		stop()
	}

	bamPath <- file.path( results.path, "BAM.Files")
	if ( ! file.exists( bamPath)) dir.create( bamPath, recursive=T)
	outfile1 <- file.path( bamPath, paste( sid, 1, "bam", sep="."))
	outfile2 <- file.path( bamPath, paste( sid, 2, "bam", sep="."))

	# allow if the pair file results already exist, to let us skip doing the alignments
	pairPath <- file.path( results.path, "PAIR.Files")
	pairFile <- file.path( pairPath, paste( sid, "ReadPairs.rda", sep="."))
	if ( ! doBowtie && file.exists( pairFile)) {
		cat( "\nSkipping Alignment step.  PAIR file already exist for: ", sid)
		return()
	}
	if ( ! doBowtie && all( file.exists( c( outfile1, outfile2)))) {
		cat( "\nSkipping Alignment step.  BAM files already exist for: ", sid)
		return()
	}

	cat( "\n\n\n-----------------------------------------------------\n\nCalling BOWTIE2 alignment on mutant sample:  ", sid)

	# build the Unix command line
	trimOptions <- paste( " --qupto", as.character( as.integer( max.reads)))
	if (trim5 > 0) trimOptions <- paste( trimOptions, " --trim5", trim5)
	if (trim3 > 0) trimOptions <- paste( trimOptions, " --trim3", trim3)

	program <- Sys.which( "bowtie2")

	# our standard scoring settings for doing "local" alignment
	options <- paste( " -q  --local  --score-min G,20,9  --very-sensitive-local  --threads 4 --ignore-quals ",
					trimOptions, " -x ")
	fileFlag <- " -U "

	# catch the statistics so we can echo them
	catchStdErrFile <- paste( sid, "bowtieStdErr.txt", sep=".")
	catchStdErrLog <- paste( "  2> ", catchStdErrFile)

	cmdline <- paste( program, options, BowtieIndex, fileFlag, file1, catchStdErrLog, " | samtools view -b - >", outfile1)
	cat( "\n\nDoing Mate 1: ", sid, "\n")
	if (verbose) cat( "\nBowtie command line: \n", cmdline, "\n")
	file.delete( catchStdErrFile)
	catch.system( cmdline)
	lines <- readLines( catchStdErrFile)
	writeLines( lines)

	cmdline <- paste( program, options, BowtieIndex, fileFlag, file2, catchStdErrLog, " | samtools view -b - >", outfile2)
	cat( "\nDoing Mate 2: ", sid, "\n")
	if (verbose) cat(  "\nBowtie command line: \n", cmdline, "\n")
	file.delete( catchStdErrFile)
	catch.system( cmdline)
	lines <- readLines( catchStdErrFile)
	writeLines( lines)

	file.delete( catchStdErrFile)
	if (verbose) cat( "\nAlignment Done:  ", sid, "\n")
	return()
}


# function to make Bowtie alignment statistics plots for both mate files of one sample
doMutantAlignStats <- function( sidSet, results.path=".") {

	bamPath <- file.path( results.path, "BAM.Files")

	# turn off any multicore for now
	ncores <- multicore.currentCoreCount()
	on.exit( multicore.setup(ncores))
	multicore.setup(1)
	
	for ( sid in sidSet) {
		# verify BAM files are found
		bamfiles <- file.path( bamPath, paste( sid, 1:2, "bam", sep="."))

		# do each mate separately
		for ( mate in 1:2) {
			thisID <- paste( sid, "_Mate", mate, sep="")
			cat( "\nPlotting alignment stats for: ", thisID)
			stats.path <- file.path( results.path, "AlignStats", thisID)
			if ( ! file.exists( stats.path)) dir.create( stats.path, recursive=TRUE)
			bamf <- bamfiles[mate]
			if ( ! file.exists( bamf)) {
				cat( "\nError:  BAM file not found:  ", bamf)
				next
			}

			# make those plots
			calcAlignStats( bamf, sampleID=thisID, statsPath=stats.path, what="SBIDMA", 
						chunkSize=500000, verbose=FALSE)
		}
	}

	# the align stats plotter closes the plot device, restore it now
	checkX11( bg="white", width=10, height=7)
	par( mai=c(0.6,0.4,0.6, 0.4))
}


# function to scan the 2 mate pair BAM files of TRIP data and tabulate the results by combining 
# the mate pairs that have the desired "Anchor - Gene" expected pairing
tabulateMutantAlignments.TRIP <- function( sid, results.path=".", makePie=TRUE, doBowtie=TRUE, 
				doPairs=TRUE, max.reads=MAX_RAW_READS, trackNoHits=TRUE) {

	# keep a subset of the most seen unaligned reads, as a debugging aid
	N_NOHIT_KEEP <- 100

	# where the BAM files will be found
	bamPath <- file.path( results.path, "BAM.Files")
	infile1 <- paste( sid, 1, "bam", sep=".")
	infile2 <- paste( sid, 2, "bam", sep=".")
	infile1 <- file.path( bamPath, infile1)
	infile2 <- file.path( bamPath, infile2)

	# where the results will be put
	pairPath <- file.path( results.path, "PAIR.Files")
	if ( ! file.exists( pairPath)) dir.create( pairPath, recursive=T)
	pairFile <- file.path( pairPath, paste( sid, "ReadPairs.rda", sep="."))

	# let's try to find/call empty adapters too
	adapter1 <- "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
	adapter2 <- "CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT"
	rcadapter1 <- myReverseComplement( adapter1)
	rcadapter2 <- myReverseComplement( adapter2)
	adapterMinScore <- 30

	# do the tabulation whenever Bowtie was run or the answer not already made
	doTabulate <- ( doBowtie || doPairs || !file.exists( pairFile))
	if ( doTabulate) {

		cat( "\n\nReading BAM alignment results for sample:  ", sid, "  ")
		# set up to iterate thru the two files
		reader1 <- bamReader( infile1)
		reader2 <- bamReader( infile2)
		refData <- getRefData( reader1)
	
		hasMore <- TRUE
		big1 <- big2 <- bigID <- vector()
		big1pos <- big2pos <- vector()
		nohitTbl1 <- nohitTbl2 <- NULL
		nReadIn <- 0
		CHUNK_SIZE <- 100000

		while (hasMore) {
			if (nReadIn >= max.reads) break
			chunk1 <- getNextChunk( reader1, n=CHUNK_SIZE)
			chunk2 <- getNextChunk( reader2, n=CHUNK_SIZE)
			nNow <- size( chunk1)
			if ( nNow < CHUNK_SIZE) hasMore <- FALSE
			if ( nNow < 1) break
			nReadIn <- nReadIn + nNow
	
			readIDs1 <- readID( chunk1)
			readIDs2 <- readID( chunk2)
			# some Sequencing cores will leave mate suffixes...  strip them off
			readIDs1 <- sub( "/1$", "", readIDs1)
			readIDs2 <- sub( "/2$", "", readIDs2)
			refIDs1 <- refID( chunk1)
			refIDs2 <- refID( chunk2)
			seqIDs1 <- refID2seqID( refIDs1, refData=refData)
			seqIDs2 <- refID2seqID( refIDs2, refData=refData)
			seqIDs1[ is.na( seqIDs1)] <- "NoHit"
			seqIDs2[ is.na( seqIDs2)] <- "NoHit"
			pos1 <- position( chunk1)
			pos2 <- position( chunk2)
			isNoHit1 <- which( seqIDs1 == "NoHit")
			if ( length( isNoHit1)) {
				seqDNA1 <- readSeq(chunk1)
				# check for empty adaptor ends
				scores <- pairwiseAlignment( seqDNA1[isNoHit1], rcadapter2, type='local', scoreOnly=TRUE)
				reLabel <- isNoHit1[ scores >= adapterMinScore]
				if (length(reLabel)) seqIDs1[ reLabel] <- "EmptyAdapter"
				if (trackNoHits) {
					smlTbl <- sort( table( seqDNA1[ isNoHit1]), decreasing=T)
					Nkeep <- min( N_NOHIT_KEEP, length(smlTbl))
					nohitTbl1 <- if (is.null(nohitTbl1)) smlTbl[ 1:Nkeep] else mergeTables( nohitTbl1, smlTbl[1:Nkeep])
				}
				# check for Poly N
				isNoHit1 <- setdiff( isNoHit1, reLabel)
				bases <- strsplit( seqDNA1[isNoHit1], split="")
				nPolyN <- sapply( bases, function(x) sum( x == "N"))
				isPolyN <- which( nPolyN > 5)
				if (length(isPolyN)) seqIDs1[ isNoHit1[isPolyN]] <- "PolyN"
			}
			isNoHit2 <- which( seqIDs2 == "NoHit")
			if ( length( isNoHit2)) {
				seqDNA2 <- readSeq(chunk2)
				# check for empty adaptor ends
				scores <- pairwiseAlignment( seqDNA2[isNoHit2], rcadapter1, type='local', scoreOnly=TRUE)
				reLabel <- isNoHit2[ scores >= adapterMinScore]
				if (length(reLabel)) seqIDs2[ reLabel] <- "EmptyAdapter"
				if (trackNoHits) {
					smlTbl <- sort( table( seqDNA2[ isNoHit2]), decreasing=T)
					Nkeep <- min( N_NOHIT_KEEP, length(smlTbl))
					nohitTbl2 <- if (is.null(nohitTbl2)) smlTbl[ 1:Nkeep] else mergeTables( nohitTbl2, smlTbl[1:Nkeep])
				}
				# check for Poly N
				isNoHit2 <- setdiff( isNoHit2, reLabel)
				bases <- strsplit( seqDNA2[isNoHit2], split="")
				nPolyN <- sapply( bases, function(x) sum( x == "N"))
				isPolyN <- which( nPolyN > 5)
				if (length(isPolyN)) seqIDs2[ isNoHit2[isPolyN]] <- "PolyN"
			}
	
			# force the mate pairs to match by their IDs
			if ( ! all( readIDs1 == readIDs2)) {
				both <- intersect( readIDs1, readIDs2)
				wh1 <- match( both, readIDs1)
				wh2 <- match( both, readIDs2)
				refIDs1 <- refIDs1[wh1]
				refIDs2 <- refIDs2[wh2]
				seqIDs1 <- seqIDs1[wh1]
				seqIDs2 <- seqIDs2[wh2]
				pos1 <- pos1[wh1]
				pos2 <- pos2[wh2]
				readIDs1 <- readIDs2 <- both
			}
			# sanity check that we have intact pairs
			if ( length( readIDs1) < (nNow/2)) stop( "Error:  Problem pairing up Mate Read IDs.")
			cat( ".")
			big1 <- c( big1, seqIDs1)
			big2 <- c( big2, seqIDs2)
			bigID <- c( bigID, readIDs1)
			big1pos <- c( big1pos, pos1)
			big2pos <- c( big2pos, pos2)
		}
	
		bamClose( reader1)
		bamClose( reader2)
		cat( "Done.\n")

		cat( "\nForce mate pairs into alphabetical order..")
		# make the pairs be in alpha order
		# we are using the fact that "Anchor" comes before all the "Rv" genes when we sort by alphabetical order
		# "Anchor" will be in 'big1' and the Rv gene will be in 'big2'
		Npairs <- length(big1)
		isAlpha <- (big1 < big2)
		isOne <- rep.int( 1, Npairs)
		isTwo <- rep.int( 2, Npairs)
		if ( any( ! isAlpha)) {
			toSwap <- which( ! isAlpha)
			isOne[toSwap] <- 2
			isTwo[toSwap] <- 1
		}
			
		# each entry should be one Anchor and one Gene
		cat( "\nAssess pairing success..")
		text1 <- text2 <- rep.int( "GENE", Npairs)
		text1[ grep( "Anchor", big1, fixed=T)] <- "ANCHOR"
		text2[ grep( "Anchor", big2, fixed=T)] <- "ANCHOR"
		text1[ grep( "NoHit", big1, fixed=T)] <- "NOHIT"
		text2[ grep( "NoHit", big2, fixed=T)] <- "NOHIT"
		text1[ grep( "EmptyAdapter", big1, fixed=T)] <- "EMPTY_ADAPT"
		text2[ grep( "EmptyAdapter", big2, fixed=T)] <- "EMPTY_ADAPT"
		text1[ grep( "PolyN", big1, fixed=T)] <- "POLY_N"
		text2[ grep( "PolyN", big2, fixed=T)] <- "POLY_N"
		bigStr <- paste( text1, text2, sep=".")

		pairDF <- data.frame( "ReadID"=bigID, "Mate1"=big1, "Mate2"=big2, "Status"=bigStr, 
				"Pos1"=big1pos, "Pos2"=big2pos, stringsAsFactors=FALSE)
		save( pairDF, file=pairFile)
		cat( "\nWrote Read Pair results:  ", sid, "\tN_Pairs: ", Npairs)

		if (trackNoHits && !is.null( nohitTbl1)) {
			nohitFile <- file.path( pairPath, paste( sid, "TopNoHits.Mate1.fasta", sep="."))
			nohitTbl1 <- sort( nohitTbl1, decreasing=T)
			if ( length(nohitTbl1) > N_NOHIT_KEEP) nohitTbl1 <- nohitTbl1[ 1:N_NOHIT_KEEP]
			NFA <- length(nohitTbl1)
			writeFasta( as.Fasta( desc=paste( "Mate1_", 1:NFA, "_Count=", nohitTbl1, sep=""),
					seq=names(nohitTbl1)), nohitFile, line.width=200)
		}
		if (trackNoHits && !is.null( nohitTbl2)) {
			nohitFile <- file.path( pairPath, paste( sid, "TopNoHits.Mate2.fasta", sep="."))
			nohitTbl2 <- sort( nohitTbl2, decreasing=T)
			if ( length(nohitTbl2) > N_NOHIT_KEEP) nohitTbl2 <- nohitTbl2[ 1:N_NOHIT_KEEP]
			NFA <- length(nohitTbl2)
			writeFasta( as.Fasta( desc=paste( "Mate2_", 1:NFA, "_Count=", nohitTbl2, sep=""),
					seq=names(nohitTbl2)), nohitFile, line.width=200)
		}

	} else {
		cat( "\nSkipping Pair Tabulation step.  PAIR file already exists for: ", sid)
	}
}


# function to summarize the TRIP pairings, and make a pie chart
summarizeMutantAlignments.TRIP <- function( sid, rowptr=0, results.path=".", makePie=TRUE) {

	# where the pair results will be found
	pairPath <- file.path( results.path, "PAIR.Files")
	pairFile <- file.path( pairPath, paste( sid, "ReadPairs.rda", sep="."))
	if ( ! file.exists( pairFile)) {
		cat( "\niRead Pair File not found: ", sid, "\t", pairFile)
		return( NULL)
	}
	cat( "\nSummarizing Read Pairs for:  ", sid)

	# read those pair results in
	load( pairFile)
	Npairs <- nrow(pairDF)

	# get all the various groups and their counts
	statusTypes <- c( "ANCHOR", "EMPTY_ADAPT", "GENE", "NOHIT", "POLY_N")
	Ntypes <- length( statusTypes)
	statusLevels <- paste( rep( statusTypes, each=Ntypes), rep( statusTypes, times=Ntypes), sep=".")
	statusFac <- factor( pairDF$Status, levels=statusLevels)
	statusSets <- tapply( 1:Npairs, statusFac, function(x) return(x))
	isDoubleAnch <- statusSets[[1]]
	isAnchEmpty <- statusSets[[2]]
	isAnchGene <- statusSets[[3]]
	isAnchNo <- statusSets[[4]]
	isAnchPolyn <- statusSets[[5]]
	isEmptyAnch <- statusSets[[6]]
	isDoubleEmpty <- statusSets[[7]]
	isEmptyGene <- statusSets[[8]]
	isEmptyNo <- statusSets[[9]]
	isEmptyPolyn <- statusSets[[10]]
	isGeneAnch <- statusSets[[11]]
	isGeneEmpty <- statusSets[[12]]
	isDoubleGene <- statusSets[[13]]
	isGeneNo <- statusSets[[14]]
	isGenePolyn <- statusSets[[15]]
	isNoAnch <- statusSets[[16]]
	isNoEmpty <- statusSets[[17]]
	isNoGene <- statusSets[[18]]
	isDoubleNo <- statusSets[[19]]
	isNoPolyn <- statusSets[[20]]
	isPolynAnch <- statusSets[[21]]
	isPolynEmpty <- statusSets[[22]]
	isPolynGene <- statusSets[[23]]
	isPolynNo <- statusSets[[24]]
	isDoublePolyn <- statusSets[[25]]

	nAnch1 <- length( isAnch1 <- c( isDoubleAnch, isAnchEmpty, isAnchGene, isAnchNo, isAnchPolyn))
	nGene1 <- length( isGene1 <- c( isGeneAnch, isGeneEmpty, isDoubleGene, isGeneNo, isGenePolyn))
	nEmpty1 <- length( isEmpty1 <- c( isEmptyAnch, isDoubleEmpty, isEmptyGene, isEmptyNo, isEmptyPolyn))
	nNoHit1 <- length( isNoHit1 <- c( isNoAnch, isNoEmpty, isNoGene, isDoubleNo, isNoPolyn))
	nPolyn1 <- length( isPolyn1 <- c( isPolynAnch, isPolynEmpty, isPolynGene, isPolynNo, isDoublePolyn))
	nAnch2 <- length( isAnch2 <- c( isDoubleAnch, isEmptyAnch, isGeneAnch, isNoAnch, isPolynAnch))
	nGene2 <- length( isGene2 <- c( isAnchGene, isEmptyGene, isDoubleGene, isNoGene, isPolynGene))
	nEmpty2 <- length( isEmpty2 <- c( isAnchEmpty, isDoubleEmpty, isGeneEmpty, isNoEmpty, isPolynEmpty))
	nNoHit2 <- length( isNoHit2 <- c( isAnchNo, isEmptyNo, isGeneNo, isDoubleNo, isPolynNo))
	nPolyn2 <- length( isPolyn2 <- c( isAnchPolyn, isEmptyPolyn, isGenePolyn, isNoPolyn, isDoublePolyn))

	# make a visual of the distribution
	piePath <- file.path( results.path, "PIE.Plots")
	if ( ! file.exists( piePath)) dir.create( piePath, recursive=T)
	pieFile <- file.path( piePath, paste( sid, "MatePair.Profile.pdf", sep="."))
	if ( makePie || ! file.exists(pieFile)) {
		cnts <- c( length(isGeneAnch), length(isAnchGene), length(isAnchNo), length( isNoAnch), length(isGeneNo), length( isNoGene),
				length(isDoubleAnch), length(isDoubleGene), length(isDoubleNo), length( isDoubleEmpty), 
				length(isAnchEmpty), length( isEmptyAnch), length( isGeneEmpty), length( isEmptyGene),
				length( isEmptyNo), length( isNoEmpty), length(isAnchPolyn), length(isPolynAnch), length( isGenePolyn),
				length( isPolynGene), length(isEmptyPolyn), length( isPolynEmpty), length( isNoPolyn),
				length( isPolynNo), length( isDoublePolyn))
		nams <- c( "Gene - Anchor", "Anchor - Gene", "Anchor - NoHit", "NoHit - Anchor", "Gene - NoHit", "NoHit - Gene", 
				"Anchor - Anchor", "Gene - Gene", "NoHit - NoHit", "EmptyAdapt - EmptyAdapt",
				"Anchor - EmptyAdapt", "EmptyAdapt - Anchor", "Gene - EmptyAdapt", "EmptyAdapt - Gene",
				"EmptyAdapt - NoHit", "NoHit - EmptyAdapt", "Anchor - PolyN", "PolyN - Anchor", "Gene - PolyN",
				"PolyN - Gene", "EmptyAdapt - PolyN", "PolyN - EmptyAdapt", "NoHit - PolyN", "PolyN - NoHit",
				"PolyN - PolyN")
		pcts <- cnts * 100 / sum(cnts)
		strPcts <- as.percent( cnts, big.value=sum(cnts))
		names( cnts) <- paste( nams, "  (", strPcts, ")", sep="")
		cols <- c( 3, 3, "gold", "gold", "orange", "orange", 
				"yellow", "pink", 2, 'purple', 
				'blue', "blue", 'dodgerblue', 'dodgerblue', 'slateblue', 'slateblue',
				rep.int( 'sandybrown',4), rep.int( 'saddlebrown',4), rep.int( 'brown',3))
		toShow <- which( pcts >= 0.75)

		pie( cnts[toShow], col=cols[toShow], clock=F, 
				main=paste( "Read Pair Profile:   ", sid, "\nTotal Reads:   ", 
				formatC( Npairs, format="d", big.mark=",")))

		# try to show each mates overall status as color bars along lower edge
		usrLims <- par( 'usr')
		xlo <- usrLims[1]
		xhi <- usrLims[2]
		xwid <- xhi - xlo
		ylo <- usrLims[3]
		yhi <- usrLims[4]
		ywid <- yhi - ylo
		ytop <- ylo + ywid * 0.05
		ymid <- (ylo + ytop) / 2
		xright <- xlo + xwid * 0.33
		rect( xlo, ylo, xright, ytop, border=1, col='red')
		text( xright, ymid, "NoHit", pos=2, offset=0.2, cex=0.9)
		xright <- xlo + (xwid*0.33) * ((nAnch1+nGene1+nEmpty1+nPolyn1)/Npairs)
		rect( xlo, ylo, xright, ytop, border=1, col='brown')
		text( xright, ymid, "PolyN", pos=2, offset=0.2, cex=0.9)
		xright <- xlo + (xwid*0.33) * ((nAnch1+nGene1+nEmpty1)/Npairs)
		rect( xlo, ylo, xright, ytop, border=1, col='purple')
		text( xright, ymid, "Adaptor", pos=2, offset=0.2, cex=0.9)
		xright <- xlo + (xwid*0.33) * ((nAnch1+nGene1)/Npairs)
		rect( xlo, ylo, xright, ytop, border=1, col='green')
		text( xright, ymid, "Gene", pos=2, offset=0.2, cex=0.9)
		xright <- xlo + (xwid*0.33) * ((nAnch1)/Npairs)
		rect( xlo, ylo, xright, ytop, border=1, col='yellow')
		text( xright, ymid, "Anchor", pos=2, offset=0.2, cex=0.9)
		xleft <- xhi - xwid * 0.33
		rect( xleft, ylo, xhi, ytop, border=1, col='red')
		text( xleft, ymid, "NoHit", pos=4, offset=0.2, cex=0.9)
		xleft <- xhi - (xwid*0.33) * ((nAnch2+nGene2+nEmpty2+nPolyn2)/Npairs)
		rect( xleft, ylo, xhi, ytop, border=1, col='brown')
		text( xleft, ymid, "PolyN", pos=4, offset=0.2, cex=0.9)
		xleft <- xhi - (xwid*0.33) * ((nAnch2+nGene2+nEmpty2)/Npairs)
		rect( xleft, ylo, xhi, ytop, border=1, col='purple')
		text( xleft, ymid, "Adaptor", pos=4, offset=0.2, cex=0.9)
		xleft <- xhi - (xwid*0.33) * ((nAnch2+nGene2)/Npairs)
		rect( xleft, ylo, xhi, ytop, border=1, col='green')
		text( xleft, ymid, "Gene", pos=4, offset=0.2, cex=0.9)
		xleft <- xhi - (xwid*0.33) * ((nAnch2)/Npairs)
		rect( xleft, ylo, xhi, ytop, border=1, col='yellow')
		text( xleft, ymid, "Anchor", pos=4, offset=0.2, cex=0.9)
		text( xlo+xwid*0.16, ytop, "Mate 1 Status", pos=3, cex=1, font=2)
		text( xhi-xwid*0.16, ytop, "Mate 2 Status", pos=3, cex=1, font=2)
		subText <- paste( "Valid Read Pairs:   ", formatC( length( c( isGeneAnch, isAnchGene)), 
				format="d", big.mark=","))
		text( 0, ylo, subText, font=2, cex=1.2, pos=3)
		dev.flush()
		dev.print( pdf, pieFile, width=14, height=9)
	}

	# final result is the counts of the GeneIDs for the good rows
	geneNames <- c( pairDF$Mate1[ isGeneAnch], pairDF$Mate2[ isAnchGene])
	out <- table( geneNames)
	cat( "\nN_BAM mate pairs considered:             ", Npairs)
	cat( "\nTotal Good 'Anchor - Gene' mate pairs:   ", sum(out))
	cat( "\nTotal Good mate pairs as percentage:    ", as.percent( sum(out), big.value=Npairs))
	cat( "\nN_Detected Genes:                        ", length(out))

	# save some overall metrics
	valid1 <- nAnch1 + nGene1
	valid2 <- nAnch2 + nGene2
	if ( rowptr > 0 && exists( "outStats")) {
		outStats[ rowptr, 'StartingReads'] <<- Npairs
		outStats[ rowptr, 'ValidMate1'] <<- valid1
		outStats[ rowptr, 'Pct_ValidMate1'] <<- round( valid1 * 100 / Npairs, digits=2)
		outStats[ rowptr, 'ValidMate2'] <<- valid2
		outStats[ rowptr, 'Pct_ValidMate2'] <<- round( valid2 * 100 / Npairs, digits=2)
		outStats[ rowptr, 'ValidPairs'] <<- sum(out)
		outStats[ rowptr, 'Pct_ValidPairs'] <<- round( sum(out) * 100 / Npairs, digits=2)
		outStats[ rowptr, 'N_Genes'] <<- length(out)
	}
	return( out)
}


# take all the results from one experiment, and write out the gene abundance data tables
writeResultTables <- function( tbl, samples, results.path, nExpectedGenes, poolGenes, experiment, 
				makeROC=FALSE, runName="") {

	# turn these raw counts of all genes into various normalized RPMHEG and Log2 transformed units
	# then write them all to disk
	outRPM <- outLog2 <- outNRPM <- outNLog2 <- tbl
	expectRPM <- expectLog2 <- expectNRPM <- expectNLog2 <- expect <- tbl[ which( rownames(tbl) %in% poolGenes), , drop=F]
	outfile <- file.path( results.path, paste( "AllGenes.RawCounts", experiment, "txt", sep="."))
	outDF <- data.frame( "GeneID"=rownames(tbl), tbl, stringsAsFactors=F)
	write.table( outDF, outfile, sep="\t", quote=F, row.names=F)
	expectfile <- file.path( results.path, paste( "PoolGenes.RawCounts", experiment, "txt", sep="."))
	expectDF <- data.frame( "GeneID"=rownames(expect), expect, stringsAsFactors=F)
	write.table( expectDF, expectfile, sep="\t", quote=F, row.names=F)

	cat( "\nNormalizing to RPM_HEG units..\n")
	nSamples <- ncol(tbl)
	nGenes <- nrow(tbl)
	if ( nrow(samples) != nSamples) stop( "Sample / matrix sizing error!")
	for ( i in 1:nSamples) {
		v <- tbl[ , i]
		expectedGenesFactor <- nExpectedGenes[i] / 100
		vnorm <- v * expectedGenesFactor * 1000000 / max( sum(v), MIN_READ_PAIRS)
		outRPM[ , i] <- vnorm
		v <- expect[ , i]
		vnorm <- v * expectedGenesFactor * 1000000 / max( sum(v), MIN_READ_PAIRS)
		expectRPM[ , i] <- vnorm
		cat( "\r", colnames(tbl)[i])
	}
	outfile <- file.path( results.path, paste( "AllGenes.RPMHEG", experiment, "txt", sep="."))
	outDF <- data.frame( "GeneID"=rownames(outRPM), round(outRPM,digits=3), stringsAsFactors=F)
	write.table( outDF, outfile, sep="\t", quote=F, row.names=F)
	expectfile <- file.path( results.path, paste( "PoolGenes.RPMHEG", experiment, "txt", sep="."))
	expectDF <- data.frame( "GeneID"=rownames(expectRPM), round(expectRPM,digits=3), stringsAsFactors=F)
	write.table( expectDF, expectfile, sep="\t", quote=F, row.names=F)

	# do the normalization based on uninduced samples, to try to compensate for changing mutant proportions over time
	outNRPM <- normalizeByUninduced( outRPM, samples=samples, nExpectedGenes)
	expectNRPM <- normalizeByUninduced( expectRPM, samples=samples, nExpectedGenes)
	outfile <- file.path( results.path, paste( "AllGenes.Normalized.RPMHEG", experiment, "txt", sep="."))
	outDF <- data.frame( "GeneID"=rownames(outNRPM), round(outNRPM,digits=3), stringsAsFactors=F)
	write.table( outDF, outfile, sep="\t", quote=F, row.names=F)
	expectfile <- file.path( results.path, paste( "PoolGenes.Normalized.RPMHEG", experiment, "txt", sep="."))
	expectDF <- data.frame( "GeneID"=rownames(expectNRPM), round(expectNRPM,digits=3), stringsAsFactors=F)
	write.table( expectDF, expectfile, sep="\t", quote=F, row.names=F)
	
	
	cat( "\nTransforming to Log2 units..\n")
	for ( i in 1:nSamples) {
		v <- tbl[ , i]
		outLog2[ , i] <- log2( outRPM[ , i] + 1)
		expectLog2[ , i] <- log2( expectRPM[ , i] + 1)
		outNLog2[ , i] <- log2( outNRPM[ , i] + 1)
		expectNLog2[ , i] <- log2( expectNRPM[ , i] + 1)
		cat( "\r", colnames(tbl)[i])
	}
	if ( makeROC) {
		makeROCimage( outLog2, poolGenes=poolGenes, text="Log2 RPMHEG", experiment=experiment, 
				results.path=results.path)
		makeROCimage( outNLog2, poolGenes=poolGenes, text="Normalized Log2 RPMHEG", experiment=experiment, 
				results.path=results.path)
	}
	outfile <- file.path( results.path, paste( "AllGenes.Log2RPMHEG", experiment, "txt", sep="."))
	outDF <- data.frame( "GeneID"=rownames(outLog2), round(outLog2,digits=5), stringsAsFactors=F)
	write.table( outDF, outfile, sep="\t", quote=F, row.names=F)
	expectfile <- file.path( results.path, paste( "PoolGenes.Log2RPMHEG", experiment, "txt", sep="."))
	expectDF <- data.frame( "GeneID"=rownames(expectLog2), round(expectLog2,digits=5), stringsAsFactors=F)
	write.table( expectDF, expectfile, sep="\t", quote=F, row.names=F)
	outfile <- file.path( results.path, paste( "AllGenes.Normalized.Log2RPMHEG", experiment, "txt", sep="."))
	outDF <- data.frame( "GeneID"=rownames(outNLog2), round(outNLog2,digits=5), stringsAsFactors=F)
	write.table( outDF, outfile, sep="\t", quote=F, row.names=F)
	expectfile <- file.path( results.path, paste( "PoolGenes.Normalized.Log2RPMHEG", experiment, "txt", sep="."))
	expectDF <- data.frame( "GeneID"=rownames(expectNLog2), round(expectNLog2,digits=5), stringsAsFactors=F)
	write.table( expectDF, expectfile, sep="\t", quote=F, row.names=F)

	# lastly, make a table with all meta data and all genes for modeling etc.
	# there is potential for confusion here, about what genes to include in the model, if the
	# samples had 2+ different mutant pool sources.
	# Try to use the file mentioned in the sample key, as long as they are all the same.
	# If 2+, use the globally defined one, with a warning
	myPoolNames <- sort( unique( samples$MutantPool))
	if ( length(myPoolNames) > 1) {
		completeGeneTbl <- read.delim( CompleteGeneFile, as.is=T)
		completeGenes <- sort( unique( completeGeneTbl$GeneID))
		cat( "\n  Warning: Found 2+ mutant pool names, using: ", CompleteGeneFile, "\n  to define the universe of genes for modeling.\n")
	} else {
		completeGenes <- getMutantPoolGenes( myPoolNames)
	}
	outCounts <- outCounts2 <- matrix( NA, nrow=ncol(expectLog2), ncol=length(completeGenes))
	colnames(outCounts) <- colnames(outCounts2) <- completeGenes
	rownames(outCounts) <- rownames(outCounts2) <- samples$SampleID

	# first the log2 RPMHEG counts
	where <- match( rownames(expectLog2), completeGenes)
	if ( any( is.na(where))) {
		cat( "\n\nSome GeneIDs not in the 'CompleteGene' set:\n")
		cat( "\nN: ", sum(is.na(where)), "  Who:  ", rownames(expectLog2)[is.na(where)])
	}
	for (i in 1:nGenes) {
		if ( is.na( where[i])) next
		outCounts[ , where[i]] <- expectLog2[ i, ]
		outCounts2[ , where[i]] <- expectNLog2[ i, ]
	}
	# then all the meta data (SampleID through Mutant Pool)
	outMeta <- samples[ , MetaDataColumns]

	# lastly, also write out the run names and the raw read counts and log2 RPMHEG values
	readCnts <- apply( tbl, MARGIN=2, sum, na.rm=T)
	out <- data.frame( "SampleID"=outMeta[ ,1], "Run"=runName, outMeta[ ,2:ncol(outMeta)], "ReadCount"=readCnts, 
				round( outCounts, digits=4), stringsAsFactors=FALSE)
	out2 <- data.frame( "SampleID"=outMeta[ ,1], "Run"=runName, outMeta[ ,2:ncol(outMeta)], "ReadCount"=readCnts, 
				round( outCounts2, digits=4), stringsAsFactors=FALSE)
	rownames(out) <- rownames(out2) <- 1:nrow(out)
	outfile <- file.path( results.path, paste( "FinalResults", experiment, "txt", sep="."))
	write.table( out, outfile, sep="\t", quote=F, row.names=F)
	outfile2 <- file.path( results.path, paste( "Normalized.FinalResults", experiment, "txt", sep="."))
	write.table( out2, outfile2, sep="\t", quote=F, row.names=F)
}


normalizeByUninduced <- function( tbl, samples, nExpectedGenes=rep.int( 100, ncol(tbl)),
				max.scale.factor=128) {

	# given a matrix of RPMHEG normalized read count data, apply a second round of normalization,
	# meant to compensate for how some mutant poolss may expand or contract overtime, separate from
	# any induction effect.  The idea is to put all time points at a similar abundance in the 
	# UnInduced samples, so Induction effect is less impacted by overall trends
	
	if ( ! all( colnames(tbl) == samples$SampleID)) stop( "Sample vs matrix sizing or naming error!")
	SAV_SAM <<- samples
	NC <- ncol(tbl)
	
	# all shared day zero samples are always in every group, and always treated like uninduced
	dayZeroSet <- which( samples$Day %in% DAY_ZERO)
	
	# know the induction status for every sample
	isInduced <- as.logical(samples$ATc_Induced)
	useForNorm <- which( ! isInduced)

	# build a set of sample keys of various facts except induction.  We process one experiment at a time, so that fact not in the key
	sDayKey <- paste( samples$Day, samples$Class, samples$SubClass1, samples$SubClass2, sep="::")
	sDayKeyFac <- factor(sDayKey)
	sDayKeyPtrs <- tapply( 1:NC, sDayKeyFac, FUN=NULL)
		
	out <- tbl
	min.scale.factor <- 1 / max.scale.factor
	min.value <- 2 ^ MIN_LOG2_RPMHEG
	for ( i in 1:nrow(tbl)) {
		vIn <- tbl[ i, ]
		# now give each sample the correct scale factor from it's set of uninduced matched sample
		allScaleFac <- rep.int(1, NC)
		for ( k in 1:NC) {
			# since the uninduced are shared by all grouping, never try to adjust them
			if ( k %in% dayZeroSet) next
			# when it's a day zero column, expand its group...
			myDaySet <-  which( sDayKeyPtrs == sDayKeyPtrs[k])
			myUnindDayPtrs <- intersect( myDaySet, useForNorm)
			# lastly, only use those with data
			myUnindDayPtrs <- intersect( myUnindDayPtrs, which(vIn >= min.value))
			if ( length(myUnindDayPtrs)) {
				# use avg of all day zero to pick a central target level for all uninduced
				globalAvg <- sqrtmean( vIn[dayZeroSet], na.rm=T)
				if (globalAvg < min.value) next
				# then use only the subset from this day subset's uninduced to make the scale term for this sample
				scaleFac <- globalAvg / sqrtmean(vIn[myUnindDayPtrs], na.rm=T)
				scaleFac[ is.na(scaleFac)] <- 1
				scaleFac[ is.nan(scaleFac)] <- 1
				scaleFac[ scaleFac > max.scale.factor] <- max.scale.factor
				scaleFac[ scaleFac < min.scale.factor] <- min.scale.factor
				allScaleFac[k] <- scaleFac
			}
		}
		vOut <- vIn * allScaleFac
		out[ i, ] <- vOut
	}
	return(out)
}


combineResults <- function( folder=SampleKeyFolder) {

	# get all the sample key files
	keyFiles <- dir( folder, patt="^SampleKey.+\\.txt$", full.name=T)
	nFiles <- length(keyFiles)
	cat( "\nN_Sample Key files: ", nFiles)

	out <- data.frame()
	NC.out <- 0
	for ( i in 1:nFiles) {
		samples <- read.delim( keyFiles[i], as.is=T)
		runName <- extractRunName( keyFiles[i])
		results.path <- createResultsFolder( keyFiles[i], create=FALSE)
		if ( is.null( results.path)) {
			cat( "\nResults folder not found: ", runName)
			next
		}
		experiments <- sort( unique( samples$Experiment))
		for (e in experiments) {
			f <- file.path( results.path, paste( "FinalResults", e, "txt", sep="."))
			if ( ! file.exists(f)) {
				cat( "\nFinal Results File not found: ", f)
				next
			}
			thisDF <- read.delim( f, as.is=T)

			# check before blindly combining
			if ( NC.out == 0) {
				NC.out <- ncol(thisDF)
			} else {
				if ( ncol(thisDF) != NC.out) {
					cat( "\nWrong size results file..  Skipping: ", basename(f))
					next
				}
			}
			out <- rbind( out, thisDF)
			cat( "\nRun: ", runName, "  Experiment: ", e, "  N_Samples: ", nrow(thisDF))
		}
	}
	nExp <- length( unique( out$Experiment))
	nSamples <- nrow(out)
	cat( "\nN_Experiments: ", nExp)
	cat( "\nN_Samples:     ", nSamples)

	outfile <- file.path( ResultsFolder, "FinalResults.All.TRIP.Log2RPMHEG.txt")
	write.table( out, outfile, sep="\t", quote=F, row.names=F)
}


# build a consistent style ID for samples from their meta data columns
getSampleLabel <- function( samples, mode=c("isolate", "group", "ID")) {

	mode <- match.arg( mode)

	# using the named columns, not the SampleID anymore
	N <- nrow(samples)
	pool <- samples$MutantPool
	nPools <- length( unique( pool))
	class <- samples$Class
	nClass <- length( unique( class))
	sub1 <- samples$SubClass1
	nSub1 <- length( unique( sub1[sub1 != ""]))
	sub2 <- samples$SubClass2
	# sub clase 2 may have decimal quantities
	sub2 <- sub( "0p", "0.", sub2, fixed=T)
	nSub2 <- length( unique( sub2[sub2 != ""]))

	# try to build a meaningful key, but not full of empty terms
	key <- rep.int( "", N)
	if (nPools > 1) key <- pool
	seps <- ifelse( nchar(key) > 0, "_", "")
	if (nClass > 1) key <- paste( key, seps, class, sep="")
	# if we have no details yet, use the class
	if ( all( key == "")) key <- class
	seps <- ifelse( nchar(key) > 0, "_", "")
	if (nSub1 > 0) key <- paste( key, seps, sub1, sep="")
	seps <- ifelse( nchar(key) > 0, "_", "")
	if (nSub2 > 0) key <- paste( key, seps, sub2, sep="")
	sampleKey <- key
	if ( mode == "group") return( sampleKey)

	induced <- ifelse( samples$Induced, "Ind", "Un")
	day <- paste( "D", samples$Day, sep="")
	repl <- paste( "r", samples$Replicate, sep="")
	nRepl <- length( unique( repl))
	txt <- paste( day, induced, sep="_")
	if ( nRepl > 1) txt <- paste( txt, repl, sep="_")
	# prepend sub class info to the text
	seps <- ifelse( nchar(sub2) > 0, "_", "")
	if (nSub2 > 1) txt <- paste( sub2, seps, txt, sep="")
	seps <- ifelse( nchar(sub1) > 0, "_", "")
	if (nSub1 > 1) txt <- paste( sub1, seps, txt, sep="")
	sampleText <- txt
	if ( mode == "isolate") return( sampleText)

	return( samples$SampleID)
}


# make a scatterplot overview that shows how many genes and reads each mutant got.  Useful for
# deciding what samples to exclude for low/unexpected counts.
plotAlignmentOverview <- function( m=outStats, samples, results.path=".") {

	require( plotrix)
	par( mai=c( 1, 1, 0.8, 0.4))

	x <- m[ ,'N_Genes']
	y <- m[,'Pct_ValidPairs']
	x[ is.na(x)] <- 0
	y[ is.na(y)] <- 0
	nPts <- nrow(m)
	sampleKey <- getSampleLabel( samples, mode="group")
	sampleText <- getSampleLabel( samples, mode="isolate")
	keys <- sort( unique( sampleKey))
	nKeys <- length( keys)
	colorSet <- rainbow( nKeys, end=0.8)
	col <- colorSet[ match( sampleKey, keys)]
	# if we have only one class & subclass, show the full IDs
	if ( nKeys < 2) sampleText <- getSampleLabel( samples, mode="ID")
	
	# set up, and size by how many to draw
	xlim <- range( x) * c( 0.3, 1.25)
	txt.cex <- 1
	if ( nPts > 6) txt.cex <- 0.9
	if ( nPts > 12) txt.cex <- 0.87
	if ( nPts > 18) txt.cex <- 0.84
	if ( nPts > 24) txt.cex <- 0.81
	if ( nPts > 30) txt.cex <- 0.78
	if ( nPts > 40) txt.cex <- 0.75
	if ( nPts > 60) txt.cex <- 0.72
	leg.cex <- 1
	if ( nKeys > 4) leg.cex <- 0.9
	if ( nKeys > 8) leg.cex <- 0.82
	if ( nKeys > 12) leg.cex <- 0.74
	if ( nKeys > 16) leg.cex <- 0.66

	plot( x, y, xlab="Number of Genes Detected per Sample", xlim=xlim, 
		ylab="Pct Valid Mate Pairs per Sample", ylim=c(0,100),
		main=paste( "Alignment Success Overview:    ", basename(results.path)),
		pch=21, cex=1.5, bg=col, col=1, font.lab=2, font.axis=2)

	# replot in random order to show colors better
	ord <- sample( length(x))
	points( x[ord], y[ord], pch=21, cex=1.5, bg=col[ord])
	thigmophobe.labels( x, y, sampleText, cex=txt.cex, font=1)

	legend( 'bottomright', keys, pch=21, pt.bg=colorSet, cex=leg.cex, pt.cex=1.3)
	outfile <- file.path( results.path, paste( "AlignmentSuccessOverview", "pdf", sep="."))
	dev.print( pdf, outfile, width=14, height=9)

	# also make a barplot of the read counts
	par( mai=c( 2.5, 1, 0.8, 0.4))

	cntM <- as.matrix( m[ , c( 'StartingReads', 'ValidPairs')])
	rownames(cntM) <- rownames(m)
	cntM[ is.na(cntM)] <- 0
	colorSet <- c( 'tan', 'green')

	barplot( t( cntM), beside=T, xlab=NA, ylab="Number of Read Pairs", main="Read Counts Overview",
		ylim=c(0,max(cntM[,1])*1.05), col=colorSet, las=3)

	lines( c( -10, 1000), rep.int( MIN_READ_PAIRS, 2), lty=2, lwd=3, col='grey50')
	text( nPts*1.5, MIN_READ_PAIRS, "Minimum Valid Pairs", col='grey50', cex=0.8, pos=3)

	legend( 'topleft', c( "Raw Read Pairs", "Valid Pairs"), fill=colorSet, cex=leg.cex)
	outfile <- file.path( results.path, paste( "ReadCountsOverview", "pdf", sep="."))
	dev.print( pdf, outfile, width=14, height=9)
}


replotAlignmentOverview <- function( sampleKeyFile) {

	# pull back in the results summary and plot it
	samples <- loadSampleKeyFile( sampleKeyFile)
	results.path <- createResultsFolder( sampleKeyFile)

	# get the file of stats
	statfile <- file.path( results.path, "AlignmentStatistics.txt")
	outStats <- read.delim( statfile, as.is=T)

	# call the plotter
	plotAlignmentOverview( outStats, samples, results.path=results.path)
}


# make an ROC image that compare expected mutant genes against unexpected genes.  Can be helpful for 
# deciding a minimum abundance threshold.
makeROCimage <- function( out, poolGenes, text, experiment, results.path=".") {

	cat( "\nMaking ROC plot..")
	isPool <- which( rownames(out) %in% poolGenes)
	isBad <- setdiff( 1:nrow(out), isPool)

	goodScores <- as.vector( out[ isPool, ])
	badScores <- as.vector( out[ isBad, ])

	textLabel <- paste( "Experiment:  ", experiment, "      Abundance Units:  ", text)
	ans <- duffy.ROC( goodScores, badScores, label=textLabel)

	plotFile <- file.path( results.path, paste( "ROC", experiment, gsub(" ",".",text), "pdf", sep="."))
	dev.print( pdf, plotFile, width=14,  height=9)
	return( ans)
}


loadSampleKeyFile <- function( sampleKeyFile) {

	if ( ! file.exists( sampleKeyFile)) stop( "'Sample key file' not found..  Tried: ", sampleKeyFile)

	allSamples <- read.delim( sampleKeyFile, as.is=T)

	# force all sample IDs to be valid R names
	checkSampleNames( allSamples$SampleID)
	
	# drop any flagged for exclusion
	if ( "Exclude" %in% colnames(allSamples)) {
		toDrop <- which( allSamples$Exclude == TRUE)
		if ( length( toDrop)) {
			cat( "\nDropping some samples as 'Excluded': ", length(toDrop))
			allSamples <- allSamples[ -toDrop, ]
		}
	}

	# force any NA to be blanks
	allSamples <- cleanMetaDataFields( allSamples)
	if ( "FailSeq" %in% colnames(allSamples)) {
		toDrop <- which( toupper(as.character(allSamples$FailSeq)) %in% c( "T", "TRUE"))
		if ( length( toDrop)) allSamples <- allSamples[ -toDrop, ]
	}
	nSamples <- nrow( allSamples)
	if ( nSamples < 1) stop( "No non-excluded samples!")

	if ( ! all( SampleKeyColumns %in% colnames(allSamples))) {
		cat( "\nSample Key file must have the following columns:  ", SampleKeyColumns)
		stop()
	}

	if ( length( isDup <- which( duplicated( allSamples$SampleID)))) {
		stop( paste( "Found duplicate SampleIDs.  Not allowed.  \nCheck: ", 
				paste( allSamples$SampleID[isDup], collapse="; ")))
	}

	cat( "\n\nSample Key: ", sampleKeyFile, "\nN_Samples:      ", nSamples, "\n")
	cat( "\nExperiment Names: ", sort( unique( allSamples$Experiment)))
	cat( "\nMutant Pools:     ", sort( unique( allSamples$MutantPool)))

	# test the existance of the files
	f1 <- file.path( allSamples$FastqPath, allSamples$File1)
	notfound <- which( ! file.exists( f1))
	if ( length(notfound)) {
		cat( "\nMate 1 files not found: ", length(notfound), "\n", allSamples$File1[notfound])
	}
	f2 <- file.path( allSamples$FastqPath, allSamples$File2)
	notfound <- which( ! file.exists( f2))
	if ( length(notfound)) {
		cat( "\nMate 2 files not found: ", length(notfound), "\n", allSamples$File2[notfound])
	}
	return( allSamples)
}


# turn the filename of the sample key into the name to attach to results, etc.
extractRunName <- function( sampleKeyFile) {

	runName <- sub( "^SampleKey.", "", basename(sampleKeyFile))
	runName <- sub( ".txt$", "", runName)
	runName <- sub( ".csv$", "", runName)
	return( runName)
}


createResultsFolder <- function( sampleKeyFile, create=TRUE) {

	runName <- extractRunName( sampleKeyFile)
	results.path <- file.path( ResultsFolder, runName)
	if ( ! file.exists( results.path)) {
		if (create) {
			dir.create( results.path, recursive=T)
		} else {
			return(NULL)
		}
	}
	return( results.path)
}


# returns the gene names that we expect to see mate pair hits to, as these are the constructs 
# in the mutant pool used for the experiment
getMutantPoolGenes <- function( myPool) {

	poolFile <- paste( "TRIP", myPool, "Primers.txt", sep=".")
	poolFile <- file.path( TRIP.PATH, "CustomMutantGenome", poolFile)
	if ( ! file.exists( poolFile)) stop( paste( "Failed to find Mutant Pool Primers file: ", poolFile))
	poolTable <- read.delim( poolFile, as.is=T)
	if ( ! ("GeneID" %in% colnames(poolTable))) stop( "Mutant Pool Primer file must have a 'GeneID' column")
	poolGenes <- poolTable$GeneID
	return( poolGenes)
}


cleanMetaDataFields <- function( tbl) {

	# force any NA to be blanks
	if ( "FailSeq" %in% colnames(tbl)) tbl$FailSeq[ is.na( tbl$FailSeq)] <- ""
	if ( "Class" %in% colnames(tbl)) tbl$Class[ is.na( tbl$Class)] <- ""
	if ( "SubClass1" %in% colnames(tbl)) tbl$SubClass1[ is.na( tbl$SubClass1)] <- ""
	if ( "SubClass2" %in% colnames(tbl)) tbl$SubClass2[ is.na( tbl$SubClass2)] <- ""
	return( tbl)
}


createClassDescriptor <- function( tbl) {

	# try to use the Class/SubClass data to make a useful descriptor for plots, etc.
	tbl <- cleanMetaDataFields( tbl)
	
	# if one field, use it...
	classes <- setdiff( unique( tbl$Class), "")
	subclass1 <- setdiff( unique( tbl$SubClass1), "")
	subclass2 <- setdiff( unique( tbl$SubClass2), "")
	descriptor <- if (length(classes) == 1) classes[1] else "MultiClass"
	if ( length(subclass1)) {
		suffix <- if ( length(subclass1) == 1) subclass1[1] else "MultiSubClass1"
		descriptor <- paste( descriptor, suffix, sep="_")
	}
	if ( length(subclass2)) {
		suffix <- if ( length(subclass2) == 1) subclass2[1] else "MultiSubClass2"
		descriptor <- paste( descriptor, suffix, sep="_")
	}
	return( descriptor)
}


# this is the low level core growth defects modeling function, to be called for one experiments set of data.
#    Note:  for the high level wrapper called by the "do.all()" function, see:  model.All.TRIP.Samples()

model.TRIP.GrowthDefects <- function( tbl, min.value=MIN_LOG2_RPMHEG, makePlots=FALSE, 
					plot.path=".", plot.prefix=NULL, asPDF=makePlots,
					normalized=FALSE) {

	# given a data frame of final results (large table with columns of meta data & Mutant gene log2RPMHEG data)
	# calculate all growth defects from 'Induced' vs 'Uninduced'
	# this routine assumes all rows go into the model!
	# if you want to subset by Class / SubClass / etc., do that in a parent function that calls this.

	if ( ! all( MetaDataColumns %in% colnames(tbl))) {
		cat( "\nMissing some required Meta Data Columns..  Invalid input..")
		stop()
	}
	# for debugging
	saveModelInput <<- tbl

	# make sure we have both Induced and Uninduced rows
	inducedRows <- which( tbl$ATc_Induced)
	if ( ! length(inducedRows)) stop( "No Induced samples in this model. Can't run LM...")
	uninducedRows <- which( ! tbl$ATc_Induced)
	if ( ! length(uninducedRows)) stop( "No Uninduced samples in this model. Can't run LM...")
	
	# we will run the same model against each gene column, after dropping out any meta data columns
	geneNames <- setdiff( colnames(tbl), ModelMetaDataColumns)
	nGenes <- length( geneNames)

	# by rule, all the 'DayZero' data will be used by both Induced and Uninduced
	dayZeroRows <- which( tbl$Day %in% DAY_ZERO)
	
	# original intent was to always have some Day Zero, and to model each later day as its own decay.
	# but if we don't not have any Day Zero, then model all days as a single linear model
	MODEL_BY_DAY <- TRUE
	if ( ! length( dayZeroRows)) {
		cat( "\n\nWarning:  no Day Zero data found!  Will try running a single model of all days together.")
		daysToModel <- "_All_Joined"
		nDays <- 1
		MODEL_BY_DAY <- FALSE
	} else {
		# and we will model each Day (past Day 1) on its own, 
		daysToModel <- setdiff( sort( unique( tbl$Day)), DAY_ZERO)
		nDays <- length( daysToModel)
	}
	
	# storage for data frames and model results in case we draw
	dfByDayIND <- dfByDayUNIND <- vector( length=nDays, mode="list")
	names(dfByDayIND) <- names(dfByDayUNIND) <- daysToModel
	modelAnsByDayIND <- modelAnsByDayUNIND <- vector( length=nDays, mode="list")
	names(modelAnsByDayIND) <- names(modelAnsByDayUNIND) <- daysToModel

	# set up results storage
	growthRate <- growthPvalue <- matrix( NA, nrow=nDays, ncol=nGenes)
	rownames(growthRate) <- rownames(growthPvalue) <- paste( "Day", daysToModel, sep="")
	colnames(growthRate) <- colnames(growthPvalue) <- geneNames

	# if we will make plots, set up for what we need
	if (makePlots) {
		if ( ! file.exists(plot.path)) dir.create( plot.path, recursive=T)
		if ( is.null( plot.prefix)) {
			plot.prefix <- createClassDescriptor(tbl)
		}
	}

	# ready to loop over all Days/Genes
	for ( j in 1:nGenes) {
		thisGene <- geneNames[j]

		# zero out storage passed to plotter
		for ( i in 1:nDays) {
			dfByDayIND[[i]] <- NULL
			dfByDayUNIND[[i]] <- NULL
			modelAnsByDayIND[[i]] <- NULL
			modelAnsByDayUNIND[[i]] <- NULL
		}
		haveDataToPlot <- FALSE

		for ( i in 1:nDays) {
			thisDay <- daysToModel[i]
			myRows <- if (MODEL_BY_DAY) which( tbl$Day == thisDay) else 1:nrow(tbl)

			# make the small data frames for both models
			wantColumns <- c( thisGene, "Day", "Class", "SubClass1", "SubClass2")
			wantRowsIND <- c( dayZeroRows, intersect( inducedRows, myRows))
			wantRowsUNIND <- c( dayZeroRows, intersect( uninducedRows, myRows))
			smlIND <- tbl[ wantRowsIND, wantColumns]
			smlUNIND <- tbl[ wantRowsUNIND, wantColumns]
			colnames(smlIND) <- colnames(smlUNIND) <- c( "Value", wantColumns[2:length(wantColumns)])

			# values below a threshold are not real -- set to NA
			# No- don't remove, just use as clipping threshold for more robust model
			smlIND$Value[ smlIND$Value < min.value] <- NA
			smlUNIND$Value[ smlUNIND$Value < min.value] <- NA

			# there may be no data values for some genes
			if (MODEL_BY_DAY) {
				if ( all( is.na( smlIND$Value[ smlIND$Day == thisDay]))) next
			} else {
				if ( all( is.na( smlIND$Value))) next
			}

			# there may be no uninduced data for this day, 
			# if so propagate the day zero data to be this day's control
			if (MODEL_BY_DAY) {
				nThisDay <- sum( smlUNIND$Day ==  thisDay)
				if ( nThisDay == 0) {
					extraData <- smlUNIND
					extraData$Day <- thisDay
					smlUNIND <- rbind( smlUNIND, extraData)
				}
				if ( all( is.na( smlUNIND$Value[ smlUNIND$Day == thisDay]))) next
			} else {
				if ( all( is.na( smlUNIND$Value))) next
			}

			# if only a single value, the modeling will bust, extend if needed
			smlIND <- forceEnoughModelData( smlIND)
			if ( is.null(smlIND)) next
			smlUNIND <- forceEnoughModelData( smlUNIND)
			if ( is.null(smlUNIND)) next

			# ready to run the model, save the values for drawing
			dfByDayIND[[i]] <- smlIND
			dfByDayUNIND[[i]] <- smlUNIND
			saveIND <<- smlIND
			saveUNIND <<- smlUNIND

			# call the linear model for both 
			modelAnsByDayIND[[i]] <- ansIND <- lm( Value ~ Day, data=smlIND)
			modelAnsByDayUNIND[[i]] <- ansUNIND <- lm( Value ~ Day, data=smlUNIND)

			# see how different they are, first by LM slope
			diffAns <- lmSlopeDifference( ansUNIND, ansIND)
			lmSlope <- diffAns$difference
			lmPval <- diffAns$p.value
			if ( is.na(lmPval) || is.nan(lmPval)) lmPval <- 1

			# also do simple Student's T test on just the day of interest
			if (MODEL_BY_DAY) {
				unY <- subset( smlUNIND, Day == thisDay)$Value
				indY <- subset( smlIND, Day == thisDay)$Value
				ttestAns <- t.test( unY, indY)
				ttestSlope <- diff( ttestAns$estimate) / thisDay
				ttestPval <- ttestAns$p.value
				if ( is.na(ttestPval) || is.nan(ttestPval)) ttestPval <- 1
			} else {
				# when we do just a single model of all days, use the 'last/largest' day
				bigDay <- max( smlIND$Day, na.rm=T)
				unY <- subset( smlUNIND, Day == bigDay)$Value
				indY <- subset( smlIND, Day == bigDay)$Value
				ttestAns <- t.test( unY, indY)
				ttestSlope <- diff( ttestAns$estimate) / bigDay
				ttestPval <- ttestAns$p.value
				if ( is.na(ttestPval) || is.nan(ttestPval)) ttestPval <- 1
			}
			
			# final answer is the average of the 2 methods
			finalGrowth <- growthRate[i,j] <- mean( c( lmSlope, ttestSlope))
			finalPvalue <- growthPvalue[i,j] <- logmean( c( lmPval, ttestPval))
			haveDataToPlot <- TRUE
			cat( "\r", j, thisGene, thisDay, nrow(smlIND), nrow(smlUNIND), lmSlope, 
					ttestSlope, lmPval, ttestPval)
		}

		# we may do a plot...
		if ( haveDataToPlot && makePlots) {
			plotOneGene.TRIP( thisGene, dfByDayIND, dfByDayUNIND, modelAnsByDayIND, modelAnsByDayUNIND, 
					plot.path=plot.path, plot.prefix=plot.prefix, asPDF=asPDF, 
					showGrowth=finalGrowth, showPvalue=finalPvalue, normalized=normalized)
		}
	}
	return( list( "GrowthDefect"=growthRate, "P_Value"=growthPvalue))
}


# linear regression fitting can break on single data points. Crude way to up-sample if we need to.
forceEnoughModelData <- function( df) {

	# make sure every day has at least 2 observations
	dayFac <- factor( df$Day)
	nGood <- tapply( !is.na(df$Value), dayFac, sum)
	if ( all( nGood > 1)) return(df)

	df$DrawIt <- TRUE
	out <- df

	for (i in 1:nlevels(dayFac)) {
		if (nGood[i] > 1) next
		thisDay <- as.numeric( levels(dayFac)[i])
		curVal <- df$Value[ df$Day == thisDay]
		goodVal <- curVal[ !is.na(curVal)][1]
		if ( is.na(goodVal)) return(NULL)
		extraRow <- df[ which( df$Value == goodVal & df$Day == thisDay)[1], ]
		extraRow$Value <- jitter(goodVal, factor=1, amount=0.5)
		extraRow$Day <- thisDay
		extraRow$DrawIt <- FALSE
		out <- rbind( out, extraRow)
	}
	return( out)
}


# make a plot image that shows one mutant gene's abundance data and model best fit lines.
# Used as hyperlink images in final HTML files.
plotOneGene.TRIP <- function( gene, dfsIND, dfsUNIND, modelsIND, modelsUNIND, 
			plot.path=".", plot.prefix="", asPDF=FALSE, 
			showGrowth=NULL, showPvalue=NULL, normalized=FALSE) {

	checkX11( bg="white", width=10, height=7)
	par( mai=c( 1,1,0.8,0.4))
	fileRoot <- paste( gene, plot.prefix, sep="_")

	# scan all the data frames for plot extents
	allDays <- allValues <- allClass <- vector()
	nModels <- length( dfsIND)
	if (nModels < 1) return()
	for (i in 1:nModels) {
		smlIND <- dfsIND[[i]]
		if ( is.null(smlIND)) next
		smlUNIND <- dfsUNIND[[i]]
		if ( is.null(smlUNIND)) next
		# small chance of extra non-drawn data to prevent breaking the model fitting
		if ( "DrawIt" %in% colnames(smlIND)) smlIND <- subset( smlIND, DrawIt)
		if ( "DrawIt" %in% colnames(smlUNIND)) smlUNIND <- subset( smlUNIND, DrawIt)
		allDays <- c( allDays, smlIND$Day)
		allDays <- c( allDays, smlUNIND$Day)
		allValues <- c( allValues, smlIND$Value)
		allValues <- c( allValues, smlUNIND$Value)
	}
	xlim <- range( allDays, na.rm=TRUE) + c( -0.5,0.5)
	ylim <- range( allValues, na.rm=TRUE)
	ylim[2] <- ylim[2] + diff(ylim)*0.35
	if ( ylim[1] > 0) ylim[1] <- ylim[1] * 0.9
	yLabel <- "Abundance   Log2(RPMHEG)"
	if (normalized) yLabel <- "Abundance   Log2(Normalized RPMHEG)"

	mainText <- paste( "Mutant TRIP Timecourse:    ", plot.prefix, "\n", gene, " -  ", gene2Product(gene))
	plot( 1, 1, type='n', main=mainText, xlim=xlim, ylim=ylim, xlab="Time  (days)", xaxt='n',
			ylab=yLabel, font.axis=2, font.lab=2)
	axis( side=1, at=sort( unique( allDays)), font=2)


	# local function to draw one "Linear Model" object...
	showOneTimeCourse <- function( df, lmAns, thisCol, thisPCH=21, pt.cex=1.5, showZero=TRUE, showLine=TRUE) {

		if ( is.null(df)) return()

		thisX <- df$Day
		thisY <- df$Value
		xFac <- factor( thisX)
		yMeans <- rep.int( 0, nlevels(xFac))
		saveDF <<- df
		
		# perhaps do not redraw the zero points (where 'Zero' includes day 1)
		toShow <- 1:length(thisX)
		if ( !showZero) toShow <- which( ! (thisX %in% DAY_ZERO))

		points( jitter(thisX,factor=0.2)[toShow], thisY[toShow], pch=thisPCH, bg=thisCol, cex=pt.cex)
		yMeans <- tapply( thisY, xFac, mean)
		mX <- as.numeric( names(yMeans))
		for ( who in which( yMeans > 0)) {
			if ( !showZero && mX[who] %in% DAY_ZERO) next
			lines( mX[who]+c(-0.2,0.2), rep.int(yMeans[who],2), col=thisCol, lwd=2)
		}

		if ( is.null(lmAns)) return()
		if ( showLine && length(yMeans) > 1) {
			abline.segment( reg=lmAns, x1=min(mX), x2=max(mX), col=thisCol, lty=2, lwd=2)
		}
	}


	# OK, all set up, draw those time courses
	colUN <- "dodgerblue"
	colIND <- "red"
	
	# for the uninduced, only show the day zero and best fit for the last model
	for ( i in 1:nModels) {
		if ( is.null( dfsUNIND[[i]])) next
		showZero <- showLine <- (i == nModels)
		showOneTimeCourse( dfsUNIND[[i]], modelsUNIND[[i]], colUN, showZero=showZero, showLine=TRUE)  # showLine)
	}

	# for the induced, never show the day zero and always show the best fit 
	for ( i in 1:nModels) {
		if ( is.null( dfsIND[[i]])) next
		showOneTimeCourse( dfsIND[[i]], modelsIND[[i]], colIND, showZero=FALSE, showLine=TRUE)
	}

	legend( "topright", c('UnInduced', 'Induced  '), col=c(colUN,colIND), lwd=5, bg='white', cex=1.1)
	if ( ! is.null(showGrowth)) {
		pValText <- formatC( showPvalue, format="e", digits=2)
		legend( 'bottom', paste( "Delta Growth = ", round(showGrowth,digits=4), 
				"    P-value = ", pValText), bg='white', cex=1.1)
	}
	dev.flush()
	
	if (asPDF) {
		pdffile <- file.path( plot.path, paste( fileRoot, "pdf", sep="."))
		dev.print( pdf, pdffile, width=14, height=9)
	}
}


# top level function to model the growth defects and make plots and growth defect summarys 
# for all experiments in one sample key file
model.All.TRIP.Samples <- function( sampleKeyFile="SampleKey.LogTRIP.txt",
				sharedDayZero=NULL, makePlots=TRUE, experiment=NULL,
				normalized=FALSE) {

	# load that sample key
	allSamples <- loadSampleKeyFile( sampleKeyFile)
	nSamples <- nrow( allSamples)

	# allow specifying just one experiment at a time, for testing, etc.
	if ( ! is.null( experiment)) {
		allSamples <- subset( allSamples, Experiment == experiment)
		nSamples <- nrow( allSamples)
	}

	# the 'Run Name' is extracted from the name of the Samples Key filename for the results folder name
	results.path <- createResultsFolder( sampleKeyFile)

	# let's allow more than one 'Experiment' per Sample Key file
	expFactor <- factor( allSamples$Experiment)
	tapply( 1:nrow(allSamples), expFactor, function(x) {

		# do the subset of samples for this experiment
		samples <- allSamples[x, ]
		nSamples <- nrow( samples)

		# get the final results table for experiment
		myExperiment <- samples$Experiment[1]
		cat( "\nModelling TRIP growth defects for experiment: ", myExperiment)
		cat( "\n  Normalization Mode:  ", normalized)

		resultsFileIn <- paste( "FinalResults", myExperiment, "txt", sep=".")
		if (normalized) resultsFileIn <- paste( "Normalized.FinalResults", myExperiment, "txt", sep=".")
		resultsFileIn <- file.path( results.path, resultsFileIn)
		if ( ! file.exists( resultsFileIn)) {
			cat( "\nResults file not found:  ", resultsFileIn)
			cat( "\nSkipping...")
			return()
		}
		tbl <- read.delim( resultsFileIn, as.is=T)
		tbl <- cleanMetaDataFields( tbl)

		# create folder for growth defect plots
		plot.path <- file.path( results.path, paste( "TRIP.GrowthDecayPlots", myExperiment, sep="."))
		if (normalized) plot.path <- file.path( results.path, paste( "TRIP.Normalized.GrowthDecayPlots", myExperiment, sep="."))
		if ( ! file.exists( plot.path)) dir.create( plot.path, recursive=TRUE)

		# default behavior is to process each sub class separately, and merge them all back together
		out <- data.frame()
		facClass <- factor(tbl$Class)
		facSub1 <- factor(tbl$SubClass1)
		facSub2 <- factor(tbl$SubClass2)
		
		# allow a common shared set of Day Zero uninduced measrurements to be shared
		# let's relax this rule to look in multiple class columns, and let the term be a vector of text strings
		commonDayZeroUNDrows <- vector()
		if ( ! is.null( sharedDayZero)) {
			# first place to check is SubClass1
			common1 <- common2 <- common3 <- vector()
			if ( any (tbl$SubClass1 %in%  sharedDayZero)) {
				common1 <- which( tbl$Day %in% DAY_ZERO & tbl$SubClass1 %in% sharedDayZero)
			}
			# second place to check is Class
			if ( any (tbl$Class %in% sharedDayZero)) {
				common2 <- which( tbl$Day %in% DAY_ZERO & tbl$Class %in% sharedDayZero)
			}
			# third place to check is SubClass2
			if ( any (tbl$SubClass2 %in% sharedDayZero)) {
				common3 <- which( tbl$Day %in% DAY_ZERO & tbl$SubClass2 %in% sharedDayZero)
			}
			commonDayZeroUNDrows <- sort( unique( c( common1, common2, common3)))
		}
		
		# change the logic here to make sure we alway allow the shared day 0, but never use those samples twice
		#by( tbl, INDICES=list( facClass, facSub1, facSub2), function(x) {
		tapply( 1:nrow(tbl), INDEX=list( facClass, facSub1, facSub2), function(xx) {

			# make sure we have data, and get the class details
			if ( ! length(xx)) return()
			smlTbl <- tbl[ xx, ]

			desc <- createClassDescriptor(smlTbl)
			cat( "\nModelling:  ", desc, "\n")
			myclass <- smlTbl$Class[1]
			mysub1 <- smlTbl$SubClass1[1]
			mysub2 <- smlTbl$SubClass2[1]
			
			# allow the addition of the common shared day zero uninduced data
			if ( !is.null(sharedDayZero) && length(commonDayZeroUNDrows)) {
				xx <- union( xx, commonDayZeroUNDrows)
				smlTbl <- tbl[ xx, ]
			}

			# there may be cases where all the rows in a SubClass are day 0/1. as when its a negative control
			# group.   If so, catch it here and do nothing
			if ( all (smlTbl$Day %in% DAY_ZERO)) return()

			# call the modeler
			ans <- model.TRIP.GrowthDefects(smlTbl, makePlots=makePlots, plot.path=plot.path,
						plot.prefix=desc, asPDF=TRUE, normalized=normalized)
			saveModelAns <<- ans

			# we get back matrices of Rate and Pvalues, linearize them
			rates <- ans$GrowthDefect
			pvals <- ans$P_Value
			days <- as.numeric( sub( "Day", "", rownames(rates)))
			genes <- colnames(rates)
			nDays <- length(days)
			nGenes <- length(genes)
			smlday <- smlgene <- smlrate <- smlpval <- vector()
			nout <- 0
			for ( i in 1:nGenes) {
				nnow <- (nout+1) : (nout+nDays)
				smlday[nnow] <- days
				smlgene[nnow] <- genes[i]
				smlrate[nnow] <- rates[,i]
				smlpval[nnow] <- pvals[,i]
				nout <- nout + nDays
			}
			plotID <- paste( smlgene, desc, sep="_")
			sml <- data.frame( "PlotID"=plotID, "Gene"=smlgene, "Product"=gene2Product(smlgene), 
					"Class"=myclass, "SubClass1"=mysub1,
					"SubClass2"=mysub2, "Day"=smlday, "Delta_Growth"=round(smlrate,digits=4),
					"Delta_Pvalue"=smlpval, stringsAsFactors=FALSE)

			# discard empty results
			drops <- which( is.na( sml$Delta_Growth) & is.na( sml$Delta_Pvalue))
			if ( length( drops)) sml <- sml[ -drops, ]

			if ( nrow(sml)) out <<- rbind( out, sml)
		})

		# when all done, we can sort and write out
		ord <- diffExpressRankOrder( -out$Delta_Growth, out$Delta_Pvalue)
		out <- out[ ord, ]
		rownames(out) <- 1:nrow(out)

		resultsFileOut <- paste( "TRIP.GrowthDefect", myExperiment, "txt", sep=".")
		if (normalized) resultsFileOut <- paste( "TRIP.Normalized.GrowthDefect", myExperiment, "txt", sep=".")
		resultsFileOut <- file.path( results.path, resultsFileOut)
		write.table( out, resultsFileOut, sep="\t", quote=F, row.names=F)

		# for the HTML, drop meta data columns with no data
		if ( all( out$SubClass2 == "")) out <- out[ , -match("SubClass2",colnames(out))]
		if ( all( out$SubClass1 == "")) out <- out[ , -match("SubClass1",colnames(out))]
		htmlFile <- sub( "txt$", "html", resultsFileOut)
		colnames(out) <- sub( "_", " ", colnames(out))
		myTitle <- "TRIP Growth Defect Results:   "
		if ( normalized) myTitle <- "TRIP Normalized Growth Defect Results:   "
		table2html( out, htmlFile, title=paste( myTitle, myExperiment),
			linkColumnNames="PlotID", linkPaths=basename(plot.path), linkExtensions=".pdf")

	})  ## done with tapply() of all experiments in the run

	cat( "\nDone.\n")
}


poolHitsSummary <- function( folder) {

	# generate a summary of valid read counts that hit the expected pool of TRIP mutant genes
	# versus false hits to other genes not in the pool
	# 'folder' is one results folder from one run of 'align.OneSample()'
	allPattern <- "^AllGenes.RawCounts.*txt$"
	allFiles <- dir( folder, pattern=allPattern, full.name=T)
	poolFiles <- sub( "/All", "/Pool", allFiles)

	experiments <- sub( "(AllGenes\\.RawCounts\\.)(.+)(\\.txt)", "\\2", basename(allFiles))
	nExp <- length( experiments)

	allCounts <- poolCounts <- expNames <- colNames <- vector()
	for ( i in 1:length( allFiles)) {
		allTbl <- read.delim( allFiles[i], as.is=T)
		poolTbl <- read.delim( poolFiles[i], as.is=T)
		if ( ncol(allTbl) != ncol(poolTbl)) {
			cat( "\nError:  column count mismatch of 'All' & 'Pool' results: ", allFiles[i])
			next
		}
		if ( any( colnames(allTbl) != colnames(poolTbl))) {
			cat( "\nError:  column names mismatch of 'All' & 'Pool' results: ", allFiles[i])
			next
		}

		# get the column sums for each, and accumulate
		myCntsAll <- apply( as.matrix( allTbl[ ,2:ncol(allTbl)]), 2, sum, na.rm=T)
		myCntsPool <- apply( as.matrix( poolTbl[ ,2:ncol(poolTbl)]), 2, sum, na.rm=T)
		allCounts <- c( allCounts, myCntsAll)
		poolCounts <- c( poolCounts, myCntsPool)
		expNames <- c( expNames, rep.int( experiments[i], length(myCntsAll)))
		colNames <- c( colNames, names(myCntsAll))
	}

	# now we can know the percent that hit the pool
	meanAll <- round( mean( allCounts, na.rm=T))
	meanPool <- round( mean( poolCounts, na.rm=T))
	pctPool <- poolCounts * 100 / allCounts
	meanPct <- round( mean( pctPool, na.rm=T), digits=2)
	sdPct <- round( sd( pctPool, na.rm=T), digits=2)

	# report results to stdout
	cat( "\nFolder:      ", folder)
	cat( "\nExperiments: ", experiments)
	cat( "\n\nOverall Pool Hit Summary:")
	cat( "\nMean 'All'  Read Count per Column:  ", prettyNum(meanAll,big.mark=','))
	cat( "\nMean 'Pool' Read Count per Column:  ", prettyNum(meanPool,big.mark=','))
	cat( "\nMean % of Reads hitting Pool Genes: ", meanPct)
	cat( "\nStandard Deviation of % Pool:       ", sdPct)
	cat( "\n\nPercentage by Experiment:\n")
	print( tapply( pctPool, factor(expNames), function(x) return( round( mean( x, na.rm=T), digits=2))))
}

