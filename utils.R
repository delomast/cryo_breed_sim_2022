# utility functions used in other script(s)

library(dplyr)
library(stringr)
library(optiSel)
library(AllocateMate)

#' create G matrix with specified base pop frequencies
#' first method of VanRaden (2008), also Endelman and Jannink (2012)
#' @param g genotypes with rows inds, cols loci, coded as 0,1,2, no missing data allowed
#' @param af base generation frequencies of the non reference allele (the allele whose dosage the number corresponds to)
#'   in the same order as the columns of `genos`
createG <- function(g, af){
	if(length(af) != ncol(g)) stop("wrong af length for number of loci")
	# P matrix in VanRaden (2008)
	p <- matrix(2*(af - 0.5), nrow = nrow(g), ncol = ncol(g), byrow = TRUE)
	return(
		# subtract 1 to convert to -1,0-1 creating M,
		# then subtract P from M to create Z
		# then ZZ' and divide by P to scale G analagous to A
		tcrossprod((g - 1) - p) / (2 * sum(af *(1 - af)))
	)
}

#' adds spaces to teh end of strings to make
#' all strings in the vector the same length (number of characters)
#' Useful for creating fixed width files, such as for input of genotypes to
#' blupf90
#' @param x the character vector
#' @param pad the character to use to pad the strings
pad_id <- function(x, pad = " "){
	if(nchar(pad) != 1) stop("pad must be only one character in length")
	nc <- nchar(x)
	m <- max(nc)
	toAdd <- m - nc
	uToAdd <- unique(toAdd[toAdd > 0])
	if(length(uToAdd) < 1) return(as.character(x))
	for(u in uToAdd){
		temp <- paste(rep(pad, u), collapse = "")
		x[toAdd == u] <- paste0(x[toAdd == u], temp)
	}
	return(x)
}


#' function to apply full-sib family testing with an AlphaSimR object
#' pulls ID's and phenotypes for specified proportion (rounded) of individuals from each
#' full-sib family
#' 
#' @param fam an AlphaSimR "Pop-class" object containing full-sib families
#' @param propTest the proportion of each family to phenotype. One value applied
#'   to all families. The actual number of indiviudals to phenotype will be rounded.
sibTestEqual <- function(fam, propTest){
	fs <- data.frame(mother = fam@mother, father = fam@father) %>% 
		count(mother, father) %>% mutate(nPheno = round(propTest * n))
	phenos <- data.frame()
	for(i in 1:nrow(fs)){
		# select id
		toPheno <- sample(fam@id[fam@mother == fs$mother[i] & fam@father == fs$father[i]], size = fs$nPheno[i])
		temp <- cbind(toPheno, fs$mother[i], fs$father[i], as.data.frame(fam@pheno[match(toPheno, fam@id),]))
		colnames(temp) <- c("id", "mother", "father", paste0("Trait_", 1:fam@nTraits))
		phenos <- phenos %>% rbind(temp)
	}
	return(phenos)
}

#' using results of `opticont` from `optiSel` package
#' calculate number of matings for each individual with proper rounding
#' distributes rounding error randomly according to weights
#' but with no individual having more than 1 randomly allocated
#' to it
#' @param ocsParent $parent componenet of optisolve result
#' @param N number of matings to perform
calcNumMatings <- function(ocsParent, N){
	# males
	males <- ocsParent %>% filter(Sex == "male") %>% 
		mutate(n = round(2 * N * oc))
	diff <- N - sum(males$n)
	if(diff != 0){
		temp <- sample(1:nrow(males), size = abs(diff), prob = males$n, replace = FALSE)
		if(diff > 0) males$n[temp] <- males$n[temp] + 1
		if(diff < 0) males$n[temp] <- males$n[temp] - 1
	}
	# females
	females <- ocsParent %>% filter(Sex == "female") %>% 
		mutate(n = round(2 * N * oc))
	diff <- N - sum(females$n)
	if(diff != 0){
		temp <- sample(1:nrow(females), size = abs(diff), prob = females$n, replace = FALSE)
		if(diff > 0) females$n[temp] <- females$n[temp] + 1
		if(diff < 0) females$n[temp] <- females$n[temp] - 1
	}
	
	return(rbind(males, females))
}

#' Calculate the maximum mean kinship corresponding to 
#' a given effective population size and starting mean
#' kinship
#' @param kBar mean kinship at time 0
#' @param Ne desired effective population size
#' @param t0 start time
#' @param t end time
#' @param L generation length (time)
ubKin <- function(kBar, Ne, t0 = 0, t = 1, L = 1){
	return(1 - ((1 - kBar)*((1 - (1 / (2*Ne)))^((t - t0)/L))))
}

#' choose matings with OCS followed by inbreeding 
#' minimization
#' @param ocsData a data frame with each row corresponding to
#'   a selection candidate. Columns are Indiv, Sex (male and female)
#'   and gebv
#' @param Gmat genomic relationship matrix (the function internally converts
#'   this to the coancestry matrix)
#' @param N the number of matings to be performed
#' @param Ne the desired minimum effective population size
runOCS <- function(ocsData, Gmat, N, Ne = 50){
	# convert to coancestry/kinship matrix
	# and making sure order/presence of individuals is correct
	Gmat <- Gmat[ocsData$Indiv,ocsData$Indiv] / 2
	# data processing
	ocsCandes <- candes(phen = ocsData, N = N * 2, kin = Gmat)
	# optimum contributions
	ocsContrib <- opticont(method = "max.gebv", cand = ocsCandes, 
												 con = list(ub.kin = ubKin(kBar = ocsCandes$mean$kin, Ne = Ne)), 
												 trace=FALSE)
	# calculate number of matings per individual from contribution proportions
	ocsMatings <- calcNumMatings(ocsParent = ocsContrib$parent, N = N)
	
	# assign crosses (limit each pair to one cross)
	# to minimize inbreeding of each family
	# This branch and bound algorithm failed frequently with ub.n=1
	# crosses <- matings(ocsMatings, Kin=Gmat, ub.n = 1)
	ocsMatings <- ocsMatings %>% mutate(ID=as.character(Indiv), 
																			SEX=ifelse(Sex == "male", "M", "F"), 
																			EBV=gebv, N_AS_PARENT=n) %>% 
		select(ID, SEX, EBV, N_AS_PARENT) %>% filter(N_AS_PARENT > 0)
	crosses <- allocate.mate.H(H = Gmat[ocsMatings$ID, ocsMatings$ID]*2, 
														 parents = ocsMatings, max_F = 1, method = "min_F")
	
	return(crosses$optimal_families)
}

#' collapse a list of "Pop" class (alphaSimR) objects into one pop
#' with just the selected individuals
#' useful for breeding individuals across generations
#' when generations are saved as different elements of a list
#' and you want to reduce memory by avoiding combining all pops
#' in there entirety
#' @param pop_list list(pop1, pop2, pop3, ..., popN)
#' @param indivs vector of individula id's that need to be in the new pop
lowMem_mergePops <- function(pop_list, indivs){
	temp_pop_list <- list()
	for(i in 1:length(pop_list)){
		temp_indivs <- indivs[indivs %in% pop_list[[i]]@id]
		if(length(temp_indivs) > 0) temp_pop_list[[length(temp_pop_list) + 1]] <- 
				pop_list[[i]][temp_indivs]
	}
	return(mergePops(temp_pop_list))
}

#' read in scrm output to load into AlphaSimR
#' for use with output that has individuals as columns, SNPs as
#' rows (scrm ... -transpose-segsites)
#' Tries to limit memory usage
#' reads in one chromosome (one file), filters by maf, and
#' selects a random subset of loci
#' maintains phase, assumes diploidy, and scrm (should) only yield biallelic loci
#' assumes no missing genotypes (these are being treated as "true" to start a
#' simulation)
#' Deals with duplicate positions by making the current position 
#' equal to the current position + `incr` if the current is less than or equal
#' to the previous position. Does this AFTER subsampling. recommend use 1 / L where
#' L is the length of the chromosome (i.e., the increment should correspond to one base pair)
#' outputs as a matrix with rows haplotypes and columns loci
#' @param path path to the input file
#' @param numLoci number of loci to randomly sample (if < 1, no subsampling is performed).
#'   If no subsampling is required, use non-transposed output and other function for quicker
#'   processing
#' @param min_maf minimum maf to keep a locus
#' @param numLines number of lines to read at one time
#' @param incr the increment used to handle duplicate positions
#' @return a matrix with rows as haplotypes (adjacent rows are individuals) and cols as loci
read_scrm_transpose_for_AlphaSimR <- function(path, numLoci, min_maf = 0.05, numLines = 20000, incr = 1e-7){
	
	# read through once to assess maf and determine which SNPs to sample from
	f <- file(path, "r") # open scrm output
	# move to end of header
	line <- readLines(f, n = 1)
	while(length(line) > 0){
		if(substr(line, 1, 8) == "position") break
		line <- readLines(f, n = 1)
	}
	posToSample <- c() # all positions with valid maf
	line <- readLines(f, n = numLines)
	numC <- length(str_split(line[1], " ")[[1]])
	lineCounter <- 0
	while(length(line) > 0){
		snps <- matrix(as.numeric(str_split(line, " ", simplify = TRUE)), ncol = numC)
		snps <- snps[,-c(1,2)] # remove position and time column
		altFreq <- rowSums(snps) / ncol(snps)
		posToSample <- c(posToSample, lineCounter + which(altFreq >= min_maf & altFreq <= (1 - min_maf)))
		lineCounter <- lineCounter + length(line)
		line <- readLines(f, n = numLines) # read next chunk
	}
	close(f)
	
	# now subsample from SNPs with valid maf
	if(numLoci >= 1) posToSample <- sort(sample(posToSample, size = numLoci, replace = FALSE))
	
	# read through file again and only save SNPs that you want
	f <- file(path, "r") # open scrm output
	# move to end of header
	line <- readLines(f, n = 1)
	while(length(line) > 0){
		if(substr(line, 1, 8) == "position") break
		line <- readLines(f, n = 1)
	}
	line <- readLines(f, n = numLines)
	saveHaplos <- matrix(nrow = 0, ncol = numC - 1)
	lineCounter <- 0
	while(length(line) > 0){
		
		tempToSample <- posToSample[posToSample <= (lineCounter + length(line))]
		if(length(tempToSample) > 0){
			snps <- matrix(as.numeric(str_split(line[tempToSample - lineCounter], " ", simplify = TRUE)), ncol = numC)
			saveHaplos <- rbind(saveHaplos, snps[,-2]) # remove time column
			posToSample <- posToSample[-(1:length(tempToSample))] # remove already sampled SNPs
		}
		lineCounter <- lineCounter + length(line)
		line <- readLines(f, n = numLines) # read next chunk
	}
	close(f)
	rm(snps) # save a bit of memory
	# find any duplicated positions
	dups <- which(saveHaplos[2:nrow(saveHaplos),1] <= saveHaplos[1:(nrow(saveHaplos) - 1), 1])
	while(length(dups) > 0){
		for(d in dups){
			saveHaplos[d,1] <- saveHaplos[(d-1),1] + incr
		}
		# find any positions that are less than or equal to previous
		# less than could happen if incr is large enough to boost one position past the following position
		dups <- which(saveHaplos[2:nrow(saveHaplos),1] <= saveHaplos[1:(nrow(saveHaplos) - 1), 1])
	}
	
	# now transpose, set column names as pos, and remove pos row
	saveHaplos <- t(saveHaplos)
	colnames(saveHaplos) <- paste0("pos_", saveHaplos[1,])
	saveHaplos <- saveHaplos[-1,]

	return(saveHaplos)
}
