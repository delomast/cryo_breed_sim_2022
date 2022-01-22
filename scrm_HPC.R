# simulating breeding program with or without
# cryopreserved sperm for previous generations
# This script looks at accuracy of and with imputation
#   over several generations
#   and uses scrm to simulate the genome
# It is meant to be called by Rscript as part of a
# slurm array job

print(Sys.time())
print("begin")

.libPaths(c(.libPaths(), "/project/oyster_gs_sim/R_packages/4.1/"))
library(AlphaSimR, lib.loc="/project/oyster_gs_sim/R_packages/4.1/")
library(tidyverse, lib.loc="/project/oyster_gs_sim/R_packages/4.1/")
# library(rrBLUP, lib.loc="/project/oyster_gs_sim/R_packages/4.1/")
library(optiSel, lib.loc="/project/oyster_gs_sim/R_packages/4.1/")
library(AllocateMate, lib.loc="/project/oyster_gs_sim/R_packages/4.1/")

# library(AlphaSimR)
# library(tidyverse)
# library(rrBLUP)
# library(optiSel)
# library(AllocateMate)


source("utils.R")

cmdArgs <- commandArgs(trailingOnly=TRUE)

# Script parameters given on the command line
#' @param randSeed random seed to set R's random number generator
#' @param iterationNumber used to set unique file output names
randSeed <- as.numeric(cmdArgs[1])
iterationNumber <- cmdArgs[2]
localTempDir <- cmdArgs[3]
useCryo <- cmdArgs[4]

set.seed(randSeed)

chrLen <- c(32650045, 65668440, 61752955, 77061148, 59691872, 98698416, 51258098, 57830854, 75944018, 104168038)

qtlPerChr <- round(1000 * (chrLen/ sum(chrLen))) # does result in 1000, no rounding error
neutralPerChr <- round(50000 * (chrLen/ sum(chrLen))) #SNP chip per chromosome # results in 50000
nChr <- 10

nFound <- 200
nOffspringPerCross <- 50
nGenerations <- 8


# simulate and genotypes from scrm (one chromosome at a time), load in,
# and set up alphaSimR simulation
print(Sys.time())
print("begin simulating founders")

haplo_list <- list()
genMap <- list()
for(i in 1:nChr){
	# call scrm to simulate genotypes
	# 1 number of haplotypes
	# 2 Ne, 
	# 3 mutation rate per base, 
	# 4 expected number of recombinations per chromosome per generation,
	# 5 output file prefix
	# 6 random seed
	# 7 chromosome length
	system2("bash", args = c("sim_oyster_chrom_scrm.sh", nFound * 2, 
													 20000, 0.00000004, 1, 
													 localTempDir, randSeed + i, chrLen[i]))
	# load in and subsample
	haplo_list[[i]] <- read_scrm_transpose_for_AlphaSimR(paste0(localTempDir, "chr.txt"),
																											 numLoci = qtlPerChr[i] + neutralPerChr[i], min_maf = 0.05, numLines = 20000, incr = 1/chrLen[i])
	genMap[[i]] <- as.numeric(gsub("^pos_", "", colnames(haplo_list[[i]])))
	# remove simulated haplotypes after loading in
	file.remove(paste0(localTempDir, "chr.txt"))
}

print(Sys.time())
print("end simulating founders")

if(any(sapply(haplo_list, nrow) != (nFound * 2))) stop("mismatch in number of founders and number of haplotypes from scrm")

# For quick desktop testing
# founderPop <- runMacs2(
# 	nInd = 100,
# 	nChr = 1,
# 	segSites = 1000, # a discovery set that you then select SNPs from
# 	Ne = 1000,
# 	bp = 7e+07,
# 	genLen = 1,
# 	mutRate = 4e-08, # higher mutation rate to reflect presumed higher rate in oysters
# 	histNe = c(1000, 1e+05),
# 	histGen = c(100, 1e+06),
# 	inbred = FALSE,
# 	split = NULL,
# 	ploidy = 2,
# 	returnCommand = FALSE,
# 	nThreads = NULL
# )
# SP$addTraitA(nQtlPerChr = 10)
# SP$addSnpChip(nSnpPerChr = 100)


founderPop <- newMapPop(genMap=genMap, haplotypes=haplo_list)
SP <- SimParam$new(founderPop)
SP$setTrackPed(isTrackPed = TRUE) # have AlphaSimR maintain pedigree records
SP$addTraitA(nQtlPerChr = qtlPerChr)
SP$setVarE(h2 = 0.3) # in the range of heritability for growth, meat yield, survival, etc
SP$setSexes("yes_sys") # at the time of breeding, all individuals will only be one sex
SP$addSnpChip(nSnpPerChr = neutralPerChr) # all non-QTL SNPs saved from simulation

pop <- list()
# pull founders from simulated pop while avoiding full and half sibs
pop[[1]] <- newPop(founderPop)
snpGen <- pullSnpGeno(pop[[1]]) # founder SNP genotypes

# save allele freqs in base pop for each panel for calculation of G
baseAlleleFreqs <- colSums(snpGen) / (2*nrow(snpGen))

# write out base pop freqs for blupf90
write.table(cbind(1:length(baseAlleleFreqs), baseAlleleFreqs), 
						paste0(localTempDir, "baseFreqs.txt"), 
						sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)
# write out parameter file for renumf90
cat("DATAFILE
", localTempDir, "f90dat.txt
TRAITS
3
FIELDS_PASSED TO OUTPUT

WEIGHT(S)

RESIDUAL_VARIANCE
2.0
EFFECT          # first fixed effect, overall mean
2 cross numer
EFFECT           # first random effect (animal)
1 cross alpha
RANDOM           ## additive effect without pedigree
animal
SNP_FILE         ## SNP marker file
", localTempDir, "f90snp.txt
(CO)VARIANCES    ## its variance component
1.0
OPTION use_yams
OPTION AlphaBeta 0.99 0.01
OPTION tunedG 0
OPTION whichG 1 # vanRaden 2008
OPTION whichfreq 0 # use freqs from file
OPTION FreqFile ", localTempDir, "baseFreqs.txt # file with frequencies (in same order as genotypes)
OPTION whichfreqScale 0 # use freqs from file
OPTION minfreq 0.0 # turning off all filters and checks
OPTION monomorphic 0
OPTION verify_parentage 0
OPTION no_quality_control
OPTION num_threads_pregs 2 # number of threads
OPTION threshold_duplicate_samples 100 # effectively ignore
OPTION high_threshold_diagonal_g 2 # effectively ignore
OPTION low_threshold_diagonal_g 0.5 # effectively ignore
", file=paste0(localTempDir, "renum.txt"), sep = "")


# initial spawning
pop[[2]] <- randCross(pop[[1]], nCrosses = nFound/2, nProgeny = nOffspringPerCross, balance = TRUE)
# save.image("testing.rda")
# load("testing.rda")
trainPhenos <- data.frame()
gebvRes <- data.frame()
saveGenGain <- data.frame()
for(gen in 1:nGenerations){
	print(Sys.time())
	print(paste("begin gen: ", gen))
	# phenotype training pop (sibs) of current generation and add to phenotype data set
	trainPhenos <- rbind(trainPhenos, sibTestEqual(fam = pop[[gen + 1]], propTest = 0.6)) # phenotype 30, select from 20
	
	# get all genotypes
	g <- pullSnpGeno(pop[[1]])
	for(j in 2:length(pop)) g <- rbind(g, pullSnpGeno(pop[[j]]))
	
	# calc GEBVs
	p <- data.frame(id = rownames(g)) %>% 
		left_join(trainPhenos %>% select(id, Trait_1) %>% rename(pheno = Trait_1), by = "id") # hard coded for first trait
	
	# write out input for blupf90
	# phenotypes
	# coding so that all phenotypes are above 100, missing is 0, and including an overall mean
	p %>% mutate(pheno = pheno + abs(min(min(pheno, na.rm = TRUE), 0)) + 100, mu = 1) %>%
		filter(!is.na(pheno)) %>% select(id, mu, pheno) %>%
		write.table(paste0(localTempDir, "f90dat.txt"), sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)
	# genotypes
	rownames(g) <- paste0(pad_id(rownames(g)), " ")
	write.table(g, paste0(localTempDir, "f90snp.txt"), sep = "", col.names = FALSE, row.names = TRUE, quote = FALSE)
	rownames(g) <- gsub(" ", "", rownames(g)) # undo padding for blupf90 input
	# estimate gebvs with airemlf90
	system2(command = "bash", args = c("run_blupf90.sh", localTempDir))
	
	# load in solutions
	sol <- read.table(paste0(localTempDir, "solutions"), row.names = NULL, skip = 1) %>%
		filter(V2 == 2) # get only animal effect
	xref <- read.table(paste0(localTempDir, "f90snp.txt_XrefID"), row.names = NULL)
	sol$levelNew <- xref$V2[match(sol$V3, xref$V1)] # append original name to solutions
	
	# NOTE: only using _current_ generation to calculate accuracy of gebvs
	comp <- data.frame(id = pop[[gen + 1]]@id, gv = gv(pop[[gen + 1]])) %>% 
		left_join(data.frame(id = as.character(sol$levelNew), gebv = sol$V4), by = "id") %>%
		left_join(p, by = "id")
	# calc accuracy of prediction and save
	gebvRes <- gebvRes %>% rbind(data.frame(genNum = gen,
																					useCryo = useCryo, 
																					# only not phenotyped for calcuation of accuracy
																					acc = cor(comp$gv[is.na(comp$pheno)], comp$gebv[is.na(comp$pheno)])))
	
	# make selections
	print(Sys.time())
	print("begin ocs")
	
	# OCS with lagrangian
	if(useCryo){
		# here we will consider a cryo strategy of cryopreserving all males that are spawned
		# candidates are all males that were previously spawned and all 
		# non-phenotyped individuals in the current generation
		
		# all males previously spawned
		allSires <- unique(SP$pedigree[,"father"])
		allSires <- allSires[allSires != 0] # remove founder placeholder
		
		# all current selection candiates
		selCands <- pop[[gen + 1]]@id
		selCands <- selCands[!(selCands %in% trainPhenos$id)]
		ocsData <- data.frame(Indiv = pop[[gen + 1]]@id, Sex = if_else(pop[[gen + 1]]@sex == "M", "male", "female")) %>%
			filter(Indiv %in% selCands) %>% 
			# and add cryo individuals
			bind_rows(data.frame(Indiv = as.character(allSires), Sex = "male")) %>%
			left_join(data.frame(Indiv = as.character(sol$levelNew), gebv = sol$V4), by = "Indiv")
		
	} else {
		# not cryopreserving, so only selection candidates are the non-phenotyped
		# individuals in the current generation
		selCands <- pop[[gen + 1]]@id
		selCands <- selCands[!(selCands %in% trainPhenos$id)]
		ocsData <- data.frame(Indiv = pop[[gen + 1]]@id, Sex = if_else(pop[[gen + 1]]@sex == "M", "male", "female")) %>%
			left_join(data.frame(Indiv = as.character(sol$levelNew), gebv = sol$V4), by = "Indiv") %>%
			filter(Indiv %in% selCands)
	}

	# make G for coancestry and inbreeding coefficients
	Amat <- createG(g, af = baseAlleleFreqs)[ocsData$Indiv,ocsData$Indiv]
	matingPlan <- runOCS(ocsData = ocsData, Gmat = Amat[ocsData$Indiv,ocsData$Indiv], 
											 N = nFound / 2, Ne = 50)
	print(Sys.time())
	print("end ocs")
	
	# Make combined pop to allow crosses across generations
	tempCombBreedingPop <- lowMem_mergePops(pop_list = pop, indivs = unique(c(matingPlan[,1], matingPlan[,2])))
	
	# save mean genetic value and mean coancestry (weighted by contributions)
	realizedContrib <- data.frame(id = c(matingPlan[,1], matingPlan[,2])) %>%
		count(id) %>% mutate(contrib = n / sum(n)) %>% select(id, contrib)
	genVal <- data.frame(id = tempCombBreedingPop@id, gv = gv(tempCombBreedingPop)) %>%
		right_join(realizedContrib, by = "id") %>% mutate(gv = gv * contrib) %>%
		pull(gv) %>% mean()
	saveGenGain <- saveGenGain %>% rbind(data.frame(genNum = gen,
																									useCryo = useCryo, 
																									genValue = genVal,
																									coAncestry = as.vector(
																										matrix(realizedContrib$contrib, nrow = 1) %*%
																											(Amat[realizedContrib$id, realizedContrib$id] / 2) %*%
																											matrix(realizedContrib$contrib, ncol = 1)
																									)
	))
	rm(Amat) # save some memory
	# create next generation
	pop[[gen + 2]] <- makeCross(tempCombBreedingPop, 
															crossPlan = as.matrix(matingPlan[,1:2]), 
															nProgeny = nOffspringPerCross)
}

# save results
save.image(paste0("/90daydata/oyster_gs_sim/cryo/cryo_", iterationNumber, "_", useCryo, ".rda"))
# for low memory use
# save(snpGen, saveGenGain, gebvRes, file = paste0("cryo_", iterationNumber, "_", useCryo, ".rda"))
