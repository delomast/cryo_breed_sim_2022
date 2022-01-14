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
library(rrBLUP, lib.loc="/project/oyster_gs_sim/R_packages/4.1/")
library(optiSel, lib.loc="/project/oyster_gs_sim/R_packages/4.1/")
library(AllocateMate, lib.loc="/project/oyster_gs_sim/R_packages/4.1/")

source("utils.R")

cmdArgs <- commandArgs(trailingOnly=TRUE)

# Script parameters given on the command line
#' @param randSeed random seed to set R's random number generator
#' @param iterationNumber used to set unique file output names
randSeed <- cmdArgs[1]
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

# initial spawning
pop[[2]] <- randCross(pop[[1]], nCrosses = nFound/2, nProgeny = nOffspringPerCross, balance = TRUE)

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
	Amat <- createG(g = g, af = baseAlleleFreqs) # G with first method of VanRaden (2008), also Endelman and Jannink (2012)
	p <- data.frame(id = rownames(Amat)) %>% 
		left_join(trainPhenos %>% select(id, Trait_1) %>% rename(pheno = Trait_1), by = "id") # hard coded for first trait
	# predict genomic breeding values
	gebv <- kin.blup(data = p, geno = "id", pheno = "pheno", K = Amat)
	# NOTE: only using _current_ generation to calculate accuracy of gebvs
	comp <- data.frame(id = pop[[gen + 1]]@id, gv = gv(pop[[gen + 1]])) %>% 
		left_join(data.frame(id = names(gebv$g), gebv = gebv$g), by = "id") %>%
		left_join(p, by = "id")
	# calc accuracy of prediction and save
	gebvRes <- gebvRes %>% rbind(data.frame(genNum = gen,
																					useCryo = useCryo, 
																					acc = cor(comp$gv[is.na(comp$pheno)], comp$gebv[is.na(comp$pheno)])))
	
	# make selections
	print(Sys.time())
	print("begin ocs")
	# OCS with lagrangian
	selCands <- comp %>% filter(is.na(pheno)) %>% pull(id)
	ocsData <- data.frame(Indiv = pop[[gen + 1]]@id, Sex = if_else(pop[[gen + 1]]@sex == "M", "male", "female")) %>%
		left_join(data.frame(Indiv = names(gebv$g), gebv = gebv$g), by = "Indiv") %>%
		filter(Indiv %in% selCands)
	matingPlan <- runOCS(ocsData = ocsData, Gmat = Amat[ocsData$Indiv,ocsData$Indiv], 
											 N = nFound / 2, Ne = 50)
	print(Sys.time())
	print("end ocs")
	
	# save mean genetic value and mean coancestry (weighted by contributions)
	realizedContrib <- data.frame(id = c(matingPlan[,1], matingPlan[,2])) %>%
		count(id) %>% mutate(contrib = n / sum(n)) %>% select(id, contrib)
	genVal <- data.frame(id = pop[[gen + 1]]@id, gv = gv(pop[[gen + 1]])) %>%
		right_join(realizedContrib, by = "id") %>% mutate(gv = gv * contrib) %>%
		pull(gv) %>% mean()
	saveGenGain <- saveGenGain %>% rbind(data.frame(genNum = gen,
																					useCryo = useCryo, 
																					genValue = genVal,
																					coAncestry = as.vector(
																						matrix(realizedContrib$contrib, nrow = 1) %*%
																							(Amat[realizedContrib$id, realizedContrib$id] / 2) %*%
																							matrix(realizedContrib$contrib, ncol = 1)
																					)))
	# create next generation
	pop[[gen + 2]] <- makeCross(pop[[gen + 1]], crossPlan = as.matrix(matingPlan[,1:2]), nProgeny = nOffspringPerCross)
	

	
}

# save results
save.image(paste0("/90daydata/oyster_gs_sim/cryo/cryo_", iterationNumber, "_", useCryo, ".rda"))
# for low memory use
# save(snpGen, saveGenGain, gebvRes, file = paste0("cryo_", iterationNumber, "_", useCryo, ".rda"))
