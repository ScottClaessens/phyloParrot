library(targets)
library(tarchetypes)
library(tidyverse)
library(future)
library(future.callr)
source("R/functionsAncestralState.R")
source("R/functionsComparisonModels.R")
source("R/functionsDAG.R")
source("R/functionsData.R")
source("R/functionsPlotPhylo.R")
source("R/functionsSignal.R")
source("R/functionsSurvivalCure.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("tidyverse", "ape", "brms", "cowplot", "dagitty", 
                            "geiger", "ggdag", "ggridges", "ggstance", "ggtree", 
                            "ggtreeExtra", "magick", "papaja", "phangorn", 
                            "phytools", "readxl", "reshape2", "rethinking", 
                            "rstan", "scales", "sdamr", "tidybayes"))
options(
  clustermq.scheduler = "slurm", 
  clustermq.template = "slurm_clustermq.tmpl"
)
plan(callr)

# cross-validation static branching (see below)
cvTargets <-
  tar_map(
    # data points to map over (25 data points with observed tool use)
    values = tibble(cv = c(83, 85, 87, 89, 102, 109, 110, 111, 114,
                           138, 140, 141, 142, 143, 144, 145, 147,
                           154, 156, 158, 160, 161, 162, 164, 172)),
    # fit cv model to data with data point edited, dynamic branching iterating over trees
    tar_target(crossValSurvCureInd, crossValSurvCureModel(d, phylogeny, iter, survCureFullCompiled, cv), pattern = map(iter), iteration = "list"),
    # combine models
    tar_target(crossValSurvCure, rstan::sflist2stanfit(lapply(crossValSurvCureInd, function(x) x))),
    # posterior samples
    tar_target(crossValPost, extract.samples(crossValSurvCure)),
    # did the cross-validation identify the tool-user?
    tar_target(crossVal, getCrossVal(d, crossValPost, cv))
  )

# full workflow
list(
  
  ######################
  # Load and plot data #
  ######################
  
  # load phylogeny
  tar_target(phylogeny_file, "data/parrotTree.nex", format = "file"),
  tar_target(phylogeny, .compressTipLabel(read.nexus(phylogeny_file))),
  # get maximum clade credibility tree
  tar_target(MCCphylogeny, mcc(phylogeny)),
  # load crowdsourcing data
  tar_target(crowdsourcing_file, "data/parrotToolUseYoutubeSurveyData.xlsx"),
  tar_target(crowdData, read_xlsx(crowdsourcing_file, sheet = 1, range = "A1:O117")),
  # load data
  tar_target(data_file, "data/parrotData.csv", format = "file"),
  tar_target(litCount_file, "data/parrotLiteratureSearch.csv", format = "file"),
  tar_target(d, loadData(data_file, litCount_file, crowdData, phylogeny)),
  # plot maximum clade credibility tree with data
  tar_target(plotPhylo1, plotPhylogeny1(d, MCCphylogeny)),
  tar_target(plotPhylo2, plotPhylogeny2(d, MCCphylogeny)),
  tar_target(plotPhylo3, plotPhylogeny3(d, MCCphylogeny)),
  tar_target(plotPhylo4, plotPhylogeny4(d, MCCphylogeny)),
  # phylogenetic tree indexes to iterate over in analyses
  tar_target(iter, getIter(100)),
  
  #######################
  # Phylogenetic signal #
  #######################
  
  # phylogenetic signal
  tar_target(phySignal1, fitPhySignal(d, phylogeny, iter, dv = "LiteratureExists_Binary"), pattern = map(iter)),
  tar_target(phySignal2, fitPhySignal(d, phylogeny, iter, dv = "ToolUse"), pattern = map(iter)),
  
  ######################
  # DAG - Causal model #
  ######################
  
  tar_target(plotDAG, drawDAG()),
  
  ####################################
  # Phylogenetic survival cure model #
  ####################################
  
  # stan files
  tar_target(survCureFull_file, "stan/survCureModelFull.stan", format = "file"),
  tar_target(survCureFullSensitivity_file, "stan/survCureModelFullSensitivity.stan", format = "file"),
  tar_target(survCureReduced_file, "stan/survCureModelReduced.stan", format = "file"),
  # compile models
  tar_target(survCureFullCompiled, stan_model(survCureFull_file, model_name = "survCureFull")),
  tar_target(survCureFullSensitivityCompiled, stan_model(survCureFullSensitivity_file, model_name = "survCureFullSensitivity")),
  tar_target(survCureReducedCompiled, stan_model(survCureReduced_file, model_name = "survCureReduced")),
  ## simulate model and recover parameters
  tar_target(survCureSimSeed, 1:100),
  tar_target(survCureSim, simSurvCureModel(survCureFullCompiled, survCureSimSeed), pattern = map(survCureSimSeed)),
  tar_target(plotSurvCureSim, plotSurvSim(survCureSim)),
  # fit models to data, dynamic branching iterating over trees
  tar_target(survCureFullInd, fitSurvCureModel(d, phylogeny, iter, survCureFullCompiled), pattern = map(iter), iteration = "list"),
  tar_target(survCureFullSensitivityInd, fitSurvCureModel(d, phylogeny, iter, survCureFullSensitivityCompiled), pattern = map(iter), iteration = "list"),
  tar_target(survCureReducedInd, fitSurvCureModel(d, phylogeny, iter, survCureReducedCompiled), pattern = map(iter), iteration = "list"),
  # combine models
  tar_target(survCureFull, rstan::sflist2stanfit(lapply(survCureFullInd, function(x) x))),
  tar_target(survCureFullSensitivity, rstan::sflist2stanfit(lapply(survCureFullSensitivityInd, function(x) x))),
  tar_target(survCureReduced, rstan::sflist2stanfit(lapply(survCureReducedInd, function(x) x))),
  # prior samples
  tar_target(priorFull, getPrior(iter, phylogeny)),
  # posterior samples
  tar_target(postFull, extract.samples(survCureFull)),
  tar_target(postFullSensitivity, extract.samples(survCureFullSensitivity)),
  tar_target(postReduced, extract.samples(survCureReduced)),
  # sum of probabilities
  tar_target(sumProbs, getSumProbs(d, postFull)),
  # trace plot
  tar_target(plotTrace, plotTraceMCMC(postFull, numTrees = length(iter), numChains = 4)),
  # main plots
  tar_target(plotSurvCure1, plotSpeciesProb(d, postFull)),
  tar_target(plotSurvCure2, plotSurvival(postFull)),
  tar_target(plotSurvCure3, plotReduced(d, postReduced)),
  tar_target(plotSurvCure4, plotSurvPred(d, postFull)),
  tar_target(plotSurvCure5, plotGP1(postFull)),
  tar_target(plotSurvCure6, plotGP2(postFull, MCCphylogeny)),
  tar_target(plotSurvCure7, plotSpeciesProb(d, priorFull, prior = TRUE)),
  tar_target(plotSurvCure8, plotSensitivity(d, postFull, postFullSensitivity)),
  tar_target(plotSurvCure9, plotROCCurve(d, postFull)),
  # calculate area under the curve
  tar_target(auc, calculateAUC(d, postFull)),
  # cross-validation
  cvTargets,
  tar_combine(cv, cvTargets[["crossVal"]]),
  
  #####################
  # Comparison models #
  #####################
  
  # pre-video-survey data
  tar_target(m1Ind, fitComparisonModel(d, phylogeny, iter, var = "LiteratureExists_Binary"), pattern = map(iter), iteration = "list"),
  tar_target(m1, rstan::sflist2stanfit(lapply(m1Ind, function(x) x@stanfit))),
  tar_target(postM1, extract.samples(m1)),
  # post-video-survey data
  tar_target(m2Ind, fitComparisonModel(d, phylogeny, iter, var = "ToolUse"), pattern = map(iter), iteration = "list"),
  tar_target(m2, rstan::sflist2stanfit(lapply(m2Ind, function(x) x@stanfit))),
  tar_target(postM2, extract.samples(m2)),
  # plot comparison
  tar_target(plotComparison, plotComparisons(postM1, postM2, postFull)),
  
  ##################################
  # Ancestral state reconstruction #
  ##################################

  # plot ancestral probabilities of tool use on MCC tree
  tar_target(plotASR1, plotASR(d, postFull, MCCphylogeny, type = 1, filename = "figures/asr1.pdf")),
  tar_target(plotASR2, plotASR(d, postFull, MCCphylogeny, type = 2, filename = "figures/asr2.pdf")),
  tar_target(plotASR3, plotASR(d, postFull, MCCphylogeny, type = 3, filename = "figures/asr3.pdf")),
  # fit ancestral state model iterating over 100 posterior phylogenies
  tar_target(fitASR1, fitASR(d, postFull, phylogeny, iter, type = 1), pattern = map(iter)),
  tar_target(fitASR2, fitASR(d, postFull, phylogeny, iter, type = 2), pattern = map(iter)),
  tar_target(fitASR3, fitASR(d, postFull, phylogeny, iter, type = 3), pattern = map(iter)),
  # most recent common ancestor of Cacatua genus
  tar_target(cacatuaSpecies, c("Cacatua_tenuirostris","Cacatua_sanguinea","Cacatua_ducorpsii","Cacatua_goffiniana",
                               "Cacatua_haematuropygia","Cacatua_sulphurea","Cacatua_galerita","Cacatua_moluccensis",
                               "Cacatua_pastinator","Cacatua_alba","Cacatua_ophthalmica")),
  tar_target(cacatuaMRCA1, getProbMRCA(fitASR1, phylogeny, cacatuaSpecies)),
  tar_target(cacatuaMRCA2, getProbMRCA(fitASR2, phylogeny, cacatuaSpecies)),
  tar_target(cacatuaMRCA3, getProbMRCA(fitASR3, phylogeny, cacatuaSpecies)),
  # most recent common ancestor of Amazona genus
  tar_target(amazonaSpecies, c("Amazona_aestiva","Amazona_albifrons","Amazona_amazonica","Amazona_auropalliata",
                               "Amazona_autumnalis","Amazona_farinosa","Amazona_leucocephala","Amazona_ochrocephala",
                               "Amazona_oratrix","Amazona_ventralis","Amazona_vinacea","Amazona_viridigenalis",
                               "Amazona_xantholora")),
  tar_target(amazonaMRCA1, getProbMRCA(fitASR1, phylogeny, amazonaSpecies)),
  tar_target(amazonaMRCA2, getProbMRCA(fitASR2, phylogeny, amazonaSpecies)),
  tar_target(amazonaMRCA3, getProbMRCA(fitASR3, phylogeny, amazonaSpecies)),
  # most recent common ancestor of Poicephalus genus
  tar_target(poicephalusSpecies, c("Poicephalus_cryptoxanthus","Poicephalus_gulielmi","Poicephalus_meyeri",
                                   "Poicephalus_robustus","Poicephalus_senegalus")),
  tar_target(poicephalusMRCA1, getProbMRCA(fitASR1, phylogeny, poicephalusSpecies)),
  tar_target(poicephalusMRCA2, getProbMRCA(fitASR2, phylogeny, poicephalusSpecies)),
  tar_target(poicephalusMRCA3, getProbMRCA(fitASR3, phylogeny, poicephalusSpecies)),
  # most recent common ancestor of Nestor genus
  tar_target(nestorSpecies, c("Nestor_notabilis","Nestor_meridionalis")),
  tar_target(nestorMRCA1, getProbMRCA(fitASR1, phylogeny, nestorSpecies)),
  tar_target(nestorMRCA2, getProbMRCA(fitASR2, phylogeny, nestorSpecies)),
  tar_target(nestorMRCA3, getProbMRCA(fitASR3, phylogeny, nestorSpecies)),
  # table of ancestral probabilities
  tar_target(tableASR, getTableASR(amazonaMRCA1, amazonaMRCA2, amazonaMRCA3,
                                   cacatuaMRCA1, cacatuaMRCA2, cacatuaMRCA3,
                                   nestorMRCA1, nestorMRCA2, nestorMRCA3,
                                   poicephalusMRCA1, poicephalusMRCA2, poicephalusMRCA3)),
  
  ###############
  ## Manuscript #
  ###############
  
  tar_render(manuscript, "manuscript.Rmd"),
  
  
  ##################
  ## Session info ##
  ##################
  
  # print session info for reproducibility
  tar_target(sessionInfo, writeLines(capture.output(sessionInfo()), "sessionInfo.txt"))
  
)
