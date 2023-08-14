# custom functions - ancestral state reconstruction

# plot ancestral probabilities of tool use on mcc tree
plotASR <- function(d, postFull, MCCphylogeny, type, filename) {
  # get post probs of tool use at tips
  probs <-
    tibble(
      species = rep(1:nrow(d), each = length(postFull$aP[,1])),
      aP1 = rep(postFull$aP[,1], times = nrow(d)),
      aP2 = rep(postFull$aP[,2], times = nrow(d)),
      bEQ = rep(postFull$bEQ, times = nrow(d)),
      k = as.vector(postFull$k),
      Feeding = rep(d$Feeding, each = length(postFull$aP[,1])),
      EQ = rep(as.numeric(scale(d$EQ)), each = length(postFull$aP[,1]))
    ) %>%
    mutate(
      aP = ifelse(Feeding == "Generalist", aP1, aP2),
      p = inv_logit(-(aP + k + bEQ*EQ))
    ) %>%
    group_by(species) %>%
    summarise(p = median(p)) %>%
    ungroup()
  # get input data for asr
  if (type == 1) {
    # literature only
    X <- data.frame(Present = d$LiteratureExists_Binary)
  } else if (type == 2) {
    # literature and videos
    X <- data.frame(Present = d$ToolUse)
  } else if (type == 3) {
    # literature, videos, and model probabilities
    X <- data.frame(Present = ifelse(d$ToolUse == 1, 1, probs$p))
  }
  X$Absent <- 1 - X$Present
  rownames(X) <- d$Name
  # run ancestral state for mcc tree
  set.seed(1000)
  asr <- ancThresh(MCCphylogeny, X, ngen = 1e+06, model = "OU")
  # abbreviated species names for labels
  X <- rownames_to_column(X, var = "Name")
  makeLabel <- function(x) paste0(substring(x, 1, 1), ". ", str_split(x, "_")[[1]][2])
  for (i in 1:length(d$Name)) {
    d$Name[i] <- makeLabel(d$Name[i])
    X$Name[i] <- makeLabel(X$Name[i])
    MCCphylogeny$tip.label[i] <- makeLabel(MCCphylogeny$tip.label[i])
  }
  # plot ancestral states
  p <-
    ggtree(MCCphylogeny) %<+% X +
    geom_tiplab(size = 2) +
    geom_tippoint(aes(colour = Present)) +
    scale_colour_gradient(low = "#E0E0E0", high = "#F8766D") +
    theme(legend.position = "none")
  out <-
    nodepie(
      data = rownames_to_column(as.data.frame(asr$ace), var = "node"),
      cols = 2:3,
      color = c("#E0E0E0", "#F8766D")
    ) %>%
    inset(p, insets = ., height = 0.04, width = 0.04, hjust = 0.5, vjust = -0.5) +
    xlim(0, 60)
  # save
  ggsave(out, filename = filename, width = 10, height = 13)
  return(out)
}

# fit ancestral state reconstruction over posterior phylogenies
fitASR <- function(d, postFull, phylogeny, iter, type) {
  # get post probs of tool use at tips
  probs <-
    tibble(
      species = rep(1:nrow(d), each = length(postFull$aP[,1])),
      aP1 = rep(postFull$aP[,1], times = nrow(d)),
      aP2 = rep(postFull$aP[,2], times = nrow(d)),
      bEQ = rep(postFull$bEQ, times = nrow(d)),
      k = as.vector(postFull$k),
      Feeding = rep(d$Feeding, each = length(postFull$aP[,1])),
      EQ = rep(as.numeric(scale(d$EQ)), each = length(postFull$aP[,1]))
    ) %>%
    mutate(
      aP = ifelse(Feeding == "Generalist", aP1, aP2),
      p = inv_logit(-(aP + k + bEQ*EQ))
    ) %>%
    group_by(species) %>%
    summarise(p = median(p)) %>%
    ungroup()
  # get input data for asr
  if (type == 1) {
    # literature only
    X <- data.frame(Present = d$LiteratureExists_Binary)
  } else if (type == 2) {
    # literature and videos
    X <- data.frame(Present = d$ToolUse)
  } else if (type == 3) {
    # literature, videos, and model probabilities
    X <- data.frame(Present = ifelse(d$ToolUse == 1, 1, probs$p))
  }
  X$Absent <- 1 - X$Present
  rownames(X) <- d$Name
  # run ancestral state for mcc tree
  set.seed(1)
  asr <- ancThresh(phylogeny[[iter]], X, ngen = 1e+06, model = "OU")
  out <- 
    asr$ace %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "node") %>%
    mutate(iter = iter)
  return(out)
}

# get posterior distribution for probability at particular ancestral node
getProbMRCA <- function(fitASR, phylogeny, species) {
  fitASR %>%
    # get most recent common ancestor of species for each phylogeny
    mutate(mrca = purrr::map(iter, function(x) getMRCA(phylogeny[[x]], tip = species))) %>%
    # filter to rows where node == mcra
    filter(node == mrca) %>%
    # extract probability that tool use is present at that node
    pull(Present)
}

# table of ancestral probabilities
getTableASR <- function(amazonaMRCA1, amazonaMRCA2, amazonaMRCA3,
                        cacatuaMRCA1, cacatuaMRCA2, cacatuaMRCA3,
                        nestorMRCA1, nestorMRCA2, nestorMRCA3,
                        poicephalusMRCA1, poicephalusMRCA2, poicephalusMRCA3) {
  # vectors
  genera <- c("\\textit{Amazona}", "\\textit{Cacatua}", "\\textit{Nestor}", "\\textit{Poicephalus}")
  datasets <- c("Pre-video-survey", "Post-video-survey", "Survival cure probabilities")
  # produce table
  tibble(
    Genus = factor(rep(genera, each = length(cacatuaMRCA1)*3), levels = genera),
    Dataset = factor(rep(datasets, each = length(cacatuaMRCA1), times = 4), levels = datasets),
    Value = c(amazonaMRCA1, amazonaMRCA2, amazonaMRCA3,
              cacatuaMRCA1, cacatuaMRCA2, cacatuaMRCA3,
              nestorMRCA1, nestorMRCA2, nestorMRCA3,
              poicephalusMRCA1, poicephalusMRCA2, poicephalusMRCA3)
  ) %>%
    group_by(Genus, Dataset) %>%
    summarise(
      Probability = 
        paste0(
          format(round(median(Value), 2), nsmall = 2),
          ", 95% CI [",
          format(round(quantile(Value, 0.025), 2), nsmall = 2),
          " ",
          format(round(quantile(Value, 0.975), 2), nsmall = 2),
          "]"
          ),
      .groups = "drop"
      ) %>%
    pivot_wider(
      names_from = "Dataset",
      values_from = "Probability"
    )
}

