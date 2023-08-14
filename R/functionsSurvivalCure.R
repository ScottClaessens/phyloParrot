# custom functions - survival cure model

# simulation
simSurvCureModel <- function(model, seed) {
  # set seed
  set.seed(seed)
  # number of species
  n <- 100
  # simulate tree
  phylo <- rtree(n)
  # scaled distance matrix from phylogeny
  distMat <- cophenetic.phylo(phylo)
  distMat <- distMat / max(distMat)
  # generate standardised EQ
  EQ <- rnorm(n)
  # generate feeding type
  f <- ifelse(rbinom(n, 1, prob = 0.5) == 1, "Generalist", "Specialist")
  # set simulation parameters
  etasq <- 2
  rhosq <- 2
  aP1 <- 1
  aP2 <- -1
  bEQ <- 1.5
  # simulate tool use
  # sigma from gaussian process
  sigma <- etasq*exp(-rhosq*distMat)
  k <- rmvnorm(1, sigma = sigma)[1,]
  # log odds probability of tool use
  p <- ifelse(f == "Generalist", aP1, aP2) + k + bEQ*EQ
  # simulate tool use
  t <- rbern(n, inv_logit(p))
  # simulate time-to-event #1 (number of publications)
  # https://inla.r-inla-download.org/r-inla.org/doc/likelihood/weibullcure.pdf
  lambda1 <- 1
  censorTime1 <- runif(n, 0, 5)
  y1 <- round(rexp(n, rate = lambda1), 0)
  # censoring if (1) not measured yet or (2) not a tool user
  censoredEvent1 <- (y1 > censorTime1) | (t == 0)
  # censored observations
  yObs1 <- y1
  yObs1[censoredEvent1] <- censorTime1[censoredEvent1]
  event1 <- as.numeric(!censoredEvent1)
  # simulate time-to-event #1 (number of videos)
  lambda2 <- 0.01
  censorTime2 <- runif(n, 0, 600)
  y2 <- rexp(n, rate = lambda2)
  # censoring if (1) not measured yet or (2) not a tool user
  censoredEvent2 <- (y2 > censorTime2) | (t == 0)
  # censored observations
  yObs2 <- y2
  yObs2[censoredEvent2] <- censorTime2[censoredEvent2]
  event2 <- as.numeric(!censoredEvent2)
  # data list for stan
  dat_list <- list(
    N        = n,
    eventLit = as.integer(event1),
    eventVid = as.integer(event2),
    numLit   = as.integer(yObs1),
    numVid   = as.integer(yObs2),
    Species  = 1:n,
    EQ       = EQ,
    Feeding  = as.integer(ifelse(f == "Generalist", 1, 2)),
    distMat  = distMat
  )
  # fit model
  m <- sampling(model, data = dat_list, pars = c("aLit", "aVid", "aP", "k", "bEQ", "etasq", "rhosq"),
                iter = 4000, chains = 1, cores = 1, seed = 2113, control = list(adapt_delta = 0.85))
  # get posterior
  post <- extract.samples(m)
  # parameters
  out <- tibble(
    seed = rep(seed, 7),
    parameter = c("lambdaLit", "lambdaVid", "aP1", "aP2", "bEQ", "etasq", "rhosq"),
    median = c(
      median(1 / exp(post$aLit)),
      median(1 / exp(post$aVid)),
      median(post$aP[,1]),
      median(post$aP[,2]),
      median(post$bEQ),
      median(post$etasq),
      median(post$rhosq)
      ),
    lower = c(
      quantile(1 / exp(post$aLit), 0.025),
      quantile(1 / exp(post$aVid), 0.025),
      quantile(post$aP[,1], 0.025),
      quantile(post$aP[,2], 0.025),
      quantile(post$bEQ, 0.025),
      quantile(post$etasq, 0.025),
      quantile(post$rhosq, 0.025)
      ),
    upper = c(
      quantile(1 / exp(post$aLit), 0.975),
      quantile(1 / exp(post$aVid), 0.975),
      quantile(post$aP[,1], 0.975),
      quantile(post$aP[,2], 0.975),
      quantile(post$bEQ, 0.975),
      quantile(post$etasq, 0.975),
      quantile(post$rhosq, 0.975)
      )
  )
  return(out)
}

# plot simulation results
plotSurvSim <- function(survCureSim) {
  out <-
    survCureSim %>%
    # parameters in linear model are inverted in posterior
    mutate(median2 = ifelse(parameter %in% c("aP1", "aP2", "bEQ"), -median, median),
           lower2  = ifelse(parameter %in% c("aP1", "aP2", "bEQ"), -upper, lower),
           upper2  = ifelse(parameter %in% c("aP1", "aP2", "bEQ"), -lower, upper)) %>%
    # plot sim results
    ggplot(aes(x = median2, xmin = lower2, xmax = upper2, y = seed)) +
    geom_pointrangeh(fatten = 1) +
    geom_vline(data = tibble(parameter = c("lambdaLit", "lambdaVid", "aP1", "aP2", "bEQ", "etasq", "rhosq"),
                             parSim = c(1, 0.01, 1, -1, 1.5, 2, 2)),
               aes(xintercept = parSim)) +
    facet_wrap(. ~ parameter, scales = "free_x") +
    xlab("Median posterior parameter") +
    theme_classic() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  # save
  ggsave(out, filename = "figures/survCureSim.pdf", width = 8.5, height = 8.5)
  return(out)
}

# fit survival cure model
fitSurvCureModel <- function(d, phylogeny, iter, model) {
  # scaled distance matrix from one phylogeny
  phy <- phylogeny[[iter]]
  distMat <- cophenetic.phylo(phy)
  distMat <- distMat / max(distMat)
  # data list for stan
  dat_list <- list(
    N        = nrow(d),
    eventLit = d$LiteratureExists_Binary,
    eventVid = ifelse(d$VideosFound > 0, 1, 0),
    numLit   = d$EX.Pre.Tool.Papers,
    numVid   = d$SearchHitsPreToolUse,
    Species  = 1:nrow(d),
    EQ       = as.numeric(scale(d$EQ)),
    Feeding  = ifelse(d$Feeding == "Generalist", 1, 2),
    distMat  = distMat
  )
  # parameters to return
  pars <- c("aLit", "aVid", "aP")
  if (model@model_name %in% c("survCureFull","survCureFullSensitivity")) {
    pars <- c(pars, "bEQ", "k", "etasq", "rhosq")
  }
  # fit model
  out <- sampling(
    model,
    data = dat_list,
    pars = pars,
    iter = 4000,
    chains = 1,
    cores = 1,
    seed = 2113,
    control = list(adapt_delta = 0.85)
  )
  return(out)
}

# get prior samples
getPrior <- function(iter, phylogeny, n = 1000) {
  # random seed
  set.seed(2113)
  # number of species
  N <- 174
  # total number of samples
  totalSamps <- n * length(iter)
  # create empty prior list to fill
  out <- list(
    aLit  = rep(NA, totalSamps),
    aVid  = rep(NA, totalSamps),
    aP    = matrix(nrow = totalSamps, ncol = 2),
    bEQ   = rep(NA, totalSamps),
    k     = matrix(nrow = totalSamps, ncol = N),
    etasq = rep(NA, totalSamps),
    rhosq = rep(NA, totalSamps)
    )
  # iterate over phylogenetic uncertainty
  for (treeIndex in 1:length(iter)) {
    # get tree and distance matrix
    phy <- phylogeny[[iter[treeIndex]]]
    distMat <- cophenetic.phylo(phy)
    distMat <- distMat / max(distMat)
    # get prior distribution
    aLit <- rnorm(n)
    aVid <- rnorm(n)
    aP <- matrix(rnorm(n*2), nrow = n, ncol = 2)
    bEQ <- rnorm(n)
    etasq <- rexp(n, 0.5)
    rhosq <- rexp(n, 0.5)
    # get k
    k <- matrix(nrow = n, ncol = N)
    for (a in 1:n) {
      sigma <- matrix(ncol = N, nrow = N)
      for (i in 1:(N-1)) {
        sigma[i,i] <- etasq[a] + 0.01
        for (j in (i+1):N) {
          sigma[i,j] <- etasq[a] * exp(-rhosq[a] * distMat[i,j])
          sigma[j,i] <- sigma[i,j]
        }
      }
      sigma[N,N] <- etasq[a] + 0.01
      k[a,] <- as.vector(t(chol(sigma)) %*% rnorm(N))
    }
    # add to list
    samples <- (((treeIndex-1)*n)+1):(treeIndex*n)
    out$aLit[samples] <- aLit
    out$aVid[samples] <- aVid
    out$aP[samples,] <- aP
    out$bEQ[samples] <- bEQ 
    out$k[samples,] <- k   
    out$etasq[samples] <- etasq
    out$rhosq[samples] <- rhosq
  }
  return(out)
}

# get sum of probabilities
getSumProbs <- function(d, postFull) {
  tibble(
    iter = rep(1:length(postFull$aP[,1]), times = nrow(d)),
    species = rep(1:nrow(d), each = length(postFull$aP[,1])),
    aP1 = rep(postFull$aP[,1], times = nrow(d)),
    aP2 = rep(postFull$aP[,2], times = nrow(d)),
    bEQ = rep(postFull$bEQ, times = nrow(d)),
    k = as.vector(postFull$k),
    Feeding = rep(d$Feeding, each = length(postFull$aP[,1])),
    EQ = rep(as.numeric(scale(d$EQ)), each = length(postFull$aP[,1])),
    ToolUse = rep(d$ToolUse, each = length(postFull$aP[,1]))
    ) %>%
    mutate(
      aP = ifelse(Feeding == "Generalist", aP1, aP2),
      p = inv_logit(-(aP + k + bEQ*EQ))
      ) %>%
    # filter to only non-tool users
    filter(ToolUse == 0) %>%
    # for each iteration, get sum of probs
    group_by(iter) %>%
    summarise(p = sum(p)) %>%
    # get posterior sum of probs
    pull(p)
}

# create traceplot for main parameters
plotTraceMCMC <- function(postFull, numTrees, numChains) {
  # number of iterations
  numIter <- dim(postFull$aLit)[1] / numTrees
  # plot
  out <-
    tibble(
      chain = rep(1:numChains, each = numIter),
      iter  = rep(1:numIter, times = numChains),
      aLit  = postFull$aLit[1:(numChains*numIter)],
      aVid  = postFull$aVid[1:(numChains*numIter)],
      aP1   = postFull$aP[1:(numChains*numIter),1],
      aP2   = postFull$aP[1:(numChains*numIter),2],
      bEQ   = postFull$bEQ[1:(numChains*numIter)],
      etasq = postFull$etasq[1:(numChains*numIter)],
      rhosq = postFull$rhosq[1:(numChains*numIter)]
      ) %>%
    pivot_longer(cols = !c(chain, iter),
                 names_to = "parameter",
                 values_to = "value") %>%
    mutate(chain = factor(chain)) %>%
    ggplot(aes(x = iter, y = value, colour = chain)) +
    geom_line(alpha = 0.3) +
    facet_wrap(~ parameter, scales = "free") +
    labs(x = "Iteration", y = NULL) +
    theme_classic() +
    theme(strip.background = element_blank(),
          legend.position = "none")
  # save plot
  ggsave(out, filename = "figures/trace.pdf", width = 6, height = 6)
  return(out)
}

# plot posterior probabilities of tool use for different species
plotSpeciesProb <- function(d, postFull, prior = FALSE) {
  # abbreviated species names for labels
  makeLabel <- function(x) paste0(substring(x, 1, 1), ". ", str_split(x, "_")[[1]][2])
  for (i in 1:length(d$Name)) d$Name[i] <- makeLabel(d$Name[i])
  # get posterior
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
    mutate(aP = ifelse(Feeding == "Generalist", aP1, aP2))
  # take into account EQ if posterior, otherwise ignore
  if (!prior) {
    probs <- mutate(probs, p = inv_logit(-(aP + k + bEQ*EQ)))
  } else {
    probs <- mutate(probs, p = inv_logit(-(aP + k)))
  }
  # summarise posterior
  probs <-
    probs %>%
    group_by(species) %>%
    summarise(med = median(p),
              low95 = quantile(p, 0.025),
              upp95 = quantile(p, 0.975),
              low50 = quantile(p, 0.250),
              upp50 = quantile(p, 0.750)) %>%
    ungroup() %>%
    # tool use variable
    mutate(ToolUse = ifelse(d$LiteratureExists_Binary == 0 & d$VideosFound == 0, "Absent",
                            ifelse(d$LiteratureExists_Binary == 1 & d$VideosFound == 0, "Present (in literature only)",
                                   ifelse(d$LiteratureExists_Binary == 0 & d$VideosFound  > 0, "Present (in videos only)", "Present (in videos and lit.)"))),
           ToolUse = factor(ToolUse, levels = c("Present (in videos only)",
                                                "Present (in literature only)",
                                                "Present (in videos and lit.)",
                                                "Absent")))
  # if posterior plot, sort by probability of tool use
  if (!prior) probs <- arrange(probs, desc(med))
  # plot
  out <-
    ggplot(probs, aes(x = fct_reorder(factor(species), med), y = med, colour = factor(ToolUse))) +
    geom_pointrange(aes(ymin = low95, ymax = upp95), size = 0.2) +
    geom_linerange(aes(ymin = low50, ymax = upp50), linewidth = 1, show.legend = FALSE) +
    coord_flip() +
    scale_x_discrete(name = NULL, labels = function(x) d$Name[parse_number(x)]) +
    scale_colour_manual(values = c("#D55E00", "#0072B2", "#CC79A7", "grey80")) +
    scale_y_continuous(name = paste0(ifelse(prior, "Prior", "Posterior"), " probability of tool use")) +
    theme_minimal() +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_text(face = "italic", size = 5),
          panel.grid = element_line(colour = "grey96"),
          legend.title = element_blank())
  # save
  filename <- paste0("figures/survCureProbs", ifelse(prior, "Prior", "Post"), ".pdf")
  ggsave(out, filename = filename, width = 7, height = 10)
  return(out)
}

# plot survival parameters
plotSurvival <- function(postFull) {
  pA <-
    tibble(par = exp(postFull$aLit)) %>%
    ggplot(aes(x = par)) +
    stat_halfeye() +
    lims(x = c(0, 100), y = 0:1) +
    labs(x = "Expected number of publications until\ntool use discovery in scientific literature",
         y = "Density") +
    theme_classic()
  pB <-
    tibble(par = exp(postFull$aVid)) %>%
    ggplot(aes(x = par)) +
    stat_halfeye() +
    lims(x = c(0, 10000), y = 0:1) +
    labs(x = "Expected number of videos until\ntool use discovery on YouTube",
         y = "") +
    theme_classic()
  # put together
  out <- plot_grid(pA, pB, nrow = 1)
  # save
  ggsave(out, filename = "figures/survCureSurvival.pdf", width = 7, height = 3)
  return(out)
}

# plot relationship between probability of tool use and 
# number of papers/videos for non-tool-users
plotReduced <- function(d, postReduced) {
  # get median posterior probabilities of tool use from reduced model
  d$probs <- 1 - apply(postReduced$aP, 2, function(x) median(inv_logit(x)))
  # plot
  out <-
    d %>%
    filter(ToolUse == 0) %>% # non-tool-users only
    ggplot(aes(x = log(EX.Pre.Tool.Papers + 1), y = log(SearchHits), colour = probs)) +
    geom_point() +
    labs(
      x = "log(Number of papers published + 1)",
      y = "log(Number of YouTube search hits)",
      colour = "Probability of\nundetected\ntool use"
    ) +
    theme_classic()
  # save
  ggsave(out, filename = "figures/survCureReduced.pdf", height = 4, width = 5)
  return(out)
}

# plot effect of predictors
plotSurvPred <- function(d, postFull) {
  # get sequence for EQ
  seqEQ <- seq(min(as.numeric(scale(d$EQ))),
               max(as.numeric(scale(d$EQ))), length.out = 100)
  # get predictions for EQ (holding feeding type at specialist)
  medianEQ <- sapply(seqEQ, function(x) median(inv_logit(-(postFull$aP[,2] + postFull$bEQ*x))))
  upper50EQ  <- sapply(seqEQ, function(x) quantile(inv_logit(-(postFull$aP[,2] + postFull$bEQ*x)), 0.250))
  lower50EQ  <- sapply(seqEQ, function(x) quantile(inv_logit(-(postFull$aP[,2] + postFull$bEQ*x)), 0.750))
  upper95EQ  <- sapply(seqEQ, function(x) quantile(inv_logit(-(postFull$aP[,2] + postFull$bEQ*x)), 0.025))
  lower95EQ  <- sapply(seqEQ, function(x) quantile(inv_logit(-(postFull$aP[,2] + postFull$bEQ*x)), 0.975))
  # get labels for plot
  d <-
    d %>%
    # tool use label
    mutate(ToolUseLabel = ifelse(d$LiteratureExists_Binary == 0 & d$VideosFound == 0, "Absent",
                                 ifelse(d$LiteratureExists_Binary == 1 & d$VideosFound == 0, "Present (in literature only)",
                                        ifelse(d$LiteratureExists_Binary == 0 & d$VideosFound > 0, "Present (in videos only)", "Present (in videos and lit.)"))),
           ToolUseLabel = factor(ToolUseLabel, levels = c("Present (in videos only)",
                                                          "Present (in literature only)",
                                                          "Present (in videos and lit.)",
                                                          "Absent")))
  # plot (feeding type)
  set.seed(2113)
  pA <-
    ggplot() +
    stat_halfeye(data = tibble(Feeding = rep(c("Generalist", "Specialist"), each = length(postFull$aP[,1])),
                               aP = c(inv_logit(-(postFull$aP[,1])), inv_logit(-(postFull$aP[,2])))),
                 aes(x = Feeding, y = aP), scale = 0.5, fill = "lightgrey", .width = c(0.50, 0.95)) +
    geom_jitter(data = d, aes(x = Feeding, y = ToolUse, colour = ToolUseLabel),
                alpha = 0.3, width = 0.3, height = 0.05) +
    scale_x_discrete(name = "Feeding strategy", expand = expansion(mult = 1)) +
    scale_y_continuous(name = "Probability of tool use", limits = c(-0.05, 1.05)) +
    scale_colour_manual(values = c("#D55E00", "#0072B2", "#CC79A7", "grey")) +
    theme_classic() +
    theme(legend.title = element_blank())
  # plot (EQ)
  pB <-
    tibble(
      EQ = seqEQ,
      median = medianEQ, 
      lower95 = lower95EQ, 
      upper95 = upper95EQ,
      lower50 = lower50EQ, 
      upper50 = upper50EQ
    ) %>%
    ggplot() +
    geom_ribbon(aes(x = EQ, y = median, ymin = lower95, ymax = upper95), fill = "lightgrey") +
    geom_ribbon(aes(x = EQ, y = median, ymin = lower50, ymax = upper50), fill = "grey") +
    geom_line(aes(x = EQ, y = median)) +
    geom_jitter(data = d, aes(x = as.numeric(scale(EQ)), y = ToolUse, colour = ToolUseLabel),
               alpha = 0.3, width = 0, height = 0.05) +
    scale_x_continuous(name = "Encephalisation quotient", 
                       labels = function(x) format(round((x * sd(d$EQ)) + mean(d$EQ), 2), nsmall = 2),
                       breaks = function(x) (seq(0.8, 2.0, by = 0.2) - mean(d$EQ)) / sd(d$EQ)) +
    scale_y_continuous(name = "", limits = c(-0.05, 1.05)) +
    scale_colour_manual(values = c("#D55E00", "#0072B2", "#CC79A7", "grey")) +
    theme_classic() +
    theme(legend.title = element_blank())
  # put together
  out <- plot_grid(pA + theme(legend.position = "none"),
                   pB, nrow = 1, rel_widths = c(0.5, 1))
  # save
  ggsave(out, filename = "figures/survCurePred.pdf", height = 4, width = 8)
  return(out)
}

# plot phylogenetic signal
plotGP1 <- function(postFull) {
  # distance vector
  dist <- seq(0, 1, length.out = 1000)
  # randomly sample 2000 iterations from prior and posterior
  set.seed(2113)
  etasqPrior <- rexp(2000, 0.5)
  rhosqPrior <- rexp(2000, 0.5)
  etasqPost  <- sample(postFull$etasq, 2000)
  rhosqPost  <- sample(postFull$rhosq, 2000)
  # prior gaussian process function
  priorFun <-
    tibble(dist) %>%
    mutate(
      type    = "Prior",
      median  = purrr::map(dist, function(x) median(etasqPrior * exp(-rhosqPrior * x) + 0.01)),
      lower95 = purrr::map(dist, function(x) quantile(etasqPrior * exp(-rhosqPrior * x) + 0.01, 0.025)),
      upper95 = purrr::map(dist, function(x) quantile(etasqPrior * exp(-rhosqPrior * x) + 0.01, 0.975)),
      lower50 = purrr::map(dist, function(x) quantile(etasqPrior * exp(-rhosqPrior * x) + 0.01, 0.250)),
      upper50 = purrr::map(dist, function(x) quantile(etasqPrior * exp(-rhosqPrior * x) + 0.01, 0.750)),
    ) %>%
    unnest(c(median, lower95, upper95, lower50, upper50))
  # posterior gaussian process function
  postFun <-
    tibble(dist) %>%
    mutate(
      type    = "Posterior",
      median  = purrr::map(dist, function(x) median(etasqPost * exp(-rhosqPost * x) + 0.01)),
      lower95 = purrr::map(dist, function(x) quantile(etasqPost * exp(-rhosqPost * x) + 0.01, 0.025)),
      upper95 = purrr::map(dist, function(x) quantile(etasqPost * exp(-rhosqPost * x) + 0.01, 0.975)),
      lower50 = purrr::map(dist, function(x) quantile(etasqPost * exp(-rhosqPost * x) + 0.01, 0.250)),
      upper50 = purrr::map(dist, function(x) quantile(etasqPost * exp(-rhosqPost * x) + 0.01, 0.750))
    ) %>%
    unnest(c(median, lower95, upper95, lower50, upper50))
  # plot both prior and posterior gaussian process function
  out <-
    ggplot(data = bind_rows(priorFun, postFun)) +
    geom_ribbon(aes(x = dist, ymin = lower95, ymax = upper95), alpha = 0.3, fill = "lightsteelblue") +
    geom_ribbon(aes(x = dist, ymin = lower50, ymax = upper50), alpha = 0.6, fill = "lightsteelblue") +
    geom_line(aes(x = dist, y = median), colour = "black") +
    facet_grid(. ~ fct_rev(type)) +
    labs(x = "Scaled phylogenetic distance", y = "Phylogenetic covariance") +
    theme_classic()
  # save
  ggsave(out, filename = "figures/survCureGaussianProcess1.pdf", height = 4, width = 7)
  return(out)
}

# plot implied correlation matrix between species from gaussian process
plotGP2 <- function(postFull, MCCphylogeny) {
  # scaled distance matrix
  distMat <- cophenetic.phylo(MCCphylogeny)
  distMat <- distMat / max(distMat)
  # calculate between-species covariance matrix
  n <- length(MCCphylogeny$tip.label)
  K <- matrix(NA, nrow = n, ncol = n)
  etasqMed <- median(postFull$etasq)
  rhosqMed <- median(postFull$rhosq)
  for (i in 1:n) {
    for (j in 1:n) {
      K[i,j] <- etasqMed * exp(-rhosqMed * distMat[i,j])
    }
  }
  diag(K) <- etasqMed + 0.01
  # convert to correlation matrix
  Rho <- round(cov2cor(K), 2)
  # lower triangle
  Rho[upper.tri(Rho)] <- NA
  # plot correlation matrix
  out <- 
    Rho %>%
    melt() %>%
    drop_na() %>%
    ggplot(aes(x = Var1, y = Var2, fill = value)) + 
    geom_tile() +
    scale_fill_gradient2(name = "Cor", low = "white", high = "red", limit = c(0,1), na.value = "white") +
    annotate("text", x = 7  , y = 13 , angle = 45, size = 3, label = "Psittaculinae") +
    annotate("text", x = 26 , y = 32 , angle = 45, size = 3, label = "Agapornithinae") +
    annotate("text", x = 48 , y = 54 , angle = 45, size = 3, label = "Platycercinae") +
    annotate("text", x = 67 , y = 73 , angle = 45, size = 3, label = "Loriinae") +
    annotate("text", x = 94 , y = 100, angle = 45, size = 3, label = "Cacatuidae") +
    annotate("text", x = 125, y = 131, angle = 45, size = 3, label = "Androglossini") +
    annotate("text", x = 152, y = 158, angle = 45, size = 3, label = "Arini") +
    theme_void() +
    theme(legend.box.margin = margin(0, 5, 0, 0))
  # save
  ggsave(out, filename = "figures/survCureGaussianProcess2.pdf", height = 6, width = 6.5)
  return(out)
}

# plot results of sensitivity analysis with different intercept prior
plotSensitivity <- function(d, postFull, postFullSensitivity) {
  # function to get probabilities
  getProbs <- function(post) {
    tibble(
      species = rep(1:nrow(d), each = length(post$aP[,1])),
      aP1 = rep(post$aP[,1], times = nrow(d)),
      aP2 = rep(post$aP[,2], times = nrow(d)),
      bEQ = rep(post$bEQ, times = nrow(d)),
      k = as.vector(post$k),
      Feeding = rep(d$Feeding, each = length(post$aP[,1])),
      EQ = rep(as.numeric(scale(d$EQ)), each = length(post$aP[,1]))
    ) %>%
      mutate(
        aP = ifelse(Feeding == "Generalist", aP1, aP2),
        p = inv_logit(-(aP + k + bEQ*EQ))
      ) %>%
      group_by(species) %>%
      summarise(p = median(p)) %>%
      ungroup() %>%
      # tool use variable
      mutate(ToolUse = ifelse(d$LiteratureExists_Binary == 0 & d$VideosFound == 0, "Absent",
                              ifelse(d$LiteratureExists_Binary == 1 & d$VideosFound == 0, "Present (in literature only)",
                                     ifelse(d$LiteratureExists_Binary == 0 & d$VideosFound  > 0, "Present (in videos only)", "Present (in videos and lit.)"))),
             ToolUse = factor(ToolUse, levels = c("Present (in videos only)",
                                                  "Present (in literature only)",
                                                  "Present (in videos and lit.)",
                                                  "Absent"))) %>%
      # rank number
      arrange(desc(p)) %>%
      mutate(
        rank = 1:nrow(d),
        species = d$Name[species]
        )
  }
  # get probs and rankings for both models
  probsFull <- getProbs(postFull) %>% rename(pFull = p, rankFull = rank)
  probsFullSensitivity <- getProbs(postFullSensitivity) %>% rename(pFullSensitivity = p, rankFullSensitivity = rank)
  probs <- full_join(probsFull, probsFullSensitivity, by = c("species", "ToolUse"))
  # rank order plot
  pA <-
    ggplot(probs, aes(x = rankFull, y = rankFullSensitivity, colour = ToolUse)) +
      geom_point(size = 2.5, alpha = 0.4) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      scale_x_continuous(
        name = "Rank from main model\nwith Normal(0, 1) intercept prior",
        breaks = c(1, 50, 100, 150)
        ) +
      scale_y_continuous(
        name = "Rank from alternative model\nwith Normal(1.78507, 2) intercept prior",
        breaks = c(1, 50, 100, 150)
        ) +
      scale_colour_manual(values = c("#D55E00", "#0072B2", "#CC79A7", "grey80")) +
      theme_classic() +
      theme(
        legend.title = element_blank(),
        legend.position = c(0.3, 0.8),
        axis.title = element_text(size = 9.5)
      )
  # probability plot
  pB <-
    ggplot(probs, aes(x = pFull, y = pFullSensitivity, colour = ToolUse)) +
    geom_point(size = 2.5, alpha = 0.4) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    scale_x_continuous(name = "Probability of tool use from main model\nwith Normal(0, 1) intercept prior") +
    scale_y_continuous(name = "Probability of tool use from alternative model\nwith Normal(1.78507, 2) intercept prior") +
    scale_colour_manual(values = c("#D55E00", "#0072B2", "#CC79A7", "grey80")) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 9.5)
    )
  # put together
  out <- 
    plot_grid(pA, pB, nrow = 2, labels = c("a","b")) + 
    theme(plot.margin = unit(c(0,0.1,0,0), "cm"))
  ggsave(out, filename = "figures/survCureSensitivity.pdf", height = 7.5, width = 4.5)
  return(out)
}

# plot ROC curve
plotROCCurve <- function(d, postFull) {
  # 1. get probabilities
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
      ungroup() %>%
      # tool use variable
      mutate(ToolUse = ifelse(d$LiteratureExists_Binary == 0 & d$VideosFound == 0, 0, 1))
  # 2. iterate over thresholds
  roc <- 
    data.frame(
      threshold = seq(0, 1, by = 0.001),
      truePositiveRate = NA,
      falsePositiveRate = NA
    )
  for (i in 1:nrow(roc)) {
    # a. calculate confusion matrix
    confusionMatrix <-
      table(
        factor(as.numeric(probs$p > roc$threshold[i]), levels = c(1, 0)),
        factor(probs$ToolUse, levels = c(1, 0))
      )
    # b. save point for ROC curve
    # true positive rate = true positives / (true positives + false negatives)
    roc$truePositiveRate[i] <- confusionMatrix[1,1] / (confusionMatrix[1,1] + confusionMatrix[2,1])
    # false positive rate = false positives / (false positives + true negatives)
    roc$falsePositiveRate[i] <- confusionMatrix[1,2] / (confusionMatrix[1,2] + confusionMatrix[2,2])
  }
  # 3. plot curve
  out <-
    ggplot() +
    geom_path(data = roc, aes(x = falsePositiveRate, y = truePositiveRate)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    labs(x = "False positive rate", y = "True positive rate") +
    theme_classic()
  # save
  ggsave(out, filename = "figures/survCureROC.pdf", width = 4, height = 4)
  return(out)
}

# calculate area under the curve
calculateAUC <- function(d, postFull) {
  # get species probabilities
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
    ungroup() %>%
    # tool use variable
    mutate(ToolUse = ifelse(d$LiteratureExists_Binary == 0 & d$VideosFound == 0, 0, 1))
  # calculate area under the curve
  roc <- 
    pROC::roc(
      response = probs$ToolUse,
      predictor = probs$p
    )
  auc <- as.numeric(roc$auc)
  return(auc)
}

# cross-validate survival cure model
crossValSurvCureModel <- function(d, phylogeny, iter, model, cv) {
  # scaled distance matrix from one phylogeny
  phylogeny <- phylogeny[[iter]]
  distMat <- cophenetic.phylo(phylogeny)
  distMat <- distMat / max(distMat)
  # set tool use to unobserved for cross validation
  eventLit <- d$LiteratureExists_Binary
  eventVid <- ifelse(d$VideosFound > 0, 1, 0)
  eventLit[cv] <- 0
  eventVid[cv] <- 0
  # data list for stan
  dat_list <- list(
    N        = nrow(d),
    eventLit = eventLit,
    eventVid = eventVid,
    numLit   = d$EX.Pre.Tool.Papers,
    numVid   = d$SearchHitsPreToolUse,
    Species  = 1:nrow(d),
    EQ       = as.numeric(scale(d$EQ)),
    Feeding  = ifelse(d$Feeding == "Generalist", 1, 2),
    distMat  = distMat
  )
  # fit model
  out <- sampling(
    model,
    data = dat_list,
    pars = c("aLit", "aVid", "aP", "bEQ", "k", "etasq", "rhosq"),
    iter = 4000,
    chains = 1,
    cores = 1,
    seed = 2113,
    control = list(adapt_delta = 0.85)
  )
  return(out)
}

# determine whether cross-validation identified tool user
getCrossVal <- function(d, crossValPost, cv) {
  # get probs from model
  probs <-
    tibble(
      species = rep(1:nrow(d), each = length(crossValPost$aP[,1])),
      aP1 = rep(crossValPost$aP[,1], times = nrow(d)),
      aP2 = rep(crossValPost$aP[,2], times = nrow(d)),
      bEQ = rep(crossValPost$bEQ, times = nrow(d)),
      k = as.vector(crossValPost$k),
      Feeding = rep(d$Feeding, each = length(crossValPost$aP[,1])),
      EQ = rep(as.numeric(scale(d$EQ)), each = length(crossValPost$aP[,1]))
    ) %>%
    mutate(
      aP = ifelse(Feeding == "Generalist", aP1, aP2),
      p = inv_logit(-(aP + k + bEQ*EQ))
    ) %>%
    # summarise posterior
    group_by(species) %>%
    summarise(p = median(p)) %>%
    ungroup() %>%
    # add tool use
    mutate(ToolUse = d$ToolUse)
  # is the target species within all other tool users?
  target <- probs$p[cv]
  threshold <-
    probs %>%
    slice(-cv) %>%
    filter(ToolUse == 1) %>%
    pull(p) %>%
    min()
  return(target >= threshold)
}
