# logistic regression with phylogenetic control
fitComparisonModel <- function(d, phylogeny, iter, var) {
  # scaled distance matrix from one phylogeny
  phylogeny <- phylogeny[[iter]]
  distMat <- cophenetic.phylo(phylogeny)
  distMat <- distMat / max(distMat)
  # reorder data
  d <- d[match(phylogeny$tip.label, d$Name),]
  # data list for stan
  dat_list <- list(
    N = nrow(d),
    Species = 1:nrow(d),
    y = pull(d, var), 
    Feeding = ifelse(d$Feeding == "Generalist", 1, 2),
    EQ = as.numeric(scale(d$EQ)),
    distMat = distMat
  )
  # fit model in rethinking
  out <- ulam(
    alist(
      y ~ bernoulli(p),
      logit(p) <- aP[Feeding] + bEQ*EQ + k[Species],
      vector[N]:k ~ multi_normal(0, SIGMA),
      matrix[N,N]:SIGMA <- cov_GPL1(distMat, etasq, rhosq, 0.01),
      aP[Feeding] ~ normal(0, 1),
      bEQ ~ normal(0, 1),
      c(rhosq, etasq) ~ exponential(0.5)
    ),
    data = dat_list,
    iter = 4000, chains = 1, cores = 1, seed = 2113, 
    control = list(adapt_delta = 0.85)
  )
  return(out)
}

# plot comparing models
plotComparisons <- function(postM1, postM2, postFull) {
  out <-
    tibble(
      Model = rep(c("Pre-survey", "Post-survey", "Survival cure"), each = length(postFull$bEQ)),
      `Effect of encephalisation quotient` = c(postM1$bEQ, postM2$bEQ, -postFull$bEQ),
      `Difference between feeding strategies` = c(
        postM1$aP[,1] - postM1$aP[,2],
        postM2$aP[,1] - postM2$aP[,2],
        postFull$aP[,2] - postFull$aP[,1]
      )
    ) %>%
    pivot_longer(cols = !Model, names_to = "parameter", values_to = "Log-odds posterior value") %>%
    mutate(Model = factor(Model, levels = c("Survival cure", "Post-survey", "Pre-survey"))) %>%
    ggplot(aes(x = `Log-odds posterior value`, y = Model)) +
    stat_halfeye(.width = c(0.50, 0.95)) +
    facet_grid(. ~ parameter) +
    xlim(c(-2, 4))
  # save
  ggsave(out, filename = "figures/comparison.pdf", width = 8, height = 4)
  return(out)
}
