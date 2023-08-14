# custom functions - phylogenetic signal

# get phylogenetic signal
fitPhySignal <- function(d, phylogeny, iter, dv) {
  # get data
  var <- factor(pull(d[,dv]))
  names(var) <- d$Name
  # fit model
  model <- fitDiscrete(phylogeny[[iter]], var, transform = "lambda")
  # extract lambda
  out <- model$opt$lambda
  return(out)
}
