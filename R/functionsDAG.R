# custom functions - dag

# draw causal model
drawDAG <- function() {
  # coordinates for plot
  dag_coords <- tibble(
    name = c("TUu", "TUoLit", "TUoVid", "numLit", "numVid", "EQ", "U", "Phy", "F"),
    x    = c(0, -1, 1, -2, 2, -1, 0, 1.5, 1),
    y    = c(1, 2, 2, 1, 1, 0.25, -0.5, -0.5, 0.25)
  )
  out <-
    dagify(TUoVid ~ TUu + numVid,
           TUoLit ~ TUu + numLit,
           TUu ~ EQ + F + U,
           EQ ~ U,
           F ~ U,
           U ~ Phy,
           latent = c("TUu", "U"),
           coords = dag_coords
           ) %>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend))
  # modify labels
  labels <- c("TUu"    = "Tool use\n(unobserved)",
              "TUoLit" = "Tool use\nobserved in\nscientific literature",
              "TUoVid" = "Tool use\nobserved on\ndigital video platform",
              "numLit" = "Number of\nscientific papers\npublished",
              "numVid" = "Number of\nvideos on\ndigital platform",
              "EQ"     = "Relative\nbrain size",
              "F"      = "Feeding\nstrategy",
              "U"      = "Unobserved\nconfounds",
              "Phy"    = "Phylogeny")
  out$data$name <- labels[out$data$name]
  out$data$to <- labels[out$data$to]
  # plot
  out <-
    out +
    geom_dag_point(
      data = function(x) dplyr::filter(x, name %in% c("Tool use\n(unobserved)", "Unobserved\nconfounds")),
      alpha = 0.5, size = 35, show.legend = FALSE, colour = "lightgrey"
    ) +
    geom_dag_text(colour = "black") +
    geom_dag_edges(
      start_cap = ggraph::circle(radius = 1.35, unit = "cm"),
      end_cap = ggraph::circle(radius = 1.35, unit = "cm")
    ) +
    ylim(c(-1, 2.3)) +
    xlim(c(-2.5, 2.5)) +
    theme_void()
  # save
  ggsave(out, filename = "figures/dag.pdf", height = 4, width = 6)
  return(out)
}
