# custom functions - plot phylogeny

# plot maximum clade credibility phylogeny with tool use data
plotPhylogeny1 <- function(d, MCCphylogeny) {
  # top ten candidate tool users from survival cure model
  candidates <- c(
    "Cacatua_ophthalmica","Poicephalus_meyeri","Guaruba_guarouba",
    "Cacatua_tenuirostris","Cacatua_ducorpsii","Poicephalus_gulielmi",
    "Poicephalus_robustus","Ognorhynchus_icterotis","Cacatua_haematuropygia",
    "Callocephalon_fimbriatum"
  )
  # mutate and rename vars
  d <-
    d %>%
    mutate(
      # total number of papers/videos
      ResearchEffort_Lit = log(Current.General.Papers + 1),
      ResearchEffort_Vid = log(SearchHits + 1),
      ToolUse_Lit = ifelse(LiteratureExists_Binary == 1, "Present in literature", "Absent"),
      ToolUse_Vid = ifelse(VideosFound > 0, "Present in videos", "Absent"),
      # add "candidate" tool users from survival cure model
      Candidates = ifelse(Name %in% candidates, "Top 10 candidate species", "Non-candidate")
      )
  # abbreviated species names for labels
  makeLabel <- function(x) paste0(substring(x, 1, 1), ". ", str_split(x, "_")[[1]][2])
  for (i in 1:length(d$Name)) {
    d$Name[i] <- makeLabel(d$Name[i])
    MCCphylogeny$tip.label[i] <- makeLabel(MCCphylogeny$tip.label[i])
  }
  # plot
  out <- 
    ggtree(MCCphylogeny, layout = "fan") %<+% d +
    geom_tiplab2(aes(colour = Candidates), size = 2.6, offset = 5) +
    geom_hilight(node = 179, fill = "#EDEDED") +
    geom_fruit(geom = geom_point, aes(colour = ToolUse_Vid, size = ResearchEffort_Vid)) +
    geom_fruit(geom = geom_point, aes(colour = ToolUse_Lit, size = ResearchEffort_Lit), offset = 0.04) +
    geom_cladelab(node = 328, offset = 15.5, fontsize = 7, label = "Cacatuidae") +
    geom_cladelab(node = 346, offset = 15.5, fontsize = 7, label = "Strigopidae") +
    geom_cladelab(node = 199, offset = 15.5, fontsize = 7, label = "\nAgapornithinae", hjust = 0.8) +
    geom_cladelab(node = 300, offset = 15.5, fontsize = 7, label = "Androglossini\n", hjust = 0.8) +
    geom_cladelab(node = 291, offset = 15.5, fontsize = 7, label = "Forpini", hjust = 1.1, vjust = -0.2) +
    geom_cladelab(node = 238, offset = 15.5, fontsize = 7, label = "Loriinae", hjust = -0.05) +
    geom_cladelab(node = 224, offset = 15.5, fontsize = 7, label = "Platycercinae", hjust = -0.1) +
    geom_cladelab(node = 323, offset = 15.5, fontsize = 7, label = "Psittacinae", hjust = -0.1) +
    geom_cladelab(node = 180, offset = 15.5, fontsize = 7, label = "Psittaculinae", hjust = 1) +
    geom_cladelab(node = 260, offset = 15.5, fontsize = 7, label = "Psittrichasiidae", hjust = 1) +
    geom_strip("A. chloropterus", "P. leucogaster", offset = 15.5, fontsize = 7, label = "Arini", hjust = 1.2) +
    scale_colour_manual(
      name = "\nTool use",
      values = c(
        "Present in literature" = "#0072B2",
        "Present in videos" = "#D55E00",
        "Top 10 candidate species" = "#F0E442",
        "Absent" = "#E0E0E0"
        )
      ) +
    scale_size(
      name = "Research effort\n(log scale)",
      labels = function(x) exp(x) - 1,
      breaks = log(c(0, 10, 100, 1000) + 1)
    ) +
    guides(colour = guide_legend(override.aes = list(size = 6), order = 1)) +
    xlim(NA, 75) +
    theme(legend.position = c(0.86, 0.055),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.box = "horizontal")
  # get clade highlight behind
  out$layers <- out$layers[c(5,1:4,6:length(out$layers))]
  # add parrot images
  out <- 
    ggdraw(out) +
    draw_image("images/parrots/Androglossini_SMALL.png"   , x = -0.01, y =  0.45, scale = 0.10) +
    draw_image("images/parrots/Arini_SMALL.png"           , x = -0.36, y =  0.20, scale = 0.14) +
    draw_image("images/parrots/Cacatuidae_SMALL.png"      , x =  0.38, y =  0.24, scale = 0.10) +
    draw_image("images/parrots/Loriinae_SMALL.png"        , x =  0.39, y = -0.24, scale = 0.10) +
    #draw_image("images/parrots/Platycercinae_SMALL.png"   , x =  0.12, y = -0.42, scale = 0.10) +
    draw_image("images/parrots/Psittacinae_SMALL.png"     , x =  0.27, y =  0.38, scale = 0.11) +
    draw_image("images/parrots/Psittaculinae_SMALL.png"   , x = -0.32, y = -0.28, scale = 0.13) +
    draw_image("images/parrots/Psittrichasiidae_SMALL.png", x = -0.47, y = -0.04, scale = 0.10) +
    draw_text("Psittaculidae", hjust = -1.5, vjust = 1.8, size = 20, colour = "#808080")
  # save
  ggsave(out, filename = "figures/phylogeny1.pdf", width = 17, height = 15)
  return(out)
}

# plot maximum clade credibility phylogeny with eq and feeding data
plotPhylogeny2 <- function(d, MCCphylogeny) {
  # abbreviated species names for labels
  makeLabel <- function(x) paste0(substring(x, 1, 1), ". ", str_split(x, "_")[[1]][2])
  for (i in 1:length(d$Name)) {
    d$Name[i] <- makeLabel(d$Name[i])
    MCCphylogeny$tip.label[i] <- makeLabel(MCCphylogeny$tip.label[i])
  }
  # plot
  out <- 
    ggtree(MCCphylogeny, layout = "fan") %<+% (d %>% rename(`Feeding Type` = Feeding)) +
    geom_tiplab2(size = 2.6, offset = 1) +
    geom_tippoint(aes(size = `EQ`, colour = `Feeding Type`)) +
    geom_cladelab(node = 328, offset = 11, fontsize = 7, label = "Cacatuidae") +
    geom_cladelab(node = 346, offset = 11, fontsize = 7, label = "Strigopidae") +
    geom_cladelab(node = 199, offset = 11, fontsize = 7, label = "\nAgapornithinae", hjust = 0.8) +
    geom_cladelab(node = 300, offset = 11, fontsize = 7, label = "Androglossini\n", hjust = 0.8) +
    geom_cladelab(node = 291, offset = 11, fontsize = 7, label = "Forpini", hjust = 1.1, vjust = -0.2) +
    geom_cladelab(node = 238, offset = 11, fontsize = 7, label = "Loriinae", hjust = -0.05) +
    geom_cladelab(node = 224, offset = 11, fontsize = 7, label = "Platycercinae", hjust = -0.1) +
    geom_cladelab(node = 323, offset = 11, fontsize = 7, label = "Psittacinae", hjust = -0.1) +
    geom_cladelab(node = 180, offset = 11, fontsize = 7, label = "Psittaculinae", hjust = 1) +
    geom_cladelab(node = 260, offset = 11, fontsize = 7, label = "Psittrichasiidae", hjust = 1) +
    geom_strip("A. chloropterus", "P. leucogaster", offset = 11, fontsize = 7, label = "Arini", hjust = 1.1) +
    scale_size(range = c(0, 6.5), breaks = c(1.0, 1.4, 1.8)) +
    scale_colour_manual(values = c("#D55E00", "#0072B2")) +
    guides(colour = guide_legend(override.aes = list(size = 6))) +
    xlim(NA, 70) +
    theme(legend.position = c(0.88, 0.05),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.box = "horizontal")
  # save
  ggsave(out, filename = "figures/phylogeny2.pdf", width = 17, height = 15)
  return(out)
}

# plot maximum clade credibility phylogeny with published papers
plotPhylogeny3 <- function(d, MCCphylogeny) {
  # abbreviated species names for labels
  makeLabel <- function(x) paste0(substring(x, 1, 1), ". ", str_split(x, "_")[[1]][2])
  for (i in 1:length(d$Name)) {
    d$Name[i] <- makeLabel(d$Name[i])
    MCCphylogeny$tip.label[i] <- makeLabel(MCCphylogeny$tip.label[i])
  }
  # plot
  out <- 
    ggtree(MCCphylogeny, layout = "fan") %<+% rename(d, `Pre-tool papers` = EX.Pre.Tool.Papers) +
    geom_tiplab2(size = 2.6, offset = 1) +
    geom_tippoint(aes(size = `Pre-tool papers`)) +
    geom_cladelab(node = 328, offset = 11, fontsize = 7, label = "Cacatuidae") +
    geom_cladelab(node = 346, offset = 11, fontsize = 7, label = "Strigopidae") +
    geom_cladelab(node = 199, offset = 11, fontsize = 7, label = "\nAgapornithinae", hjust = 0.8) +
    geom_cladelab(node = 300, offset = 11, fontsize = 7, label = "Androglossini\n", hjust = 0.8) +
    geom_cladelab(node = 291, offset = 11, fontsize = 7, label = "Forpini", hjust = 1.1, vjust = -0.2) +
    geom_cladelab(node = 238, offset = 11, fontsize = 7, label = "Loriinae", hjust = -0.05) +
    geom_cladelab(node = 224, offset = 11, fontsize = 7, label = "Platycercinae", hjust = -0.1) +
    geom_cladelab(node = 323, offset = 11, fontsize = 7, label = "Psittacinae", hjust = -0.1) +
    geom_cladelab(node = 180, offset = 11, fontsize = 7, label = "Psittaculinae", hjust = 1) +
    geom_cladelab(node = 260, offset = 11, fontsize = 7, label = "Psittrichasiidae", hjust = 1) +
    geom_strip("A. chloropterus", "P. leucogaster", offset = 11, fontsize = 7, label = "Arini", hjust = 1.1) +
    scale_size(range = c(0, 6.5), trans = "log1p", breaks = c(1, 10, 100)) +
    xlim(NA, 70) +
    theme(legend.position = c(0.88, 0.05),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.box = "horizontal")
  # save
  ggsave(out, filename = "figures/phylogeny3.pdf", width = 17, height = 15)
  return(out)
}

# plot maximum clade credibility phylogeny with search hits
plotPhylogeny4 <- function(d, MCCphylogeny) {
  # abbreviated species names for labels
  makeLabel <- function(x) paste0(substring(x, 1, 1), ". ", str_split(x, "_")[[1]][2])
  for (i in 1:length(d$Name)) {
    d$Name[i] <- makeLabel(d$Name[i])
    MCCphylogeny$tip.label[i] <- makeLabel(MCCphylogeny$tip.label[i])
  }
  # plot
  out <- 
    ggtree(MCCphylogeny, layout = "fan") %<+% rename(d, `Pre-tool search hits` = SearchHitsPreToolUse) +
    geom_tiplab2(size = 2.6, offset = 1) +
    geom_tippoint(aes(size = `Pre-tool search hits`)) +
    geom_cladelab(node = 328, offset = 11, fontsize = 7, label = "Cacatuidae") +
    geom_cladelab(node = 346, offset = 11, fontsize = 7, label = "Strigopidae") +
    geom_cladelab(node = 199, offset = 11, fontsize = 7, label = "\nAgapornithinae", hjust = 0.8) +
    geom_cladelab(node = 300, offset = 11, fontsize = 7, label = "Androglossini\n", hjust = 0.8) +
    geom_cladelab(node = 291, offset = 11, fontsize = 7, label = "Forpini", hjust = 1.1, vjust = -0.2) +
    geom_cladelab(node = 238, offset = 11, fontsize = 7, label = "Loriinae", hjust = -0.05) +
    geom_cladelab(node = 224, offset = 11, fontsize = 7, label = "Platycercinae", hjust = -0.1) +
    geom_cladelab(node = 323, offset = 11, fontsize = 7, label = "Psittacinae", hjust = -0.1) +
    geom_cladelab(node = 180, offset = 11, fontsize = 7, label = "Psittaculinae", hjust = 1) +
    geom_cladelab(node = 260, offset = 11, fontsize = 7, label = "Psittrichasiidae", hjust = 1) +
    geom_strip("A. chloropterus", "P. leucogaster", offset = 11, fontsize = 7, label = "Arini", hjust = 1.1) +
    scale_size(range = c(0, 6.5), trans = "log", breaks = c(100, 1000, 10000, 100000), labels = scales::comma) +
    xlim(NA, 70) +
    theme(legend.position = c(0.88, 0.07),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.box = "horizontal")
  # save
  ggsave(out, filename = "figures/phylogeny4.pdf", width = 17, height = 15)
  return(out)
}
