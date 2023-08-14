# custom functions - load data

# convert search hits variable to "search hits before tool use discovered"
convertSearchHits <- function(dateFirstVideoPublished, searchHits) {
  # estimate how many search hits there were when the first video of tool use was published
  # fit model, assuming 1 search hit at youtube inception and exponential growth
  dateSearchHitsRecorded <- as.Date("2020-03-01")
  dateYouTubeInception <- as.Date("2005-02-14")
  d <- data.frame(
    daysSinceYouTubeInception = c(0, as.numeric(dateSearchHitsRecorded - dateYouTubeInception)),
    numSearchHits = c(1, searchHits)
  )
  m <- lm(numSearchHits ~ daysSinceYouTubeInception, data = d)
  # make prediction for date video published
  newdata <- data.frame(daysSinceYouTubeInception = as.numeric(dateFirstVideoPublished - dateYouTubeInception))
  pred <- predict(m, newdata = newdata)[[1]]
  out <- round(pred, 0)
  return(out)
}

# load main data
loadData <- function(data_file, litCount_file, crowdData, phylogeny) {
  # load lit count data for matching
  litCount <- read.csv(litCount_file)
  # edit crowdsourcing data
  crowdData <- 
    crowdData %>%
    # mutate species name for matching
    mutate(Species = str_replace(crowdData$Species, " ", "_")) %>%
    # for each species, keep earliest video of true tool use
    group_by(Species) %>%
    filter(ToolUse == "Y") %>%
    filter(PublishingDate == min(PublishingDate)) %>%
    dplyr::select(Species, PublishingDate)
  # load main data
  out <- 
    read_csv(data_file, show_col_types = FALSE) %>% 
    mutate(
      # match names with phylogeny
      Name = str_replace(Name, " ", "_"),
      # tool use
      ToolUse = ifelse(VideosFound == 0 & LiteratureExists_Binary == 0, 0, 1)
      ) %>% 
    # keep only species in phylogeny
    filter(Name %in% phylogeny[[1]]$tip.label) %>%
    # join lit count data
    left_join(litCount, by = "Name") %>%
    # join publishing date
    left_join(crowdData, by = c("Name" = "Species")) %>%
    # estimate search hits before tool use
    mutate(
      SearchHitsPreToolUse = ifelse(
        VideosFound > 0,
        # if tool use found on YouTube, estimate search hits when tool use found
        unlist(map2(as.Date(PublishingDate), SearchHits, convertSearchHits)),
        # else, use total search hits
        SearchHits
      )
    )
  # match dataset order to phylogeny tip labels
  out <- out[match(phylogeny[[1]]$tip.label, out$Name),]
  return(out)
}

# get indexes to iterate over
getIter <- function(n) {
  set.seed(2113)
  sample(1:1000, n)
}
