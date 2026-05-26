#### MANUAL CODING
## 1. LABELS


####----

library(tidyverse)


comments <- read_csv("data/comments/blahblahblah.csv")

  
## ead comments data for classification




##----

set.seed(666)

sample_labels <- comments |> # draw a sample of 200 obs for manual labelling from full dataset
  mutate(
    depth = NA, 
    feeling = NA,
    breadth = NA,
    pol_op = NA,
    valence = NA,
    contr = NA) 

sample_labels_reli <- sample_labels |> 
  sample_n(n = 200) |> # both code 200 obs, then we calculate reliability
  mutate(coder = NA) |> 
  relocate(coder, .before = post_id) # post_id ändern so dass coder in erster spalte

writexl::write_xlsx(sample_labels, "data/sample_labels.xlsx") # write excel for full labelling
writexl::write_xlsx(sample_labels_reli, "data/sample_labels_reli.xlsx") # write excel for icr test


## 2. ICR
labels_j <- readxl::read_excel("data/sample_labels_relij_coded.xlsx")
labels_l <- readxl::read_excel("data/sample_labels_reliy_coded.xlsx")

labels_jl <- labels_j |> 
  rbind(labels_l)

icr <- labels_jl |> 
  tidycomm::test_icr(unit_var = post_id, 
                     coder_var = coder, 
                     labels, 
                     kripp_alpha = TRUE,
                     cohens_kappa = TRUE,
                     lotus = TRUE, 
                     s_lotus = TRUE, 
                     na.omit = TRUE)

icr
writexl::write_xlsx(icr, "data/icr.xlsx") # write excel for icr table

rm(list = setdiff(ls(), c("bluesky_pp"))) # clean env


## 2. FINAL LABELS
bluesky_samplelabelled <- readxl::read_excel("data/sample_labels_coded.xlsx")
labelled_observations <- bluesky_samplelabelled

bluesky_samplelabelled |> # check if everything is coded
  count(labels)

training_data <- bluesky_samplelabelled |> 
  filter(!is.na(labels)) # remove NAs

training_data |> 
  count(labels)

writexl::write_xlsx(training_data, "data/training_data.xlsx") # use this to train classifier
rm(list = setdiff(ls(), c("bluesky_pp", "labelled_observations"))) # clean env


## 3. JOIN MANUALLY LABELLED DATA WITH REST OF DATASET

bluesky_pp <- bluesky_pp |> 
  left_join(labelled_observations |>  
              select(post_id, labels),
    by = "post_id"
  ) |> 
  relocate(labels, .after = "text")

write.csv(bluesky_pp, "data/bluesky_pp.csv", row.names = FALSE)


## END