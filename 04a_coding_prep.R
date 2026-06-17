#### MANUAL CODING

library(tidyverse)

comments <- readxl::read_xlsx("data/comments/data_posts_VM1.xlsx")

set.seed(666)

comments <- comments |>
  mutate(
    coder = NA,
    pers_exp = NA, 
    emot_exp = NA,
    pol_opin = NA,
    breadth = NA,
    valence = NA,
    contr = NA) |> 
  relocate(c(pers_exp:contr), .after = raw) |> 
  relocate(coder, .before = post_number)

comments_reli <- comments |> 
  slice_sample(n = 100) # both code 100 obs, then we calculate reliability


writexl::write_xlsx(comments, "data/comments.xlsx")
writexl::write_xlsx(comments_reli, "data/comments_reli.xlsx") # write excel for icr test

## 2. ICR
labels_j <- readxl::read_excel("data/comments_relij_coded.xlsx")
labels_l <- readxl::read_excel("data/comments_reliy_coded.xlsx")

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