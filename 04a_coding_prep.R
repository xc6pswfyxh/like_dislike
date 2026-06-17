#### MANUAL CODING

library(tidyverse)

## def reply chain function
get_reply_chain <- function(df, topic_name, post_num) {
  chain <- c()
  current <- post_num
  visited <- c()
  
  while (!is.na(current) && !current %in% visited) {
    visited <- c(visited, current)
    row <- df |> filter(topic == topic_name, post_number == current)
    if (nrow(row) == 0) break
    chain <- c(chain, row$raw[1])
    current <- row$reply_to_post_number[1]
  }
  
  paste(rev(chain), collapse = " >>> ")
}


## data pp
comments1 <- read_csv("data/comments/data_posts_VM1.csv") |> 
  filter(!str_detect(topic, 
    regex("netiquette | studie", ignore_case = TRUE)))

comments1 <- comments1 |>
  rowwise() |>
  mutate(reply_chain = if_else(
    !is.na(reply_to_post_number),
    get_reply_chain(comments1, topic, reply_to_post_number),
    NA_character_
  )) |>
  ungroup() |> 
  filter(created_at > as.Date("2017-08-01")) 


comments2 <- read_csv("data/comments/data_posts_VM2.csv") |> 
  mutate(created_at = as_datetime(created_at)) |> 
  filter(!str_detect(topic, 
    regex("netiquette | studie", ignore_case = TRUE)))

comments2 <- comments2 |>
  rowwise() |>
  mutate(reply_chain = if_else(
    !is.na(reply_to_post_number),
    get_reply_chain(comments2, topic, reply_to_post_number),
    NA_character_
  )) |>
  ungroup() |> 
  filter(created_at > as.Date("2017-08-01")) 


comments3 <- read_csv("data/comments/data_posts_VM3.csv") |> 
  mutate(created_at = as_datetime(created_at)) |> 
  filter(!str_detect(topic, 
    regex("netiquette | studie", ignore_case = TRUE)))

comments3 <- comments3 |>
  rowwise() |>
  mutate(reply_chain = if_else(
    !is.na(reply_to_post_number),
    get_reply_chain(comments3, topic, reply_to_post_number),
    NA_character_
  )) |>
  ungroup() |> 
  filter(created_at > as.Date("2017-08-01")) 


## joining all comments
comments <- comments1 |> mutate(vm = "vm1") |> 
  bind_rows(comments2 |> mutate(vm = "vm2")) |>
  bind_rows(comments3 |> mutate(vm = "vm3"))

## prepare coding file 
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


# weighted sample
set.seed(666)
comments_reli <- comments |> 
  group_by(vm) |> 
  slice_sample(n = 50) |> # 50 per group (vm)
  ungroup() |> 
  select(-c(vm))


write_csv(comments, "data/comments/comments.csv", na = "")
write_csv(comments_reli, "data/comments/comments_reli.csv", na = "") # write csv for icr test


####-------------------------------------------------------------------------------------------------


## 2. ICR
coding_j <- readxl::read_excel("data/comments/comments_relij_coded.xlsx")
coding_y <- readxl::read_excel("data/comments/comments_reliy_coded.xlsx")

coding_j <- coding_j |> head(5)
coding_y <- coding_y |> head(5)


coding_comb <- coding_j |> 
  rbind(coding_y)

icr <- coding_comb |> 
  tidycomm::test_icr(unit_var = post_number, 
                     coder_var = coder, 
                     c(pers_exp:contr), 
                     kripp_alpha = TRUE,
                     cohens_kappa = TRUE,
                     lotus = TRUE, 
                     s_lotus = TRUE, 
                     na.omit = TRUE)

icr
writexl::write_xlsx(icr, "data/icr.xlsx") # write excel for icr table


####-------------------------------------------------------------------------------------------------


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