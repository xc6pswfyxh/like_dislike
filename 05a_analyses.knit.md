---
title: "Analyses"
bibliography: 
  - bibliography/references.bib

format:
  html:
    toc: true
    toc-depth: 2
    code-fold: true
    code-summary: "Show code"
---





In what follows, you can find the analyses for the results reported in the paper. To see the underlying code, click on the button "code".

# Setup
## Packages

Let's first load the packages we require for our analyses. Note that some need to be installed manually from github.



::: {.cell}

```{.r .cell-code}
# install github packages
# devtools::install_github("tdienlin/td@v0.0.2.6")  # uncomment to install
# devtools::install_github("yrosseel/lavaan") 
# remotes::install_github("stan-dev/cmdstanr")

# define packages
packages <- 
  c(
    "brms", 
    "cowplot", 
    "devtools", 
    "english", 
    "faoutlier", 
    "GGally",
    "kableExtra", 
    "knitr", 
    "lavaan", 
    "magrittr",
    "marginaleffects",
    "MVN",
    "pwr", 
    "quanteda", 
    "semTools", 
    "tidyverse",
    "cmdstanr",
    "td"
    )

# load packages
lapply(
  packages
  # , install.packages
  , library
  , character.only = TRUE
  , quietly = TRUE
  )

# If not wanting to compute all transformations etc., simply load workspace
# load("data/workspace.RData")
```
:::



## Custom Functions

We now build some custom functions required for the analyses. Note that here you can find the syntax for the main analyses wrapped in helper functions, which we will then call below.



::: {.cell}

```{.r .cell-code}
# define function to silence result output
# needed for brms
hush <- 
  function(
    code
    ){
    sink("/dev/null") # use /dev/null in UNIX
    tmp = code
    sink()
    return(tmp)
    }

# plot distributions of items
plot_dist <-
  function(
    item,
    n_rows = 1
    ) {
    ggplot(
      gather(
        select(
          d, 
          starts_with(
            item
            )
          )
        ), 
      mapping = aes(
        x = value
        )
      ) +
      geom_bar() +
      facet_wrap(
        ~key, 
        nrow = n_rows
        ) + 
      theme_bw()
  }

get_desc <- 
  function(
    item,
    name
  ) {
    desc <- 
      select(
        d, 
        starts_with(
          item
          )
        ) %>% 
      apply(
        1, 
        mean, 
        na.rm = TRUE
        ) %>% 
      as.data.frame() %>% 
      summarise(
        m = mean(., na.rm = TRUE), 
        sd = sd(., na.rm = TRUE)
        )
    assign(
      paste0(
        "des_",
        name
      ),
      desc,
      envir = .GlobalEnv
    )
  }

# Function to display diagnostics of fitted hurdle models
model_diagnostics <- 
  function(
    model
  ){
  plot_grid(
    pp_check(
      model,
      type = "dens_overlay",
      nsamples = 100
    ),
    pp_check(
      model,
      type = "loo_pit_qq",
      nsamples = 100
    ),
    pp_check(
      model,
      type = "stat",
      stat = "median",
      nsamples = 100
    ),
    pp_check(
      model,
      type = "stat",
      stat = "mean",
      nsamples = 100
    ),
    labels = c("Density overlay", "LOO-PIT QQ", "Predicted medians", "Predicted means"),
    ncol = 2,
    label_size = 8,
    hjust = 0,
    vjust = 0,
    label_x = 0,
    label_y = 0.93
    )
  }

# wrapper function that runs relevant analyses for hurdle models
run_hurdles <- 
  function(
    object,
    name,
    outcome,
    predictor,
    x_lab,
    y_lab,
    plot = FALSE,
    ...
    ){
    
    # define define general model
    formula_mdl <- 
      formula(
        paste0(
          outcome,
          " ~ 1 + age + male + edu +",
          predictor
        )
    )
    
    # define hurdle model
    formula_hrdl <- 
      formula(
        paste0(
          "hu ~ 1 + age + male + edu +",
          predictor
        )
    )
   
    # fit model
    fit <- 
      hush(
        brms::brm(
          bf(
            formula_mdl, 
            formula_hrdl
            ), 
          data = object, 
          family = hurdle_gamma(),
          silent = TRUE,
          refresh = 0,
          iter = 1000,
          backend = "rstan"
          )
      )
    
    # export fit
    assign(
      paste0(
        "fit_",
        name),
      fit,
      envir = .GlobalEnv
    )
    
    # plot model
    if(
      isTRUE(plot
             )) {
      plot(fit, ask = FALSE)
  }
    
    # run diagnostics
    fit_diagnostics <- 
      model_diagnostics(fit)
    plot(fit_diagnostics)
    
    # print summary
    fit_summary <- 
      summary(fit)
    print(fit_summary)

    # calculate slopes
    slopes <-
      cbind(
        Outcome = outcome,
        avg_slopes(fit)
      )
    
    # print slopes
    cat("\nResults of marginal effects:\n")
    print(slopes)
    
    # export as object to environment
    assign(
      paste0(
        "slopes_",
        name
        ),
      slopes,
      envir = .GlobalEnv
    )
    
    # make figure
    fig_res <- 
      conditional_effects(
        fit
        , ask = FALSE
      )
    
    # make figure
    fig <- 
      plot(
        fig_res,
        plot = FALSE
      )[[predictor]] +
      labs(
        x = x_lab,
        y = y_lab
        )
    
    # plot figure
    print(fig)
    
    # export as object to environment
    assign(
      paste0(
        "fig_",
        name
      ),
      fig,
      envir = .GlobalEnv
    )
}
```
:::



## Data-Wrangling

Let's next tidy up and prepare the data, so we can analyze them.



::: {.cell}

```{.r .cell-code}
# load data
# please note that in a prior step several different data sets were merged using information such as case_taken or IP-addresses. To guarantee the privacy of the participants, these data were deleted after merging.
d_raw <- read_csv("data/data_raw.csv")

# recode variables
d <- 
  d_raw %>% 
  dplyr::rename(
    "male" = "SD01", 
    "age" = "SD02_01", 
    "state" = "SD03", 
    "edu_fac" = "SD04",
    "TR02_01" = "TR01_01",
    "TR02_02" = "TR01_05",
    "TR02_03" = "TR01_09",
    "TR01_01" = "TR01_02",
    "TR01_02" = "TR01_03",
    "TR01_03" = "TR01_04",
    "TR01_04" = "TR01_06",
    "TR01_05" = "TR01_07",
    "TR01_06" = "TR01_08",
    "TR01_07" = "TR01_10",
    "TR01_08" = "TR01_11",
    "TR01_09" = "TR01_12"
    ) %>% 
  mutate(
    
    #general
    version = factor(
      VE01_01, 
      labels = c(
        "Control", 
        "Like", 
        "Like & dislike"
        )
      ),
    version_lkdslk = factor(
      version, 
      levels = c(
        "Like & dislike", 
        "Like", 
        "Control"
        )
      ),
    version_no = recode(
      version, 
      "Control" = 1, 
      "Like" = 2, 
      "Like & dislike" = 3
      ),
    male = replace_na(
      male, 
      "[3]") %>% 
      dplyr::recode(
        "männlich [1]" = 1, 
        "weiblich [2]" = 0
        ),
    age_grp = cut(
      age, 
      breaks = c(0, 29, 39, 49, 59, 100), 
      labels = c("<30", "30-39", "40-49", "50-59", ">59")
      ),
    
    # contrasts
    like = recode(
      version, 
      "Control" = 0, 
      "Like" = 1, 
      "Like & dislike" = 0
      ),
    likedislike = recode(
      version, 
      "Control" = 0, 
      "Like" = 0, 
      "Like & dislike" = 1
      ),
    control = recode(
      version, 
      "Control" = 1, 
      "Like" = 0, 
      "Like & dislike" = 0
      ),
    
    # recode
    edu = recode(
      edu_fac, 
      "Hauptschulabschluss/Volksschulabschluss" = 1, 
      "Noch Schüler" = 1,
      "Realschulabschluss (Mittlere Reife)" = 1,
      "Abschluss Polytechnische Oberschule 10. Klasse (vor 1965: 8. Klasse)" = 1,
      "Anderer Schulabschluss:" = 1,
      "Fachhochschulreife (Abschluss einer Fachoberschule)" = 2,
      "Abitur, allgemeine oder fachgebundene Hochschulreife (Gymnasium bzw. EOS)" = 2,
                 "Hochschulabschluss" = 3),
    PC01_03 = 8 - PC01_03, 
    SE01_05 = 8 - SE01_05, 
    SE01_06 = 8 - SE01_06,
    
    # behavioral data
    words = replace_na(
      words, 
      0
      ),
    self_dis = words + (reactions * 2),
    time_read_MIN = time_read / 60,
    time_sum_min_t2 = TIME_SUM_t2 / 60,
    
    # take logarithm
    words_log = log1p(words),
    self_dis_log = log1p(self_dis)
  )

# variable labels
var_names <- c(
  "Privacy concerns", 
  "Gratifications general", 
  "Gratifications specific",
  "Privacy deliberation", 
  "Self-efficacy", 
  "Trust general", 
  "Trust specific",
  "Words log"
  )

var_names_breaks <- c(
  "Privacy\nconcerns", 
  "Gratifications\ngeneral", 
  "Gratifications\nspecific",
  "Privacy\ndeliberation", 
  "Self-\nefficacy", 
  "Trust\ngeneral", 
  "Trust\nspecific",
  "Words (log)"
  )

# Extract sample characteristics and descriptives.
# sample descriptives t1
n_t1 <- nrow(d)
age_t1_m <- mean(d$age, na.rm = TRUE)
male_t1_m <- mean(d$male, na.rm = TRUE)
college_t1_m <- table(d$edu)[3] / n_t1

# sample descriptives t2
# participants
n_t2 <- 
  filter(
    d, 
    !is.na(
      id_survey_t2
      )
    ) %>% 
  nrow()

# age
age_t2_m <- 
  filter(
    d, 
    !is.na(
      id_survey_t2
      )
    ) %>% 
  select(
    age
    ) %>% 
  mean()

# gender
male_t2_m <- 
  filter(
    d, 
    !is.na(
      id_survey_t2
      )
    ) %>% 
  select(
    male
    ) %>% 
  mean(
    na.rm = TRUE
    )

# education
college_t2_m <- 
  filter(
    d, 
    !is.na(
      id_survey_t2
      )
    ) %>% 
  select(
    edu
    ) %>% 
  table() %>% 
  .[3] / n_t2

# descriptives users
n_users <- 
  filter(
    d, 
    !is.na(
      post_count
      )
    ) %>% 
  nrow()

# characteristics of posts
n_comments <- 
  sum(
    d$post_count, 
    na.rm = TRUE
    )

n_words <- 
  sum(
    d$words, 
    na.rm = TRUE
    )

n_time <- 
  sum(
    d$time_read_MIN, 
    na.rm = TRUE
    )

n_posts <- 
  sum(
    d$post_count, 
    na.rm = TRUE
    )

# filter unmatched cases, use only completes
d <- 
  filter(
    d, 
    !is.na(
      post_count
      ) & !is.na(
        id_survey_t2
        )
    ) 

n_matched <- nrow(d)

# save data file with all participants to compare results
d_all <- d
```
:::



# Filter Participants

I filtered participants who answered the questionnaire in less than three minutes, which I considered to be unreasonably fast.



::: {.cell}

```{.r .cell-code}
# filter speeders
time_crit <- 3 # minimum time on survey

# count number of speeders
n_speeding <- 
  nrow(
    filter(
      d, 
      time_sum_min_t2 < time_crit
      )
    )

# Deletion of fast respondents
d <- 
  filter(
    d, 
    time_sum_min_t2 >= time_crit
    )
```
:::



I inspected the data manually for cases with obvious response patterns. The following cases show extreme response patterns (alongside fast response times), and were hence removed.



::: {.cell}

```{.r .cell-code}
# identify response patterns
resp_pattern <- c(
  "ANIEVLK9F2SW", 
  "BN4MAOWZO7W2"  # clear response pattern
  )

n_resp_pattern <- 
  length(
    resp_pattern
    )

d %>% 
  filter(
    case_token %in% resp_pattern
    ) %>% 
  select(
    case_token, 
    GR01_01:SE01_06, 
    topics_entered:reactions, 
    -SO01_01, 
    TIME_SUM_t1, 
    TIME_SUM_t2
    ) %>% 
  kable() %>% 
  kable_styling("striped") %>% 
  scroll_box(width = "100%")
```

::: {.cell-output-display}
`````{=html}
<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; "><table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> case_token </th>
   <th style="text-align:right;"> GR01_01 </th>
   <th style="text-align:right;"> GR01_02 </th>
   <th style="text-align:right;"> GR01_03 </th>
   <th style="text-align:right;"> GR01_04 </th>
   <th style="text-align:right;"> GR01_05 </th>
   <th style="text-align:right;"> GR01_06 </th>
   <th style="text-align:right;"> GR01_07 </th>
   <th style="text-align:right;"> GR01_08 </th>
   <th style="text-align:right;"> GR01_09 </th>
   <th style="text-align:right;"> GR01_10 </th>
   <th style="text-align:right;"> GR01_11 </th>
   <th style="text-align:right;"> GR01_12 </th>
   <th style="text-align:right;"> GR01_13 </th>
   <th style="text-align:right;"> GR01_14 </th>
   <th style="text-align:right;"> GR01_15 </th>
   <th style="text-align:right;"> GR02_01 </th>
   <th style="text-align:right;"> GR02_02 </th>
   <th style="text-align:right;"> GR02_03 </th>
   <th style="text-align:right;"> GR02_04 </th>
   <th style="text-align:right;"> GR02_05 </th>
   <th style="text-align:right;"> PC01_01 </th>
   <th style="text-align:right;"> PC01_02 </th>
   <th style="text-align:right;"> PC01_03 </th>
   <th style="text-align:right;"> PC01_04 </th>
   <th style="text-align:right;"> PC01_05 </th>
   <th style="text-align:right;"> PC01_06 </th>
   <th style="text-align:right;"> PC01_07 </th>
   <th style="text-align:right;"> TR02_01 </th>
   <th style="text-align:right;"> TR01_01 </th>
   <th style="text-align:right;"> TR01_02 </th>
   <th style="text-align:right;"> TR01_03 </th>
   <th style="text-align:right;"> TR02_02 </th>
   <th style="text-align:right;"> TR01_04 </th>
   <th style="text-align:right;"> TR01_05 </th>
   <th style="text-align:right;"> TR01_06 </th>
   <th style="text-align:right;"> TR02_03 </th>
   <th style="text-align:right;"> TR01_07 </th>
   <th style="text-align:right;"> TR01_08 </th>
   <th style="text-align:right;"> TR01_09 </th>
   <th style="text-align:right;"> PD01_01 </th>
   <th style="text-align:right;"> PD01_02 </th>
   <th style="text-align:right;"> PD01_03 </th>
   <th style="text-align:right;"> PD01_04 </th>
   <th style="text-align:right;"> PD01_05 </th>
   <th style="text-align:right;"> SE01_01 </th>
   <th style="text-align:right;"> SE01_02 </th>
   <th style="text-align:right;"> SE01_03 </th>
   <th style="text-align:right;"> SE01_04 </th>
   <th style="text-align:right;"> SE01_05 </th>
   <th style="text-align:right;"> SE01_06 </th>
   <th style="text-align:right;"> topics_entered </th>
   <th style="text-align:right;"> posts_read_count </th>
   <th style="text-align:right;"> time_read </th>
   <th style="text-align:right;"> topic_count </th>
   <th style="text-align:right;"> post_count </th>
   <th style="text-align:right;"> words </th>
   <th style="text-align:right;"> reactions </th>
   <th style="text-align:right;"> TIME_SUM_t1 </th>
   <th style="text-align:right;"> TIME_SUM_t2 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> BN4MAOWZO7W2 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 27 </td>
   <td style="text-align:right;"> 211 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 37 </td>
   <td style="text-align:right;"> 246 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANIEVLK9F2SW </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 151 </td>
   <td style="text-align:right;"> 1518 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 73 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 66 </td>
   <td style="text-align:right;"> 469 </td>
  </tr>
</tbody>
</table></div>

`````
:::

```{.r .cell-code}
d %<>% 
  filter(
    !case_token %in% resp_pattern
    )

# sample descriptives final data set
n_final <- nrow(d)
age_final_m <- mean(d$age)
male_final_m <- mean(d$male, na.rm = TRUE)
college_final_m <- table(d$edu)[[3]] / n_final
```
:::



# Measures

Let's next inspect all measures. I list all items, their distributions and CFAs.

## Privacy concerns
### Items

Using the participation platform I had ...

 1. ... concerns about what happens to my data.
 2. ... concerns about disclosing information about myself.
 3. ... no concerns. (reversed)
 4. ... concerns that others could discover my real identity (i.e. my last and first name).
 5. ... concerns that information about myself could fall into wrong hands.
 6. ... concerns that others could discover what my political views are.
 7. ... concerns about my privacy.

### Distributions



::: {.cell}

```{.r .cell-code}
# plot distribution
plot_dist(
  item = "PC01"
)
```

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-6-1.png){width=768}
:::

```{.r .cell-code}
# extract results
get_desc(
  item = "PC01",
  name = "pricon"
)
```
:::



### CFA



::: {.cell}

```{.r .cell-code}
name <- "pricon"
model <- "
pri_con =~ PC01_01 + PC01_02 + PC01_03 + PC01_04 + PC01_05 + PC01_06 + PC01_07
"
fit <- 
  lavaan::sem(
    model = model, 
    data = d, 
    estimator = "MLR", 
    missing = "ML"
    )
```
:::



Model fit:



::: {.cell}

```{.r .cell-code}
facval <- 
  fit_tab(
    fit, 
    # reliability = TRUE, 
    # scaled = TRUE
    ) %T>%
  print()
```

::: {.cell-output .cell-output-stdout}

```
  chisq df  pvalue   cfi   tli  rmsea   srmr
1  31.9 14 0.00414 0.995 0.993 0.0478 0.0104
```


:::
:::



Factor loadings:



::: {.cell}

```{.r .cell-code}
inspect(
  fit,
  what = "std") %$% 
    lambda
```

::: {.cell-output .cell-output-stdout}

```
        pri_cn
PC01_01  0.931
PC01_02  0.901
PC01_03  0.547
PC01_04  0.890
PC01_05  0.908
PC01_06  0.795
PC01_07  0.928
```


:::
:::



Shows that PC01_03 doesn't load well. As it's an inverted item that's not surprising. Also from a theoretic perspective it's suboptimal, because it doesn't explicitly focus on privacy, but just concerns in general. Will be deleted.

### CFA 2



::: {.cell}

```{.r .cell-code}
model <- "
pri_con =~ PC01_01 + PC01_02 + PC01_04 + PC01_05 + PC01_06 + PC01_07
"
fit <-
  lavaan::sem(
    model = model,
    data = d,
    estimator = "MLR",
    missing = "ML"
    )
```
:::



Model fit:



::: {.cell}

```{.r .cell-code}
facval <- 
  fit_tab(
    fit, 
    # reliability = TRUE, 
    # scaled = TRUE
    ) %T>% 
  print()
```

::: {.cell-output .cell-output-stdout}

```
  chisq df  pvalue   cfi   tli  rmsea    srmr
1  25.8  9 0.00217 0.995 0.992 0.0578 0.00963
```


:::

```{.r .cell-code}
assign(
  paste0(
    name, 
    "_facval"
    ), 
  facval
  )
```
:::



Factor loadings:



::: {.cell}

```{.r .cell-code}
inspect(
  fit,
  what = "std"
  ) %$% 
  lambda
```

::: {.cell-output .cell-output-stdout}

```
        pri_cn
PC01_01  0.931
PC01_02  0.901
PC01_04  0.891
PC01_05  0.909
PC01_06  0.796
PC01_07  0.926
```


:::
:::



Updated version shows good fit.

## Gratifications general
### Items

Using the participation platform ...
 
 1.	... had many benefits for me.
 2.	... has paid off for me.
 3.	... was worthwhile.
 4.	... was fun.
 5.	... has brought me further regarding content.

### Distributions



::: {.cell}

```{.r .cell-code}
plot_dist(
  item = "GR02"
)
```

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-13-1.png){width=768}
:::

```{.r .cell-code}
get_desc(
  item = "GR02",
  name = "gratsgen"
)
```
:::



### CFA



::: {.cell}

```{.r .cell-code}
name <- "gratsgen"
model <- "
  grats_gen =~ GR02_01 + GR02_02 + GR02_03 + GR02_04 + GR02_05
"
fit <- 
  lavaan::sem(
    model = model, 
    data = d, 
    estimator = "MLR", 
    missing = "ML"
    )
```
:::



Model fit:



::: {.cell}

```{.r .cell-code}
facval <- 
  fit_tab(
    fit, 
    # reliability = TRUE, 
    # scaled = TRUE
    ) %T>%
  print()
```

::: {.cell-output .cell-output-stdout}

```
  chisq df   pvalue   cfi   tli rmsea   srmr
1  53.1  5 3.22e-10 0.979 0.958 0.131 0.0193
```


:::

```{.r .cell-code}
assign(
  paste0(
    name, 
    "_facval"
    ), 
  facval
  )
```
:::



Factor loadings:



::: {.cell}

```{.r .cell-code}
inspect(
  fit,
  what = "std"
  ) %$% 
  lambda
```

::: {.cell-output .cell-output-stdout}

```
        grts_g
GR02_01  0.853
GR02_02  0.908
GR02_03  0.851
GR02_04  0.837
GR02_05  0.848
```


:::
:::



Model fit was good.

## Gratifications specific
### Items

Using the participation platform it has been possible for me ...

_Information_

 1.	... to learn things I would not otherwise have noticed.
 2.	... to hear the opinion of others.
 3.	... to learn how other people tick.

_Relevance_

 4.	... to react to a subject that is very dear to me.
 5.	... to react to a subject that is important to me.
 6.	... to react to a subject that I am affected by.

_Political participation_

 7.	... to engage politically.
 8.	... to discuss political issues.
 9.	... to pursue my political interest.

_Idealism_

 10. ... to try to improve society.
 11. ... to advocate to something meaningful.
 12. ... to serve a good purpose.

_Extrinsic benefits_

 13. ... to do the responsible persons a favor.
 14. ... to soothe my guilty consciences.
 15. ... to fulfil my civic duty.

### Distributions



::: {.cell}

```{.r .cell-code}
plot_dist(
  item = "GR01",
  n_rows = 2
)
```

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-17-1.png){width=768}
:::

```{.r .cell-code}
get_desc(
  item = "GR01",
  name = "gratsspec"
)
```
:::



### CFA



::: {.cell}

```{.r .cell-code}
name <- "gratsspec"
model <- "
  grats_inf =~ GR01_01 + GR01_02 + GR01_03 
  grats_rel =~ GR01_04 + GR01_05 + GR01_06 
  grats_par =~ GR01_07 + GR01_08 + GR01_09
  grats_ide =~ GR01_10 + GR01_11 + GR01_12 
  grats_ext =~ GR01_13 + GR01_14 + GR01_15
  grats_spec =~ grats_inf + grats_rel + grats_par + grats_ide + grats_ext
"
fit <- lavaan::sem(
  model = model, 
  data = d, 
  estimator = "MLR", 
  missing = "ML"
  )
```
:::



Model fit:



::: {.cell}

```{.r .cell-code}
facval <- 
  fit_tab(
    fit, 
    # reliability = TRUE, 
    # scaled = TRUE
    ) %T>%
  print()
```

::: {.cell-output .cell-output-stdout}

```
  chisq df pvalue   cfi   tli  rmsea   srmr
1   441 85      0 0.934 0.918 0.0866 0.0527
```


:::

```{.r .cell-code}
assign(
  paste0(
    name, 
    "_facval"
    ), 
  facval
  )
```
:::



Factor loadings:



::: {.cell}

```{.r .cell-code}
inspect(
  fit,
  what = "std"
  ) %$% 
  lambda
```

::: {.cell-output .cell-output-stdout}

```
        grts_n grts_r grts_p grts_d grts_x grts_s
GR01_01  0.688  0.000  0.000  0.000  0.000      0
GR01_02  0.819  0.000  0.000  0.000  0.000      0
GR01_03  0.851  0.000  0.000  0.000  0.000      0
GR01_04  0.000  0.891  0.000  0.000  0.000      0
GR01_05  0.000  0.852  0.000  0.000  0.000      0
GR01_06  0.000  0.704  0.000  0.000  0.000      0
GR01_07  0.000  0.000  0.826  0.000  0.000      0
GR01_08  0.000  0.000  0.811  0.000  0.000      0
GR01_09  0.000  0.000  0.816  0.000  0.000      0
GR01_10  0.000  0.000  0.000  0.796  0.000      0
GR01_11  0.000  0.000  0.000  0.882  0.000      0
GR01_12  0.000  0.000  0.000  0.762  0.000      0
GR01_13  0.000  0.000  0.000  0.000  0.519      0
GR01_14  0.000  0.000  0.000  0.000  0.513      0
GR01_15  0.000  0.000  0.000  0.000  0.848      0
```


:::
:::



Model fit was okay.

## Privacy deliberation
### Items

Using the participation platform ...

 1. ... I have considered whether I could be disadvantaged by writing a comment.
 2. ... I have considered whether I could be advantaged by writing a comment.
 3. ... I have weighed up the advantages and disadvantages of writing a comment.
 4. ... I have thought about consequences of a possible comment.
 5. ... I have considered whether I should write a comment or not.

### Distributions



::: {.cell}

```{.r .cell-code}
plot_dist(
  item = "PD01"
)
```

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-21-1.png){width=768}
:::

```{.r .cell-code}
get_desc(
  item = "PD01",
  name = "pridel"
)
```
:::



### CFA



::: {.cell}

```{.r .cell-code}
name <- "pridel"
model <- "
  pri_delib =~ PD01_01 + PD01_02 + PD01_03 + PD01_04 + PD01_05
"
fit <- lavaan::sem(
  model = model, 
  data = d, 
  estimator = "MLR", 
  missing = "ML"
  )
```
:::



Model fit:



::: {.cell}

```{.r .cell-code}
facval <- 
  fit_tab(
    fit, 
    # reliability = TRUE, 
    # scaled = TRUE
    ) %T>%
  print()
```

::: {.cell-output .cell-output-stdout}

```
  chisq df   pvalue   cfi   tli  rmsea   srmr
1  27.4  5 4.71e-05 0.979 0.958 0.0896 0.0235
```


:::

```{.r .cell-code}
assign(
  paste0(
    name, 
    "_facval"), 
  facval
  )
```
:::



Factor loadings:



::: {.cell}

```{.r .cell-code}
inspect(
  fit,
  what = "std"
  ) %$%
  lambda
```

::: {.cell-output .cell-output-stdout}

```
        pr_dlb
PD01_01  0.849
PD01_02  0.653
PD01_03  0.691
PD01_04  0.752
PD01_05  0.656
```


:::
:::



Model fit was very good (which is nice, given it's a newly designed scale).

## Trust general
### Items

 1. The other users seemed trustworthy.  
 2. The operators of the participation platform seemed trustworthy.  
 3. The website seemed trustworthy.  

### Distributions



::: {.cell}

```{.r .cell-code}
plot_dist(
  item = "TR02"
)
```

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-25-1.png){width=768}
:::

```{.r .cell-code}
get_desc(
  item = "TR02",
  name = "trustgen"
)
```
:::



### CFA

Note that I constrained two loadings to be able to report model fit. I used operators and website, as these are theoretically closer.



::: {.cell}

```{.r .cell-code}
name <- "trustgen"
model <- "
trust =~ TR02_01 + a*TR02_02 + a*TR02_03
"
fit <- 
  lavaan::sem(
    model = model, 
    data = d, 
    estimator = "MLR", 
    missing = "ML"
    )
```
:::



Model fit:



::: {.cell}

```{.r .cell-code}
facval <- 
  fit_tab(
    fit, 
    # reliability = TRUE, 
    # scaled = TRUE
    ) %T>%
  print()
```

::: {.cell-output .cell-output-stdout}

```
  chisq df pvalue   cfi   tli  rmsea   srmr
1  1.86  1  0.173 0.999 0.997 0.0392 0.0125
```


:::

```{.r .cell-code}
assign(
  paste0(
    name, 
    "_facval"
    ), 
  facval
  )
```
:::



Factor loadings:



::: {.cell}

```{.r .cell-code}
inspect(
  fit,
  what = "std"
  ) %$%
  lambda
```

::: {.cell-output .cell-output-stdout}

```
        trust
TR02_01 0.660
TR02_02 0.923
TR02_03 0.898
```


:::
:::



Note that I constrained Items 5 and Item 9 to be equal. Explanation: First, they are theoretically related. Second, not constraining would yield to just-identified model, for which model fit cannot be interpreted meaningfully.

## Trust specific
### Items

_Community_

 1.	The comments of other users were useful.
 2.	The other users had good intentions.
 3.	I could rely on the statements of other users.

_Provider_

 4.	The operators of the participation platform have done a good job.
 5.	It was important to the operators that the users are satisfied with the participation platform.
 6.	I could rely on the statements of the operators of the participation platform.

_Information System_

 7.	The website worked well.
 8.	I had the impression that my data was necessary for the use of the website.
 9.	I found the website useful.

### Distributions



::: {.cell}

```{.r .cell-code}
plot_dist(
  item = "TR01",
  n_rows = 2
)
```

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-29-1.png){width=768}
:::

```{.r .cell-code}
get_desc(
  item = "TR01",
  name = "trustspec"
)
```
:::



### CFA



::: {.cell}

```{.r .cell-code}
name <- "trustspec"
model <- "
  trust_community =~ TR01_01 + TR01_02 + TR01_03
  trust_provider =~ TR01_04 + TR01_05 + TR01_06
  trust_system =~ TR01_07 + TR01_08 + TR01_09
  
  trust_spec =~ trust_community + trust_provider + trust_system
"
fit <- lavaan::sem(
  model = model, 
  data = d, 
  estimator = "MLR", 
  missing = "ML"
  )
```
:::



Model fit:



::: {.cell}

```{.r .cell-code}
facval <- 
  fit_tab(
    fit, 
    # reliability = TRUE, 
    # scaled = TRUE
    ) %T>% 
  print()
```

::: {.cell-output .cell-output-stdout}

```
  chisq df pvalue   cfi   tli  rmsea   srmr
1   153 24      0 0.959 0.939 0.0981 0.0351
```


:::

```{.r .cell-code}
assign(
  paste0(
    name, 
    "_facval"
    ), 
  facval
  )
```
:::



Factor loadings:



::: {.cell}

```{.r .cell-code}
inspect(
  fit,
  what = "std") %$%
  lambda
```

::: {.cell-output .cell-output-stdout}

```
        trst_c trst_p trst_sy trst_sp
TR01_01  0.814  0.000   0.000       0
TR01_02  0.765  0.000   0.000       0
TR01_03  0.822  0.000   0.000       0
TR01_04  0.000  0.884   0.000       0
TR01_05  0.000  0.779   0.000       0
TR01_06  0.000  0.793   0.000       0
TR01_07  0.000  0.000   0.690       0
TR01_08  0.000  0.000   0.662       0
TR01_09  0.000  0.000   0.817       0
```


:::
:::



Model fit was good (there's a Highwood case on level 1, which is not too problematic in this case.)

## Self-efficacy
### Items

 1. In principle, I felt able to write a comment.
 2. I felt technically competent enough to write a comment.
 3. In terms of the topic, I felt competent enough to express my opinion.
 4. I found it easy to express my opinion regarding the topic.
 5. I found it complicated to write a comment. (reversed)
 6. I was overburdened to write a comment. (reversed)

### Distributions



::: {.cell}

```{.r .cell-code}
plot_dist(
  item = "SE01"
)
```

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-33-1.png){width=768}
:::

```{.r .cell-code}
get_desc(
  item = "SE01",
  name = "selfeff"
)
```
:::



### CFA



::: {.cell}

```{.r .cell-code}
name <- "self-eff"
model <- "
  self_eff =~ SE01_01 + SE01_02 + SE01_03 + SE01_04 + SE01_05 + SE01_06
"
fit <- 
  lavaan::sem(
    model = model, 
    data = d, 
    estimator = "MLR", 
    missing = "ML"
    )
```
:::



Model fit:



::: {.cell}

```{.r .cell-code}
facval <- 
  fit_tab(
    fit, 
    # reliability = TRUE, 
    # scaled = TRUE
    ) %T>%
  print()
```

::: {.cell-output .cell-output-stdout}

```
  chisq df pvalue   cfi   tli rmsea   srmr
1   208  9      0 0.873 0.789 0.199 0.0671
```


:::

```{.r .cell-code}
assign(
  paste0(
    name, 
    "_facval"
    ), 
  facval
  )
```
:::



Factor loadings:



::: {.cell}

```{.r .cell-code}
inspect(
  fit,
  what = "std") %$%
  lambda
```

::: {.cell-output .cell-output-stdout}

```
        slf_ff
SE01_01  0.842
SE01_02  0.689
SE01_03  0.764
SE01_04  0.801
SE01_05  0.528
SE01_06  0.633
```


:::
:::



Shows significant misfit. I will delete inverted items, while allowing covariations between Items 1 and 2 (tech-oriented) and Items 3 and 4 (topic-oriented). 

### CFA 2



::: {.cell}

```{.r .cell-code}
name <- "selfeff"
model <- "
  self_eff_pos =~ SE01_01 + SE01_02 + SE01_03 + SE01_04
  SE01_01 ~~ x*SE01_02
  SE01_03 ~~ x*SE01_04
"
fit <- 
  lavaan::sem(
    model = model, 
    data = d, 
    estimator = "MLR", 
    missing = "ML"
    )
```
:::



Model fit:



::: {.cell}

```{.r .cell-code}
facval <- 
  fit_tab(
    fit, 
    # reliability = TRUE, 
    # scaled = TRUE
    ) %T>% 
  print()
```

::: {.cell-output .cell-output-stdout}

```
  chisq df  pvalue   cfi   tli rmsea   srmr
1  10.5  1 0.00117 0.991 0.945 0.131 0.0136
```


:::

```{.r .cell-code}
assign(
  paste0(
    name, 
    "_facval"
    ), 
  facval
  )
```
:::



Factor loadings:



::: {.cell}

```{.r .cell-code}
inspect(
  fit, 
  what = "std") %$%
  lambda
```

::: {.cell-output .cell-output-stdout}

```
        slf_f_
SE01_01  0.828
SE01_02  0.675
SE01_03  0.779
SE01_04  0.787
```


:::
:::



Adapted version shows better and adequate fit.

## Communication
### Density plot



::: {.cell}

```{.r .cell-code}
ggplot(
  gather(
    select(
      d, 
      words
      )
    ), 
  mapping = aes(
    x = value
    )
  ) +
  geom_density() +
  facet_wrap(
    ~key, 
    nrow = 1
    ) + 
  theme_bw()
```

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-40-1.png){width=768}
:::
:::



We can see that Communication is severely skewed. Let's calculate how many participants actually commented.



::: {.cell}

```{.r .cell-code}
no_words_perc <- 
  d %>% 
  select(words) %>% 
  table() %>% 
  prop.table() %>% 
  .[1] %>% 
  unlist()

words_m <- 
  mean(
    d$words, 
    na.rm = TRUE
    )
```
:::



0.58 percent wrote a comment. The average number of communicated words was 76.533.

Will hence be log-scaled for SEMs.

## Communication Logged



::: {.cell}

```{.r .cell-code}
ggplot(
  gather(
    select(
      d, 
      words_log
      )
    ), 
  mapping = aes(
    x = value
    )
  ) +
  geom_density() +
  facet_wrap(
    ~key, 
    nrow = 1
    ) + 
  theme_bw()
```

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-42-1.png){width=768}
:::
:::



## Baseline model

In what follows, please find the results of all variables combined in one model. This model will be used to extract factor scores.



::: {.cell}

```{.r .cell-code}
model_baseline <- "
  pri_con =~ PC01_01 + PC01_02 + PC01_04 + PC01_05 + PC01_06 + PC01_07
  grats_gen =~ GR02_01 + GR02_02 + GR02_03 + GR02_04 + GR02_05
  grats_inf =~ GR01_01 + GR01_02 + GR01_03 
  grats_rel =~ GR01_04 + GR01_05 + GR01_06 
  grats_par =~ GR01_07 + GR01_08 + GR01_09
  grats_ide =~ GR01_10 + GR01_11 + GR01_12 
  grats_ext =~ GR01_13 + GR01_14 + GR01_15
  grats_spec =~ grats_inf + grats_rel + grats_par + grats_ide + grats_ext
  pri_delib =~ PD01_01 + PD01_02 + PD01_03 + PD01_04 + PD01_05
  trust_gen =~ TR02_01 + TR02_02 + TR02_03
  trust_community =~ TR01_01 + TR01_02 + TR01_03
  trust_provider =~ TR01_04 + TR01_05 + TR01_06
  trust_system =~ TR01_07 + TR01_08 + TR01_09
  trust_spec =~ trust_community + trust_provider + trust_system
  self_eff =~ SE01_01 + SE01_02 + SE01_03 + SE01_04
    SE01_01 ~~ x*SE01_02
    SE01_03 ~~ x*SE01_04
  Words_log =~ words_log
  
  Words_log ~~ a1*pri_con + b1*grats_gen + c1*pri_delib + d1*self_eff + e1*trust_spec + f1*trust_gen + g1*grats_spec
"

fit_baseline <- 
  lavaan::sem(
    model_baseline, 
    data = d, 
    missing = "ML"
    )

summary(
  fit_baseline, 
  standardized = TRUE, 
  fit.measures = TRUE
  )
```

::: {.cell-output .cell-output-stdout}

```
lavaan 0.6-20 ended normally after 169 iterations

  Estimator                                         ML
  Optimization method                           NLMINB
  Number of model parameters                       181
  Number of equality constraints                     1

  Number of observations                           559
  Number of missing patterns                         3

Model Test User Model:
                                                      
  Test statistic                              3210.778
  Degrees of freedom                              1044
  P-value (Chi-square)                           0.000

Model Test Baseline Model:

  Test statistic                             22544.020
  Degrees of freedom                              1128
  P-value                                        0.000

User Model versus Baseline Model:

  Comparative Fit Index (CFI)                    0.899
  Tucker-Lewis Index (TLI)                       0.891
                                                      
  Robust Comparative Fit Index (CFI)             0.899
  Robust Tucker-Lewis Index (TLI)                0.891

Loglikelihood and Information Criteria:

  Loglikelihood user model (H0)             -37699.172
  Loglikelihood unrestricted model (H1)     -36093.783
                                                      
  Akaike (AIC)                               75758.344
  Bayesian (BIC)                             76537.051
  Sample-size adjusted Bayesian (SABIC)      75965.644

Root Mean Square Error of Approximation:

  RMSEA                                          0.061
  90 Percent confidence interval - lower         0.059
  90 Percent confidence interval - upper         0.063
  P-value H_0: RMSEA <= 0.050                    0.000
  P-value H_0: RMSEA >= 0.080                    0.000
                                                      
  Robust RMSEA                                   0.061
  90 Percent confidence interval - lower         0.059
  90 Percent confidence interval - upper         0.063
  P-value H_0: Robust RMSEA <= 0.050             0.000
  P-value H_0: Robust RMSEA >= 0.080             0.000

Standardized Root Mean Square Residual:

  SRMR                                           0.065

Parameter Estimates:

  Standard errors                             Standard
  Information                                 Observed
  Observed information based on                Hessian

Latent Variables:
                     Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
  pri_con =~                                                              
    PC01_01             1.000                               1.602    0.929
    PC01_02             0.994    0.027   36.762    0.000    1.592    0.901
    PC01_04             0.977    0.027   35.602    0.000    1.565    0.892
    PC01_05             1.002    0.026   38.002    0.000    1.605    0.910
    PC01_06             0.855    0.032   26.995    0.000    1.369    0.798
    PC01_07             0.996    0.025   40.248    0.000    1.595    0.925
  grats_gen =~                                                            
    GR02_01             1.000                               1.128    0.841
    GR02_02             1.120    0.040   27.925    0.000    1.264    0.893
    GR02_03             1.032    0.040   25.952    0.000    1.164    0.870
    GR02_04             0.990    0.040   25.023    0.000    1.118    0.850
    GR02_05             1.073    0.042   25.275    0.000    1.210    0.845
  grats_inf =~                                                            
    GR01_01             1.000                               0.978    0.696
    GR01_02             1.018    0.062   16.518    0.000    0.996    0.816
    GR01_03             1.114    0.066   16.803    0.000    1.089    0.847
  grats_rel =~                                                            
    GR01_04             1.000                               1.175    0.891
    GR01_05             0.943    0.034   27.545    0.000    1.108    0.857
    GR01_06             0.878    0.046   19.086    0.000    1.032    0.698
  grats_par =~                                                            
    GR01_07             1.000                               1.192    0.819
    GR01_08             0.938    0.042   22.086    0.000    1.119    0.816
    GR01_09             0.961    0.043   22.316    0.000    1.146    0.819
  grats_ide =~                                                            
    GR01_10             1.000                               1.148    0.791
    GR01_11             1.000    0.043   23.251    0.000    1.148    0.885
    GR01_12             0.928    0.048   19.451    0.000    1.065    0.764
  grats_ext =~                                                            
    GR01_13             1.000                               0.850    0.518
    GR01_14             1.003    0.108    9.311    0.000    0.853    0.509
    GR01_15             1.528    0.144   10.596    0.000    1.299    0.851
  grats_spec =~                                                           
    grats_inf           1.000                               0.843    0.843
    grats_rel           1.314    0.089   14.744    0.000    0.922    0.922
    grats_par           1.379    0.097   14.187    0.000    0.954    0.954
    grats_ide           1.306    0.095   13.801    0.000    0.938    0.938
    grats_ext           0.805    0.088    9.162    0.000    0.781    0.781
  pri_delib =~                                                            
    PD01_01             1.000                               1.494    0.866
    PD01_02             0.676    0.041   16.543    0.000    1.009    0.657
    PD01_03             0.704    0.043   16.324    0.000    1.052    0.676
    PD01_04             0.848    0.044   19.139    0.000    1.267    0.743
    PD01_05             0.718    0.045   15.904    0.000    1.073    0.648
  trust_gen =~                                                            
    TR02_01             1.000                               0.814    0.706
    TR02_02             1.253    0.064   19.587    0.000    1.020    0.886
    TR02_03             1.355    0.068   19.999    0.000    1.103    0.914
  trust_community =~                                                      
    TR01_01             1.000                               1.024    0.808
    TR01_02             0.824    0.043   18.982    0.000    0.844    0.768
    TR01_03             0.931    0.044   21.029    0.000    0.953    0.826
  trust_provider =~                                                       
    TR01_04             1.000                               1.042    0.869
    TR01_05             0.866    0.038   22.891    0.000    0.902    0.779
    TR01_06             0.864    0.036   24.282    0.000    0.900    0.812
  trust_system =~                                                         
    TR01_07             1.000                               0.812    0.690
    TR01_08             1.037    0.069   14.917    0.000    0.841    0.649
    TR01_09             1.376    0.073   18.761    0.000    1.117    0.831
  trust_spec =~                                                           
    trust_communty      1.000                               0.853    0.853
    trust_provider      1.166    0.062   18.798    0.000    0.977    0.977
    trust_system        0.959    0.061   15.790    0.000    1.032    1.032
  self_eff =~                                                             
    SE01_01             1.000                               1.135    0.821
    SE01_02             0.809    0.046   17.415    0.000    0.918    0.680
    SE01_03             0.923    0.048   19.106    0.000    1.048    0.781
    SE01_04             0.941    0.047   19.947    0.000    1.067    0.791
  Words_log =~                                                            
    words_log           1.000                               2.311    1.000

Covariances:
                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
 .SE01_01 ~~                                                            
   .SE01_02    (x)    0.115    0.029    3.993    0.000    0.115    0.148
 .SE01_03 ~~                                                            
   .SE01_04    (x)    0.115    0.029    3.993    0.000    0.115    0.166
  pri_con ~~                                                            
    Words_log (a1)   -0.559    0.161   -3.460    0.001   -0.151   -0.151
  grats_gen ~~                                                          
    Words_log (b1)    0.299    0.115    2.594    0.009    0.115    0.115
  pri_delib ~~                                                          
    Words_log (c1)   -0.653    0.160   -4.090    0.000   -0.189   -0.189
  self_eff ~~                                                           
    Words_log (d1)    1.016    0.132    7.711    0.000    0.387    0.387
  trust_spec ~~                                                         
    Words_log (e1)    0.358    0.091    3.920    0.000    0.177    0.177
  trust_gen ~~                                                          
    Words_log (f1)    0.323    0.086    3.767    0.000    0.171    0.171
  grats_spec ~~                                                         
    Words_log (g1)    0.426    0.089    4.771    0.000    0.224    0.224
  pri_con ~~                                                            
    grats_gen        -0.284    0.082   -3.465    0.001   -0.157   -0.157
    grats_spc        -0.109    0.060   -1.825    0.068   -0.083   -0.083
    pri_delib         1.356    0.131   10.332    0.000    0.567    0.567
    trust_gen        -0.548    0.068   -8.083    0.000   -0.420   -0.420
    trust_spc        -0.414    0.068   -6.120    0.000   -0.296   -0.296
    self_eff         -0.384    0.088   -4.352    0.000   -0.211   -0.211
  grats_gen ~~                                                          
    grats_spc         0.733    0.071   10.351    0.000    0.787    0.787
    pri_delib        -0.073    0.080   -0.918    0.359   -0.043   -0.043
    trust_gen         0.561    0.056    9.943    0.000    0.610    0.610
    trust_spc         0.748    0.067   11.174    0.000    0.759    0.759
    self_eff          0.464    0.066    7.002    0.000    0.362    0.362
  grats_spec ~~                                                         
    pri_delib         0.011    0.059    0.189    0.850    0.009    0.009
    trust_gen         0.441    0.049    9.002    0.000    0.657    0.657
    trust_spc         0.562    0.058    9.619    0.000    0.780    0.780
    self_eff          0.496    0.059    8.428    0.000    0.530    0.530
  pri_delib ~~                                                          
    trust_gen        -0.309    0.062   -4.974    0.000   -0.254   -0.254
    trust_spc        -0.136    0.063   -2.150    0.032   -0.104   -0.104
    self_eff         -0.335    0.087   -3.865    0.000   -0.197   -0.197
  trust_gen ~~                                                          
    trust_spc         0.664    0.060   11.079    0.000    0.934    0.934
    self_eff          0.486    0.055    8.838    0.000    0.526    0.526
  trust_spec ~~                                                         
    self_eff          0.534    0.059    9.072    0.000    0.539    0.539

Intercepts:
                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
   .PC01_01           3.293    0.073   45.160    0.000    3.293    1.910
   .PC01_02           3.327    0.075   44.525    0.000    3.327    1.883
   .PC01_04           3.222    0.074   43.394    0.000    3.222    1.835
   .PC01_05           3.263    0.075   43.748    0.000    3.263    1.850
   .PC01_06           3.004    0.073   41.409    0.000    3.004    1.751
   .PC01_07           3.224    0.073   44.188    0.000    3.224    1.869
   .GR02_01           4.281    0.057   75.413    0.000    4.281    3.190
   .GR02_02           4.596    0.060   76.742    0.000    4.596    3.246
   .GR02_03           5.131    0.057   90.628    0.000    5.131    3.833
   .GR02_04           5.089    0.056   91.559    0.000    5.089    3.873
   .GR02_05           4.692    0.061   77.447    0.000    4.692    3.276
   .GR01_01           4.878    0.059   82.009    0.000    4.878    3.469
   .GR01_02           5.436    0.052  105.391    0.000    5.436    4.458
   .GR01_03           5.283    0.054   97.076    0.000    5.283    4.106
   .GR01_04           4.925    0.056   88.265    0.000    4.925    3.733
   .GR01_05           5.086    0.055   93.032    0.000    5.086    3.935
   .GR01_06           4.660    0.063   74.538    0.000    4.660    3.153
   .GR01_07           4.682    0.062   76.077    0.000    4.682    3.218
   .GR01_08           5.066    0.058   87.375    0.000    5.066    3.696
   .GR01_09           4.841    0.059   81.782    0.000    4.841    3.459
   .GR01_10           4.547    0.061   74.048    0.000    4.547    3.132
   .GR01_11           4.964    0.055   90.449    0.000    4.964    3.826
   .GR01_12           4.760    0.059   80.678    0.000    4.760    3.412
   .GR01_13           4.079    0.069   58.781    0.000    4.079    2.486
   .GR01_14           3.039    0.071   42.918    0.000    3.039    1.815
   .GR01_15           4.410    0.065   68.283    0.000    4.410    2.888
   .PD01_01           3.658    0.073   50.136    0.000    3.658    2.121
   .PD01_02           3.352    0.065   51.628    0.000    3.352    2.184
   .PD01_03           4.191    0.066   63.662    0.000    4.191    2.693
   .PD01_04           4.081    0.072   56.578    0.000    4.081    2.393
   .PD01_05           4.351    0.070   62.149    0.000    4.351    2.629
   .TR02_01           4.846    0.049   99.399    0.000    4.846    4.204
   .TR02_02           5.383    0.049  110.608    0.000    5.383    4.678
   .TR02_03           5.390    0.051  105.542    0.000    5.390    4.464
   .TR01_01           4.764    0.054   88.924    0.000    4.764    3.761
   .TR01_02           4.844    0.046  104.196    0.000    4.844    4.407
   .TR01_03           4.615    0.049   94.568    0.000    4.615    4.000
   .TR01_04           5.403    0.051  106.479    0.000    5.403    4.504
   .TR01_05           5.200    0.049  106.180    0.000    5.200    4.491
   .TR01_06           5.129    0.047  109.405    0.000    5.129    4.627
   .TR01_07           5.725    0.050  114.939    0.000    5.725    4.864
   .TR01_08           4.834    0.055   88.152    0.000    4.834    3.728
   .TR01_09           5.179    0.057   91.089    0.000    5.179    3.853
   .SE01_01           5.277    0.058   90.256    0.000    5.277    3.820
   .SE01_02           5.523    0.057   96.667    0.000    5.523    4.092
   .SE01_03           5.224    0.057   92.028    0.000    5.224    3.895
   .SE01_04           5.137    0.057   89.906    0.000    5.137    3.805
   .words_log         1.834    0.098   18.765    0.000    1.834    0.794

Variances:
                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
   .PC01_01           0.406    0.032   12.792    0.000    0.406    0.137
   .PC01_02           0.587    0.042   14.012    0.000    0.587    0.188
   .PC01_04           0.631    0.044   14.284    0.000    0.631    0.205
   .PC01_05           0.534    0.039   13.741    0.000    0.534    0.172
   .PC01_06           1.066    0.068   15.606    0.000    1.066    0.362
   .PC01_07           0.430    0.033   13.040    0.000    0.430    0.145
   .GR02_01           0.528    0.038   13.980    0.000    0.528    0.293
   .GR02_02           0.407    0.033   12.355    0.000    0.407    0.203
   .GR02_03           0.436    0.033   13.158    0.000    0.436    0.243
   .GR02_04           0.478    0.035   13.822    0.000    0.478    0.277
   .GR02_05           0.587    0.042   14.018    0.000    0.587    0.286
   .GR01_01           1.021    0.072   14.111    0.000    1.021    0.516
   .GR01_02           0.496    0.042   11.733    0.000    0.496    0.334
   .GR01_03           0.468    0.045   10.366    0.000    0.468    0.283
   .GR01_04           0.359    0.034   10.458    0.000    0.359    0.206
   .GR01_05           0.443    0.036   12.232    0.000    0.443    0.265
   .GR01_06           1.121    0.074   15.181    0.000    1.121    0.513
   .GR01_07           0.696    0.052   13.418    0.000    0.696    0.329
   .GR01_08           0.628    0.047   13.436    0.000    0.628    0.334
   .GR01_09           0.646    0.049   13.309    0.000    0.646    0.330
   .GR01_10           0.790    0.056   14.010    0.000    0.790    0.374
   .GR01_11           0.366    0.034   10.622    0.000    0.366    0.217
   .GR01_12           0.811    0.056   14.458    0.000    0.811    0.417
   .GR01_13           1.969    0.131   14.980    0.000    1.969    0.731
   .GR01_14           2.076    0.138   15.094    0.000    2.076    0.741
   .GR01_15           0.643    0.097    6.622    0.000    0.643    0.276
   .PD01_01           0.745    0.082    9.135    0.000    0.745    0.250
   .PD01_02           1.338    0.090   14.870    0.000    1.338    0.568
   .PD01_03           1.317    0.092   14.323    0.000    1.317    0.543
   .PD01_04           1.303    0.096   13.632    0.000    1.303    0.448
   .PD01_05           1.588    0.107   14.850    0.000    1.588    0.580
   .TR02_01           0.666    0.044   15.225    0.000    0.666    0.501
   .TR02_02           0.284    0.024   12.036    0.000    0.284    0.215
   .TR02_03           0.241    0.024   10.204    0.000    0.241    0.165
   .TR01_01           0.556    0.045   12.395    0.000    0.556    0.346
   .TR01_02           0.496    0.037   13.434    0.000    0.496    0.411
   .TR01_03           0.423    0.036   11.680    0.000    0.423    0.317
   .TR01_04           0.353    0.028   12.654    0.000    0.353    0.245
   .TR01_05           0.527    0.036   14.562    0.000    0.527    0.393
   .TR01_06           0.418    0.030   14.162    0.000    0.418    0.341
   .TR01_07           0.727    0.045   16.067    0.000    0.727    0.524
   .TR01_08           0.973    0.060   16.122    0.000    0.973    0.579
   .TR01_09           0.560    0.045   12.449    0.000    0.560    0.310
   .SE01_01           0.621    0.053   11.716    0.000    0.621    0.325
   .SE01_02           0.979    0.066   14.732    0.000    0.979    0.537
   .SE01_03           0.701    0.056   12.593    0.000    0.701    0.390
   .SE01_04           0.684    0.054   12.760    0.000    0.684    0.375
   .words_log         0.000                               0.000    0.000
    pri_con           2.567    0.177   14.473    0.000    1.000    1.000
    grats_gen         1.273    0.105   12.120    0.000    1.000    1.000
   .grats_inf         0.277    0.038    7.350    0.000    0.290    0.290
   .grats_rel         0.207    0.033    6.248    0.000    0.150    0.150
   .grats_par         0.129    0.031    4.095    0.000    0.091    0.091
   .grats_ide         0.159    0.031    5.120    0.000    0.120    0.120
   .grats_ext         0.282    0.054    5.183    0.000    0.390    0.390
    grats_spec        0.680    0.092    7.415    0.000    1.000    1.000
    pri_delib         2.231    0.185   12.033    0.000    1.000    1.000
    trust_gen         0.663    0.071    9.319    0.000    1.000    1.000
   .trust_communty    0.286    0.036    7.946    0.000    0.272    0.272
   .trust_provider    0.050    0.019    2.591    0.010    0.046    0.046
   .trust_system     -0.043    0.016   -2.739    0.006   -0.065   -0.065
    trust_spec        0.763    0.082    9.256    0.000    1.000    1.000
    self_eff          1.288    0.116   11.135    0.000    1.000    1.000
    Words_log         5.341    0.319   16.718    0.000    1.000    1.000
```


:::
:::



# Descriptive analyses

I first report the factor validity of all variables combined.



::: {.cell}

```{.r .cell-code}
# extract model factor scores / predicted values for items & calc means
d_fs <- 
  lavPredict(
    fit_baseline, 
    type = "ov"
    ) %>% 
  as.data.frame() %>% 
  mutate(
    version = d$version, 
    grats_gen_fs = rowMeans(select(., starts_with("GR02"))),
    grats_spec_fs = rowMeans(select(., starts_with("GR01"))), 
    pri_con_fs = rowMeans(select(., starts_with("PC01"))),
    trust_gen_fs = rowMeans(select(., starts_with("TR02"))),
    trust_spec_fs = rowMeans(select(., starts_with("TR01"))),
    pri_del_fs = rowMeans(select(., starts_with("PD01"))),
    self_eff_fs = rowMeans(select(., starts_with("SE01")))) %>%
  select(
    version, 
    pri_con_fs, 
    grats_gen_fs, 
    grats_spec_fs, 
    pri_del_fs, 
    self_eff_fs, 
    trust_gen_fs, 
    trust_spec_fs, 
    words_log
    )

# combine d with d factor scores
d %<>% 
  cbind(
    select(
      d_fs, 
      -version, 
      -words_log
      )
    )

# rename for plotting
d_fs %<>% 
  set_names(
    c(
      "version", 
      var_names_breaks
      )
    )

# means of model predicted values
des <- 
  rbind(
    des_pricon,
    des_gratsgen,
    des_gratsspec,
    des_pridel,
    des_trustgen,
    des_trustspec,
    des_selfeff
    )

facval_tab <- 
  rbind(
    pricon_facval,
    gratsgen_facval,
    gratsspec_facval,
    pridel_facval,
    selfeff_facval,
    trustgen_facval,
    trustspec_facval
    ) %>%
      cbind(des[-c(8), ], .) %>%
      set_rownames(var_names[-c(7)])

facval_tab %>% 
  kable() %>% 
  kableExtra::kable_styling("striped") %>% 
  kableExtra::scroll_box(width = "100%")
```

::: {.cell-output-display}
`````{=html}
<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; "><table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">  </th>
   <th style="text-align:right;"> m </th>
   <th style="text-align:right;"> sd </th>
   <th style="text-align:right;"> chisq </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> pvalue </th>
   <th style="text-align:right;"> cfi </th>
   <th style="text-align:right;"> tli </th>
   <th style="text-align:right;"> rmsea </th>
   <th style="text-align:right;"> srmr </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Privacy concerns </td>
   <td style="text-align:right;"> 3.21 </td>
   <td style="text-align:right;"> 1.514 </td>
   <td style="text-align:right;"> 25.84 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 0.002 </td>
   <td style="text-align:right;"> 0.995 </td>
   <td style="text-align:right;"> 0.992 </td>
   <td style="text-align:right;"> 0.058 </td>
   <td style="text-align:right;"> 0.010 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Gratifications general </td>
   <td style="text-align:right;"> 4.76 </td>
   <td style="text-align:right;"> 1.219 </td>
   <td style="text-align:right;"> 53.09 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.979 </td>
   <td style="text-align:right;"> 0.958 </td>
   <td style="text-align:right;"> 0.131 </td>
   <td style="text-align:right;"> 0.019 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Gratifications specific </td>
   <td style="text-align:right;"> 4.71 </td>
   <td style="text-align:right;"> 1.019 </td>
   <td style="text-align:right;"> 441.39 </td>
   <td style="text-align:right;"> 85 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.934 </td>
   <td style="text-align:right;"> 0.918 </td>
   <td style="text-align:right;"> 0.087 </td>
   <td style="text-align:right;"> 0.053 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Privacy deliberation </td>
   <td style="text-align:right;"> 3.93 </td>
   <td style="text-align:right;"> 1.285 </td>
   <td style="text-align:right;"> 27.43 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.979 </td>
   <td style="text-align:right;"> 0.958 </td>
   <td style="text-align:right;"> 0.090 </td>
   <td style="text-align:right;"> 0.024 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Self-efficacy </td>
   <td style="text-align:right;"> 5.21 </td>
   <td style="text-align:right;"> 1.039 </td>
   <td style="text-align:right;"> 10.53 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.001 </td>
   <td style="text-align:right;"> 0.991 </td>
   <td style="text-align:right;"> 0.945 </td>
   <td style="text-align:right;"> 0.131 </td>
   <td style="text-align:right;"> 0.014 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Trust general </td>
   <td style="text-align:right;"> 5.08 </td>
   <td style="text-align:right;"> 0.942 </td>
   <td style="text-align:right;"> 1.86 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.173 </td>
   <td style="text-align:right;"> 0.999 </td>
   <td style="text-align:right;"> 0.997 </td>
   <td style="text-align:right;"> 0.039 </td>
   <td style="text-align:right;"> 0.012 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Words log </td>
   <td style="text-align:right;"> 5.25 </td>
   <td style="text-align:right;"> 1.118 </td>
   <td style="text-align:right;"> 153.12 </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.959 </td>
   <td style="text-align:right;"> 0.939 </td>
   <td style="text-align:right;"> 0.098 </td>
   <td style="text-align:right;"> 0.035 </td>
  </tr>
</tbody>
</table></div>

`````
:::
:::



In what follows, I report zero-order correlations, distributions, and scatterplots of the variables' factor scores.



::: {.cell}

```{.r .cell-code}
# corr_plot <- 
#   ggpairs(
#     select(d_fs, -version),
#     upper = list(
#       continuous = cor_plot
#     ),
#     lower = list(
#       continuous = wrap(
#         td::scat_plot, 
#         coords = c(1, 7, 0, 7)
#       )
#     )
#   ) + 
#   theme_bw()

# print(corr_plot)

# ggsave(
#   "figures/results/cor_plot.png",
#   plot = corr_plot,
#   bg = "white",
#   width = 8,
#   height = 8
# )
```
:::



# Power analyses

In what follows, I report power analyses for our study. Please note that I conduct a rudimentary power-analysis, assuming bivariate correlations. (At the time I was not yet aware of the existence of power analyses for multivariate structural equation models.)







I first estimate the sample size necessary to find small effects in 95% of all cases.



::: {.cell}

```{.r .cell-code}
# estimate pwr-samplesize
pwr.r.test(
  r = r_sesoi, 
  sig.level = alpha, 
  power = power_desired, 
  alternative = "greater"
  ) %T>%
  print %$% 
  n %>% 
  round(0) ->
  n_desired
```

::: {.cell-output .cell-output-stdout}

```

     approximate correlation power calculation (arctangh transformation) 

              n = 1077
              r = 0.1
      sig.level = 0.05
          power = 0.95
    alternative = greater
```


:::
:::



I then compute the power we have achieved with our final sample size to detect small effects.



::: {.cell}

```{.r .cell-code}
# compute pwr-achieved
pwr.r.test(
  n = n_final, 
  r = r_sesoi, 
  sig.level = alpha, 
  alternative = "greater"
  ) %T>%
  print %$% 
  power %>% 
  round(2) ->
  power_achieved
```

::: {.cell-output .cell-output-stdout}

```

     approximate correlation power calculation (arctangh transformation) 

              n = 559
              r = 0.1
      sig.level = 0.05
          power = 0.765
    alternative = greater
```


:::
:::



I finally compute what effect size we are likely to find in 95% of all cases given our final sample size.



::: {.cell}

```{.r .cell-code}
# estimate pwr-sensitivity
pwr.r.test(
  n = n_final, 
  power = power_desired, 
  sig.level = alpha, 
  alternative = "greater"
  ) %T>%
  print %$% 
  r %>% 
  round(2) ->
  r_sensitive
```

::: {.cell-output .cell-output-stdout}

```

     approximate correlation power calculation (arctangh transformation) 

              n = 559
              r = 0.138
      sig.level = 0.05
          power = 0.95
    alternative = greater
```


:::
:::



# Assumptions
## Multivariate normal distribution



::: {.cell}

```{.r .cell-code}
# create subset of data with all items that were used
items_used <- 
  c(
    "words", 
    "GR02_01", "GR02_02", "GR02_03", "GR02_04", "GR02_05", 
    "PC01_01", "PC01_02", "PC01_04", "PC01_05", "PC01_06", "PC01_07", 
    "TR01_01", "TR01_02", "TR01_03", "TR01_04", "TR01_05", "TR01_06", "TR01_07", "TR01_08", "TR01_09",
    "TR02_01", "TR02_02", "TR02_03", 
    "PD01_01", "PD01_02", "PD01_03", "PD01_04", "PD01_05", 
    "SE01_01", "SE01_02", "SE01_03", "SE01_04", 
    "male", "age", "edu"
    )
d_sub <- d[, items_used]

# test multivariate normal distribution
mvn_result <- mvn(d_sub, mvn_test = "mardia")
mvn_result$multivariateNormality
```

::: {.cell-output .cell-output-stdout}

```
NULL
```


:::
:::



Shows that multivariate normal distribution is violated. I hence use maximum likelihood estimation with robust standard errors and a Satorra-Bentler scaled test statistic.

## Influential cases

In what follows I test for influential cases in the baseline model, to detect potentially corrupt data (e.g., people who provided response patterns). Specifically, I compute Cook's distance.

[Note: The following lines stopped working after some time, potentially due to changes in package. I can hence not recreate them here, but refer to results obtained before.]



::: {.cell}

```{.r .cell-code}
# cooks_dis <- 
#   gCD(
#     d, 
#     model_baseline
#     ) %T>% 
#   plot()
```
:::



The following ten cases have a particularly strong influence on the baseline model.



::: {.cell}

```{.r .cell-code}
# infl_cases <- 
#   invisible(
#     rownames(
#       print(
#         cooks_dis
#         )
#       )
#     )
```
:::



Let's inspect these cases.



::: {.cell}

```{.r .cell-code}
# infl_cases_tokens <- 
#   d[infl_cases, "case_token"] %>% 
#   as_vector()

# d %>% 
#   filter(
#     case_token %in% infl_cases_tokens
#     ) %>% 
#   select(
#     case_token, 
#     GR01_01:SE01_06, 
#     topics_entered:reactions, 
#     -SO01_01, 
#     TIME_SUM_t1, 
#     TIME_SUM_t2
#     ) %>% 
#   kable() %>% 
#   kable_styling("striped") %>% 
#   scroll_box(width = "100%")
```
:::



These data do not reveal potential cases of response patterns. Indeed, answer times suggest that respondents were diligent.

# Results
## Preregistered
### Privacy calculus



::: {.cell}

```{.r .cell-code}
model <- "
  pri_con =~ PC01_01 + PC01_02 + PC01_04 + PC01_05 + PC01_06 + PC01_07 
  grats_gen =~ GR02_01 + GR02_02 + GR02_03 + GR02_04 + GR02_05
  pri_delib =~ PD01_01 + PD01_02 + PD01_03 + PD01_04 + PD01_05
  self_eff =~ SE01_01 + SE01_02 + SE01_03 + SE01_04
  SE01_01 ~~ x*SE01_02
  SE01_03 ~~ x*SE01_04
  trust_community =~ TR01_01 + TR01_02 + TR01_03
  trust_provider =~ TR01_04 + TR01_05 + TR01_06
  trust_system =~ TR01_07 + TR01_08 + TR01_09
  
  trust_spec =~ trust_community + trust_provider + trust_system

words_log ~ a1*pri_con + b1*grats_gen + c1*pri_delib + d1*self_eff + e1*trust_spec

# Covariates
words_log + GR02_01 + GR02_02 + GR02_03 + GR02_04 + GR02_05 + PC01_01 + PC01_02 + PC01_04 + PC01_05 + PC01_06 + PC01_07 + TR01_01 + TR01_02 + TR01_03 + TR01_04 + TR01_05 + TR01_06 + TR01_07 + TR01_08 + TR01_09 + PD01_01 + PD01_02 + PD01_03 + PD01_04 + PD01_05 + SE01_01 + SE01_02 + SE01_03 + SE01_04 ~ male + age + edu

# Covariances
male ~~ age + edu
age ~~ edu
"

fit_prereg <- 
  lavaan::sem(
    model, 
    data = d, 
    estimator = "MLR", 
    missing = "ML"
    )

summary(
  fit_prereg, 
  fit = TRUE, 
  std = TRUE
  )
```

::: {.cell-output .cell-output-stdout}

```
lavaan 0.6-20 ended normally after 343 iterations

  Estimator                                         ML
  Optimization method                           NLMINB
  Number of model parameters                       208
  Number of equality constraints                     1

  Number of observations                           559
  Number of missing patterns                         4

Model Test User Model:
                                              Standard      Scaled
  Test Statistic                              1220.157     934.980
  Degrees of freedom                               387         387
  P-value (Chi-square)                           0.000       0.000
  Scaling correction factor                                  1.305
    Yuan-Bentler correction (Mplus variant)                       

Model Test Baseline Model:

  Test statistic                             13379.428   10044.844
  Degrees of freedom                               528         528
  P-value                                        0.000       0.000
  Scaling correction factor                                  1.332

User Model versus Baseline Model:

  Comparative Fit Index (CFI)                    0.935       0.942
  Tucker-Lewis Index (TLI)                       0.912       0.921
                                                                  
  Robust Comparative Fit Index (CFI)                         0.944
  Robust Tucker-Lewis Index (TLI)                            0.924

Loglikelihood and Information Criteria:

  Loglikelihood user model (H0)             -27300.107  -27300.107
  Scaling correction factor                                  1.242
      for the MLR correction                                      
  Loglikelihood unrestricted model (H1)     -26690.028  -26690.028
  Scaling correction factor                                  1.285
      for the MLR correction                                      
                                                                  
  Akaike (AIC)                               55014.214   55014.214
  Bayesian (BIC)                             55909.727   55909.727
  Sample-size adjusted Bayesian (SABIC)      55252.609   55252.609

Root Mean Square Error of Approximation:

  RMSEA                                          0.062       0.050
  90 Percent confidence interval - lower         0.058       0.047
  90 Percent confidence interval - upper         0.066       0.054
  P-value H_0: RMSEA <= 0.050                    0.000       0.434
  P-value H_0: RMSEA >= 0.080                    0.000       0.000
                                                                  
  Robust RMSEA                                               0.057
  90 Percent confidence interval - lower                     0.052
  90 Percent confidence interval - upper                     0.062
  P-value H_0: Robust RMSEA <= 0.050                         0.007
  P-value H_0: Robust RMSEA >= 0.080                         0.000

Standardized Root Mean Square Residual:

  SRMR                                           0.049       0.049

Parameter Estimates:

  Standard errors                             Sandwich
  Information bread                           Observed
  Observed information based on                Hessian

Latent Variables:
                     Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
  pri_con =~                                                              
    PC01_01             1.000                               1.596    0.926
    PC01_02             0.990    0.027   36.315    0.000    1.581    0.895
    PC01_04             0.972    0.027   35.734    0.000    1.551    0.884
    PC01_05             1.002    0.024   42.571    0.000    1.600    0.907
    PC01_06             0.855    0.038   22.731    0.000    1.365    0.796
    PC01_07             0.995    0.023   43.861    0.000    1.588    0.920
  grats_gen =~                                                            
    GR02_01             1.000                               1.129    0.841
    GR02_02             1.120    0.033   33.538    0.000    1.265    0.893
    GR02_03             1.025    0.048   21.437    0.000    1.158    0.865
    GR02_04             0.988    0.048   20.466    0.000    1.115    0.849
    GR02_05             1.076    0.040   27.087    0.000    1.215    0.848
  pri_delib =~                                                            
    PD01_01             1.000                               1.472    0.853
    PD01_02             0.669    0.048   13.902    0.000    0.985    0.642
    PD01_03             0.709    0.055   12.939    0.000    1.043    0.670
    PD01_04             0.843    0.047   17.856    0.000    1.240    0.727
    PD01_05             0.717    0.050   14.304    0.000    1.056    0.638
  self_eff =~                                                             
    SE01_01             1.000                               1.116    0.808
    SE01_02             0.816    0.057   14.291    0.000    0.911    0.675
    SE01_03             0.934    0.046   20.343    0.000    1.043    0.778
    SE01_04             0.953    0.043   22.351    0.000    1.063    0.788
  trust_community =~                                                      
    TR01_01             1.000                               1.026    0.810
    TR01_02             0.814    0.053   15.343    0.000    0.835    0.760
    TR01_03             0.917    0.048   19.282    0.000    0.941    0.815
  trust_provider =~                                                       
    TR01_04             1.000                               1.057    0.881
    TR01_05             0.857    0.040   21.572    0.000    0.906    0.782
    TR01_06             0.833    0.040   20.611    0.000    0.881    0.794
  trust_system =~                                                         
    TR01_07             1.000                               0.800    0.680
    TR01_08             1.058    0.081   13.016    0.000    0.846    0.653
    TR01_09             1.398    0.086   16.276    0.000    1.119    0.832
  trust_spec =~                                                           
    trust_communty      1.000                               0.847    0.847
    trust_provider      1.160    0.080   14.545    0.000    0.954    0.954
    trust_system        0.975    0.079   12.345    0.000    1.059    1.059

Regressions:
                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
  words_log ~                                                           
    pri_con   (a1)   -0.035    0.079   -0.447    0.655   -0.057   -0.025
    grats_gen (b1)   -0.010    0.147   -0.068    0.946   -0.011   -0.005
    pri_delib (c1)   -0.164    0.093   -1.764    0.078   -0.242   -0.105
    self_eff  (d1)    0.770    0.142    5.424    0.000    0.860    0.372
    trust_spc (e1)   -0.087    0.240   -0.363    0.716   -0.076   -0.033
    male              0.020    0.199    0.098    0.922    0.020    0.004
    age               0.005    0.006    0.805    0.421    0.005    0.033
    edu               0.230    0.117    1.967    0.049    0.230    0.084
  GR02_01 ~                                                             
    male             -0.127    0.115   -1.098    0.272   -0.127   -0.047
    age               0.000    0.004    0.089    0.929    0.000    0.004
    edu               0.006    0.068    0.082    0.935    0.006    0.004
  GR02_02 ~                                                             
    male             -0.068    0.120   -0.565    0.572   -0.068   -0.024
    age               0.006    0.004    1.538    0.124    0.006    0.068
    edu              -0.078    0.071   -1.107    0.268   -0.078   -0.047
  GR02_03 ~                                                             
    male             -0.027    0.116   -0.230    0.818   -0.027   -0.010
    age               0.001    0.004    0.303    0.762    0.001    0.013
    edu              -0.080    0.067   -1.199    0.231   -0.080   -0.050
  GR02_04 ~                                                             
    male              0.027    0.113    0.239    0.811    0.027    0.010
    age               0.005    0.004    1.296    0.195    0.005    0.057
    edu              -0.070    0.067   -1.035    0.301   -0.070   -0.045
  GR02_05 ~                                                             
    male             -0.141    0.123   -1.141    0.254   -0.141   -0.049
    age              -0.004    0.004   -0.878    0.380   -0.004   -0.039
    edu               0.014    0.072    0.193    0.847    0.014    0.008
  PC01_01 ~                                                             
    male             -0.184    0.151   -1.218    0.223   -0.184   -0.053
    age              -0.004    0.005   -0.830    0.407   -0.004   -0.037
    edu               0.114    0.087    1.309    0.190    0.114    0.056
  PC01_02 ~                                                             
    male             -0.304    0.154   -1.977    0.048   -0.304   -0.086
    age              -0.008    0.005   -1.674    0.094   -0.008   -0.072
    edu               0.051    0.089    0.577    0.564    0.051    0.025
  PC01_04 ~                                                             
    male             -0.227    0.152   -1.487    0.137   -0.227   -0.065
    age              -0.010    0.005   -1.990    0.047   -0.010   -0.086
    edu               0.117    0.089    1.320    0.187    0.117    0.056
  PC01_05 ~                                                             
    male             -0.100    0.154   -0.650    0.516   -0.100   -0.028
    age              -0.006    0.005   -1.174    0.240   -0.006   -0.051
    edu               0.094    0.090    1.049    0.294    0.094    0.045
  PC01_06 ~                                                             
    male             -0.110    0.150   -0.734    0.463   -0.110   -0.032
    age              -0.005    0.005   -1.065    0.287   -0.005   -0.046
    edu               0.047    0.087    0.540    0.589    0.047    0.023
  PC01_07 ~                                                             
    male             -0.176    0.150   -1.172    0.241   -0.176   -0.051
    age              -0.007    0.005   -1.347    0.178   -0.007   -0.059
    edu               0.086    0.087    0.987    0.324    0.086    0.042
  TR01_01 ~                                                             
    male             -0.296    0.108   -2.742    0.006   -0.296   -0.117
    age              -0.004    0.004   -1.101    0.271   -0.004   -0.049
    edu               0.005    0.061    0.076    0.940    0.005    0.003
  TR01_02 ~                                                             
    male             -0.139    0.095   -1.464    0.143   -0.139   -0.063
    age              -0.002    0.003   -0.557    0.577   -0.002   -0.025
    edu               0.021    0.053    0.385    0.700    0.021    0.016
  TR01_03 ~                                                             
    male             -0.133    0.099   -1.343    0.179   -0.133   -0.058
    age              -0.004    0.003   -1.201    0.230   -0.004   -0.054
    edu              -0.006    0.060   -0.098    0.922   -0.006   -0.004
  TR01_04 ~                                                             
    male             -0.086    0.104   -0.825    0.409   -0.086   -0.036
    age               0.000    0.003    0.114    0.909    0.000    0.005
    edu              -0.053    0.058   -0.903    0.366   -0.053   -0.037
  TR01_05 ~                                                             
    male             -0.043    0.099   -0.429    0.668   -0.043   -0.018
    age               0.001    0.003    0.357    0.721    0.001    0.016
    edu               0.014    0.058    0.241    0.810    0.014    0.010
  TR01_06 ~                                                             
    male              0.048    0.096    0.498    0.619    0.048    0.021
    age              -0.004    0.003   -1.235    0.217   -0.004   -0.053
    edu               0.021    0.056    0.370    0.712    0.021    0.016
  TR01_07 ~                                                             
    male              0.093    0.100    0.927    0.354    0.093    0.039
    age              -0.004    0.003   -1.166    0.244   -0.004   -0.050
    edu              -0.058    0.058   -0.995    0.320   -0.058   -0.042
  TR01_08 ~                                                             
    male              0.028    0.112    0.255    0.799    0.028    0.011
    age               0.003    0.004    0.831    0.406    0.003    0.036
    edu              -0.096    0.065   -1.473    0.141   -0.096   -0.062
  TR01_09 ~                                                             
    male             -0.120    0.115   -1.039    0.299   -0.120   -0.045
    age              -0.002    0.004   -0.401    0.688   -0.002   -0.018
    edu              -0.148    0.068   -2.182    0.029   -0.148   -0.092
  PD01_01 ~                                                             
    male             -0.176    0.148   -1.194    0.232   -0.176   -0.051
    age              -0.015    0.005   -3.274    0.001   -0.015   -0.137
    edu              -0.027    0.085   -0.322    0.748   -0.027   -0.013
  PD01_02 ~                                                             
    male             -0.118    0.131   -0.900    0.368   -0.118   -0.038
    age              -0.014    0.004   -3.439    0.001   -0.014   -0.141
    edu               0.030    0.077    0.384    0.701    0.030    0.016
  PD01_03 ~                                                             
    male             -0.321    0.132   -2.425    0.015   -0.321   -0.103
    age              -0.004    0.004   -1.024    0.306   -0.004   -0.044
    edu               0.065    0.080    0.811    0.418    0.065    0.035
  PD01_04 ~                                                             
    male             -0.411    0.144   -2.846    0.004   -0.411   -0.120
    age              -0.009    0.005   -1.904    0.057   -0.009   -0.082
    edu               0.102    0.085    1.206    0.228    0.102    0.051
  PD01_05 ~                                                             
    male             -0.205    0.142   -1.441    0.150   -0.205   -0.062
    age              -0.012    0.004   -2.698    0.007   -0.012   -0.111
    edu              -0.002    0.084   -0.018    0.985   -0.002   -0.001
  SE01_01 ~                                                             
    male              0.121    0.118    1.020    0.308    0.121    0.044
    age               0.000    0.004    0.015    0.988    0.000    0.001
    edu               0.207    0.068    3.046    0.002    0.207    0.126
  SE01_02 ~                                                             
    male              0.060    0.112    0.537    0.591    0.060    0.022
    age              -0.013    0.004   -3.590    0.000   -0.013   -0.151
    edu               0.194    0.066    2.938    0.003    0.194    0.121
  SE01_03 ~                                                             
    male              0.195    0.114    1.705    0.088    0.195    0.073
    age               0.001    0.004    0.254    0.800    0.001    0.011
    edu               0.138    0.067    2.065    0.039    0.138    0.087
  SE01_04 ~                                                             
    male              0.054    0.115    0.473    0.636    0.054    0.020
    age               0.007    0.004    2.055    0.040    0.007    0.086
    edu               0.122    0.066    1.837    0.066    0.122    0.076

Covariances:
                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
 .SE01_01 ~~                                                            
   .SE01_02    (x)    0.106    0.044    2.400    0.016    0.106    0.140
 .SE01_03 ~~                                                            
   .SE01_04    (x)    0.106    0.044    2.400    0.016    0.106    0.156
  male ~~                                                               
    age               0.757    0.328    2.312    0.021    0.757    0.097
    edu               0.052    0.018    2.966    0.003    0.052    0.125
  age ~~                                                                
    edu              -1.018    0.547   -1.862    0.063   -1.018   -0.078
  pri_con ~~                                                            
    grats_gen        -0.279    0.095   -2.924    0.003   -0.155   -0.155
    pri_delib         1.320    0.130   10.121    0.000    0.562    0.562
    self_eff         -0.383    0.091   -4.210    0.000   -0.215   -0.215
    trust_spec       -0.411    0.069   -5.932    0.000   -0.296   -0.296
  grats_gen ~~                                                          
    pri_delib        -0.068    0.102   -0.669    0.503   -0.041   -0.041
    self_eff          0.466    0.067    6.999    0.000    0.370    0.370
    trust_spec        0.746    0.080    9.308    0.000    0.760    0.760
  pri_delib ~~                                                          
    self_eff         -0.323    0.094   -3.437    0.001   -0.197   -0.197
    trust_spec       -0.140    0.079   -1.768    0.077   -0.109   -0.109
  self_eff ~~                                                           
    trust_spec        0.532    0.059    8.938    0.000    0.548    0.548

Intercepts:
                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
   .PC01_01           3.360    0.292   11.525    0.000    3.360    1.949
   .PC01_02           3.760    0.304   12.368    0.000    3.760    2.128
   .PC01_04           3.562    0.297   11.993    0.000    3.562    2.029
   .PC01_05           3.405    0.304   11.200    0.000    3.405    1.931
   .PC01_06           3.207    0.288   11.128    0.000    3.207    1.870
   .PC01_07           3.452    0.294   11.749    0.000    3.452    2.001
   .GR02_01           4.318    0.224   19.261    0.000    4.318    3.217
   .GR02_02           4.489    0.244   18.372    0.000    4.489    3.170
   .GR02_03           5.239    0.221   23.652    0.000    5.239    3.914
   .GR02_04           4.983    0.221   22.508    0.000    4.983    3.792
   .GR02_05           4.902    0.254   19.325    0.000    4.902    3.422
   .PD01_01           4.495    0.290   15.524    0.000    4.495    2.605
   .PD01_02           4.000    0.248   16.125    0.000    4.000    2.605
   .PD01_03           4.432    0.269   16.447    0.000    4.432    2.847
   .PD01_04           4.507    0.295   15.294    0.000    4.507    2.643
   .PD01_05           4.999    0.276   18.098    0.000    4.999    3.020
   .SE01_01           4.833    0.249   19.394    0.000    4.833    3.500
   .SE01_02           5.739    0.234   24.499    0.000    5.739    4.252
   .SE01_03           4.828    0.226   21.365    0.000    4.828    3.602
   .SE01_04           4.542    0.234   19.415    0.000    4.542    3.366
   .TR01_01           5.084    0.219   23.244    0.000    5.084    4.014
   .TR01_02           4.956    0.189   26.210    0.000    4.956    4.509
   .TR01_03           4.877    0.200   24.387    0.000    4.877    4.227
   .TR01_04           5.524    0.204   27.020    0.000    5.524    4.605
   .TR01_05           5.142    0.197   26.057    0.000    5.142    4.441
   .TR01_06           5.240    0.189   27.798    0.000    5.240    4.728
   .TR01_07           5.960    0.193   30.944    0.000    5.960    5.063
   .TR01_08           4.860    0.218   22.299    0.000    4.860    3.749
   .TR01_09           5.584    0.231   24.200    0.000    5.584    4.154
   .words_log         1.169    0.380    3.079    0.002    1.169    0.506
    male              0.493    0.021   23.296    0.000    0.493    0.986
    age              46.132    0.658   70.151    0.000   46.132    2.967
    edu               1.852    0.036   51.966    0.000    1.852    2.198

Variances:
                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
   .PC01_01           0.403    0.050    8.075    0.000    0.403    0.135
   .PC01_02           0.579    0.103    5.652    0.000    0.579    0.186
   .PC01_04           0.626    0.077    8.134    0.000    0.626    0.203
   .PC01_05           0.533    0.064    8.359    0.000    0.533    0.171
   .PC01_06           1.067    0.116    9.239    0.000    1.067    0.363
   .PC01_07           0.430    0.065    6.570    0.000    0.430    0.144
   .GR02_01           0.523    0.053    9.777    0.000    0.523    0.290
   .GR02_02           0.389    0.040    9.823    0.000    0.389    0.194
   .GR02_03           0.446    0.073    6.073    0.000    0.446    0.249
   .GR02_04           0.473    0.048    9.904    0.000    0.473    0.274
   .GR02_05           0.567    0.062    9.160    0.000    0.567    0.277
   .PD01_01           0.742    0.110    6.717    0.000    0.742    0.249
   .PD01_02           1.332    0.127   10.455    0.000    1.332    0.565
   .PD01_03           1.301    0.128   10.186    0.000    1.301    0.537
   .PD01_04           1.297    0.147    8.844    0.000    1.297    0.446
   .PD01_05           1.576    0.128   12.340    0.000    1.576    0.575
   .SE01_01           0.625    0.087    7.177    0.000    0.625    0.328
   .SE01_02           0.917    0.118    7.772    0.000    0.917    0.504
   .SE01_03           0.684    0.095    7.212    0.000    0.684    0.381
   .SE01_04           0.667    0.078    8.544    0.000    0.667    0.366
   .TR01_01           0.524    0.067    7.792    0.000    0.524    0.327
   .TR01_02           0.505    0.057    8.907    0.000    0.505    0.418
   .TR01_03           0.437    0.045    9.653    0.000    0.437    0.329
   .TR01_04           0.317    0.034    9.461    0.000    0.317    0.220
   .TR01_05           0.520    0.053    9.727    0.000    0.520    0.388
   .TR01_06           0.449    0.042   10.628    0.000    0.449    0.365
   .TR01_07           0.739    0.057   13.027    0.000    0.739    0.533
   .TR01_08           0.955    0.079   12.047    0.000    0.955    0.568
   .TR01_09           0.534    0.061    8.705    0.000    0.534    0.295
   .words_log         4.459    0.215   20.710    0.000    4.459    0.835
    male              0.250    0.000  817.190    0.000    0.250    1.000
    age             241.746    9.947   24.303    0.000  241.746    1.000
    edu               0.710    0.021   34.591    0.000    0.710    1.000
    pri_con           2.549    0.144   17.662    0.000    1.000    1.000
    grats_gen         1.275    0.114   11.180    0.000    1.000    1.000
    pri_delib         2.167    0.157   13.796    0.000    1.000    1.000
    self_eff          1.245    0.113   10.991    0.000    1.000    1.000
   .trust_communty    0.297    0.045    6.672    0.000    0.282    0.282
   .trust_provider    0.100    0.031    3.213    0.001    0.089    0.089
   .trust_system     -0.078    0.019   -4.092    0.000   -0.122   -0.122
    trust_spec        0.756    0.095    7.920    0.000    1.000    1.000
```


:::

```{.r .cell-code}
rsquare_fit_prereg <- 
  inspect(
    fit_prereg, 
    what = "rsquare"
    )["words"]
```
:::



Results show that there's only one significant predictor of Communication, being self-efficacy. The other predictors are in the direction as planned, albeit not significant. Trust, however, shows the inverse relation as effect, that is, more trust, less communication.

### Effects of popularity cues
#### Confidence intervals

The easiest way to assess the effect of the experimental manipulation on the variables is by visualizing their means. If 83% confidence intervals don't overlap, the variables differ significantly across the conditions. One can quickly see that there aren't any major effects.



::: {.cell}

```{.r .cell-code}
# violin plot
fig_fs_m <- 
  ggplot(
    gather(
      d_fs, 
      variable, 
      value, 
      -version) %>% 
      mutate(
        variable = factor(
          variable, 
          levels = var_names_breaks
          )
        ),
    aes(
      x = version, 
      y = value, 
      fill = version
      )
    ) +
  geom_violin(
    trim = TRUE
    ) +
  stat_summary(
    fun.y = mean, 
    geom = "point"
    ) +
  stat_summary(
    fun.data = mean_se, 
    fun.args = list(mult = 1.39), 
    geom = "errorbar", 
    width = .5
    ) + 
  facet_wrap(
    ~ variable, 
    nrow = 2
    ) + 
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(hjust = .5),
    panel.spacing = unit(.9, "lines"),
    text = element_text(size = 12),
    legend.position="none",
    legend.title = element_blank()) +
  coord_cartesian(
    ylim = c(0, 7)
    ) +
  scale_fill_brewer(
    palette = "Greys"
    )
ggsave(
  "figures/results/violin_plot.png", 
  width = 8, 
  height = 6
  )
fig_fs_m
```

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-55-1.png){width=768}
:::
:::



#### SEM

In what follows, I also report explicit statistical tests of the differences between the conditions using contrasts. 

**Like & Like-Dislike vs. Control**



::: {.cell}

```{.r .cell-code}
model <- "
  pri_con =~ PC01_01 + PC01_02 + PC01_04 + PC01_05 + PC01_06 + PC01_07
  grats_gen =~ GR02_01 + GR02_02 + GR02_03 + GR02_04 + GR02_05
  pri_delib =~ PD01_01 + PD01_02 + PD01_03 + PD01_04 + PD01_05
  self_eff =~ SE01_01 + SE01_02 + SE01_03 + SE01_04
  SE01_01 ~~ x*SE01_02
  SE01_03 ~~ x*SE01_04
  trust_community =~ TR01_01 + TR01_02 + TR01_03
  trust_provider =~ TR01_04 + TR01_05 + TR01_06
  trust_system =~ TR01_07 + TR01_08 + TR01_09
  trust_spec =~ trust_community + trust_provider + trust_system

  pri_con + grats_gen + pri_delib + self_eff + trust_spec ~ like + likedislike
  words ~ a*pri_con + b*grats_gen + c*pri_delib + d*self_eff + e*trust_spec + f*like + g*likedislike

# Covariates
words + GR02_01 + GR02_02 + GR02_03 + GR02_04 + GR02_05 + PC01_01 + PC01_02 + PC01_04 + PC01_05 + PC01_06 + PC01_07 + TR01_01 + TR01_02 + TR01_03 + TR01_04 + TR01_05 + TR01_06 + TR01_07 + TR01_08 + TR01_09 + PD01_01 + PD01_02 + PD01_03 + PD01_04 + PD01_05 + SE01_01 + SE01_02 + SE01_03 + SE01_04 ~ male + age + edu
"

fit_lik_ctrl <- 
  lavaan::sem(
    model = model, 
    data = d, 
    estimator = "MLR", 
    missing = "ML", 
    fixed.x = FALSE
    )

summary(
  fit_lik_ctrl, 
  fit = TRUE, 
  std = TRUE
  )
```

::: {.cell-output .cell-output-stdout}

```
lavaan 0.6-20 ended normally after 595 iterations

  Estimator                                         ML
  Optimization method                           NLMINB
  Number of model parameters                       221
  Number of equality constraints                     1

  Number of observations                           559
  Number of missing patterns                         4

Model Test User Model:
                                              Standard      Scaled
  Test Statistic                              2019.407    1485.951
  Degrees of freedom                               445         445
  P-value (Chi-square)                           0.000       0.000
  Scaling correction factor                                  1.359
    Yuan-Bentler correction (Mplus variant)                       

Model Test Baseline Model:

  Test statistic                             13359.444   10346.817
  Degrees of freedom                               585         585
  P-value                                        0.000       0.000
  Scaling correction factor                                  1.291

User Model versus Baseline Model:

  Comparative Fit Index (CFI)                    0.877       0.893
  Tucker-Lewis Index (TLI)                       0.838       0.860
                                                                  
  Robust Comparative Fit Index (CFI)                         0.888
  Robust Tucker-Lewis Index (TLI)                            0.853

Loglikelihood and Information Criteria:

  Loglikelihood user model (H0)             -30981.342  -30981.342
  Scaling correction factor                                  1.201
      for the MLR correction                                      
  Loglikelihood unrestricted model (H1)     -29971.639  -29971.639
  Scaling correction factor                                  1.308
      for the MLR correction                                      
                                                                  
  Akaike (AIC)                               62402.685   62402.685
  Bayesian (BIC)                             63354.438   63354.438
  Sample-size adjusted Bayesian (SABIC)      62656.052   62656.052

Root Mean Square Error of Approximation:

  RMSEA                                          0.080       0.065
  90 Percent confidence interval - lower         0.076       0.062
  90 Percent confidence interval - upper         0.083       0.068
  P-value H_0: RMSEA <= 0.050                    0.000       0.000
  P-value H_0: RMSEA >= 0.080                    0.422       0.000
                                                                  
  Robust RMSEA                                               0.075
  90 Percent confidence interval - lower                     0.071
  90 Percent confidence interval - upper                     0.080
  P-value H_0: Robust RMSEA <= 0.050                         0.000
  P-value H_0: Robust RMSEA >= 0.080                         0.036

Standardized Root Mean Square Residual:

  SRMR                                           0.188       0.188

Parameter Estimates:

  Standard errors                             Sandwich
  Information bread                           Observed
  Observed information based on                Hessian

Latent Variables:
                     Estimate   Std.Err  z-value      P(>|z|)   Std.lv   Std.all
  pri_con =~                                                                    
    PC01_01              1.000                                    1.599    0.927
    PC01_02              0.988    0.027       36.235    0.000     1.580    0.895
    PC01_04              0.969    0.027       35.679    0.000     1.550    0.883
    PC01_05              0.999    0.023       42.633    0.000     1.598    0.906
    PC01_06              0.851    0.038       22.657    0.000     1.361    0.793
    PC01_07              0.994    0.023       43.524    0.000     1.590    0.922
  grats_gen =~                                                                  
    GR02_01              1.000                                    1.145    0.853
    GR02_02              1.120    0.033       33.637    0.000     1.283    0.906
    GR02_03              0.992    0.046       21.721    0.000     1.136    0.849
    GR02_04              0.958    0.046       20.688    0.000     1.098    0.835
    GR02_05              1.064    0.039       27.133    0.000     1.219    0.851
  pri_delib =~                                                                  
    PD01_01              1.000                                    1.443    0.836
    PD01_02              0.679    0.051       13.436    0.000     0.979    0.638
    PD01_03              0.737    0.058       12.626    0.000     1.063    0.683
    PD01_04              0.873    0.049       17.786    0.000     1.260    0.739
    PD01_05              0.740    0.052       14.183    0.000     1.067    0.645
  self_eff =~                                                                   
    SE01_01              1.000                                    1.126    0.815
    SE01_02              0.805    0.060       13.435    0.000     0.907    0.671
    SE01_03              0.918    0.051       17.902    0.000     1.033    0.772
    SE01_04              0.939    0.044       21.178    0.000     1.058    0.785
  trust_community =~                                                            
    TR01_01              1.000                                    1.022    0.807
    TR01_02              0.819    0.054       15.111    0.000     0.837    0.762
    TR01_03              0.922    0.048       19.178    0.000     0.943    0.817
  trust_provider =~                                                             
    TR01_04              1.000                                    1.058    0.882
    TR01_05              0.854    0.040       21.372    0.000     0.903    0.780
    TR01_06              0.834    0.041       20.194    0.000     0.882    0.796
  trust_system =~                                                               
    TR01_07              1.000                                    0.810    0.688
    TR01_08              1.060    0.084       12.677    0.000     0.859    0.662
    TR01_09              1.348    0.084       16.058    0.000     1.093    0.813
  trust_spec =~                                                                 
    trust_communty       1.000                                    0.838    0.838
    trust_provider       1.185    0.083       14.342    0.000     0.961    0.961
    trust_system         1.009    0.081       12.402    0.000     1.067    1.067

Regressions:
                   Estimate   Std.Err  z-value      P(>|z|)   Std.lv   Std.all
  pri_con ~                                                                   
    like               0.032    0.163        0.197    0.844     0.020    0.010
    likedislik         0.177    0.169        1.043    0.297     0.111    0.051
  grats_gen ~                                                                 
    like              -0.167    0.120       -1.388    0.165    -0.146   -0.069
    likedislik        -0.221    0.122       -1.811    0.070    -0.193   -0.090
  pri_delib ~                                                                 
    like               0.002    0.159        0.011    0.992     0.001    0.001
    likedislik        -0.084    0.169       -0.499    0.617    -0.058   -0.027
  self_eff ~                                                                  
    like              -0.088    0.126       -0.693    0.488    -0.078   -0.037
    likedislik        -0.066    0.131       -0.501    0.616    -0.058   -0.027
  trust_spec ~                                                                
    like              -0.154    0.091       -1.702    0.089    -0.180   -0.086
    likedislik        -0.151    0.094       -1.594    0.111    -0.176   -0.082
  words ~                                                                     
    pri_con    (a)     0.620    6.851        0.090    0.928     0.991    0.004
    grats_gen  (b)    11.167   13.072        0.854    0.393    12.790    0.051
    pri_delib  (c)   -18.340    7.948       -2.307    0.021   -26.458   -0.106
    self_eff   (d)    52.682   14.678        3.589    0.000    59.340    0.238
    trust_spec (e)   -10.062   14.918       -0.674    0.500    -8.624   -0.035
    like       (f)   -29.288   26.879       -1.090    0.276   -29.288   -0.056
    likedislik (g)   -40.080   27.169       -1.475    0.140   -40.080   -0.075
    male              -8.863   23.351       -0.380    0.704    -8.863   -0.018
    age                0.339    0.558        0.608    0.543     0.339    0.021
    edu               12.181   15.227        0.800    0.424    12.181    0.041
  GR02_01 ~                                                                   
    male              -0.130    0.116       -1.129    0.259    -0.130   -0.049
    age                0.000    0.004        0.079    0.937     0.000    0.003
    edu                0.004    0.069        0.051    0.959     0.004    0.002
  GR02_02 ~                                                                   
    male              -0.072    0.120       -0.600    0.548    -0.072   -0.025
    age                0.006    0.004        1.529    0.126     0.006    0.067
    edu               -0.081    0.071       -1.142    0.253    -0.081   -0.048
  GR02_03 ~                                                                   
    male              -0.030    0.115       -0.263    0.793    -0.030   -0.011
    age                0.001    0.004        0.292    0.771     0.001    0.013
    edu               -0.082    0.067       -1.230    0.219    -0.082   -0.052
  GR02_04 ~                                                                   
    male               0.023    0.113        0.207    0.836     0.023    0.009
    age                0.005    0.004        1.287    0.198     0.005    0.056
    edu               -0.072    0.067       -1.063    0.288    -0.072   -0.046
  GR02_05 ~                                                                   
    male              -0.145    0.123       -1.172    0.241    -0.145   -0.051
    age               -0.004    0.004       -0.890    0.374    -0.004   -0.040
    edu                0.012    0.073        0.162    0.871     0.012    0.007
  PC01_01 ~                                                                   
    male              -0.179    0.151       -1.191    0.234    -0.179   -0.052
    age               -0.004    0.005       -0.795    0.427    -0.004   -0.035
    edu                0.119    0.087        1.361    0.173     0.119    0.058
  PC01_02 ~                                                                   
    male              -0.300    0.154       -1.950    0.051    -0.300   -0.085
    age               -0.008    0.005       -1.643    0.100    -0.008   -0.070
    edu                0.056    0.089        0.624    0.533     0.056    0.026
  PC01_04 ~                                                                   
    male              -0.223    0.152       -1.463    0.144    -0.223   -0.063
    age               -0.009    0.005       -1.959    0.050    -0.009   -0.084
    edu                0.121    0.089        1.370    0.171     0.121    0.058
  PC01_05 ~                                                                   
    male              -0.096    0.154       -0.622    0.534    -0.096   -0.027
    age               -0.006    0.005       -1.137    0.256    -0.006   -0.049
    edu                0.099    0.090        1.097    0.272     0.099    0.047
  PC01_06 ~                                                                   
    male              -0.106    0.150       -0.711    0.477    -0.106   -0.031
    age               -0.005    0.005       -1.034    0.301    -0.005   -0.045
    edu                0.050    0.087        0.582    0.560     0.050    0.025
  PC01_07 ~                                                                   
    male              -0.172    0.150       -1.146    0.252    -0.172   -0.050
    age               -0.006    0.005       -1.313    0.189    -0.006   -0.057
    edu                0.090    0.086        1.039    0.299     0.090    0.044
  TR01_01 ~                                                                   
    male              -0.298    0.108       -2.751    0.006    -0.298   -0.118
    age               -0.004    0.004       -1.091    0.275    -0.004   -0.048
    edu                0.004    0.061        0.066    0.947     0.004    0.003
  TR01_02 ~                                                                   
    male              -0.140    0.096       -1.467    0.142    -0.140   -0.064
    age               -0.002    0.003       -0.547    0.585    -0.002   -0.025
    edu                0.020    0.053        0.377    0.706     0.020    0.015
  TR01_03 ~                                                                   
    male              -0.134    0.100       -1.349    0.177    -0.134   -0.058
    age               -0.004    0.003       -1.190    0.234    -0.004   -0.054
    edu               -0.006    0.060       -0.107    0.915    -0.006   -0.005
  TR01_04 ~                                                                   
    male              -0.088    0.104       -0.847    0.397    -0.088   -0.037
    age                0.000    0.003        0.127    0.899     0.000    0.006
    edu               -0.053    0.058       -0.922    0.357    -0.053   -0.037
  TR01_05 ~                                                                   
    male              -0.044    0.099       -0.446    0.655    -0.044   -0.019
    age                0.001    0.003        0.371    0.711     0.001    0.016
    edu                0.013    0.058        0.233    0.816     0.013    0.010
  TR01_06 ~                                                                   
    male               0.046    0.096        0.479    0.632     0.046    0.021
    age               -0.004    0.003       -1.229    0.219    -0.004   -0.052
    edu                0.020    0.056        0.360    0.719     0.020    0.015
  TR01_07 ~                                                                   
    male               0.090    0.100        0.903    0.367     0.090    0.038
    age               -0.004    0.003       -1.156    0.248    -0.004   -0.049
    edu               -0.058    0.058       -1.003    0.316    -0.058   -0.042
  TR01_08 ~                                                                   
    male               0.027    0.112        0.237    0.812     0.027    0.010
    age                0.003    0.004        0.841    0.400     0.003    0.036
    edu               -0.096    0.065       -1.479    0.139    -0.096   -0.063
  TR01_09 ~                                                                   
    male              -0.122    0.115       -1.060    0.289    -0.122   -0.045
    age               -0.002    0.004       -0.389    0.697    -0.002   -0.018
    edu               -0.148    0.067       -2.201    0.028    -0.148   -0.093
  PD01_01 ~                                                                   
    male              -0.178    0.148       -1.208    0.227    -0.178   -0.052
    age               -0.015    0.005       -3.294    0.001    -0.015   -0.138
    edu               -0.030    0.085       -0.352    0.725    -0.030   -0.015
  PD01_02 ~                                                                   
    male              -0.119    0.131       -0.909    0.363    -0.119   -0.039
    age               -0.014    0.004       -3.453    0.001    -0.014   -0.142
    edu                0.028    0.077        0.362    0.718     0.028    0.015
  PD01_03 ~                                                                   
    male              -0.322    0.132       -2.435    0.015    -0.322   -0.104
    age               -0.004    0.004       -1.043    0.297    -0.004   -0.045
    edu                0.063    0.080        0.788    0.431     0.063    0.034
  PD01_04 ~                                                                   
    male              -0.413    0.144       -2.861    0.004    -0.413   -0.121
    age               -0.009    0.005       -1.927    0.054    -0.009   -0.083
    edu                0.100    0.085        1.183    0.237     0.100    0.050
  PD01_05 ~                                                                   
    male              -0.206    0.142       -1.453    0.146    -0.206   -0.062
    age               -0.012    0.004       -2.714    0.007    -0.012   -0.112
    edu               -0.003    0.084       -0.041    0.968    -0.003   -0.002
  SE01_01 ~                                                                   
    male               0.123    0.118        1.043    0.297     0.123    0.045
    age                0.000    0.004        0.052    0.959     0.000    0.002
    edu                0.210    0.068        3.097    0.002     0.210    0.128
  SE01_02 ~                                                                   
    male               0.062    0.112        0.557    0.577     0.062    0.023
    age               -0.013    0.004       -3.542    0.000    -0.013   -0.149
    edu                0.197    0.066        2.974    0.003     0.197    0.123
  SE01_03 ~                                                                   
    male               0.198    0.115        1.722    0.085     0.198    0.074
    age                0.001    0.004        0.289    0.773     0.001    0.013
    edu                0.142    0.067        2.110    0.035     0.142    0.089
  SE01_04 ~                                                                   
    male               0.057    0.115        0.495    0.621     0.057    0.021
    age                0.008    0.004        2.088    0.037     0.008    0.087
    edu                0.125    0.067        1.884    0.060     0.125    0.078

Covariances:
                   Estimate   Std.Err  z-value      P(>|z|)   Std.lv   Std.all
 .SE01_01 ~~                                                                  
   .SE01_02    (x)     0.110    0.045        2.455    0.014     0.110    0.147
 .SE01_03 ~~                                                                  
   .SE01_04    (x)     0.110    0.045        2.455    0.014     0.110    0.161
  like ~~                                                                     
    likedislik        -0.109    0.007      -16.483    0.000    -0.109   -0.494
    male               0.006    0.010        0.596    0.551     0.006    0.025
    age                0.364    0.313        1.161    0.246     0.364    0.049
    edu                0.016    0.017        0.916    0.360     0.016    0.039
  likedislike ~~                                                              
    male              -0.009    0.010       -0.942    0.346    -0.009   -0.040
    age               -0.321    0.306       -1.047    0.295    -0.321   -0.044
    edu               -0.019    0.017       -1.169    0.243    -0.019   -0.050
  male ~~                                                                     
    age                0.758    0.328        2.312    0.021     0.758    0.097
    edu                0.052    0.018        2.964    0.003     0.052    0.125
  age ~~                                                                      
    edu               -1.018    0.547       -1.862    0.063    -1.018   -0.078

Intercepts:
                   Estimate   Std.Err  z-value      P(>|z|)   Std.lv   Std.all
   .PC01_01            3.275    0.298       11.003    0.000     3.275    1.899
   .PC01_02            3.676    0.307       11.968    0.000     3.676    2.080
   .PC01_04            3.480    0.300       11.601    0.000     3.480    1.982
   .PC01_05            3.320    0.306       10.833    0.000     3.320    1.882
   .PC01_06            3.134    0.289       10.855    0.000     3.134    1.828
   .PC01_07            3.367    0.299       11.276    0.000     3.367    1.952
   .GR02_01            4.453    0.231       19.294    0.000     4.453    3.318
   .GR02_02            4.640    0.253       18.336    0.000     4.640    3.277
   .GR02_03            5.373    0.230       23.328    0.000     5.373    4.014
   .GR02_04            5.113    0.230       22.193    0.000     5.113    3.890
   .GR02_05            5.046    0.261       19.336    0.000     5.046    3.522
   .PD01_01            4.532    0.308       14.699    0.000     4.532    2.627
   .PD01_02            4.025    0.260       15.495    0.000     4.025    2.622
   .PD01_03            4.459    0.286       15.583    0.000     4.459    2.865
   .PD01_04            4.539    0.312       14.529    0.000     4.539    2.662
   .PD01_05            5.026    0.291       17.286    0.000     5.026    3.037
   .SE01_01            4.871    0.265       18.408    0.000     4.871    3.526
   .SE01_02            5.769    0.245       23.522    0.000     5.769    4.271
   .SE01_03            4.863    0.241       20.181    0.000     4.863    3.631
   .SE01_04            4.578    0.250       18.348    0.000     4.578    3.395
   .TR01_01            5.185    0.228       22.713    0.000     5.185    4.094
   .TR01_02            5.039    0.194       25.967    0.000     5.039    4.584
   .TR01_03            4.971    0.206       24.072    0.000     4.971    4.308
   .TR01_04            5.644    0.213       26.486    0.000     5.644    4.705
   .TR01_05            5.245    0.204       25.677    0.000     5.245    4.529
   .TR01_06            5.340    0.196       27.213    0.000     5.340    4.818
   .TR01_07            6.062    0.201       30.123    0.000     6.062    5.150
   .TR01_08            4.968    0.229       21.729    0.000     4.968    3.832
   .TR01_09            5.721    0.244       23.490    0.000     5.721    4.256
   .words             68.057   30.401        2.239    0.025    68.057    0.273
    like               0.347    0.020       17.237    0.000     0.347    0.729
    likedislike        0.315    0.020       16.027    0.000     0.315    0.678
    male               0.493    0.021       23.297    0.000     0.493    0.986
    age               46.132    0.658       70.151    0.000    46.132    2.967
    edu                1.852    0.036       51.966    0.000     1.852    2.198

Variances:
                   Estimate   Std.Err  z-value      P(>|z|)   Std.lv   Std.all
   .PC01_01            0.395    0.050        7.916    0.000     0.395    0.133
   .PC01_02            0.580    0.103        5.618    0.000     0.580    0.186
   .PC01_04            0.632    0.078        8.072    0.000     0.632    0.205
   .PC01_05            0.539    0.065        8.277    0.000     0.539    0.173
   .PC01_06            1.078    0.117        9.225    0.000     1.078    0.366
   .PC01_07            0.424    0.064        6.634    0.000     0.424    0.142
   .GR02_01            0.486    0.052        9.401    0.000     0.486    0.270
   .GR02_02            0.344    0.036        9.592    0.000     0.344    0.172
   .GR02_03            0.496    0.074        6.721    0.000     0.496    0.277
   .GR02_04            0.513    0.051       10.118    0.000     0.513    0.297
   .GR02_05            0.558    0.061        9.156    0.000     0.558    0.272
   .PD01_01            0.828    0.129        6.408    0.000     0.828    0.278
   .PD01_02            1.344    0.130       10.342    0.000     1.344    0.570
   .PD01_03            1.259    0.131        9.606    0.000     1.259    0.520
   .PD01_04            1.248    0.147        8.468    0.000     1.248    0.429
   .PD01_05            1.553    0.129       12.052    0.000     1.553    0.567
   .SE01_01            0.602    0.093        6.452    0.000     0.602    0.315
   .SE01_02            0.927    0.121        7.672    0.000     0.927    0.508
   .SE01_03            0.698    0.105        6.678    0.000     0.698    0.389
   .SE01_04            0.673    0.080        8.429    0.000     0.673    0.370
   .TR01_01            0.532    0.070        7.549    0.000     0.532    0.331
   .TR01_02            0.501    0.057        8.724    0.000     0.501    0.415
   .TR01_03            0.434    0.045        9.648    0.000     0.434    0.326
   .TR01_04            0.316    0.034        9.218    0.000     0.316    0.220
   .TR01_05            0.524    0.055        9.610    0.000     0.524    0.391
   .TR01_06            0.447    0.043       10.458    0.000     0.447    0.364
   .TR01_07            0.723    0.056       12.840    0.000     0.723    0.522
   .TR01_08            0.934    0.079       11.819    0.000     0.934    0.556
   .TR01_09            0.592    0.067        8.808    0.000     0.592    0.328
   .words          57409.646    0.006 10017820.230    0.000 57409.646    0.921
   .pri_con            2.551    0.144       17.751    0.000     0.998    0.998
   .grats_gen          1.303    0.113       11.534    0.000     0.993    0.993
   .pri_delib          2.080    0.167       12.433    0.000     0.999    0.999
   .self_eff           1.267    0.115       10.996    0.000     0.999    0.999
   .trust_communty     0.311    0.048        6.500    0.000     0.297    0.297
   .trust_provider     0.087    0.038        2.255    0.024     0.077    0.077
   .trust_system      -0.091    0.026       -3.509    0.000    -0.138   -0.138
   .trust_spec         0.729    0.093        7.834    0.000     0.993    0.993
    like               0.227    0.006       36.792    0.000     0.227    1.000
    likedislike        0.216    0.007       29.655    0.000     0.216    1.000
    male               0.250    0.000      817.431    0.000     0.250    1.000
    age              241.746    9.947       24.303    0.000   241.746    1.000
    edu                0.710    0.021       34.591    0.000     0.710    1.000
```


:::
:::



No significant effects of popularity cues on privacy calculus.

**Like-Dislike & Control vs. Like**



::: {.cell}

```{.r .cell-code}
model <- "
  pri_con =~ PC01_01 + PC01_02 + PC01_04 + PC01_05 + PC01_06 + PC01_07
  grats_gen =~ GR02_01 + GR02_02 + GR02_03 + GR02_04 + GR02_05
  pri_delib =~ PD01_01 + PD01_02 + PD01_03 + PD01_04 + PD01_05
  self_eff =~ SE01_01 + SE01_02 + SE01_03 + SE01_04
  SE01_01 ~~ x*SE01_02
  SE01_03 ~~ x*SE01_04
  trust_community =~ TR01_01 + TR01_02 + TR01_03
  trust_provider =~ TR01_04 + TR01_05 + TR01_06
  trust_system =~ TR01_07 + TR01_08 + TR01_09
  trust_spec =~ trust_community + trust_provider + trust_system

  pri_con + grats_gen + pri_delib + self_eff + trust_spec + words ~ likedislike + control

# Covariates
words + GR02_01 + GR02_02 + GR02_03 + GR02_04 + GR02_05 + PC01_01 + PC01_02 + PC01_04 + PC01_05 + PC01_06 + PC01_07 + TR01_01 + TR01_02 + TR01_03 + TR01_04 + TR01_05 + TR01_06 + TR01_07 + TR01_08 + TR01_09 + PD01_01 + PD01_02 + PD01_03 + PD01_04 + PD01_05 + SE01_01 + SE01_02 + SE01_03 + SE01_04 ~ male + age + edu
"
fit_lik_ctrl <- 
  lavaan::sem(
    model = model, 
    data = d, 
    estimator = "MLR", 
    missing = "ML")

summary(
  fit_lik_ctrl, 
  fit = TRUE, 
  std = TRUE
  )
```

::: {.cell-output .cell-output-stdout}

```
lavaan 0.6-20 ended normally after 557 iterations

  Estimator                                         ML
  Optimization method                           NLMINB
  Number of model parameters                       211
  Number of equality constraints                     1

                                                  Used       Total
  Number of observations                           558         559
  Number of missing patterns                         3            

Model Test User Model:
                                              Standard      Scaled
  Test Statistic                              1254.996     923.368
  Degrees of freedom                               435         435
  P-value (Chi-square)                           0.000       0.000
  Scaling correction factor                                  1.359
    Yuan-Bentler correction (Mplus variant)                       

Model Test Baseline Model:

  Test statistic                             13325.570   10323.740
  Degrees of freedom                               585         585
  P-value                                        0.000       0.000
  Scaling correction factor                                  1.291

User Model versus Baseline Model:

  Comparative Fit Index (CFI)                    0.936       0.950
  Tucker-Lewis Index (TLI)                       0.913       0.933
                                                                  
  Robust Comparative Fit Index (CFI)                         0.948
  Robust Tucker-Lewis Index (TLI)                            0.930

Loglikelihood and Information Criteria:

  Loglikelihood user model (H0)             -26475.947  -26475.947
  Scaling correction factor                                  1.246
      for the MLR correction                                      
  Loglikelihood unrestricted model (H1)     -25848.449  -25848.449
  Scaling correction factor                                  1.324
      for the MLR correction                                      
                                                                  
  Akaike (AIC)                               53371.893   53371.893
  Bayesian (BIC)                             54280.008   54280.008
  Sample-size adjusted Bayesian (SABIC)      53613.368   53613.368

Root Mean Square Error of Approximation:

  RMSEA                                          0.058       0.045
  90 Percent confidence interval - lower         0.054       0.041
  90 Percent confidence interval - upper         0.062       0.048
  P-value H_0: RMSEA <= 0.050                    0.000       0.993
  P-value H_0: RMSEA >= 0.080                    0.000       0.000
                                                                  
  Robust RMSEA                                               0.052
  90 Percent confidence interval - lower                     0.047
  90 Percent confidence interval - upper                     0.057
  P-value H_0: Robust RMSEA <= 0.050                         0.250
  P-value H_0: Robust RMSEA >= 0.080                         0.000

Standardized Root Mean Square Residual:

  SRMR                                           0.046       0.046

Parameter Estimates:

  Standard errors                             Sandwich
  Information bread                           Observed
  Observed information based on                Hessian

Latent Variables:
                     Estimate   Std.Err  z-value     P(>|z|)   Std.lv   Std.all
  pri_con =~                                                                   
    PC01_01              1.000                                   1.595    0.926
    PC01_02              0.990    0.027      36.192    0.000     1.579    0.894
    PC01_04              0.972    0.027      35.662    0.000     1.551    0.884
    PC01_05              1.002    0.024      42.459    0.000     1.599    0.907
    PC01_06              0.855    0.038      22.667    0.000     1.363    0.795
    PC01_07              0.995    0.023      43.750    0.000     1.587    0.920
  grats_gen =~                                                                 
    GR02_01              1.000                                   1.130    0.841
    GR02_02              1.120    0.033      33.549    0.000     1.265    0.893
    GR02_03              1.025    0.048      21.463    0.000     1.158    0.865
    GR02_04              0.988    0.048      20.449    0.000     1.116    0.849
    GR02_05              1.075    0.040      27.055    0.000     1.215    0.847
  pri_delib =~                                                                 
    PD01_01              1.000                                   1.478    0.856
    PD01_02              0.665    0.048      13.894    0.000     0.983    0.640
    PD01_03              0.703    0.055      12.854    0.000     1.038    0.667
    PD01_04              0.841    0.047      17.753    0.000     1.243    0.729
    PD01_05              0.714    0.050      14.264    0.000     1.056    0.637
  self_eff =~                                                                  
    SE01_01              1.000                                   1.111    0.805
    SE01_02              0.816    0.058      14.127    0.000     0.906    0.672
    SE01_03              0.936    0.046      20.258    0.000     1.040    0.776
    SE01_04              0.960    0.044      21.890    0.000     1.066    0.790
  trust_community =~                                                           
    TR01_01              1.000                                   1.028    0.811
    TR01_02              0.811    0.053      15.321    0.000     0.834    0.759
    TR01_03              0.914    0.047      19.282    0.000     0.940    0.815
  trust_provider =~                                                            
    TR01_04              1.000                                   1.059    0.882
    TR01_05              0.854    0.040      21.557    0.000     0.904    0.782
    TR01_06              0.830    0.040      20.584    0.000     0.879    0.794
  trust_system =~                                                              
    TR01_07              1.000                                   0.799    0.679
    TR01_08              1.060    0.082      12.941    0.000     0.846    0.653
    TR01_09              1.403    0.087      16.210    0.000     1.120    0.833
  trust_spec =~                                                                
    trust_communty       1.000                                   0.847    0.847
    trust_provider       1.160    0.080      14.542    0.000     0.954    0.954
    trust_system         0.970    0.079      12.353    0.000     1.058    1.058

Regressions:
                   Estimate   Std.Err  z-value     P(>|z|)   Std.lv   Std.all
  pri_con ~                                                                  
    likedislike        0.158    0.174       0.904    0.366     0.099    0.046
    control           -0.031    0.162      -0.191    0.849    -0.019   -0.009
  grats_gen ~                                                                
    likedislike       -0.045    0.125      -0.361    0.718    -0.040   -0.018
    control            0.164    0.119       1.386    0.166     0.145    0.069
  pri_delib ~                                                                
    likedislike       -0.083    0.166      -0.501    0.617    -0.056   -0.026
    control            0.004    0.162       0.028    0.978     0.003    0.001
  self_eff ~                                                                 
    likedislike        0.007    0.127       0.051    0.959     0.006    0.003
    control            0.083    0.125       0.666    0.506     0.075    0.035
  trust_spec ~                                                               
    likedislike       -0.002    0.095      -0.016    0.987    -0.002   -0.001
    control            0.159    0.092       1.728    0.084     0.183    0.087
  words ~                                                                    
    likedislike       -9.593   19.548      -0.491    0.624    -9.593   -0.018
    control           32.972   28.926       1.140    0.254    32.972    0.062
    male              -8.817   23.401      -0.377    0.706    -8.817   -0.018
    age                0.336    0.558       0.603    0.547     0.336    0.021
    edu               11.696   15.254       0.767    0.443    11.696    0.039
  GR02_01 ~                                                                  
    male              -0.130    0.116      -1.124    0.261    -0.130   -0.048
    age                0.000    0.004       0.083    0.934     0.000    0.004
    edu                0.003    0.069       0.048    0.962     0.003    0.002
  GR02_02 ~                                                                  
    male              -0.071    0.120      -0.592    0.554    -0.071   -0.025
    age                0.006    0.004       1.535    0.125     0.006    0.068
    edu               -0.082    0.071      -1.156    0.248    -0.082   -0.049
  GR02_03 ~                                                                  
    male              -0.029    0.116      -0.250    0.802    -0.029   -0.011
    age                0.001    0.004       0.301    0.763     0.001    0.013
    edu               -0.084    0.067      -1.264    0.206    -0.084   -0.053
  GR02_04 ~                                                                  
    male               0.025    0.113       0.219    0.827     0.025    0.009
    age                0.005    0.004       1.297    0.195     0.005    0.057
    edu               -0.074    0.067      -1.096    0.273    -0.074   -0.047
  GR02_05 ~                                                                  
    male              -0.144    0.124      -1.165    0.244    -0.144   -0.050
    age               -0.004    0.004      -0.883    0.377    -0.004   -0.039
    edu                0.011    0.073       0.147    0.883     0.011    0.006
  PC01_01 ~                                                                  
    male              -0.177    0.150      -1.175    0.240    -0.177   -0.051
    age               -0.004    0.005      -0.781    0.435    -0.004   -0.034
    edu                0.114    0.087       1.308    0.191     0.114    0.056
  PC01_02 ~                                                                  
    male              -0.297    0.153      -1.937    0.053    -0.297   -0.084
    age               -0.008    0.005      -1.629    0.103    -0.008   -0.070
    edu                0.051    0.089       0.571    0.568     0.051    0.024
  PC01_04 ~                                                                  
    male              -0.220    0.152      -1.448    0.148    -0.220   -0.063
    age               -0.009    0.005      -1.946    0.052    -0.009   -0.084
    edu                0.117    0.089       1.320    0.187     0.117    0.056
  PC01_05 ~                                                                  
    male              -0.093    0.154      -0.605    0.545    -0.093   -0.026
    age               -0.006    0.005      -1.123    0.261    -0.006   -0.049
    edu                0.094    0.090       1.046    0.295     0.094    0.045
  PC01_06 ~                                                                  
    male              -0.104    0.150      -0.695    0.487    -0.104   -0.030
    age               -0.005    0.005      -1.022    0.307    -0.005   -0.044
    edu                0.046    0.087       0.535    0.593     0.046    0.023
  PC01_07 ~                                                                  
    male              -0.169    0.150      -1.130    0.259    -0.169   -0.049
    age               -0.006    0.005      -1.299    0.194    -0.006   -0.056
    edu                0.085    0.086       0.987    0.324     0.085    0.042
  TR01_01 ~                                                                  
    male              -0.299    0.109      -2.757    0.006    -0.299   -0.118
    age               -0.004    0.004      -1.095    0.274    -0.004   -0.048
    edu                0.005    0.061       0.076    0.939     0.005    0.003
  TR01_02 ~                                                                  
    male              -0.142    0.095      -1.489    0.137    -0.142   -0.065
    age               -0.002    0.003      -0.557    0.577    -0.002   -0.025
    edu                0.023    0.053       0.426    0.670     0.023    0.017
  TR01_03 ~                                                                  
    male              -0.136    0.099      -1.374    0.169    -0.136   -0.059
    age               -0.004    0.003      -1.202    0.229    -0.004   -0.054
    edu               -0.003    0.060      -0.055    0.956    -0.003   -0.002
  TR01_04 ~                                                                  
    male              -0.089    0.104      -0.857    0.392    -0.089   -0.037
    age                0.000    0.003       0.120    0.904     0.000    0.005
    edu               -0.052    0.058      -0.899    0.369    -0.052   -0.037
  TR01_05 ~                                                                  
    male              -0.047    0.099      -0.473    0.636    -0.047   -0.020
    age                0.001    0.003       0.356    0.722     0.001    0.015
    edu                0.017    0.058       0.301    0.763     0.017    0.013
  TR01_06 ~                                                                  
    male               0.044    0.095       0.457    0.648     0.044    0.020
    age               -0.004    0.003      -1.247    0.212    -0.004   -0.053
    edu                0.024    0.056       0.434    0.664     0.024    0.018
  TR01_07 ~                                                                  
    male               0.090    0.100       0.894    0.371     0.090    0.038
    age               -0.004    0.003      -1.165    0.244    -0.004   -0.050
    edu               -0.056    0.058      -0.963    0.336    -0.056   -0.040
  TR01_08 ~                                                                  
    male               0.025    0.112       0.224    0.822     0.025    0.010
    age                0.003    0.004       0.832    0.405     0.003    0.036
    edu               -0.094    0.065      -1.442    0.149    -0.094   -0.061
  TR01_09 ~                                                                  
    male              -0.124    0.115      -1.071    0.284    -0.124   -0.046
    age               -0.002    0.004      -0.396    0.692    -0.002   -0.018
    edu               -0.147    0.067      -2.176    0.030    -0.147   -0.092
  PD01_01 ~                                                                  
    male              -0.179    0.148      -1.210    0.226    -0.179   -0.052
    age               -0.015    0.005      -3.294    0.001    -0.015   -0.138
    edu               -0.028    0.085      -0.335    0.738    -0.028   -0.014
  PD01_02 ~                                                                  
    male              -0.121    0.132      -0.915    0.360    -0.121   -0.039
    age               -0.014    0.004      -3.455    0.001    -0.014   -0.142
    edu                0.030    0.077       0.387    0.699     0.030    0.016
  PD01_03 ~                                                                  
    male              -0.323    0.133      -2.434    0.015    -0.323   -0.104
    age               -0.004    0.004      -1.041    0.298    -0.004   -0.045
    edu                0.063    0.080       0.790    0.430     0.063    0.034
  PD01_04 ~                                                                  
    male              -0.414    0.145      -2.861    0.004    -0.414   -0.121
    age               -0.009    0.005      -1.925    0.054    -0.009   -0.083
    edu                0.101    0.085       1.189    0.234     0.101    0.050
  PD01_05 ~                                                                  
    male              -0.206    0.142      -1.450    0.147    -0.206   -0.062
    age               -0.012    0.004      -2.710    0.007    -0.012   -0.112
    edu               -0.004    0.084      -0.043    0.966    -0.004   -0.002
  SE01_01 ~                                                                  
    male               0.118    0.118       1.003    0.316     0.118    0.043
    age                0.000    0.004       0.014    0.989     0.000    0.001
    edu                0.211    0.068       3.104    0.002     0.211    0.129
  SE01_02 ~                                                                  
    male               0.058    0.112       0.520    0.603     0.058    0.022
    age               -0.013    0.004      -3.578    0.000    -0.013   -0.151
    edu                0.198    0.066       2.986    0.003     0.198    0.124
  SE01_03 ~                                                                  
    male               0.193    0.115       1.686    0.092     0.193    0.072
    age                0.001    0.004       0.251    0.801     0.001    0.011
    edu                0.142    0.067       2.120    0.034     0.142    0.090
  SE01_04 ~                                                                  
    male               0.052    0.115       0.453    0.651     0.052    0.019
    age                0.007    0.004       2.050    0.040     0.007    0.086
    edu                0.126    0.067       1.897    0.058     0.126    0.079

Covariances:
                   Estimate   Std.Err  z-value     P(>|z|)   Std.lv   Std.all
 .SE01_01 ~~                                                                 
   .SE01_02    (x)     0.105    0.044       2.407    0.016     0.105    0.138
 .SE01_03 ~~                                                                 
   .SE01_04    (x)     0.105    0.044       2.407    0.016     0.105    0.156
 .pri_con ~~                                                                 
   .grats_gen         -0.277    0.096      -2.902    0.004    -0.155   -0.155
   .pri_delib          1.331    0.130      10.213    0.000     0.566    0.566
   .self_eff          -0.373    0.092      -4.055    0.000    -0.211   -0.211
   .trust_spec        -0.405    0.070      -5.824    0.000    -0.293   -0.293
   .words            -35.118   18.522      -1.896    0.058   -22.044   -0.088
 .grats_gen ~~                                                               
   .pri_delib         -0.073    0.102      -0.708    0.479    -0.044   -0.044
   .self_eff           0.465    0.067       6.914    0.000     0.372    0.372
   .trust_spec         0.744    0.080       9.326    0.000     0.761    0.761
   .words             29.451   11.287       2.609    0.009    26.151    0.105
 .pri_delib ~~                                                               
   .self_eff          -0.326    0.095      -3.416    0.001    -0.199   -0.199
   .trust_spec        -0.145    0.080      -1.819    0.069    -0.113   -0.113
   .words            -52.100   15.972      -3.262    0.001   -35.266   -0.141
 .self_eff ~~                                                                
   .trust_spec         0.526    0.060       8.774    0.000     0.546    0.546
   .words             70.815   17.043       4.155    0.000    63.779    0.256
 .trust_spec ~~                                                              
   .words             26.191    6.881       3.806    0.000    30.181    0.121

Intercepts:
                   Estimate   Std.Err  z-value     P(>|z|)   Std.lv   Std.all
   .PC01_01            3.310    0.309      10.718    0.000     3.310    1.921
   .PC01_02            3.711    0.322      11.523    0.000     3.711    2.102
   .PC01_04            3.514    0.316      11.132    0.000     3.514    2.003
   .PC01_05            3.355    0.326      10.306    0.000     3.355    1.904
   .PC01_06            3.165    0.302      10.477    0.000     3.165    1.846
   .PC01_07            3.402    0.313      10.860    0.000     3.402    1.974
   .GR02_01            4.284    0.244      17.527    0.000     4.284    3.189
   .GR02_02            4.453    0.265      16.820    0.000     4.453    3.142
   .GR02_03            5.208    0.242      21.495    0.000     5.208    3.890
   .GR02_04            4.954    0.238      20.836    0.000     4.954    3.768
   .GR02_05            4.867    0.270      17.994    0.000     4.867    3.395
   .PD01_01            4.528    0.306      14.787    0.000     4.528    2.622
   .PD01_02            4.020    0.257      15.611    0.000     4.020    2.616
   .PD01_03            4.456    0.276      16.136    0.000     4.456    2.860
   .PD01_04            4.536    0.302      14.998    0.000     4.536    2.657
   .PD01_05            5.025    0.284      17.674    0.000     5.025    3.033
   .SE01_01            4.793    0.270      17.755    0.000     4.793    3.474
   .SE01_02            5.705    0.250      22.800    0.000     5.705    4.231
   .SE01_03            4.790    0.242      19.806    0.000     4.790    3.575
   .SE01_04            4.503    0.248      18.129    0.000     4.503    3.339
   .TR01_01            5.030    0.227      22.161    0.000     5.030    3.968
   .TR01_02            4.909    0.196      25.023    0.000     4.909    4.466
   .TR01_03            4.823    0.207      23.259    0.000     4.823    4.181
   .TR01_04            5.460    0.217      25.199    0.000     5.460    4.549
   .TR01_05            5.082    0.205      24.848    0.000     5.082    4.395
   .TR01_06            5.181    0.197      26.315    0.000     5.181    4.682
   .TR01_07            5.904    0.201      29.369    0.000     5.904    5.016
   .TR01_08            4.801    0.229      20.969    0.000     4.801    3.703
   .TR01_09            5.509    0.249      22.145    0.000     5.509    4.096
   .words             35.618   39.123       0.910    0.363    35.618    0.142

Variances:
                   Estimate   Std.Err  z-value     P(>|z|)   Std.lv   Std.all
   .PC01_01            0.404    0.050       8.070    0.000     0.404    0.136
   .PC01_02            0.581    0.102       5.666    0.000     0.581    0.186
   .PC01_04            0.627    0.077       8.133    0.000     0.627    0.204
   .PC01_05            0.534    0.064       8.368    0.000     0.534    0.172
   .PC01_06            1.069    0.116       9.242    0.000     1.069    0.364
   .PC01_07            0.431    0.065       6.580    0.000     0.431    0.145
   .GR02_01            0.524    0.054       9.788    0.000     0.524    0.290
   .GR02_02            0.391    0.040       9.826    0.000     0.391    0.195
   .GR02_03            0.445    0.073       6.069    0.000     0.445    0.248
   .GR02_04            0.472    0.048       9.876    0.000     0.472    0.273
   .GR02_05            0.570    0.062       9.191    0.000     0.570    0.277
   .PD01_01            0.730    0.111       6.594    0.000     0.730    0.245
   .PD01_02            1.339    0.127      10.505    0.000     1.339    0.567
   .PD01_03            1.315    0.129      10.230    0.000     1.315    0.542
   .PD01_04            1.295    0.147       8.829    0.000     1.295    0.444
   .PD01_05            1.582    0.128      12.342    0.000     1.582    0.577
   .SE01_01            0.631    0.088       7.138    0.000     0.631    0.332
   .SE01_02            0.921    0.120       7.680    0.000     0.921    0.507
   .SE01_03            0.687    0.096       7.164    0.000     0.687    0.382
   .SE01_04            0.658    0.078       8.461    0.000     0.658    0.362
   .TR01_01            0.522    0.067       7.745    0.000     0.522    0.325
   .TR01_02            0.506    0.057       8.908    0.000     0.506    0.419
   .TR01_03            0.438    0.045       9.650    0.000     0.438    0.329
   .TR01_04            0.316    0.034       9.416    0.000     0.316    0.220
   .TR01_05            0.519    0.053       9.733    0.000     0.519    0.388
   .TR01_06            0.448    0.042      10.621    0.000     0.448    0.366
   .TR01_07            0.741    0.057      13.002    0.000     0.741    0.535
   .TR01_08            0.956    0.079      12.044    0.000     0.956    0.569
   .TR01_09            0.533    0.061       8.675    0.000     0.533    0.294
   .words          62215.808    0.042 1467939.554    0.000 62215.808    0.993
   .pri_con            2.538    0.144      17.580    0.000     0.997    0.997
   .grats_gen          1.268    0.113      11.182    0.000     0.994    0.994
   .pri_delib          2.183    0.157      13.864    0.000     0.999    0.999
   .self_eff           1.233    0.117      10.546    0.000     0.999    0.999
   .trust_communty     0.298    0.045       6.656    0.000     0.282    0.282
   .trust_provider     0.100    0.031       3.198    0.001     0.089    0.089
   .trust_system      -0.077    0.019      -4.050    0.000    -0.120   -0.120
   .trust_spec         0.753    0.095       7.912    0.000     0.992    0.992
```


:::
:::



No significant effects of popularity cues on privacy calculus.

## Exploratory analyses

The bivariate relations showed several strong correlations. 
People who trusted the providers more also experienced more gratifications (_r_ = .79).
In multiple regression, we analze the potential effect of one variable on the outcome while holding all others constant. 
In this case, I tested whether trust increases communication while holding constant gratifications, privacy concerns, privacy deliberations, and self-efficacy. 
Given the close theoretical and empirical interrelations between the predictors, this creates an unlikely and artifical scenario. 
If people say loose trust in a provider, they will likely experience more concerns and less benefits. 
But it is not even necessary to control for these mediating factors. 
When trying to analyze the causal effect, it is necessary to control for counfounding variables, but importantly not for mediators [@rohrerThinkingClearlyCorrelations2018].
Confunding variables impact both the independent variable and the outcome. 
Sociodemographic variables are often ideal candidates for confounders.
For example, men are generally less concerned about their privacy and post more online.
Hence, controlling for gender hence helps isolate the actual causal effect. 
As a result, I reran the analyses controlling for confounding variables that are not mediators.

### Experimental Factors
#### Words

I first compare likes and likes plus dislikes to the control condition.



::: {.cell}

```{.r .cell-code}
run_hurdles(
  object = d,
  name = "wrds_ctrl",
  outcome = "words", 
  predictor = "version",
  y_lab = "Words",
  x_lab = "Experimental conditions"
  )
```

::: {.cell-output .cell-output-stdout}

```
Some PIT values larger than 1! Largest:  1 
Rounding PIT > 1 to 1.
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/words-1.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: words ~ 1 + age + male + edu + version 
         hu ~ 1 + age + male + edu + version
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
                       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept                  5.61      0.37     4.85     6.36 1.00     2843     1521
hu_Intercept               1.05      0.38     0.31     1.80 1.00     3436     1864
age                       -0.00      0.01    -0.01     0.01 1.00     4656     1359
male                      -0.18      0.16    -0.49     0.14 1.00     2885     1635
edu                        0.01      0.10    -0.18     0.21 1.00     3409     1490
versionLike               -0.39      0.19    -0.76    -0.01 1.00     2723     1527
versionLike&dislike       -0.46      0.20    -0.88    -0.07 1.00     2265     1472
hu_age                    -0.01      0.01    -0.02     0.00 1.00     4426     1387
hu_male                   -0.06      0.17    -0.40     0.29 1.00     3578     1839
hu_edu                    -0.20      0.11    -0.41     0.00 1.00     2971     1652
hu_versionLike             0.01      0.21    -0.39     0.42 1.00     2614     1517
hu_versionLike&dislike     0.24      0.22    -0.21     0.68 1.00     2760     1571

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape     0.70      0.06     0.59     0.81 1.00     3850     1695

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
  Outcome    term                 contrast estimate conf.low conf.high
1   words     age                    dY/dX    0.258   -0.772      1.32
2   words     edu                    dY/dX   10.224   -7.466     28.19
3   words    male                    1 - 0  -11.053  -42.378     17.85
4   words version Like & dislike - Control  -47.041  -90.329     -9.05
5   words version           Like - Control  -34.671  -77.846      3.20
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/words-2.png){width=768}
:::
:::



Let's find out by how much people in like and like plus dislike condition wrote fewer words.



::: {.cell}

```{.r .cell-code}
words_m_ctrl <- 
  fig_wrds_ctrl$data %>% 
  filter(effect1__ == "Control") %>% 
  select(estimate__)

words_m_lk <- 
  fig_wrds_ctrl$data %>% 
  filter(effect1__ == "Like") %>% 
  select(estimate__)

words_m_lkdslk <- 
  fig_wrds_ctrl$data %>% 
  filter(effect1__ == "Like & dislike") %>% 
  select(estimate__)

effect_lk_wrds <- words_m_ctrl - words_m_lk
effect_lkdslk_wrds <- words_m_ctrl - words_m_lkdslk

effect_lk_prcnt <- 
  (effect_lk_wrds / words_m_ctrl) * 100 %>% 
  round()
effect_lkdlk_prcnt <- 
  (effect_lkdslk_wrds / words_m_ctrl) * 100 %>% 
  round()
```
:::



Compared to the control condition, participants wrote 44.93 percent fewer words on the platform with likes and dislikes.

I then compare likes and likes plus dislikes to the control condition.



::: {.cell}

```{.r .cell-code}
run_hurdles(
  object = d, 
  name = "wrds_lkdslk",
  outcome = "words", 
  predictor = "version_lkdslk",
  y_lab = "Words",
  x_lab = "Experimental conditions"
  )
```

::: {.cell-output .cell-output-stdout}

```
Some PIT values larger than 1! Largest:  1 
Rounding PIT > 1 to 1.
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/words-lkdslk-1.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: words ~ 1 + age + male + edu + version_lkdslk 
         hu ~ 1 + age + male + edu + version_lkdslk
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
                         Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept                    5.16      0.36     4.46     5.90 1.00     3617     1477
hu_Intercept                 1.29      0.38     0.56     2.03 1.00     4549     1349
age                         -0.00      0.01    -0.01     0.01 1.01     5123     1674
male                        -0.17      0.17    -0.51     0.15 1.00     3900     1665
edu                          0.01      0.10    -0.18     0.21 1.00     2884     1342
version_lkdslkLike           0.07      0.20    -0.33     0.47 1.00     2100     1680
version_lkdslkControl        0.46      0.20     0.06     0.86 1.00     2674     1742
hu_age                      -0.01      0.01    -0.02     0.00 1.01     5089     1501
hu_male                     -0.06      0.17    -0.41     0.27 1.00     4213     1110
hu_edu                      -0.21      0.10    -0.41    -0.00 1.00     4717     1538
hu_version_lkdslkLike       -0.21      0.22    -0.63     0.23 1.00     2739     1644
hu_version_lkdslkControl    -0.23      0.23    -0.69     0.20 1.00     2694     1475

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape     0.70      0.06     0.59     0.81 1.00     3822     1266

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
  Outcome           term                 contrast estimate conf.low conf.high
1   words            age                    dY/dX    0.243   -0.785      1.26
2   words            edu                    dY/dX   10.246   -7.131     29.19
3   words           male                    1 - 0  -10.490  -42.798     16.65
4   words version_lkdslk Control - Like & dislike   47.062    8.369     91.03
5   words version_lkdslk    Like - Like & dislike   12.601  -17.879     45.08
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/words-lkdslk-2.png){width=768}
:::
:::



#### Privacy Concerns

I first compare likes and likes plus dislikes to the control condition.



::: {.cell}

```{.r .cell-code}
run_hurdles(
  object = d,
  name = "pricon_ctrl",
  outcome = "pri_con_fs", 
  predictor = "version",
  y_lab = "Privacy concerns",
  x_lab = "Experimental conditions"
  )
```

::: {.cell-output-display}
![](05a_analyses_files/figure-html/pricon-ctrl-1.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: pri_con_fs ~ 1 + age + male + edu + version 
         hu ~ 1 + age + male + edu + version
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
                       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept                  1.21      0.09     1.04     1.39 1.00     1793     1570
hu_Intercept              -9.94      6.29   -24.45     0.58 1.00     1109      804
age                       -0.00      0.00    -0.00     0.00 1.01     1850     1323
male                      -0.05      0.04    -0.14     0.03 1.00     1827     1478
edu                        0.03      0.03    -0.02     0.08 1.00     2462     1542
versionLike                0.01      0.05    -0.09     0.11 1.00     1138     1363
versionLike&dislike        0.06      0.05    -0.05     0.16 1.00     1197     1312
hu_age                     0.00      0.08    -0.15     0.17 1.00     1742     1258
hu_male                   -0.18      3.06    -6.70     6.01 1.00     1623     1276
hu_edu                    -0.21      1.72    -3.73     3.14 1.00     1452     1109
hu_versionLike             0.02      4.50    -9.31     9.34 1.00      984      847
hu_versionLike&dislike    -0.35      4.85   -10.13     9.97 1.00      968      796

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape     3.95      0.22     3.53     4.41 1.00     1811     1518

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
     Outcome    term                 contrast estimate conf.low conf.high
1 pri_con_fs     age                    dY/dX -0.00603  -0.0150   0.00251
2 pri_con_fs     edu                    dY/dX  0.08441  -0.0723   0.25200
3 pri_con_fs    male                    1 - 0 -0.17139  -0.4463   0.08962
4 pri_con_fs version Like & dislike - Control  0.18857  -0.1520   0.53437
5 pri_con_fs version           Like - Control  0.03483  -0.2940   0.35976
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/pricon-ctrl-2.png){width=768}
:::
:::



I then compare the likes and control conditions to the likes plus dislikes condition.



::: {.cell}

```{.r .cell-code}
run_hurdles(
  object = d, 
  name = "pricon_lkdslk",
  outcome = "pri_con_fs", 
  predictor = "version_lkdslk",
  y_lab = "Privacy concerns",
  x_lab = "Experimental conditions"
  )
```

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-59-1.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: pri_con_fs ~ 1 + age + male + edu + version_lkdslk 
         hu ~ 1 + age + male + edu + version_lkdslk
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
                         Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept                    1.27      0.09     1.10     1.45 1.00     1957     1744
hu_Intercept               -10.40      6.27   -24.09    -0.25 1.00     1439     1060
age                         -0.00      0.00    -0.00     0.00 1.01     2036     1456
male                        -0.06      0.04    -0.14     0.03 1.00     2101     1341
edu                          0.03      0.02    -0.02     0.08 1.00     2018     1182
version_lkdslkLike          -0.05      0.05    -0.15     0.06 1.00     1757     1596
version_lkdslkControl       -0.06      0.06    -0.17     0.05 1.00     1735     1491
hu_age                       0.00      0.08    -0.17     0.17 1.00     1657     1490
hu_male                     -0.11      3.18    -6.77     6.61 1.00     1711     1063
hu_edu                      -0.09      1.71    -3.74     3.30 1.00     1779     1406
hu_version_lkdslkLike        0.19      4.64    -8.86    10.19 1.00      847      707
hu_version_lkdslkControl     0.33      4.36    -8.07    10.44 1.00      969      795

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape     3.94      0.23     3.53     4.39 1.01     1891     1410

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
     Outcome           term                 contrast estimate conf.low conf.high
1 pri_con_fs            age                    dY/dX -0.00615  -0.0147    0.0024
2 pri_con_fs            edu                    dY/dX  0.08941  -0.0761    0.2508
3 pri_con_fs           male                    1 - 0 -0.17608  -0.4534    0.0843
4 pri_con_fs version_lkdslk Control - Like & dislike -0.18595  -0.5515    0.1710
5 pri_con_fs version_lkdslk    Like - Like & dislike -0.14961  -0.5158    0.1897
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-59-2.png){width=768}
:::
:::



#### Gratifications General

I first compare likes and likes plus dislikes to the control condition.



::: {.cell}

```{.r .cell-code}
run_hurdles(
  object = d,
  name = "grats_ctrl",
  outcome = "grats_gen_fs", 
  predictor = "version",
  y_lab = "Gratifications general",
  x_lab = "Experimental conditions"
  )
```

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-60-1.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: grats_gen_fs ~ 1 + age + male + edu + version 
         hu ~ 1 + age + male + edu + version
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
                       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept                  1.60      0.05     1.51     1.70 1.00     2213     1640
hu_Intercept              -9.99      6.25   -24.16     0.24 1.01     1152     1100
age                        0.00      0.00    -0.00     0.00 1.00     2020     1526
male                      -0.02      0.02    -0.06     0.03 1.00     2163     1388
edu                       -0.01      0.01    -0.04     0.02 1.01     2432     1480
versionLike               -0.04      0.03    -0.09     0.01 1.00     2054     1544
versionLike&dislike       -0.04      0.03    -0.10     0.01 1.00     1923     1566
hu_age                     0.00      0.08    -0.15     0.17 1.01     1386     1117
hu_male                   -0.17      3.16    -6.85     5.95 1.00     1593     1189
hu_edu                    -0.19      1.79    -4.01     3.04 1.01     1462     1034
hu_versionLike            -0.19      4.52   -10.08     9.44 1.01      856      685
hu_versionLike&dislike    -0.49      4.90   -11.34     9.48 1.01      850      660

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape    13.73      0.83    12.15    15.34 1.00     2057     1222

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
       Outcome    term                 contrast estimate conf.low conf.high
1 grats_gen_fs     age                    dY/dX  0.00134 -0.00591   0.00832
2 grats_gen_fs     edu                    dY/dX -0.04882 -0.18226   0.08287
3 grats_gen_fs    male                    1 - 0 -0.07841 -0.28216   0.14365
4 grats_gen_fs version Like & dislike - Control -0.20768 -0.47726   0.05443
5 grats_gen_fs version           Like - Control -0.17352 -0.42764   0.08316
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-60-2.png){width=768}
:::
:::



I then compare the likes and control conditions to the likes plus dislikes condition.



::: {.cell}

```{.r .cell-code}
run_hurdles(
  object = d, 
  name = "grats_lkdslk",
  outcome = "grats_gen_fs", 
  predictor = "version_lkdslk",
  y_lab = "Gratifications",
  x_lab = "Experimental conditions"
  )
```

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-61-1.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: grats_gen_fs ~ 1 + age + male + edu + version_lkdslk 
         hu ~ 1 + age + male + edu + version_lkdslk
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
                         Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept                    1.56      0.05     1.46     1.65 1.00     2200     1396
hu_Intercept               -10.03      6.00   -23.07    -0.23 1.00     1714     1366
age                          0.00      0.00    -0.00     0.00 1.00     2091     1586
male                        -0.01      0.02    -0.06     0.03 1.00     2582     1152
edu                         -0.01      0.01    -0.04     0.02 1.00     2559     1218
version_lkdslkLike           0.01      0.03    -0.05     0.07 1.00     1621     1397
version_lkdslkControl        0.04      0.03    -0.02     0.10 1.00     1460     1307
hu_age                       0.01      0.08    -0.15     0.18 1.00     2462     1498
hu_male                     -0.12      3.09    -6.25     6.30 1.00     1827     1102
hu_edu                      -0.25      1.69    -4.00     2.81 1.00     1655     1089
hu_version_lkdslkLike        0.10      4.37    -8.85     8.99 1.01     1257      934
hu_version_lkdslkControl     0.12      4.49    -9.55     8.99 1.00     1180     1033

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape    13.72      0.79    12.20    15.27 1.00     2332     1356

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
       Outcome           term                 contrast estimate conf.low conf.high
1 grats_gen_fs            age                    dY/dX  0.00135 -0.00542   0.00831
2 grats_gen_fs            edu                    dY/dX -0.04862 -0.18907   0.09542
3 grats_gen_fs           male                    1 - 0 -0.06818 -0.28451   0.14724
4 grats_gen_fs version_lkdslk Control - Like & dislike  0.20348 -0.08725   0.48402
5 grats_gen_fs version_lkdslk    Like - Like & dislike  0.03449 -0.24401   0.31864
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-61-2.png){width=768}
:::
:::



#### Privacy Deliberation

I first compare likes and likes plus dislikes to the control condition.



::: {.cell}

```{.r .cell-code}
run_hurdles(
  object = d,
  name = "pridel_ctrl",
  outcome = "pri_del_fs", 
  predictor = "version",
  y_lab = "Privacy deliberation",
  x_lab = "Experimental conditions"
  )
```

::: {.cell-output .cell-output-stdout}

```
Some PIT values larger than 1! Largest:  1 
Rounding PIT > 1 to 1.
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-62-1.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: pri_del_fs ~ 1 + age + male + edu + version 
         hu ~ 1 + age + male + edu + version
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
                       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept                  1.50      0.05     1.40     1.61 1.00     2297     1794
hu_Intercept             -10.16      6.49   -24.49    -0.07 1.00     1454     1124
age                       -0.00      0.00    -0.00    -0.00 1.00     2171     1253
male                      -0.05      0.02    -0.10    -0.00 1.01     2294     1475
edu                        0.00      0.01    -0.03     0.03 1.00     2470     1540
versionLike               -0.00      0.03    -0.06     0.06 1.00     2115     1628
versionLike&dislike       -0.01      0.03    -0.07     0.05 1.00     2105     1369
hu_age                     0.00      0.08    -0.15     0.16 1.00     2155     1542
hu_male                    0.03      3.23    -6.27     7.09 1.00     1492      981
hu_edu                    -0.21      1.71    -3.92     3.02 1.00     1758     1290
hu_versionLike             0.22      4.59    -9.01    10.30 1.00     1138      897
hu_versionLike&dislike     0.05      4.52    -8.55    10.23 1.00     1124      886

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape    11.37      0.66    10.02    12.71 1.01     2314     1391

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
     Outcome    term                 contrast estimate conf.low conf.high
1 pri_del_fs     age                    dY/dX -0.00936  -0.0157  -0.00299
2 pri_del_fs     edu                    dY/dX  0.01353  -0.1038   0.12805
3 pri_del_fs    male                    1 - 0 -0.20687  -0.3962  -0.01175
4 pri_del_fs version Like & dislike - Control -0.04518  -0.2783   0.19111
5 pri_del_fs version           Like - Control -0.00777  -0.2460   0.23641
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-62-2.png){width=768}
:::
:::



I then compare the likes and control conditions to the likes plus dislikes condition.



::: {.cell}

```{.r .cell-code}
run_hurdles(
  object = d, 
  name = "pridel_lkdslk",
  outcome = "pri_del_fs", 
  predictor = "version_lkdslk",
  y_lab = "Privacy deliberation",
  x_lab = "Experimental conditions"
  )
```

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-63-1.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: pri_del_fs ~ 1 + age + male + edu + version_lkdslk 
         hu ~ 1 + age + male + edu + version_lkdslk
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
                         Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept                    1.49      0.05     1.39     1.59 1.00     2306     1740
hu_Intercept               -10.35      6.16   -24.13    -0.21 1.00     1293      998
age                         -0.00      0.00    -0.00    -0.00 1.00     2078     1038
male                        -0.05      0.03    -0.10    -0.00 1.00     2280     1383
edu                          0.00      0.02    -0.03     0.03 1.00     2451     1373
version_lkdslkLike           0.01      0.03    -0.05     0.07 1.00     1848     1557
version_lkdslkControl        0.01      0.03    -0.05     0.07 1.00     1711     1521
hu_age                       0.00      0.08    -0.15     0.17 1.00     1849     1254
hu_male                     -0.09      3.21    -6.86     6.76 1.00     1503      899
hu_edu                      -0.17      1.77    -4.06     3.20 1.01     1602     1115
hu_version_lkdslkLike        0.35      4.52    -8.48    10.11 1.01      916      893
hu_version_lkdslkControl     0.31      4.64    -9.06     9.78 1.01      894      816

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape    11.41      0.68    10.15    12.78 1.00     2056     1470

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
     Outcome           term                 contrast estimate conf.low conf.high
1 pri_del_fs            age                    dY/dX -0.00942  -0.0157  -0.00308
2 pri_del_fs            edu                    dY/dX  0.00926  -0.1090   0.13157
3 pri_del_fs           male                    1 - 0 -0.20888  -0.4076  -0.00539
4 pri_del_fs version_lkdslk Control - Like & dislike  0.04111  -0.1951   0.27915
5 pri_del_fs version_lkdslk    Like - Like & dislike  0.03740  -0.2113   0.29158
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-63-2.png){width=768}
:::
:::



#### Trust general

I first compare likes and likes plus dislikes to the control condition.



::: {.cell}

```{.r .cell-code}
run_hurdles(
  object = d,
  name = "trust_ctrl",
  outcome = "trust_gen_fs", 
  predictor = "version",
  y_lab = "Trust general",
  x_lab = "Experimental conditions"
  )
```

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-64-1.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: trust_gen_fs ~ 1 + age + male + edu + version 
         hu ~ 1 + age + male + edu + version
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
                       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept                  1.69      0.03     1.63     1.76 1.00     2224     1589
hu_Intercept             -10.15      6.41   -24.55    -0.02 1.00     1176      949
age                       -0.00      0.00    -0.00     0.00 1.00     2167     1427
male                       0.00      0.02    -0.03     0.04 1.00     1932     1227
edu                        0.00      0.01    -0.02     0.02 1.00     2346     1405
versionLike               -0.03      0.02    -0.07     0.01 1.00     1674     1459
versionLike&dislike       -0.02      0.02    -0.06     0.01 1.00     1566     1425
hu_age                     0.01      0.09    -0.16     0.19 1.00     1386     1175
hu_male                   -0.22      3.11    -7.09     5.82 1.00      741      712
hu_edu                    -0.25      1.75    -4.05     3.02 1.00     1308     1044
hu_versionLike             0.03      4.32    -8.66     8.46 1.00     1021     1052
hu_versionLike&dislike    -0.25      4.64   -10.25     8.92 1.01      664      770

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape    27.48      1.69    24.25    30.68 1.00     2065     1683

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
       Outcome    term                 contrast estimate conf.low conf.high
1 trust_gen_fs     age                    dY/dX -0.00295 -0.00891   0.00242
2 trust_gen_fs     edu                    dY/dX  0.00396 -0.11039   0.11216
3 trust_gen_fs    male                    1 - 0  0.01070 -0.15953   0.19304
4 trust_gen_fs version Like & dislike - Control -0.12539 -0.33124   0.07588
5 trust_gen_fs version           Like - Control -0.14677 -0.35347   0.05887
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-64-2.png){width=768}
:::
:::



I then compare the likes and control conditions to the likes plus dislikes condition.



::: {.cell}

```{.r .cell-code}
run_hurdles(
  object = d, 
  name = "trust_lkdslk",
  outcome = "trust_gen_fs", 
  predictor = "version_lkdslk",
  y_lab = "Trust general",
  x_lab = "Experimental conditions"
  )
```

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-65-1.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: trust_gen_fs ~ 1 + age + male + edu + version_lkdslk 
         hu ~ 1 + age + male + edu + version_lkdslk
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
                         Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept                    1.67      0.03     1.60     1.73 1.00     2096     1596
hu_Intercept               -10.34      6.29   -24.97    -0.12 1.00     1279     1123
age                         -0.00      0.00    -0.00     0.00 1.00     2019     1329
male                         0.00      0.02    -0.03     0.03 1.00     2018     1638
edu                         -0.00      0.01    -0.02     0.02 1.00     2466     1432
version_lkdslkLike          -0.00      0.02    -0.04     0.03 1.00     1554     1443
version_lkdslkControl        0.02      0.02    -0.01     0.07 1.00     1640     1340
hu_age                       0.01      0.08    -0.16     0.18 1.00     1421     1082
hu_male                     -0.24      3.37    -7.22     6.61 1.00     1332      936
hu_edu                      -0.18      1.79    -4.06     3.01 1.00     1337      836
hu_version_lkdslkLike        0.21      4.44    -8.54     9.60 1.01      830      859
hu_version_lkdslkControl     0.36      4.44    -8.77    10.47 1.00     1006      925

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape    27.48      1.66    24.33    30.74 1.00     1679     1110

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
       Outcome           term                 contrast estimate conf.low conf.high
1 trust_gen_fs            age                    dY/dX -0.00301  -0.0088   0.00277
2 trust_gen_fs            edu                    dY/dX -0.00161  -0.1011   0.10916
3 trust_gen_fs           male                    1 - 0  0.01412  -0.1510   0.17935
4 trust_gen_fs version_lkdslk Control - Like & dislike  0.12205  -0.0796   0.34126
5 trust_gen_fs version_lkdslk    Like - Like & dislike -0.02307  -0.2291   0.17767
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-65-2.png){width=768}
:::
:::



#### Self-Efficacy

I first compare likes and likes plus dislikes to the control condition.



::: {.cell}

```{.r .cell-code}
run_hurdles(
  object = d,
  name = "selfeff_ctrl",
  outcome = "self_eff_fs", 
  predictor = "version",
  y_lab = "Self-efficacy",
  x_lab = "Experimental conditions"
  )
```

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-66-1.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: self_eff_fs ~ 1 + age + male + edu + version 
         hu ~ 1 + age + male + edu + version
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
                       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept                  1.62      0.03     1.56     1.69 1.00     2104     1740
hu_Intercept              -9.84      6.34   -24.23     0.22 1.00     1250     1023
age                       -0.00      0.00    -0.00     0.00 1.00     2210     1465
male                       0.01      0.02    -0.02     0.04 1.00     2138     1375
edu                        0.02      0.01     0.01     0.04 1.01     2759     1524
versionLike               -0.02      0.02    -0.06     0.02 1.00     1473     1549
versionLike&dislike       -0.02      0.02    -0.06     0.02 1.00     1370     1405
hu_age                     0.00      0.09    -0.16     0.18 1.00     1492     1276
hu_male                   -0.08      3.09    -6.45     6.65 1.00     1326     1112
hu_edu                    -0.23      1.74    -4.09     2.72 1.00     1277      975
hu_versionLike            -0.22      4.58   -10.32     9.15 1.00     1093      896
hu_versionLike&dislike    -0.36      4.69   -10.25     8.75 1.00      994      828

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape    27.92      1.66    24.78    31.32 1.00     1980     1401

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
      Outcome    term                 contrast  estimate conf.low conf.high
1 self_eff_fs     age                    dY/dX -0.000121 -0.00599   0.00534
2 self_eff_fs     edu                    dY/dX  0.132233  0.02890   0.24105
3 self_eff_fs    male                    1 - 0  0.077955 -0.08572   0.23549
4 self_eff_fs version Like & dislike - Control -0.091942 -0.29888   0.12772
5 self_eff_fs version           Like - Control -0.088026 -0.30151   0.10745
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-66-2.png){width=768}
:::
:::



I then compare the likes and control conditions to the likes plus dislikes condition.



::: {.cell}

```{.r .cell-code}
run_hurdles(
  object = d, 
  name = "selfeff_lkdslk",
  outcome = "self_eff_fs", 
  predictor = "version_lkdslk",
  y_lab = "Self-efficacy",
  x_lab = "Experimental conditions"
  )
```

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-67-1.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: self_eff_fs ~ 1 + age + male + edu + version_lkdslk 
         hu ~ 1 + age + male + edu + version_lkdslk
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
                         Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept                    1.61      0.03     1.54     1.67 1.00     2202     1645
hu_Intercept               -10.37      6.32   -25.25    -0.56 1.00      940      634
age                         -0.00      0.00    -0.00     0.00 1.00     1957     1401
male                         0.01      0.02    -0.02     0.05 1.00     2266     1480
edu                          0.02      0.01     0.00     0.04 1.01     2354     1256
version_lkdslkLike           0.00      0.02    -0.04     0.04 1.00     1689     1598
version_lkdslkControl        0.02      0.02    -0.02     0.06 1.00     1503     1490
hu_age                       0.00      0.08    -0.15     0.19 1.00     1334      949
hu_male                     -0.11      3.09    -6.49     5.98 1.00     1469      996
hu_edu                      -0.19      1.69    -3.93     2.93 1.00     1362     1227
hu_version_lkdslkLike        0.25      4.78    -9.42    10.58 1.00      664      388
hu_version_lkdslkControl     0.44      4.75    -9.09    10.82 1.00      717      382

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape    27.94      1.66    24.81    31.21 1.00     2019     1385

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
      Outcome           term                 contrast  estimate conf.low conf.high
1 self_eff_fs            age                    dY/dX -9.91e-05 -0.00606   0.00556
2 self_eff_fs            edu                    dY/dX  1.32e-01  0.02024   0.24161
3 self_eff_fs           male                    1 - 0  7.93e-02 -0.09587   0.25617
4 self_eff_fs version_lkdslk Control - Like & dislike  9.64e-02 -0.11753   0.30120
5 self_eff_fs version_lkdslk    Like - Like & dislike  1.86e-03 -0.20942   0.21240
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-67-2.png){width=768}
:::
:::



### Words
#### Privacy Concerns



::: {.cell}

```{.r .cell-code}
run_hurdles(
  object = d, 
  name = "wrds_pricon",
  outcome = "words", 
  predictor = "pri_con_fs",
  y_lab = "Words",
  x_lab = "Privacy concerns"
  )
```

::: {.cell-output .cell-output-stdout}

```
Some PIT values larger than 1! Largest:  1 
Rounding PIT > 1 to 1.
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-68-1.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: words ~ 1 + age + male + edu + pri_con_fs 
         hu ~ 1 + age + male + edu + pri_con_fs
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept         5.64      0.40     4.89     6.47 1.00     4174     1741
hu_Intercept      0.53      0.41    -0.27     1.35 1.00     3733     1624
age              -0.00      0.01    -0.01     0.01 1.00     4583     1462
male             -0.07      0.16    -0.40     0.25 1.01     4482     1490
edu              -0.01      0.09    -0.19     0.17 1.00     4300     1284
pri_con_fs       -0.10      0.06    -0.21     0.02 1.00     3952     1489
hu_age           -0.01      0.01    -0.02     0.00 1.00     4139     1414
hu_male          -0.04      0.17    -0.38     0.28 1.00     3866     1619
hu_edu           -0.23      0.10    -0.44    -0.02 1.00     3052     1421
hu_pri_con_fs     0.19      0.06     0.08     0.30 1.01     4240     1205

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape     0.69      0.05     0.59     0.80 1.00     3569     1186

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
  Outcome       term contrast estimate conf.low conf.high
1   words        age    dY/dX    0.196   -0.894      1.18
2   words        edu    dY/dX    9.061   -7.606     26.02
3   words       male    1 - 0   -3.850  -33.832     25.48
4   words pri_con_fs    dY/dX  -15.644  -27.360     -5.41
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-68-2.png){width=768}
:::
:::



Let's explore how number of communicated words differ for people very much concerned versus not concerned at all.



::: {.cell}

```{.r .cell-code}
# let's extract effects for 1
wrds_pricon_1 <- 
  fig_wrds_pricon$data %>% 
  filter(pri_con_fs == min(pri_con_fs)) %>% 
  pull(estimate__) %>% 
  round(0)

# let's extract effects for 7
wrds_pricon_7 <- 
  fig_wrds_pricon$data %>% 
  filter(pri_con_fs == max(pri_con_fs)) %>% 
  pull(estimate__) %>% 
  round(0)
```
:::



People who were not concerned at all communicated on average 115 words, whereas people strongly concerned communicated on average 33 words.

#### Gratifications General



::: {.cell}

```{.r .cell-code}
run_hurdles(
  object = d, 
  name = "wrds_grats",
  outcome = "words", 
  predictor = "grats_gen_fs",
  y_lab = "Words",
  x_lab = "Gratifications general"
  )
```

::: {.cell-output .cell-output-stdout}

```
Some PIT values larger than 1! Largest:  1 
Rounding PIT > 1 to 1.
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-70-1.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: words ~ 1 + age + male + edu + grats_gen_fs 
         hu ~ 1 + age + male + edu + grats_gen_fs
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
                Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept           4.27      0.47     3.39     5.19 1.00     3398     1657
hu_Intercept        2.27      0.54     1.22     3.34 1.00     3857     1122
age                -0.00      0.01    -0.01     0.01 1.00     3707     1716
male               -0.08      0.16    -0.39     0.22 1.00     3399     1512
edu                 0.03      0.09    -0.14     0.21 1.00     4560     1479
grats_gen_fs        0.22      0.07     0.08     0.36 1.00     3132     1677
hu_age             -0.01      0.01    -0.02     0.00 1.00     3834     1440
hu_male            -0.08      0.17    -0.43     0.24 1.00     3432     1309
hu_edu             -0.22      0.10    -0.43    -0.02 1.00     3952     1524
hu_grats_gen_fs    -0.22      0.08    -0.39    -0.07 1.00     3647     1480

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape     0.71      0.06     0.60     0.82 1.01     4200     1636

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
  Outcome         term contrast estimate conf.low conf.high
1   words          age    dY/dX   0.0772   -0.995      1.11
2   words          edu    dY/dX  12.1441   -4.511     29.08
3   words grats_gen_fs    dY/dX  26.8528   13.383     43.19
4   words         male    1 - 0  -2.1608  -32.251     26.05
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-70-2.png){width=768}
:::
:::



#### Privacy Deliberation



::: {.cell}

```{.r .cell-code}
run_hurdles(
  object = d, 
  name = "wrds_pridel",
  outcome = "words", 
  predictor = "pri_del_fs",
  y_lab = "Words",
  x_lab = "Privacy deliberation"
  )
```

::: {.cell-output .cell-output-stdout}

```
Some PIT values larger than 1! Largest:  1 
Rounding PIT > 1 to 1.
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-71-1.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: words ~ 1 + age + male + edu + pri_del_fs 
         hu ~ 1 + age + male + edu + pri_del_fs
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept         6.50      0.49     5.55     7.50 1.00     3651     1552
hu_Intercept     -0.15      0.50    -1.12     0.82 1.00     4066     1739
age              -0.01      0.01    -0.02     0.01 1.00     3600     1655
male             -0.08      0.16    -0.38     0.23 1.00     3782     1520
edu              -0.02      0.09    -0.21     0.16 1.00     3999     1453
pri_del_fs       -0.26      0.07    -0.41    -0.12 1.00     3475     1709
hu_age           -0.01      0.01    -0.02     0.00 1.00     5239     1689
hu_male          -0.01      0.18    -0.36     0.33 1.00     3433     1588
hu_edu           -0.22      0.10    -0.42    -0.02 1.00     4095     1500
hu_pri_del_fs     0.31      0.08     0.15     0.47 1.00     3723     1417

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape     0.71      0.06     0.61     0.83 1.00     4791     1332

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
  Outcome       term contrast estimate conf.low conf.high
1   words        age    dY/dX   -0.162    -1.19     0.837
2   words        edu    dY/dX    7.478    -8.77    24.583
3   words       male    1 - 0   -5.784   -34.74    23.740
4   words pri_del_fs    dY/dX  -33.422   -51.55   -19.205
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-71-2.png){width=768}
:::
:::



#### Trust General



::: {.cell}

```{.r .cell-code}
run_hurdles(
  object = d, 
  name = "wrds_trust",
  outcome = "words", 
  predictor = "trust_gen_fs",
  y_lab = "Words",
  x_lab = "Trust general"
  )
```

::: {.cell-output .cell-output-stdout}

```
Some PIT values larger than 1! Largest:  1 
Rounding PIT > 1 to 1.
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-72-1.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: words ~ 1 + age + male + edu + trust_gen_fs 
         hu ~ 1 + age + male + edu + trust_gen_fs
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
                Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept           4.24      0.60     3.10     5.44 1.00     3590     1798
hu_Intercept        3.10      0.63     1.89     4.34 1.00     4184     1533
age                -0.00      0.01    -0.01     0.01 1.00     4685     1605
male               -0.07      0.17    -0.39     0.25 1.00     3482     1514
edu                -0.01      0.10    -0.20     0.18 1.00     3924     1738
trust_gen_fs        0.21      0.10     0.01     0.41 1.00     3433     1612
hu_age             -0.01      0.01    -0.02     0.00 1.01     3962     1478
hu_male            -0.06      0.18    -0.42     0.30 1.01     4001     1366
hu_edu             -0.21      0.11    -0.42    -0.00 1.00     3729     1681
hu_trust_gen_fs    -0.36      0.09    -0.54    -0.18 1.00     3719     1484

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape     0.70      0.06     0.59     0.81 1.00     3331     1713

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
  Outcome         term contrast estimate conf.low conf.high
1   words          age    dY/dX    0.265   -0.836      1.28
2   words          edu    dY/dX    8.338   -9.218     26.66
3   words         male    1 - 0   -3.136  -33.737     28.36
4   words trust_gen_fs    dY/dX   31.799   13.922     52.97
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-72-2.png){width=768}
:::
:::



#### Self-Efficacy



::: {.cell}

```{.r .cell-code}
run_hurdles(
  object = d, 
  name = "wrds_selfeff",
  outcome = "words", 
  predictor = "self_eff_fs",
  y_lab = "Words",
  x_lab = "Self-efficacy"
  )
```

::: {.cell-output .cell-output-stdout}

```
Some PIT values larger than 1! Largest:  1 
Rounding PIT > 1 to 1.
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-73-1.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: words ~ 1 + age + male + edu + self_eff_fs 
         hu ~ 1 + age + male + edu + self_eff_fs
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept          2.17      0.47     1.28     3.05 1.00     3335     1514
hu_Intercept       5.79      0.68     4.46     7.11 1.00     4031     1346
age               -0.00      0.01    -0.01     0.01 1.00     3794     1489
male              -0.15      0.15    -0.44     0.14 1.01     5059     1534
edu               -0.06      0.08    -0.23     0.10 1.01     2935     1386
self_eff_fs        0.58      0.07     0.43     0.72 1.00     3030     1604
hu_age            -0.01      0.01    -0.02     0.00 1.00     3388     1590
hu_male           -0.01      0.19    -0.37     0.37 1.00     5339     1706
hu_edu            -0.13      0.11    -0.33     0.10 1.00     4355     1459
hu_self_eff_fs    -0.88      0.11    -1.09    -0.67 1.00     4023     1460

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape     0.80      0.06     0.68     0.93 1.00     4652     1617

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
  Outcome        term contrast estimate conf.low conf.high
1   words         age    dY/dX   0.0257   -0.875      0.89
2   words         edu    dY/dX  -0.6898  -15.538     13.76
3   words        male    1 - 0 -11.5440  -38.094     13.43
4   words self_eff_fs    dY/dX  72.8927   55.437     95.20
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-73-2.png){width=768}
:::
:::



The effects is markedly exponential. Let's inspect the changes across values in self-efficacy.



::: {.cell}

```{.r .cell-code}
effects_wrds_selfeff <- 
  slopes(
    fit_wrds_selfeff, 
    newdata = datagrid(
      self_eff_fs = c(1, 6, 7)
      ), 
   grid_type = "counterfactual"
    ) %>%
  filter(
    term == "self_eff_fs"
    ) %T>%
  print()
```

::: {.cell-output .cell-output-stdout}

```

 self_eff_fs Estimate    2.5 %  97.5 %
           1    0.251   0.0909   0.647
           6  111.285  80.6053 155.076
           7  219.493 139.9139 345.562

Term: self_eff_fs
Type: response
Comparison: dY/dX
```


:::

```{.r .cell-code}
# let's extract effects for 1
effects_wrds_selfeff_12 <- 
  effects_wrds_selfeff %>% 
  filter(
    self_eff_fs == 1 &
      term == "self_eff_fs") %>% 
  select(estimate)

# let's extract effects for 7
effects_wrds_selfeff_67 <- 
  effects_wrds_selfeff %>% 
  filter(
    self_eff_fs == 6 &
      term == "self_eff_fs") %>% 
  select(estimate)
```
:::



Whereas a change in self-efficacy from 1 to 2 led to an increase of 0 words, a change from 6 to 7 led to an increase of 111 words.

Let's explore  how number of communicated words differ for people very much self-efficacious versus not self-efficacious at all.



::: {.cell}

```{.r .cell-code}
# let's extract effects for 1
wrds_selfeff_1 <- 
  fig_wrds_selfeff$data %>% 
  filter(self_eff_fs == min(self_eff_fs)) %>% 
  pull(estimate__) %>% 
  round(0)

# let's extract effects for 7
wrds_selfeff_7 <- 
  fig_wrds_selfeff$data %>% 
  filter(self_eff_fs == max(self_eff_fs)) %>% 
  pull(estimate__) %>% 
  round(0)
```
:::



People who reported no self-efficacy at all communicated on average 1 words, whereas people with very high self-efficacy communicated on average 254 words.

## Visualization



::: {.cell}

```{.r .cell-code}
fig_results <- 
  plot_grid(
    fig_wrds_pricon, 
    fig_wrds_grats, 
    fig_wrds_pridel, 
    fig_wrds_trust, 
    fig_wrds_selfeff,
    fig_wrds_ctrl
    ) %T>% 
  print()
```

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-76-1.png){width=768}
:::

```{.r .cell-code}
# Save visualizations
ggsave("figures/results/effects.pdf")
ggsave("figures/results/effects.png")
```
:::



## Table

Make table of slopes



::: {.cell}

```{.r .cell-code}
tab_slopes <- 
  rbind(
    slopes_wrds_pricon,
    slopes_wrds_grats,
    slopes_wrds_pridel,
    slopes_wrds_trust,
    slopes_wrds_selfeff,
    slopes_wrds_ctrl,
    slopes_wrds_lkdslk,
    slopes_pricon_ctrl,
    slopes_pricon_lkdslk,
    slopes_grats_ctrl,
    slopes_grats_lkdslk,
    slopes_pridel_ctrl,
    slopes_pridel_lkdslk,
    slopes_trust_ctrl,
    slopes_trust_lkdslk,
    slopes_selfeff_ctrl,
    slopes_selfeff_lkdslk
    ) %>% 
  as.data.frame() %>% 
  mutate(
    term = ifelse(
      .$contrast == "Like - Control", 
      "Like - Control",
      ifelse(
        .$contrast == "Like & dislike - Control", 
        "Like & dislike - Control",
        ifelse(
          .$contrast == "Like - Like & dislike", 
          "Like - Like & dislike", 
          .$term
          )
        )
      ),
    Predictor = recode(
      term, 
      "Like - Control" = "Like vs. control",
      "Like & dislike - Control" = "Like & dislike vs. control",
      "Like - Like & dislike" = "Like vs. like & dislike",
      "pri_con_fs" = "Privacy concerns",
      "pri_del_fs" = "Privacy deliberation",
      "grats_gen_fs" = "Expected gratifications",
      "trust_gen_fs" = "Trust",
      "self_eff_fs" = "Self-efficacy"
    ),
    Outcome = recode(
      Outcome,
      "pri_con_fs" = "Privacy concerns",
      "pri_del_fs" = "Privacy deliberation",
      "grats_gen_fs" = "Expected gratifications",
      "trust_gen_fs" = "Trust",
      "self_eff_fs" = "Self-efficacy"
    )
  ) %>% 
  filter(
    Predictor %in% c(
      "Like vs. control",
      "Like & dislike vs. control",
      "Like vs. like & dislike",
      "Privacy concerns",
      "Privacy deliberation",
      "Expected gratifications",
      "Trust",
      "Self-efficacy")
  ) %>% 
  mutate(
    Predictor = factor(Predictor, levels = c(
      "Self-efficacy",
      "Trust",
      "Privacy deliberation",
      "Expected gratifications",
      "Privacy concerns",
      "Like vs. like & dislike",
      "Like & dislike vs. control",
      "Like vs. control"
      )
    ),
    Estimate = 
      ifelse(
        Outcome == "words",
        round(estimate),
        estimate
      ),
    LL = 
      ifelse(
        Outcome == "words",
        round(conf.low),
        conf.low
      ),
    UL = 
      ifelse(
        Outcome == "words",
        round(conf.high),
        conf.high
      ),
  ) %>% 
  select(
    Outcome,
    Predictor, 
    Estimate,
    LL, 
    UL
    ) 

tab_slopes %>% 
  kable() %>% 
  kable_styling("striped") %>% 
  scroll_box(width = "100%")
```

::: {.cell-output-display}
`````{=html}
<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; "><table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Outcome </th>
   <th style="text-align:left;"> Predictor </th>
   <th style="text-align:right;"> Estimate </th>
   <th style="text-align:right;"> LL </th>
   <th style="text-align:right;"> UL </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> words </td>
   <td style="text-align:left;"> Privacy concerns </td>
   <td style="text-align:right;"> -16.000 </td>
   <td style="text-align:right;"> -27.000 </td>
   <td style="text-align:right;"> -5.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> words </td>
   <td style="text-align:left;"> Expected gratifications </td>
   <td style="text-align:right;"> 27.000 </td>
   <td style="text-align:right;"> 13.000 </td>
   <td style="text-align:right;"> 43.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> words </td>
   <td style="text-align:left;"> Privacy deliberation </td>
   <td style="text-align:right;"> -33.000 </td>
   <td style="text-align:right;"> -52.000 </td>
   <td style="text-align:right;"> -19.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> words </td>
   <td style="text-align:left;"> Trust </td>
   <td style="text-align:right;"> 32.000 </td>
   <td style="text-align:right;"> 14.000 </td>
   <td style="text-align:right;"> 53.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> words </td>
   <td style="text-align:left;"> Self-efficacy </td>
   <td style="text-align:right;"> 73.000 </td>
   <td style="text-align:right;"> 55.000 </td>
   <td style="text-align:right;"> 95.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> words </td>
   <td style="text-align:left;"> Like &amp; dislike vs. control </td>
   <td style="text-align:right;"> -47.000 </td>
   <td style="text-align:right;"> -90.000 </td>
   <td style="text-align:right;"> -9.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> words </td>
   <td style="text-align:left;"> Like vs. control </td>
   <td style="text-align:right;"> -35.000 </td>
   <td style="text-align:right;"> -78.000 </td>
   <td style="text-align:right;"> 3.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> words </td>
   <td style="text-align:left;"> Like vs. like &amp; dislike </td>
   <td style="text-align:right;"> 13.000 </td>
   <td style="text-align:right;"> -18.000 </td>
   <td style="text-align:right;"> 45.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Privacy concerns </td>
   <td style="text-align:left;"> Like &amp; dislike vs. control </td>
   <td style="text-align:right;"> 0.189 </td>
   <td style="text-align:right;"> -0.152 </td>
   <td style="text-align:right;"> 0.534 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Privacy concerns </td>
   <td style="text-align:left;"> Like vs. control </td>
   <td style="text-align:right;"> 0.035 </td>
   <td style="text-align:right;"> -0.294 </td>
   <td style="text-align:right;"> 0.360 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Privacy concerns </td>
   <td style="text-align:left;"> Like vs. like &amp; dislike </td>
   <td style="text-align:right;"> -0.150 </td>
   <td style="text-align:right;"> -0.516 </td>
   <td style="text-align:right;"> 0.190 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Expected gratifications </td>
   <td style="text-align:left;"> Like &amp; dislike vs. control </td>
   <td style="text-align:right;"> -0.208 </td>
   <td style="text-align:right;"> -0.477 </td>
   <td style="text-align:right;"> 0.054 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Expected gratifications </td>
   <td style="text-align:left;"> Like vs. control </td>
   <td style="text-align:right;"> -0.174 </td>
   <td style="text-align:right;"> -0.428 </td>
   <td style="text-align:right;"> 0.083 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Expected gratifications </td>
   <td style="text-align:left;"> Like vs. like &amp; dislike </td>
   <td style="text-align:right;"> 0.034 </td>
   <td style="text-align:right;"> -0.244 </td>
   <td style="text-align:right;"> 0.319 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Privacy deliberation </td>
   <td style="text-align:left;"> Like &amp; dislike vs. control </td>
   <td style="text-align:right;"> -0.045 </td>
   <td style="text-align:right;"> -0.278 </td>
   <td style="text-align:right;"> 0.191 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Privacy deliberation </td>
   <td style="text-align:left;"> Like vs. control </td>
   <td style="text-align:right;"> -0.008 </td>
   <td style="text-align:right;"> -0.246 </td>
   <td style="text-align:right;"> 0.236 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Privacy deliberation </td>
   <td style="text-align:left;"> Like vs. like &amp; dislike </td>
   <td style="text-align:right;"> 0.037 </td>
   <td style="text-align:right;"> -0.211 </td>
   <td style="text-align:right;"> 0.292 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Trust </td>
   <td style="text-align:left;"> Like &amp; dislike vs. control </td>
   <td style="text-align:right;"> -0.125 </td>
   <td style="text-align:right;"> -0.331 </td>
   <td style="text-align:right;"> 0.076 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Trust </td>
   <td style="text-align:left;"> Like vs. control </td>
   <td style="text-align:right;"> -0.147 </td>
   <td style="text-align:right;"> -0.353 </td>
   <td style="text-align:right;"> 0.059 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Trust </td>
   <td style="text-align:left;"> Like vs. like &amp; dislike </td>
   <td style="text-align:right;"> -0.023 </td>
   <td style="text-align:right;"> -0.229 </td>
   <td style="text-align:right;"> 0.178 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Self-efficacy </td>
   <td style="text-align:left;"> Like &amp; dislike vs. control </td>
   <td style="text-align:right;"> -0.092 </td>
   <td style="text-align:right;"> -0.299 </td>
   <td style="text-align:right;"> 0.128 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Self-efficacy </td>
   <td style="text-align:left;"> Like vs. control </td>
   <td style="text-align:right;"> -0.088 </td>
   <td style="text-align:right;"> -0.302 </td>
   <td style="text-align:right;"> 0.107 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Self-efficacy </td>
   <td style="text-align:left;"> Like vs. like &amp; dislike </td>
   <td style="text-align:right;"> 0.002 </td>
   <td style="text-align:right;"> -0.209 </td>
   <td style="text-align:right;"> 0.212 </td>
  </tr>
</tbody>
</table></div>

`````
:::
:::





### Self Disclosure

#### Self Disclosure Variable


::: {.cell}

```{.r .cell-code}
# self-disclosure index: quantity (log words), breadth and depth, z-standardized and averaged
# LLM codes per comment from 07a, aggregated per person
comments_coded <- read_csv("data/comments/final/comments_classified_final.csv")

sd_components <- 
  comments_coded %>% 
  filter(!is.na(id)) %>% 
  group_by(id) %>% 
  summarise(
    breadth_mean = mean(breadth, na.rm = TRUE),   # 0-3, mean topics per comment
    depth_mean   = mean(depth, na.rm = TRUE)      # 0-3, mean of pers_exp + emot_exp + pol_opin
    )

d <- 
  d %>% 
  select(-any_of(c("breadth_mean", "depth_mean"))) %>% 
  left_join(sd_components, by = "id") %>% 
  mutate(
    # non-commenters disclosed nothing
    across(c(breadth_mean, depth_mean), ~ replace_na(.x, 0)),
    z_quantity = as.numeric(scale(words_log)),
    z_breadth  = as.numeric(scale(breadth_mean)),
    z_depth    = as.numeric(scale(depth_mean)),
    self_disc_raw = (z_quantity + z_breadth + z_depth) / 3,
    # shift so that "no disclosure" is exactly 0 (needed for the hurdle model)
    self_disc = self_disc_raw - self_disc_raw[words == 0][1],
    self_disc = if_else(words == 0, 0, self_disc)
    )

# checks everything
stopifnot(all(d$self_disc[d$words > 0] > 0))
table(d$self_disc == 0, d$words == 0)
```

::: {.cell-output .cell-output-stdout}

```
       
        FALSE TRUE
  FALSE   235    0
  TRUE      0  324
```


:::

```{.r .cell-code}
summary(d$self_disc[d$words > 0])
```

::: {.cell-output .cell-output-stdout}

```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.0999  1.0496  1.6330  1.6123  2.1463  3.9570 
```


:::

```{.r .cell-code}
cor(d %>% filter(words > 0) %>% select(words_log, breadth_mean, depth_mean))
```

::: {.cell-output .cell-output-stdout}

```
             words_log breadth_mean depth_mean
words_log        1.000        0.564      0.497
breadth_mean     0.564        1.000      0.597
depth_mean       0.497        0.597      1.000
```


:::
:::

::: {.cell}

```{.r .cell-code}
   run_hurdles(object = d, name = "self_disc_ctrl", outcome = "self_disc",
               predictor = "version", y_lab = "Self Disclosure", x_lab = "Experimental conditions")
```

::: {.cell-output .cell-output-stdout}

```
Some PIT values larger than 1! Largest:  1 
Rounding PIT > 1 to 1.
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-79-1.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: self_disc ~ 1 + age + male + edu + version 
         hu ~ 1 + age + male + edu + version
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
                       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept                  0.63      0.16     0.31     0.95 1.00     2312     1880
hu_Intercept               1.06      0.38     0.30     1.80 1.00     3396     1870
age                       -0.00      0.00    -0.01     0.00 1.00     2499     1398
male                       0.01      0.08    -0.14     0.16 1.00     2386     1663
edu                       -0.02      0.04    -0.10     0.07 1.00     2075     1803
versionLike               -0.06      0.09    -0.23     0.11 1.00     2124     1606
versionLike&dislike       -0.11      0.10    -0.30     0.08 1.00     1949     1536
hu_age                    -0.01      0.01    -0.02     0.00 1.00     3632     1425
hu_male                   -0.07      0.19    -0.44     0.29 1.00     2128     1141
hu_edu                    -0.21      0.11    -0.42     0.00 1.00     2577     1662
hu_versionLike             0.02      0.20    -0.37     0.41 1.00     1837     1763
hu_versionLike&dislike     0.23      0.21    -0.19     0.65 1.00     1811     1291

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape     3.14      0.27     2.63     3.70 1.00     2567     1675

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
    Outcome    term                 contrast estimate conf.low conf.high
1 self_disc     age                    dY/dX  0.00239 -0.00339   0.00791
2 self_disc     edu                    dY/dX  0.06737 -0.03228   0.17000
3 self_disc    male                    1 - 0  0.03192 -0.13738   0.20683
4 self_disc version Like & dislike - Control -0.16464 -0.37315   0.04368
5 self_disc version           Like - Control -0.05237 -0.26705   0.14868
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-79-2.png){width=768}
:::

```{.r .cell-code}
   run_hurdles(object = d, name = "self_disc_lkdslk", outcome = "self_disc",
               predictor = "version_lkdslk", y_lab = "Self Disclosure", x_lab = "Experimental conditions")
```

::: {.cell-output .cell-output-stdout}

```
Some PIT values larger than 1! Largest:  1 
Rounding PIT > 1 to 1.
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-79-3.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: self_disc ~ 1 + age + male + edu + version_lkdslk 
         hu ~ 1 + age + male + edu + version_lkdslk
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
                         Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept                    0.52      0.16     0.20     0.84 1.00     2326     1892
hu_Intercept                 1.30      0.36     0.59     1.98 1.00     3080     1855
age                         -0.00      0.00    -0.01     0.00 1.01     2620     1296
male                         0.01      0.07    -0.13     0.16 1.00     2399     1316
edu                         -0.02      0.04    -0.10     0.07 1.00     2383     1605
version_lkdslkLike           0.05      0.09    -0.13     0.23 1.00     1810     1535
version_lkdslkControl        0.11      0.09    -0.06     0.29 1.00     1731     1523
hu_age                      -0.01      0.01    -0.02     0.00 1.00     3935     1459
hu_male                     -0.06      0.18    -0.42     0.27 1.00     2189     1330
hu_edu                      -0.21      0.10    -0.41    -0.01 1.00     2446     1384
hu_version_lkdslkLike       -0.21      0.22    -0.65     0.20 1.00     1931     1710
hu_version_lkdslkControl    -0.24      0.21    -0.64     0.20 1.00     1954     1664

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape     3.15      0.29     2.61     3.72 1.00     2253     1596

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
    Outcome           term                 contrast estimate conf.low conf.high
1 self_disc            age                    dY/dX  0.00233 -0.00278    0.0079
2 self_disc            edu                    dY/dX  0.06527 -0.03033    0.1644
3 self_disc           male                    1 - 0  0.03230 -0.13330    0.2116
4 self_disc version_lkdslk Control - Like & dislike  0.16348 -0.03324    0.3653
5 self_disc version_lkdslk    Like - Like & dislike  0.10970 -0.08825    0.3116
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-79-4.png){width=768}
:::
:::



#### Privacy Concerns



::: {.cell}

```{.r .cell-code}
run_hurdles(
  object = d, 
  name = "self_disc_pricon",
  outcome = "self_disc", 
  predictor = "pri_con_fs",
  y_lab = "Self Disclosure",
  x_lab = "Privacy concerns"
  )
```

::: {.cell-output .cell-output-stdout}

```
Some PIT values larger than 1! Largest:  1 
Rounding PIT > 1 to 1.
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-80-1.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: self_disc ~ 1 + age + male + edu + pri_con_fs 
         hu ~ 1 + age + male + edu + pri_con_fs
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept         0.62      0.17     0.28     0.95 1.00     2345     1675
hu_Intercept      0.53      0.40    -0.26     1.34 1.00     3097     1327
age              -0.00      0.00    -0.01     0.00 1.00     2305     1391
male              0.02      0.07    -0.13     0.16 1.00     1935     1401
edu              -0.02      0.04    -0.10     0.07 1.00     2446     1568
pri_con_fs       -0.02      0.03    -0.07     0.03 1.00     2092     1424
hu_age           -0.01      0.01    -0.02     0.00 1.00     4304     1498
hu_male          -0.03      0.17    -0.37     0.32 1.00     2235     1666
hu_edu           -0.23      0.11    -0.44    -0.03 1.00     2554     1606
hu_pri_con_fs     0.19      0.06     0.07     0.30 1.00     2673     1758

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape     3.15      0.29     2.63     3.78 1.00     2194     1457

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
    Outcome       term contrast estimate conf.low conf.high
1 self_disc        age    dY/dX  0.00212  -0.0032    0.0074
2 self_disc        edu    dY/dX  0.07462  -0.0181    0.1698
3 self_disc       male    1 - 0  0.02648  -0.1348    0.1902
4 self_disc pri_con_fs    dY/dX -0.08208  -0.1377   -0.0289
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-80-2.png){width=768}
:::
:::



Let's explore how number of communicated words differ for people very much concerned versus not concerned at all.



::: {.cell}

```{.r .cell-code}
# let's extract effects for 1
self_disc_pricon_1 <- 
  fig_self_disc_pricon$data %>% 
  filter(pri_con_fs == min(pri_con_fs)) %>% 
  pull(estimate__) %>% 
  round(0)

# let's extract effects for 7
self_dics_pricon_7 <- 
  fig_self_disc_pricon$data %>% 
  filter(pri_con_fs == max(pri_con_fs)) %>% 
  pull(estimate__) %>% 
  round(0)
```
:::



People who were not concerned at all communicated on average 115 words, whereas people strongly concerned communicated on average 33 words.

#### Gratifications General



::: {.cell}

```{.r .cell-code}
run_hurdles(
  object = d, 
  name = "self_disc_grats",
  outcome = "self_disc", 
  predictor = "grats_gen_fs",
  y_lab = "Self Disclosure",
  x_lab = "Gratifications general"
  )
```

::: {.cell-output .cell-output-stdout}

```
Some PIT values larger than 1! Largest:  1 
Rounding PIT > 1 to 1.
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-82-1.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: self_disc ~ 1 + age + male + edu + grats_gen_fs 
         hu ~ 1 + age + male + edu + grats_gen_fs
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
                Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept           0.70      0.24     0.24     1.18 1.00     2339     1855
hu_Intercept        2.28      0.54     1.25     3.38 1.00     2643     1377
age                -0.00      0.00    -0.01     0.00 1.00     2321     1294
male                0.02      0.08    -0.13     0.16 1.00     2519     1401
edu                -0.02      0.04    -0.11     0.07 1.00     2335     1663
grats_gen_fs       -0.03      0.04    -0.09     0.04 1.00     2738     1597
hu_age             -0.01      0.01    -0.02     0.00 1.00     4575     1730
hu_male            -0.08      0.18    -0.44     0.26 1.00     2673     1651
hu_edu             -0.22      0.10    -0.43    -0.02 1.00     2344     1696
hu_grats_gen_fs    -0.23      0.08    -0.38    -0.08 1.00     2603     1395

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape     3.14      0.28     2.62     3.75 1.00     2400     1410

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
    Outcome         term contrast estimate conf.low conf.high
1 self_disc          age    dY/dX  0.00242 -0.00298   0.00766
2 self_disc          edu    dY/dX  0.07100 -0.02871   0.16869
3 self_disc grats_gen_fs    dY/dX  0.07020 -0.00969   0.14065
4 self_disc         male    1 - 0  0.04674 -0.13116   0.21075
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-82-2.png){width=768}
:::
:::



#### Privacy Deliberation



::: {.cell}

```{.r .cell-code}
run_hurdles(
  object = d, 
  name = "self_disc_pridel",
  outcome = "self_disc", 
  predictor = "pri_del_fs",
  y_lab = "Self Disclosure",
  x_lab = "Privacy deliberation"
  )
```

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-83-1.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: self_disc ~ 1 + age + male + edu + pri_del_fs 
         hu ~ 1 + age + male + edu + pri_del_fs
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept         0.86      0.21     0.45     1.27 1.00     2667     1747
hu_Intercept     -0.16      0.49    -1.16     0.78 1.00     3221     1753
age              -0.00      0.00    -0.01     0.00 1.00     2402     1451
male              0.01      0.07    -0.14     0.16 1.00     2595     1541
edu              -0.03      0.04    -0.11     0.06 1.00     2700     1566
pri_del_fs       -0.07      0.03    -0.13    -0.01 1.00     2824     1633
hu_age           -0.01      0.01    -0.02     0.00 1.00     3792     1538
hu_male          -0.01      0.18    -0.35     0.34 1.00     2368     1680
hu_edu           -0.22      0.10    -0.42    -0.02 1.00     2980     1603
hu_pri_del_fs     0.30      0.08     0.15     0.47 1.00     2610     1635

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape     3.20      0.28     2.68     3.76 1.00     2250     1631

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
    Outcome       term contrast estimate conf.low conf.high
1 self_disc        age    dY/dX  0.00125 -0.00398    0.0064
2 self_disc        edu    dY/dX  0.06443 -0.03394    0.1633
3 self_disc       male    1 - 0  0.00572 -0.15429    0.1765
4 self_disc pri_del_fs    dY/dX -0.15967 -0.23809   -0.0855
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-83-2.png){width=768}
:::
:::



#### Trust General



::: {.cell}

```{.r .cell-code}
run_hurdles(
  object = d, 
  name = "self_disc_trust",
  outcome = "self_disc", 
  predictor = "trust_gen_fs",
  y_lab = "Self Disclosure",
  x_lab = "Trust general"
  )
```

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-84-1.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: self_disc ~ 1 + age + male + edu + trust_gen_fs 
         hu ~ 1 + age + male + edu + trust_gen_fs
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
                Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept           0.48      0.28    -0.07     1.01 1.00     2497     1450
hu_Intercept        3.10      0.63     1.91     4.34 1.00     2553     1844
age                -0.00      0.00    -0.01     0.00 1.00     2189     1228
male                0.02      0.07    -0.12     0.16 1.00     1877     1514
edu                -0.02      0.04    -0.10     0.06 1.00     2440     1726
trust_gen_fs        0.02      0.04    -0.07     0.10 1.00     2795     1533
hu_age             -0.01      0.01    -0.02     0.00 1.00     5495     1588
hu_male            -0.06      0.17    -0.39     0.28 1.00     2702     1525
hu_edu             -0.22      0.10    -0.42    -0.01 1.00     3150     1521
hu_trust_gen_fs    -0.36      0.10    -0.55    -0.17 1.00     2654     1567

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape     3.15      0.29     2.63     3.73 1.00     2724     1739

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
    Outcome         term contrast estimate conf.low conf.high
1 self_disc          age    dY/dX  0.00315 -0.00214   0.00829
2 self_disc          edu    dY/dX  0.06645 -0.02674   0.15695
3 self_disc         male    1 - 0  0.03747 -0.12724   0.19493
4 self_disc trust_gen_fs    dY/dX  0.14314  0.06100   0.23590
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-84-2.png){width=768}
:::
:::



#### Self-Efficacy



::: {.cell}

```{.r .cell-code}
run_hurdles(
  object = d, 
  name = "self_disc_selfeff",
  outcome = "self_disc", 
  predictor = "self_eff_fs",
  y_lab = "Self Disclosure",
  x_lab = "Self-efficacy"
  )
```

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-85-1.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: self_disc ~ 1 + age + male + edu + self_eff_fs 
         hu ~ 1 + age + male + edu + self_eff_fs
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept         -0.33      0.28    -0.87     0.22 1.00     2103     1600
hu_Intercept       5.76      0.70     4.41     7.19 1.00     2102     1382
age               -0.00      0.00    -0.01     0.00 1.00     2039     1043
male               0.01      0.07    -0.13     0.16 1.00     1894     1565
edu               -0.04      0.05    -0.13     0.05 1.00     1777     1077
self_eff_fs        0.16      0.04     0.08     0.24 1.00     2051     1562
hu_age            -0.01      0.01    -0.02     0.00 1.00     2734     1502
hu_male           -0.01      0.19    -0.38     0.36 1.00     1994     1617
hu_edu            -0.12      0.11    -0.34     0.10 1.00     2484     1529
hu_self_eff_fs    -0.88      0.11    -1.10    -0.67 1.00     1854     1414

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape     3.33      0.29     2.77     3.91 1.00     1657     1498

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
    Outcome        term contrast estimate conf.low conf.high
1 self_disc         age    dY/dX  0.00278 -0.00196   0.00788
2 self_disc         edu    dY/dX  0.01524 -0.07957   0.10647
3 self_disc        male    1 - 0  0.01031 -0.14329   0.15823
4 self_disc self_eff_fs    dY/dX  0.39069  0.30710   0.47238
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-85-2.png){width=768}
:::
:::



The effects is markedly exponential. Let's inspect the changes across values in self-efficacy.



::: {.cell}

```{.r .cell-code}
effects_self_disc_selfeff <- 
  slopes(
    fit_self_disc_selfeff, 
    newdata = datagrid(
      self_eff_fs = c(1, 6, 7)
      ), 
   grid_type = "counterfactual"
    ) %>%
  filter(
    term == "self_eff_fs"
    ) %T>%
  print()
```

::: {.cell-output .cell-output-stdout}

```

 self_eff_fs Estimate   2.5 % 97.5 %
           1   0.0119 0.00497  0.027
           6   0.5063 0.37311  0.655
           7   0.5432 0.37997  0.765

Term: self_eff_fs
Type: response
Comparison: dY/dX
```


:::

```{.r .cell-code}
# let's extract effects for 1
effects_self_disc_selfeff_12 <- 
  effects_self_disc_selfeff %>% 
  filter(
    self_eff_fs == 1 &
      term == "self_eff_fs") %>% 
  select(estimate)

# let's extract effects for 7
effects_self_disc_selfeff_67 <- 
  effects_self_disc_selfeff %>% 
  filter(
    self_eff_fs == 6 &
      term == "self_eff_fs") %>% 
  select(estimate)
```
:::



Whereas a change in self-efficacy from 1 to 2 led to an increase of 0 words, a change from 6 to 7 led to an increase of 1 words.

Let's explore  how number of communicated words differ for people very much self-efficacious versus not self-efficacious at all.



::: {.cell}

```{.r .cell-code}
# let's extract effects for 1
self_disc_selfeff_1 <- 
  fig_self_disc_selfeff$data %>% 
  filter(self_eff_fs == min(self_eff_fs)) %>% 
  pull(estimate__) %>% 
  round(0)

# let's extract effects for 7
self_disc_selfeff_7 <- 
  fig_self_disc_selfeff$data %>% 
  filter(self_eff_fs == max(self_eff_fs)) %>% 
  pull(estimate__) %>% 
  round(0)
```
:::



People who reported no self-efficacy at all communicated on average 1 words, whereas people with very high self-efficacy communicated on average 254 words.

## Visualization



::: {.cell}

```{.r .cell-code}
fig_results <- 
  plot_grid(
    fig_self_disc_pricon, 
    fig_self_disc_grats, 
    fig_self_disc_pridel, 
    fig_self_disc_trust, 
    fig_self_disc_selfeff,
    fig_self_disc_ctrl
    ) %T>% 
  print()
```

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-88-1.png){width=768}
:::

```{.r .cell-code}
# Save visualizations
ggsave("figures/results/self_disc/effects.pdf")
ggsave("figures/results/self_disc/effects.png")
```
:::



## Table

Make table of slopes



::: {.cell}

```{.r .cell-code}
tab_slopes_self_disc <- 
  rbind(
    slopes_self_disc_pricon,
    slopes_self_disc_grats,
    slopes_self_disc_pridel,
    slopes_self_disc_trust,
    slopes_self_disc_selfeff,
    slopes_self_disc_ctrl,
    slopes_self_disc_lkdslk,
    slopes_pricon_ctrl,
    slopes_pricon_lkdslk,
    slopes_grats_ctrl,
    slopes_grats_lkdslk,
    slopes_pridel_ctrl,
    slopes_pridel_lkdslk,
    slopes_trust_ctrl,
    slopes_trust_lkdslk,
    slopes_selfeff_ctrl,
    slopes_selfeff_lkdslk
    ) %>% 
  as.data.frame() %>% 
  mutate(
    term = ifelse(
      .$contrast == "Like - Control", 
      "Like - Control",
      ifelse(
        .$contrast == "Like & dislike - Control", 
        "Like & dislike - Control",
        ifelse(
          .$contrast == "Like - Like & dislike", 
          "Like - Like & dislike", 
          .$term
          )
        )
      ),
    Predictor = recode(
      term, 
      "Like - Control" = "Like vs. control",
      "Like & dislike - Control" = "Like & dislike vs. control",
      "Like - Like & dislike" = "Like vs. like & dislike",
      "pri_con_fs" = "Privacy concerns",
      "pri_del_fs" = "Privacy deliberation",
      "grats_gen_fs" = "Expected gratifications",
      "trust_gen_fs" = "Trust",
      "self_eff_fs" = "Self-efficacy"
    ),
    Outcome = recode(
      Outcome,
      "pri_con_fs" = "Privacy concerns",
      "pri_del_fs" = "Privacy deliberation",
      "grats_gen_fs" = "Expected gratifications",
      "trust_gen_fs" = "Trust",
      "self_eff_fs" = "Self-efficacy"
    )
  ) %>% 
  filter(
    Predictor %in% c(
      "Like vs. control",
      "Like & dislike vs. control",
      "Like vs. like & dislike",
      "Privacy concerns",
      "Privacy deliberation",
      "Expected gratifications",
      "Trust",
      "Self-efficacy")
  ) %>% 
  mutate(
    Predictor = factor(Predictor, levels = c(
      "Self-efficacy",
      "Trust",
      "Privacy deliberation",
      "Expected gratifications",
      "Privacy concerns",
      "Like vs. like & dislike",
      "Like & dislike vs. control",
      "Like vs. control"
      )
    ),
    Estimate = 
      ifelse(
        Outcome == "words",
        round(estimate),
        estimate
      ),
    LL = 
      ifelse(
        Outcome == "words",
        round(conf.low),
        conf.low
      ),
    UL = 
      ifelse(
        Outcome == "words",
        round(conf.high),
        conf.high
      ),
  ) %>% 
  select(
    Outcome,
    Predictor, 
    Estimate,
    LL, 
    UL
    ) 

tab_slopes_self_disc %>% 
  kable() %>% 
  kable_styling("striped") %>% 
  scroll_box(width = "100%")
```

::: {.cell-output-display}
`````{=html}
<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; "><table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Outcome </th>
   <th style="text-align:left;"> Predictor </th>
   <th style="text-align:right;"> Estimate </th>
   <th style="text-align:right;"> LL </th>
   <th style="text-align:right;"> UL </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> self_disc </td>
   <td style="text-align:left;"> Privacy concerns </td>
   <td style="text-align:right;"> -0.082 </td>
   <td style="text-align:right;"> -0.138 </td>
   <td style="text-align:right;"> -0.029 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> self_disc </td>
   <td style="text-align:left;"> Expected gratifications </td>
   <td style="text-align:right;"> 0.070 </td>
   <td style="text-align:right;"> -0.010 </td>
   <td style="text-align:right;"> 0.141 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> self_disc </td>
   <td style="text-align:left;"> Privacy deliberation </td>
   <td style="text-align:right;"> -0.160 </td>
   <td style="text-align:right;"> -0.238 </td>
   <td style="text-align:right;"> -0.085 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> self_disc </td>
   <td style="text-align:left;"> Trust </td>
   <td style="text-align:right;"> 0.143 </td>
   <td style="text-align:right;"> 0.061 </td>
   <td style="text-align:right;"> 0.236 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> self_disc </td>
   <td style="text-align:left;"> Self-efficacy </td>
   <td style="text-align:right;"> 0.391 </td>
   <td style="text-align:right;"> 0.307 </td>
   <td style="text-align:right;"> 0.472 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> self_disc </td>
   <td style="text-align:left;"> Like &amp; dislike vs. control </td>
   <td style="text-align:right;"> -0.165 </td>
   <td style="text-align:right;"> -0.373 </td>
   <td style="text-align:right;"> 0.044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> self_disc </td>
   <td style="text-align:left;"> Like vs. control </td>
   <td style="text-align:right;"> -0.052 </td>
   <td style="text-align:right;"> -0.267 </td>
   <td style="text-align:right;"> 0.149 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> self_disc </td>
   <td style="text-align:left;"> Like vs. like &amp; dislike </td>
   <td style="text-align:right;"> 0.110 </td>
   <td style="text-align:right;"> -0.088 </td>
   <td style="text-align:right;"> 0.312 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Privacy concerns </td>
   <td style="text-align:left;"> Like &amp; dislike vs. control </td>
   <td style="text-align:right;"> 0.189 </td>
   <td style="text-align:right;"> -0.152 </td>
   <td style="text-align:right;"> 0.534 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Privacy concerns </td>
   <td style="text-align:left;"> Like vs. control </td>
   <td style="text-align:right;"> 0.035 </td>
   <td style="text-align:right;"> -0.294 </td>
   <td style="text-align:right;"> 0.360 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Privacy concerns </td>
   <td style="text-align:left;"> Like vs. like &amp; dislike </td>
   <td style="text-align:right;"> -0.150 </td>
   <td style="text-align:right;"> -0.516 </td>
   <td style="text-align:right;"> 0.190 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Expected gratifications </td>
   <td style="text-align:left;"> Like &amp; dislike vs. control </td>
   <td style="text-align:right;"> -0.208 </td>
   <td style="text-align:right;"> -0.477 </td>
   <td style="text-align:right;"> 0.054 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Expected gratifications </td>
   <td style="text-align:left;"> Like vs. control </td>
   <td style="text-align:right;"> -0.174 </td>
   <td style="text-align:right;"> -0.428 </td>
   <td style="text-align:right;"> 0.083 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Expected gratifications </td>
   <td style="text-align:left;"> Like vs. like &amp; dislike </td>
   <td style="text-align:right;"> 0.034 </td>
   <td style="text-align:right;"> -0.244 </td>
   <td style="text-align:right;"> 0.319 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Privacy deliberation </td>
   <td style="text-align:left;"> Like &amp; dislike vs. control </td>
   <td style="text-align:right;"> -0.045 </td>
   <td style="text-align:right;"> -0.278 </td>
   <td style="text-align:right;"> 0.191 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Privacy deliberation </td>
   <td style="text-align:left;"> Like vs. control </td>
   <td style="text-align:right;"> -0.008 </td>
   <td style="text-align:right;"> -0.246 </td>
   <td style="text-align:right;"> 0.236 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Privacy deliberation </td>
   <td style="text-align:left;"> Like vs. like &amp; dislike </td>
   <td style="text-align:right;"> 0.037 </td>
   <td style="text-align:right;"> -0.211 </td>
   <td style="text-align:right;"> 0.292 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Trust </td>
   <td style="text-align:left;"> Like &amp; dislike vs. control </td>
   <td style="text-align:right;"> -0.125 </td>
   <td style="text-align:right;"> -0.331 </td>
   <td style="text-align:right;"> 0.076 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Trust </td>
   <td style="text-align:left;"> Like vs. control </td>
   <td style="text-align:right;"> -0.147 </td>
   <td style="text-align:right;"> -0.353 </td>
   <td style="text-align:right;"> 0.059 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Trust </td>
   <td style="text-align:left;"> Like vs. like &amp; dislike </td>
   <td style="text-align:right;"> -0.023 </td>
   <td style="text-align:right;"> -0.229 </td>
   <td style="text-align:right;"> 0.178 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Self-efficacy </td>
   <td style="text-align:left;"> Like &amp; dislike vs. control </td>
   <td style="text-align:right;"> -0.092 </td>
   <td style="text-align:right;"> -0.299 </td>
   <td style="text-align:right;"> 0.128 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Self-efficacy </td>
   <td style="text-align:left;"> Like vs. control </td>
   <td style="text-align:right;"> -0.088 </td>
   <td style="text-align:right;"> -0.302 </td>
   <td style="text-align:right;"> 0.107 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Self-efficacy </td>
   <td style="text-align:left;"> Like vs. like &amp; dislike </td>
   <td style="text-align:right;"> 0.002 </td>
   <td style="text-align:right;"> -0.209 </td>
   <td style="text-align:right;"> 0.212 </td>
  </tr>
</tbody>
</table></div>

`````
:::
:::





### Self Disclosure without quantity



::: {.cell}

```{.r .cell-code}
# self-disclosure without quantity: breadth + depth only (both 0-3, so no standardization needed)
d <- d %>% 
  mutate(self_disc2 = breadth_mean + depth_mean)   # 0-6, non-commenters = 0

# check: commenters whose comments contain no topic and no disclosure also get 0
d %>% 
  summarise(
    non_commenters        = sum(words == 0),
    commenters            = sum(words > 0),
    commenters_with_zero  = sum(words > 0 & self_disc2 == 0)
  )
```

::: {.cell-output .cell-output-stdout}

```
  non_commenters commenters commenters_with_zero
1            324        235                   35
```


:::

```{.r .cell-code}
summary(d$self_disc2[d$self_disc2 > 0])
```

::: {.cell-output .cell-output-stdout}

```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   0.50    1.67    2.00    2.38    3.00    6.00 
```


:::
:::

::: {.cell}

```{.r .cell-code}
# hurdle models for self_disc2, same specification as above
predictors_sd2 <- c(
  pricon  = "pri_con_fs",
  grats   = "grats_gen_fs",
  pridel  = "pri_del_fs",
  trust   = "trust_gen_fs",
  selfeff = "self_eff_fs",
  ctrl    = "version",
  lkdslk  = "version_lkdslk"
)
x_labs_sd2 <- c(
  pricon = "Privacy concerns", grats = "Gratifications general", pridel = "Privacy deliberation",
  trust = "Trust general", selfeff = "Self-efficacy",
  ctrl = "Experimental conditions", lkdslk = "Experimental conditions"
)

for (k in names(predictors_sd2)) {
  run_hurdles(
    object    = d,
    name      = paste0("self_disc2_", k),
    outcome   = "self_disc2",
    predictor = predictors_sd2[[k]],
    y_lab     = "Self-disclosure (breadth + depth)",
    x_lab     = x_labs_sd2[[k]]
  )
}
```

::: {.cell-output .cell-output-stdout}

```
Some PIT values larger than 1! Largest:  1 
Rounding PIT > 1 to 1.
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-91-1.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: self_disc2 ~ 1 + age + male + edu + pri_con_fs 
         hu ~ 1 + age + male + edu + pri_con_fs
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept         0.92      0.15     0.63     1.23 1.00     2424     1823
hu_Intercept      0.71      0.40    -0.05     1.49 1.00     2690     1457
age               0.00      0.00    -0.00     0.01 1.00     2240     1481
male              0.04      0.07    -0.10     0.17 1.01     1743     1435
edu              -0.06      0.04    -0.13     0.02 1.00     2201     1583
pri_con_fs        0.00      0.02    -0.04     0.05 1.00     2206     1421
hu_age           -0.00      0.01    -0.02     0.01 1.01     3113     1625
hu_male          -0.06      0.18    -0.40     0.29 1.00     2218     1671
hu_edu           -0.26      0.11    -0.46    -0.05 1.00     2330     1531
hu_pri_con_fs     0.19      0.06     0.08     0.31 1.00     2123     1695

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape     4.67      0.45     3.81     5.57 1.00     1953     1486

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
     Outcome       term contrast estimate conf.low conf.high
1 self_disc2        age    dY/dX  0.00321 -0.00374    0.0108
2 self_disc2        edu    dY/dX  0.08386 -0.04750    0.2074
3 self_disc2       male    1 - 0  0.06212 -0.15850    0.2854
4 self_disc2 pri_con_fs    dY/dX -0.10036 -0.17067   -0.0292
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-91-2.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
Some PIT values larger than 1! Largest:  1 
Rounding PIT > 1 to 1.
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-91-3.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: self_disc2 ~ 1 + age + male + edu + grats_gen_fs 
         hu ~ 1 + age + male + edu + grats_gen_fs
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
                Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept           1.20      0.22     0.75     1.62 1.00     2326     1761
hu_Intercept        2.34      0.54     1.31     3.40 1.00     2089     1755
age                 0.00      0.00    -0.00     0.00 1.00     2143     1545
male                0.03      0.07    -0.10     0.16 1.00     2211     1543
edu                -0.06      0.04    -0.14     0.01 1.00     2471     1603
grats_gen_fs       -0.05      0.03    -0.11     0.01 1.00     2557     1383
hu_age             -0.01      0.01    -0.02     0.01 1.00     3151     1597
hu_male            -0.10      0.18    -0.46     0.25 1.00     1884     1397
hu_edu             -0.25      0.11    -0.46    -0.04 1.01     2169     1604
hu_grats_gen_fs    -0.20      0.08    -0.37    -0.05 1.00     2173     1530

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape     4.73      0.47     3.87     5.71 1.00     2532     1378

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
     Outcome         term contrast estimate conf.low conf.high
1 self_disc2          age    dY/dX   0.0035 -0.00375    0.0109
2 self_disc2          edu    dY/dX   0.0815 -0.05316    0.2079
3 self_disc2 grats_gen_fs    dY/dX   0.0642 -0.03422    0.1654
4 self_disc2         male    1 - 0   0.0783 -0.13467    0.3037
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-91-4.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
Some PIT values larger than 1! Largest:  1 
Rounding PIT > 1 to 1.
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-91-5.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: self_disc2 ~ 1 + age + male + edu + pri_del_fs 
         hu ~ 1 + age + male + edu + pri_del_fs
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept         1.11      0.18     0.76     1.45 1.00     2185     1777
hu_Intercept     -0.06      0.51    -1.08     0.93 1.00     2511     1354
age               0.00      0.00    -0.00     0.00 1.00     2132     1596
male              0.03      0.07    -0.11     0.15 1.00     2400     1305
edu              -0.06      0.04    -0.14     0.01 1.00     2258     1609
pri_del_fs       -0.04      0.03    -0.10     0.01 1.00     2463     1749
hu_age           -0.00      0.01    -0.01     0.01 1.00     3113     1750
hu_male          -0.03      0.18    -0.39     0.33 1.00     2183     1607
hu_edu           -0.24      0.11    -0.46    -0.02 1.00     2262     1519
hu_pri_del_fs     0.33      0.08     0.16     0.49 1.00     2620     1614

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape     4.73      0.45     3.89     5.68 1.00     2665     1670

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
     Outcome       term contrast estimate conf.low conf.high
1 self_disc2        age    dY/dX  0.00216 -0.00497    0.0094
2 self_disc2        edu    dY/dX  0.07241 -0.05788    0.2005
3 self_disc2       male    1 - 0  0.04104 -0.19043    0.2591
4 self_disc2 pri_del_fs    dY/dX -0.20838 -0.30891   -0.1088
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-91-6.png){width=768}
:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-91-7.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: self_disc2 ~ 1 + age + male + edu + trust_gen_fs 
         hu ~ 1 + age + male + edu + trust_gen_fs
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
                Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept           1.24      0.25     0.75     1.75 1.00     1921     1521
hu_Intercept        3.48      0.69     2.16     4.84 1.00     1943     1495
age                 0.00      0.00    -0.00     0.01 1.00     2035     1525
male                0.03      0.07    -0.11     0.16 1.00     2267     1453
edu                -0.06      0.04    -0.14     0.02 1.00     2371     1575
trust_gen_fs       -0.05      0.04    -0.13     0.02 1.00     2107     1597
hu_age             -0.01      0.01    -0.02     0.00 1.00     3587     1410
hu_male            -0.09      0.18    -0.45     0.28 1.00     2087     1488
hu_edu             -0.24      0.11    -0.45    -0.03 1.00     2226     1354
hu_trust_gen_fs    -0.40      0.11    -0.61    -0.19 1.00     2193     1421

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape     4.71      0.47     3.87     5.66 1.00     2381     1214

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
     Outcome         term contrast estimate conf.low conf.high
1 self_disc2          age    dY/dX  0.00411 -0.00276    0.0112
2 self_disc2          edu    dY/dX  0.07736 -0.05778    0.1999
3 self_disc2         male    1 - 0  0.06893 -0.14064    0.2990
4 self_disc2 trust_gen_fs    dY/dX  0.16156  0.04553    0.2776
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-91-8.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
Some PIT values larger than 1! Largest:  1 
Rounding PIT > 1 to 1.
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-91-9.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: self_disc2 ~ 1 + age + male + edu + self_eff_fs 
         hu ~ 1 + age + male + edu + self_eff_fs
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept          0.58      0.26     0.07     1.08 1.00     1999     1855
hu_Intercept       6.75      0.76     5.29     8.26 1.00     1989     1468
age                0.00      0.00    -0.00     0.01 1.00     2446     1462
male               0.03      0.07    -0.10     0.16 1.00     2108     1266
edu               -0.07      0.04    -0.15     0.01 1.00     2302     1609
self_eff_fs        0.06      0.04    -0.01     0.14 1.00     1954     1389
hu_age            -0.01      0.01    -0.02     0.00 1.00     3173     1314
hu_male           -0.03      0.20    -0.41     0.35 1.00     1951     1788
hu_edu            -0.15      0.12    -0.37     0.08 1.00     2346     1071
hu_self_eff_fs    -1.02      0.12    -1.26    -0.79 1.00     1897     1559

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape     4.73      0.46     3.87     5.73 1.00     2255     1708

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
     Outcome        term contrast estimate conf.low conf.high
1 self_disc2         age    dY/dX  0.00415 -0.00234    0.0106
2 self_disc2         edu    dY/dX  0.01191 -0.11276    0.1309
3 self_disc2        male    1 - 0  0.04111 -0.17000    0.2501
4 self_disc2 self_eff_fs    dY/dX  0.50897  0.40746    0.6148
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-91-10.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
Some PIT values larger than 1! Largest:  1 
Rounding PIT > 1 to 1.
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-91-11.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: self_disc2 ~ 1 + age + male + edu + version 
         hu ~ 1 + age + male + edu + version
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
                       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept                  1.00      0.14     0.71     1.28 1.00     2136     1777
hu_Intercept               1.26      0.38     0.50     2.01 1.00     2347     1524
age                        0.00      0.00    -0.00     0.00 1.00     2089     1380
male                       0.03      0.07    -0.10     0.17 1.00     2210     1394
edu                       -0.06      0.04    -0.14     0.02 1.00     2349     1316
versionLike               -0.07      0.08    -0.22     0.08 1.00     1746     1630
versionLike&dislike       -0.14      0.08    -0.30     0.02 1.00     1688     1392
hu_age                    -0.01      0.01    -0.02     0.01 1.00     3662     1545
hu_male                   -0.08      0.18    -0.42     0.26 1.00     1974     1564
hu_edu                    -0.23      0.11    -0.44    -0.02 1.00     2171     1520
hu_versionLike             0.01      0.21    -0.40     0.43 1.00     1978     1620
hu_versionLike&dislike     0.18      0.22    -0.28     0.61 1.00     1832     1418

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape     4.70      0.46     3.86     5.68 1.00     2418     1547

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
     Outcome    term                 contrast estimate conf.low conf.high
1 self_disc2     age                    dY/dX  0.00348 -0.00386    0.0109
2 self_disc2     edu                    dY/dX  0.07573 -0.06291    0.2021
3 self_disc2    male                    1 - 0  0.06968 -0.14941    0.2894
4 self_disc2 version Like & dislike - Control -0.20994 -0.47815    0.0609
5 self_disc2 version           Like - Control -0.06993 -0.34885    0.1914
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-91-12.png){width=768}
:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-91-13.png){width=768}
:::

::: {.cell-output .cell-output-stdout}

```
 Family: hurdle_gamma 
  Links: mu = log; hu = logit 
Formula: self_disc2 ~ 1 + age + male + edu + version_lkdslk 
         hu ~ 1 + age + male + edu + version_lkdslk
   Data: object (Number of observations: 558) 
  Draws: 4 chains, each with iter = 1000; warmup = 500; thin = 1;
         total post-warmup draws = 2000

Regression Coefficients:
                         Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept                    0.86      0.14     0.57     1.13 1.00     2151     1850
hu_Intercept                 1.44      0.40     0.68     2.21 1.00     2657     1800
age                          0.00      0.00    -0.00     0.00 1.01     2297     1697
male                         0.03      0.07    -0.11     0.17 1.00     1731     1022
edu                         -0.06      0.04    -0.14     0.03 1.00     2479     1318
version_lkdslkLike           0.07      0.08    -0.09     0.23 1.00     1738     1500
version_lkdslkControl        0.14      0.08    -0.02     0.30 1.00     1740     1637
hu_age                      -0.01      0.01    -0.02     0.01 1.00     3393     1595
hu_male                     -0.08      0.18    -0.46     0.26 1.00     1933     1500
hu_edu                      -0.23      0.11    -0.45    -0.01 1.00     2203     1372
hu_version_lkdslkLike       -0.16      0.21    -0.58     0.27 1.00     1424     1089
hu_version_lkdslkControl    -0.17      0.22    -0.62     0.26 1.00     1518     1289

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape     4.70      0.47     3.81     5.63 1.00     2185     1377

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```


:::

::: {.cell-output .cell-output-stdout}

```

Results of marginal effects:
     Outcome           term                 contrast estimate conf.low conf.high
1 self_disc2            age                    dY/dX  0.00363 -0.00374    0.0108
2 self_disc2            edu                    dY/dX  0.07821 -0.05972    0.2048
3 self_disc2           male                    1 - 0  0.06562 -0.13917    0.3048
4 self_disc2 version_lkdslk Control - Like & dislike  0.20804 -0.05327    0.4731
5 self_disc2 version_lkdslk    Like - Like & dislike  0.14144 -0.12441    0.3896
```


:::

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-91-14.png){width=768}
:::
:::

::: {.cell}

```{.r .cell-code}
model_labels <- c(
  pricon = "Privacy concerns", grats = "Expected gratifications", pridel = "Privacy deliberation",
  trust = "Trust", selfeff = "Self-efficacy", ctrl = "Conditions", lkdslk = "Conditions"
)

tab_slopes_self_disc2 <- 
  map_dfr(names(predictors_sd2), function(k) {
    get(paste0("slopes_self_disc2_", k)) %>% 
      as_tibble() %>% 
      filter(term == predictors_sd2[[k]]) %>%                       # only the predictor, no controls
      mutate(Predictor = if_else(contrast == "dY/dX", model_labels[[k]], contrast))
  }) %>% 
  distinct(Predictor, .keep_all = TRUE) %>% 
  mutate(sig = if_else(conf.low > 0 | conf.high < 0, "*", "")) %>% 
  select(Predictor, Estimate = estimate, LL = conf.low, UL = conf.high, sig)

tab_slopes_self_disc2 %>% 
  filter(Predictor != "Control - Like & dislike") %>%   # duplicate of "Like & dislike - Control"
  kable(digits = 2) %>% 
  kable_styling("striped")
```

::: {.cell-output-display}
`````{=html}
<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Predictor </th>
   <th style="text-align:right;"> Estimate </th>
   <th style="text-align:right;"> LL </th>
   <th style="text-align:right;"> UL </th>
   <th style="text-align:left;"> sig </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Privacy concerns </td>
   <td style="text-align:right;"> -0.10 </td>
   <td style="text-align:right;"> -0.17 </td>
   <td style="text-align:right;"> -0.03 </td>
   <td style="text-align:left;"> * </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Expected gratifications </td>
   <td style="text-align:right;"> 0.06 </td>
   <td style="text-align:right;"> -0.03 </td>
   <td style="text-align:right;"> 0.17 </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Privacy deliberation </td>
   <td style="text-align:right;"> -0.21 </td>
   <td style="text-align:right;"> -0.31 </td>
   <td style="text-align:right;"> -0.11 </td>
   <td style="text-align:left;"> * </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Trust </td>
   <td style="text-align:right;"> 0.16 </td>
   <td style="text-align:right;"> 0.05 </td>
   <td style="text-align:right;"> 0.28 </td>
   <td style="text-align:left;"> * </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Self-efficacy </td>
   <td style="text-align:right;"> 0.51 </td>
   <td style="text-align:right;"> 0.41 </td>
   <td style="text-align:right;"> 0.61 </td>
   <td style="text-align:left;"> * </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Like &amp; dislike - Control </td>
   <td style="text-align:right;"> -0.21 </td>
   <td style="text-align:right;"> -0.48 </td>
   <td style="text-align:right;"> 0.06 </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Like - Control </td>
   <td style="text-align:right;"> -0.07 </td>
   <td style="text-align:right;"> -0.35 </td>
   <td style="text-align:right;"> 0.19 </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Like - Like &amp; dislike </td>
   <td style="text-align:right;"> 0.14 </td>
   <td style="text-align:right;"> -0.12 </td>
   <td style="text-align:right;"> 0.39 </td>
   <td style="text-align:left;">  </td>
  </tr>
</tbody>
</table>

`````
:::
:::

::: {.cell}

```{.r .cell-code}
# decompose the effects: hurdle part vs. conditional (gamma) part
get_parts <- function(k) {
  fit  <- get(paste0("fit_self_disc2_", k))
  pred <- predictors_sd2[[k]]
  part <- function(dpar, label) {
    avg_slopes(fit, variables = pred, dpar = dpar) %>%
      as_tibble() %>%
      mutate(part = label,
             Predictor = if_else(contrast == "dY/dX", model_labels[[k]], contrast))
  }
  bind_rows(
    avg_slopes(fit, variables = pred) %>% as_tibble() %>%
      mutate(part = "Combined",
             Predictor = if_else(contrast == "dY/dX", model_labels[[k]], contrast)),
    part("hu", "Hurdle") %>%                                   # hu = P(zero) -> flip to P(disclosing)
      mutate(across(c(estimate, conf.low, conf.high), ~ -.x)) %>%
      rename(conf.low = conf.high, conf.high = conf.low),
    part("mu", "Conditional")
  )
}

tab_parts_sd2 <- 
  map_dfr(c("pricon", "grats", "pridel", "trust", "selfeff", "ctrl", "lkdslk"), get_parts) %>%
  filter(Predictor != "Control - Like & dislike") %>%          # duplicate of Like & dislike - Control
  distinct(Predictor, part, .keep_all = TRUE) %>%
  mutate(
    sig  = if_else(conf.low > 0 | conf.high < 0, "*", ""),
    cell = sprintf("%.2f [%.2f, %.2f]%s", estimate, conf.low, conf.high, sig)
  ) %>%
  select(Predictor, part, cell) %>%
  pivot_wider(names_from = part, values_from = cell) %>%
  select(Predictor, Combined, Hurdle, Conditional)

tab_parts_sd2 %>%
  kable(caption = "Average marginal effects on self-disclosure (breadth + depth). Hurdle = change in probability of disclosing anything; Conditional = change in self-disclosure among those who disclose. * = 95% interval excludes zero.") %>%
  kable_styling("striped")
```

::: {.cell-output-display}
`````{=html}
<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<caption>Average marginal effects on self-disclosure (breadth + depth). Hurdle = change in probability of disclosing anything; Conditional = change in self-disclosure among those who disclose. * = 95% interval excludes zero.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Predictor </th>
   <th style="text-align:left;"> Combined </th>
   <th style="text-align:left;"> Hurdle </th>
   <th style="text-align:left;"> Conditional </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Privacy concerns </td>
   <td style="text-align:left;"> -0.10 [-0.17, -0.03]* </td>
   <td style="text-align:left;"> -0.04 [-0.07, -0.02]* </td>
   <td style="text-align:left;"> 0.00 [-0.10, 0.12] </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Expected gratifications </td>
   <td style="text-align:left;"> 0.06 [-0.03, 0.17] </td>
   <td style="text-align:left;"> 0.05 [0.01, 0.08]* </td>
   <td style="text-align:left;"> -0.13 [-0.28, 0.03] </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Privacy deliberation </td>
   <td style="text-align:left;"> -0.21 [-0.31, -0.11]* </td>
   <td style="text-align:left;"> -0.07 [-0.10, -0.04]* </td>
   <td style="text-align:left;"> -0.10 [-0.25, 0.03] </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Trust </td>
   <td style="text-align:left;"> 0.16 [0.05, 0.28]* </td>
   <td style="text-align:left;"> 0.09 [0.04, 0.13]* </td>
   <td style="text-align:left;"> -0.13 [-0.32, 0.04] </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Self-efficacy </td>
   <td style="text-align:left;"> 0.51 [0.41, 0.61]* </td>
   <td style="text-align:left;"> 0.19 [0.16, 0.22]* </td>
   <td style="text-align:left;"> 0.14 [-0.03, 0.31] </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Like &amp; dislike - Control </td>
   <td style="text-align:left;"> -0.21 [-0.48, 0.06] </td>
   <td style="text-align:left;"> -0.04 [-0.14, 0.06] </td>
   <td style="text-align:left;"> -0.35 [-0.72, 0.04] </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Like - Control </td>
   <td style="text-align:left;"> -0.07 [-0.35, 0.19] </td>
   <td style="text-align:left;"> -0.00 [-0.10, 0.09] </td>
   <td style="text-align:left;"> -0.18 [-0.54, 0.20] </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Like - Like &amp; dislike </td>
   <td style="text-align:left;"> 0.14 [-0.12, 0.39] </td>
   <td style="text-align:left;"> 0.03 [-0.06, 0.13] </td>
   <td style="text-align:left;"> 0.17 [-0.20, 0.54] </td>
  </tr>
</tbody>
</table>

`````
:::
:::

::: {.cell}

```{.r .cell-code}
# numeric version of the decomposition for plotting
parts_num_sd2 <- 
  map_dfr(c("pricon", "grats", "pridel", "trust", "selfeff", "ctrl", "lkdslk"), get_parts) %>%
  filter(Predictor != "Control - Like & dislike") %>%
  distinct(Predictor, part, .keep_all = TRUE) %>%
  mutate(
    part      = factor(part, levels = c("Combined", "Hurdle", "Conditional"),
                       labels = c("Combined", "Hurdle\n(probability of disclosing)", 
                                  "Conditional\n(among those who disclose)")),
    Predictor = factor(Predictor, levels = rev(c(
      "Privacy concerns", "Expected gratifications", "Privacy deliberation", "Trust", 
      "Self-efficacy", "Like - Control", "Like & dislike - Control", "Like - Like & dislike"))),
    sig       = conf.low > 0 | conf.high < 0
  )

fig_parts_sd2 <- 
  ggplot(parts_num_sd2, aes(x = estimate, y = Predictor)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high, shape = sig), size = .4) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  facet_wrap(~ part, scales = "free_x") +
  theme_bw() +
  labs(x = "Average marginal effect (95% credible interval)", y = NULL,
       caption = "Filled points: interval excludes zero. Outcome: self-disclosure (breadth + depth, 0-6).")

print(fig_parts_sd2)
```

::: {.cell-output-display}
![](05a_analyses_files/figure-html/unnamed-chunk-94-1.png){width=960}
:::

```{.r .cell-code}
ggsave("figures/results/self_disc/effects_sd2_parts.png", fig_parts_sd2, width = 10, height = 5)
```
:::
