![R](https://img.shields.io/badge/R-4.x-276DC3?style=flat&logo=r&logoColor=white)
![Python](https://img.shields.io/badge/python-3.10+-blue)
![Quarto](https://img.shields.io/badge/Quarto-1.9-39729E?style=flat)

# Replication materials: *"Likes or Dislikes, Gratifications or Concerns?"*

This repository contains the replication materials for the paper "Likes or Dislikes, Gratifications or Concerns? A One-Week Field Experiment Analyzing the Effects of Popularity Cues on the Privacy Calculus in Online Communication".

> **Design:** Preregistered one-week field experiment (*N* = 590). Participants discussed a political topic on an online platform with like buttons, like and dislike buttons, or no buttons.

All analyses are written as Quarto/R Markdown documents and are also rendered as a website (`docs/`). The data are included in `data/`. The preregistration is on [OSF](https://osf.io/a6tzc/?view_only=5d0ef9fe5e1745878cd1b19273cdf859).

## Reproduce

Install [R](https://www.r-project.org/), [Quarto](https://quarto.org/), and Python 3.10+. Then:

```bash
# Python (only for the 04*.ipynb notebooks)
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

```r
# R packages
install.packages(c("brms", "english", "kableExtra", "knitr", "lavaan",
                   "readr", "readxl", "remotes", "tidycomm", "tidyverse",
                   "writexl"))
remotes::install_github("crsh/papaja")
remotes::install_github("tdienlin/td")
```

```bash
quarto render                                              # whole site
quarto render 01b_manuscript_rev.qmd --to apaquarto-pdf    # manuscript PDF (run alone, last)
```

The LLM classification (`04b`, `04c`) used a local model via [Ollama](https://ollama.com/); saved outputs are in the notebooks, so re-running is optional.

## Workflow

```
03_coding_manual_icr.qmd      ──►  manual coding + intercoder reliability
04a – 04c (notebooks)         ──►  classification of comments
05a – 05c (analyses)          ──►  results reported in the paper, additional, revision
01b_manuscript_rev.qmd        ──►  manuscript (HTML + APA PDF)
```

## Structure

| File or folder | Content |
|---|---|
| `01a_manuscript_original.rmd`, `01b_manuscript_rev.qmd` | Manuscript (first submission, revision) |
| `02_items.rmd` | Questionnaire items |
| `03_coding_manual_icr.qmd` | Comment coding and intercoder reliability |
| `04a`–`04c_*.ipynb` | Comment classification (pretrained, LLM, final) |
| `05a`–`05c_*.qmd` | Analyses |
| `08_preregistration_changes.rmd` | Deviations from the preregistration |
| `appendix_*.qmd` | Appendices (codebook, model comparison, preregistration note, system prompt, website) |
| `data/` | Data and coded comments |
| `figures/`, `images/`, `codebook/`, `bibliography/` | Figures, study materials, codebook PDF, references |
| `docs/` | Rendered website |
| `_quarto.yml`, `styles.css`, `_extensions/` | Site configuration and APA template |
