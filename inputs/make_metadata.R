library(readr)
library(dplyr)
library(janitor)

setwd("~/github/2022-txome-minus-gnome/")

doryteuthis <- read_csv("inputs/make_metadata/Doryteuthis-pealeii-sra-run-info.csv") %>%
  clean_names() %>%
  filter(bio_project == "PRJNA641326") %>%
  mutate(sample_description = tolower(sample_name),
         sample_description = gsub("doryteuthis pealeii ", "", sample_description),
         sample_description = gsub("- ", "", sample_description),
         sample_description = gsub(" ", "_", sample_description)) %>%
  mutate(gx_accession = "GCA_023376005.1")

mus_sample_description <- read_csv("inputs/make_metadata/Mus-musculus-sample-descriptions.csv")
mus <- read_csv("inputs/make_metadata/Mus-musculus-adar-sra-run-info.csv") %>%
  clean_names() %>%
  left_join(mus_sample_description) %>%
  mutate(gx_accession = "GCA_000001635.9")

metadata <- bind_rows(mus, doryteuthis) %>%
  select(tx_run_accession = run, scientific_name, tx_study_accession = bio_project, sample_name, sample_description, gx_accession) 

write_tsv(metadata, "inputs/metadata.tsv")
