library(readr)
library(dplyr)
library(janitor)
library(purrr)

setwd("~/github/2022-tx-not-in-gx/")

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

octopus <- Sys.glob("inputs/make_metadata/Octopus-bimaculoides-PRJNA*-sra-run-info.csv") %>%
  map_dfr(read_csv) %>%
  clean_names() %>%
  mutate(gx_accession = "GCF_001194135.1",
         library_name = gsub("Obimac_", "", library_name),
         library_name = gsub("Sample_CR-19_Ova", "ova", library_name),
         library_name = gsub("GSM6213764", "gut", library_name),
         library_name = gsub("GSM6213768", "testis", library_name),
         library_name = gsub("GSM6213767", "skin", library_name),
         library_name = gsub("GSM6213766", "muscle", library_name),
         library_name = gsub("GSM6213765", "kidney", library_name),
         sample_description = library_name)

metadata <- bind_rows(mus, doryteuthis, octopus) %>%
  select(tx_run_accession = run, scientific_name, tx_study_accession = bio_project, sample_name, sample_description, gx_accession) 

write_tsv(metadata, "inputs/metadata.tsv")
