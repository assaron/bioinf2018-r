## ----message=FALSE-------------------------------------------------------
# source("requirements.R")
source("functions.R")


## ------------------------------------------------------------------------
es <- read.gct("data/lgg.gct")

## ------------------------------------------------------------------------
es

## ------------------------------------------------------------------------
str(exprs(es))
exprs(es)[1:10, 1:2]

## ------------------------------------------------------------------------
str(fData(es))

## ------------------------------------------------------------------------
str(pData(es), list.len=5)

## ------------------------------------------------------------------------
clinical <- pData(es)
str(clinical, list.len=10)

## ---- message=FALSE------------------------------------------------------
library(magrittr)

mtcars %>%
  subset(hp > 100) %>%
  aggregate(. ~ cyl, data = ., FUN = . %>% mean %>% round(2)) %>%
  transform(kpl = mpg %>% multiply_by(0.4251))

## ---- message=FALSE------------------------------------------------------
library(dplyr)

starwars %>% 
    filter(species == "Droid") %>%
    select(name, ends_with("color"))

## ------------------------------------------------------------------------
table(clinical$sample_type)

## ------------------------------------------------------------------------
clinical %>% 
    filter(sample_type == "Primary solid Tumor") %>%
    str(list.len=5)

## ------------------------------------------------------------------------
clinical <- clinical %>%
    filter(sample_type == "Primary solid Tumor")
str(clinical, list.len=10)

## ------------------------------------------------------------------------
clinical %<>%
    filter(sample_type == "Primary solid Tumor")

## ------------------------------------------------------------------------
clinical %<>%
    mutate_at(vars(starts_with("days_")), funs(as.numeric))
str(clinical, list.len=10)

## ------------------------------------------------------------------------
table(clinical$vital_status)

## ------------------------------------------------------------------------

clinical %<>%
    mutate(years_to_event = case_when(
        vital_status == "alive" ~ days_to_last_followup / 365,
        vital_status == "dead"  ~ days_to_death / 365)) %>%
    mutate(event_status = vital_status == "dead")


## ------------------------------------------------------------------------
summary(clinical$years_to_event)
table(clinical$event_status)

## ------------------------------------------------------------------------
library(survival)
clinical %>%
    survfit(Surv(years_to_event, event_status) ~ 1, data=.)

## ----surv1, fig.height=4, fig.width=6,dev='svg', message=FALSE-----------
library(GGally)
clinical %>%
    survfit(Surv(years_to_event, event_status) ~ 1, data=.) %>%
    ggsurv()

## ----surv-gender, fig.height=3.5, fig.width=7,dev='svg'------------------
clinical %>%
    survfit(Surv(years_to_event, event_status) ~ gender, data=.) %>%
    ggsurv()

## ----surv-histology, fig.height=3.5, fig.width=7,dev='svg'---------------
clinical %>%
    survfit(Surv(years_to_event, event_status) ~ histological_type, data=.) %>%
    ggsurv()

## ------------------------------------------------------------------------
clinical %>%
    survdiff(Surv(years_to_event, event_status) ~ histological_type, data=.)

## ------------------------------------------------------------------------
clinical %>%
    survdiff(Surv(years_to_event, event_status) ~ (histological_type == "astrocytoma"),
             data=.)

## ----gdcRes-tip,cache=TRUE,cache.path="survival_cache/"------------------
gdcRes <- GenomicDataCommons::cases() %>%
    GenomicDataCommons::filter( ~ project.project_id == 'TCGA-LGG') %>%
    GenomicDataCommons::select(c("diagnoses.days_to_last_follow_up",
             "diagnoses.days_to_death",
             "diagnoses.morphology",
             "diagnoses.vital_status")) %>%
    GenomicDataCommons::results()
str(gdcRes, list.len = 3)

## ----gdcRes-full,cache=TRUE,cache.path="survival_cache/"-----------------
gdcRes <- GenomicDataCommons::cases() %>%
    GenomicDataCommons::filter( ~ project.project_id == 'TCGA-LGG') %>%
    GenomicDataCommons::select(c("diagnoses.days_to_last_follow_up",
             "diagnoses.days_to_death",
             "diagnoses.morphology",
             "diagnoses.vital_status")) %>%
    GenomicDataCommons::results_all()
str(gdcRes, list.len = 3)

## ------------------------------------------------------------------------
gdcClinical <- bind_rows(gdcRes$diagnoses)
head(gdcClinical)

## ------------------------------------------------------------------------
library(readr)
icdo3 <- read_tsv("./data/icdo3_histology.tsv")
icdo3

## ------------------------------------------------------------------------
icdo3_histology <- icdo3 %>%
    select(morphology=`Histology/Behavior`, morphology_desc=`Histology/Behavior Description`) %>%
    unique
icdo3_histology

## ------------------------------------------------------------------------

gdcClinical <- gdcClinical %>%
    left_join(icdo3_histology, by=c("morphology"))
head(gdcClinical)

## ------------------------------------------------------------------------
table(gdcClinical$vital_status)
gdcClinical %<>%
    mutate(years_to_event = case_when(
        vital_status == "alive" ~ days_to_last_follow_up / 365,
        vital_status == "dead"  ~ days_to_death / 365)) %>%
    mutate(event_status = vital_status == "dead")

## ----surv-morphology, fig.height=3.5, fig.width=7,dev='svg'--------------
gdcClinical %>%
    survfit(Surv(years_to_event, event_status) ~ morphology_desc, data=.) %>%
    ggsurv()

