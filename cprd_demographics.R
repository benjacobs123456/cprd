###########################
# load packages
###########################

library(tidyverse)
library(gridExtra)
library(Hmisc)
library(survival)
library(survminer)
library(MASS)

setwd("/data/Wolfson-UKBB-Dobson/CPRD/bj/outputs/")

###########################
# define functions 
###########################

source("../scripts/cprd_functions.R")

###########################
# read in cleaned data 
###########################

cohorts = readRDS("./datasets/cprd_aurum_cleaned_exposures_patient.rds")

#################################
# descriptive stats
#################################

# now make demographics tables
make_demographics(input_data = cohorts$primary, outfile = "whole_cohort")
make_demographics(input_data = cohorts$hes_ms, outfile = "hes_ms")
make_demographics(input_data = cohorts$post_97, outfile = "post_97")
make_demographics(input_data = cohorts$post_06, outfile = "post_06")
make_demographics(input_data = cohorts$ten_y_antecedent, outfile = "ten_year")

