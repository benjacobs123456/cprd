###########################
# load packages
###########################

library(tidyverse)
library(gridExtra)
library(Hmisc)
library(survival)
library(survminer)
library(MASS)
library(survey)

setwd("/data/Wolfson-UKBB-Dobson/CPRD/bj/outputs/")

###########################
# define functions 
###########################

source("../scripts/cprd_functions.R")

###########################
# read in cleaned data 
###########################

cohorts = readRDS("./datasets/cprd_aurum_cleaned_exposures_patient.rds")

vars =c("bmi_z","smoking_status","alcohol_status","ethnicity_pc_hes_hybrid","e2019_imd_5")

missingness_df = data.frame()
for(i in c(1:length(cohorts))){
  cohort = cohorts[[i]]
  cohort_name = names(cohorts)[i]
for(var in vars){
  tbl = get_prop_miss(cohort,x=var) %>% dplyr::select(1,6)
  tbl = tbl[tbl[,1]==T,]
  tbl = data.frame("variable" = var,"missingness" = tbl$pastable, "cohort" = cohort_name)
  missingness_df <<- bind_rows(missingness_df,tbl)
  }
}

missingness_df = missingness_df %>% pivot_wider(id_cols = variable,names_from = cohort, values_from= missingness)
write_csv(missingness_df,"./tables/missing_data.csv")

#################################
# colliding & ethnicity 
#################################

post_1997_cohort = cohorts$post_97

# define some useful new / transformed predictors for missingness

# find out who has nonmissing ethnicity & smoking data 
miss_predictors = c("gender",
                    "age_at_registration",
                    "numeric_year_of_reg",
                    "fu_time",
                    "yob",
                    "MS_status",
                    "ethnicity_pc_hes_hybrid",
                    "e2019_imd_5")

p1=check_missingness_predictors(input_data = post_1997_cohort,
                                missing_col = "bmi_category_simple",
                                miss_pred_cols = miss_predictors)+ggtitle("BMI")
p2=check_missingness_predictors(input_data = post_1997_cohort,
                                missing_col = "smoking_status",
                                miss_pred_cols = miss_predictors)+ggtitle("Smoking")
p3=check_missingness_predictors(input_data = post_1997_cohort,
                                missing_col = "alcohol_status",
                                miss_pred_cols = miss_predictors)+ggtitle("Alcohol")

miss_predictors2 = c("gender",
                    "age_at_registration",
                    "numeric_year_of_reg",
                    "fu_time",
                    "yob",
                    "MS_status",
                    "e2019_imd_5")

p4=check_missingness_predictors(input_data = post_1997_cohort,
                                missing_col = "ethnicity_pc_hes_hybrid",
                                miss_pred_cols = miss_predictors2)+ggtitle("Ethnicity")


plot_fx(grid.arrange(p1,p2,p3,p4,nrow=2),
        "./figs/missingness_plots.png",
        width=12,height=10)

# create vars for IPW model
post_1997_cohort$sq_age_at_reg = post_1997_cohort$age_at_registration^2
post_1997_cohort$cub_age_at_reg = post_1997_cohort$age_at_registration^3
post_1997_cohort$sqrt_yor = sqrt(post_1997_cohort$numeric_year_of_reg)
post_1997_cohort$sq_fu_time = post_1997_cohort$fu_time^2
post_1997_cohort$sq_yob = post_1997_cohort$yob^2
post_1997_cohort$sq_years_dat_pre_index = sqrt(post_1997_cohort$years_data_before_index)
post_1997_cohort$numeric_imd = as.numeric(post_1997_cohort$e2019_imd_5)

# create multimissing var 
post_1997_cohort = post_1997_cohort %>% 
  mutate(any_missing = ifelse(
    is.na(bmi_z) | is.na(smoking_status) | is.na(ethnicity_pc_hes_hybrid),
    NA,
    "nonmissing"
  ))

check_missingness_predictors(input_data = post_1997_cohort,
                            missing_col = "any_missing",
                            miss_pred_cols = miss_predictors2)

# derive IPW weights 
ipw_res_step = get_weights_stepwise(input_data = post_1997_cohort %>% 
                                      filter(!is.na(e2019_imd_5)),
                                    missing_col = "any_missing",
                                    miss_pred_cols = c("sq_years_dat_pre_index",
                                                       "age_at_registration",
                                                       "numeric_year_of_reg",
                                                       "gender",
                                                       "numeric_imd"))

plot_from_model_output(ipw_res_step$missingness_model)
plot_missingness_prob(ipw_res_step$ipw_results)

# IPW interaction
model_data = ipw_res_step$ipw_results
model_data$bmi_category_simple = factor(model_data$bmi_category_simple,ordered=F)

vars = c("bmi_category_simple","smoking_status","im_status","vd_def_status","head_injury_status")
overall_res_list = data.frame()
for(var in vars){
model_formula = paste0("MS_status ~ index_age + ",var," * ethnicity_pc_hes_hybrid")

# print parameters
message("variable: ",var)
  
# now do weighted regression
design = survey::svydesign(id =~1,
                           weights =~ ipw,
                           data = model_data)

# regression
formula = model_formula
ethnicity_interaction_model = survey::svyglm(formula, design = design,family = quasibinomial(link="logit"))

# compare interaction and main effects model 
aov = anova(ethnicity_interaction_model,test="Chisq")
pval = aov[[length(aov)]]$p
message("P value for interaction: ",pval)

res = summary(ethnicity_interaction_model)$coefficients %>%
  data.frame() %>%
  mutate(pval_int = pval)
overall_res_list <<- bind_rows(overall_res_list,res)
}

overall_res_list


# repeat for just main effects 

model_data = ipw_res_step$ipw_results
model_data$MS_status
vars = c("bmi_category_simple","smoking_status","im_status","vd_def_status","head_injury_status","bmi_z","bmi_category","smoking_status_detailed")
overall_res_list = data.frame()
for(var in vars){
  model_formula = paste0("MS_status ~ index_age + ",var)
  
  # print parameters
  message("variable: ",var)
  
  # now do weighted regression
  design = survey::svydesign(id =~1,
                             weights =~ ipw,
                             data = model_data)
  
  # regression
  formula = model_formula
  ethnicity_interaction_model = survey::svyglm(formula, design = design,family = quasibinomial(link="logit"))
  

  res = summary(ethnicity_interaction_model)$coefficients %>%
    data.frame()
  overall_res_list <<- bind_rows(overall_res_list,res)
}

overall_res_list = overall_res_list %>%
  filter(!grepl("Intercept",rownames(overall_res_list)))
overall_res_list = overall_res_list %>%
  filter(!grepl("index_age",rownames(overall_res_list)))
overall_res_list = overall_res_list %>%
  mutate(pval = simplify_pval(`Pr...t..`))


overall_res_list$or = exp(overall_res_list$Estimate)
overall_res_list$lower_ci = exp(overall_res_list$Estimate - 1.96* overall_res_list$Std..Error)
overall_res_list$upper_ci = exp(overall_res_list$Estimate + 1.96* overall_res_list$Std..Error)
overall_res_list$var = rownames(overall_res_list)
overall_res_list$var = c(
  "Overweight/obese",
  "Ever-smoker",
  "IM",
  "Vit D def",
  "Head injury",
  "BMI Z-score",
  "Underweight",
  "Overweight",
  "Obese",
  "Morbidly obese",
  "Current smoker",
  "Ex-smoker"
)

overall_res_list$var = factor(overall_res_list$var,
                              levels = c(
                                "BMI Z-score",
                                "Overweight/obese",
                                "Underweight",
                                "Overweight",
                                "Obese",
                                "Morbidly obese",
                                "Ever-smoker",
                                "Current smoker",
                                "Ex-smoker",
                                "IM",
                                "Vit D def",
                                "Head injury"),
                              ordered=T)

ggplot(overall_res_list,
       aes(or,var))+
  geom_point()+
  theme_minimal()+
  geom_errorbarh(mapping=aes(xmin = lower_ci,xmax = upper_ci,y=var),height=0.1)+
  scale_x_continuous(limits=c(0.5,5))+
  geom_vline(xintercept = 1,alpha=0.5)+
  labs(x="Odds Ratio for MS (IPW-adjusted)",y="Variable")

