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

patient = cohorts$primary
post_1997_cohort = cohorts$post_97
post_2006_cohort = cohorts$post_06
hes_ms_cohort = cohorts$hes_ms

#################################
# case-control study
#################################

# counts 
overall_res = data.frame()
get_counts = function(x){
  message(x)
  rf = levels(factor(patient[[x]]))[2]
  non_rf = levels(factor(patient[[x]]))[1]
  
  res = patient %>%
    filter(!is.na(.data[[x]])) %>%
    dplyr::count(MS_status,.data[[x]]) %>% 
    pivot_wider(id_cols = c(x),names_from = c(MS_status),values_from=n) %>%
    replace_na(list(Case = 0,Control = 0)) %>%
    mutate(odds = Case/Control) %>%
    mutate(Exposure = x) %>%
    pivot_wider(id_cols = Exposure,names_from=x,values_from = c(Control,Case,odds)) %>%
    mutate(OR = .data[[paste0("odds_",rf)]] / .data[[paste0("odds_",non_rf)]]  )
  
  colnames(res) = c("Exposure","Control_unexposed","Control_exposed","Case_unexposed","Case_exposed","Odds_unexposed","Odds_exposed","OR")
  res$exposed_level = rf
  
  miss = patient %>%
    group_by(MS_status) %>%
    dplyr::count("var" = is.na(.data[[x]])) %>%
    mutate(percent = 100 * n / sum(n)) %>%
    filter(var == F) %>%
    mutate(nonmiss = paste0(n," (",round(percent,2),"%)")) 
  
  res = res %>%
    mutate(summ = paste0(
      "Control/unexposed=",Control_unexposed,
      "; Control/exposed=",Control_exposed,
      "; Case/unexposed=",Case_unexposed,
      "; Case/exposed=",Case_exposed,
      "; Controls non-missing data=",miss[miss$MS_status=="Control",'nonmiss'],
      "; Cases non-missing data=",miss[miss$MS_status=="Case",'nonmiss']))
      
    
  overall_res <<- bind_rows(overall_res,res)
}


vars = c("bmi_category_simple",
         "smoking_status",
         "im_status",
         "vd_def_status",
         "alcohol_status",
         "head_injury_status")
sapply(vars,get_counts)
write_csv(overall_res,"./tables/raw_exposure_counts_all.csv")


# first test predictors 
check_regression_predictors = function(predictor){
dat = patient %>% 
  mutate(rounded_predictor = round(.data[[predictor]],0)) %>%
  group_by(rounded_predictor) %>% 
  count(MS_status) %>%
  pivot_wider(id_cols = rounded_predictor,values_from = n,names_from = MS_status) %>%
  mutate(odds = Case/Control) %>%
  mutate(logodds = log(odds))
p=ggplot(dat,
       aes(rounded_predictor,logodds))+
  geom_point()+
  theme_minimal()+
  labs(y="Log(Odds) of MS",x = predictor)+
  geom_smooth(method="lm",col="red",alpha=0.3,se=F)
return(p)
}

# make basic univariable models 
make_univariable_case_control_model = function(
  var_to_test = "age_at_registration"){
  
  # print parameters
  message("variable: ",var_to_test)
  
  # subset to non-missing obs 
  model_data = patient %>% filter(!is.na(.data[[var_to_test]]))
  
  # if var is ordered, unorder is 
  model_data[[var_to_test]] = if(is.factor(model_data[[var_to_test]])){
    factor(model_data[[var_to_test]],ordered=F)
  } else{
    model_data[[var_to_test]]
  }
  
  message("make model")
  
  uni_model  = glm(
    data = model_data,
    MS_status ~ model_data[[var_to_test]],
    family=binomial(link="logit")) %>% 
    extract_coefs_from_model_obj() %>%
    mutate(model = "univariable") %>% 
    filter(var != "(Intercept)") 
  uni_model$var = str_remove(uni_model$var,pattern="model_data\\[\\[var_to_test\\]\\]")
  uni_model$variable = var_to_test
  return(uni_model)
}

covars_to_test = list("age_at_registration","yob","index_age","numeric_imd","numeric_year_of_reg",
                      "gender",
                      "e2011_urban_rural",
                      "e2019_imd_5",
                      "ethnicity_pc_hes_hybrid")

univariable_model_results = lapply(covars_to_test,
       make_univariable_case_control_model)

univariable_model_results = do.call("bind_rows",univariable_model_results)
univariable_model_results$simple_p = simplify_pval(univariable_model_results$pval)
p=ggplot(univariable_model_results,
       aes(or,paste0(variable,var),label=simple_p))+
  geom_point()+
  geom_errorbarh(mapping = aes(xmin = lower_ci,xmax = upper_ci,y=paste0(variable,var)),height=0.2)+
  theme_minimal()+
  scale_x_log10(limits=c(0.3,3))+
  geom_text(mapping = aes(x=2,y =paste0(variable,var) ))+
  annotate("text",y=14.4,x=2,label="P value")+
  geom_vline(xintercept=1,alpha=0.3)+
  labs(y="Variable",x="Odds ratio for MS")+ggtitle("A")

p1=check_regression_predictors("age_at_registration")+ggtitle("B")
p2=check_regression_predictors("yob")+ggtitle("C")
p3=check_regression_predictors("index_age")+ggtitle("D")
p4=check_regression_predictors("numeric_imd")+ggtitle("E")
p5=check_regression_predictors("numeric_year_of_reg")+ggtitle("F")

plot_fx(grid.arrange(p,p1,p2,p3,p4,p5,nrow=2),
        "./figs/covariate_selection.png",
        width=14,height=10)


# create model function 
make_case_control_model = function(
  data_for_model = patient,
  var_to_test = "bmi_category"){
  
  # print parameters
  message("variable: ",var_to_test)
  
  # subset to non-missing obs 
  model_data = data_for_model %>% filter(!is.na(.data[[var_to_test]]))
  
  # if var is ordered, unorder is 
  model_data[[var_to_test]] = if(is.factor(model_data[[var_to_test]])){
    factor(model_data[[var_to_test]],ordered=F)
  } else{
    model_data[[var_to_test]]
  }
  
  uni_model  = glm(
    data = model_data,
    MS_status ~ model_data[[var_to_test]],
    family=binomial(link="logit")) %>% 
    extract_coefs_from_model_obj() %>%
    mutate(model = "univariable")
  
  agesex_model  = glm(
    data = model_data,
    MS_status ~ index_age + gender + model_data[[var_to_test]],
    family=binomial(link="logit")) %>% 
    extract_coefs_from_model_obj() %>%
    mutate(model = "age + sex")
  
  agesex_ageatreg_model  = glm(
    data = model_data,
    MS_status ~ age_at_registration + index_age + gender + model_data[[var_to_test]],
    family=binomial(link="logit")) %>% 
    extract_coefs_from_model_obj() %>%
    mutate(model = "age + sex + age_at_reg")
  
  agesex_ethnic_model  = glm(
    data = model_data,
    MS_status ~ index_age + gender +
      ethnicity_pc_hes_hybrid + model_data[[var_to_test]],
    family=binomial(link="logit")) %>% 
    extract_coefs_from_model_obj() %>%
    mutate(model = "age + sex + ethnicity")
  
  agesex_imd_model  = glm(
    data = model_data,
    MS_status ~ index_age + gender + ethnicity_pc_hes_hybrid +
      numeric_imd + model_data[[var_to_test]],
    family=binomial(link="logit")) %>% 
    extract_coefs_from_model_obj()%>%
    mutate(model = "age + sex + IMD")
  
  covar_names = c("(Intercept)",
                  "index_age",
                  "age_at_registration",
                  "genderMale",
                  "ethnicity_pc_hes_hybridBlack",
                  "ethnicity_pc_hes_hybridMixed/Other",
                  "ethnicity_pc_hes_hybridAsian",
                  "numeric_imd") 
  
  overall_model_df = bind_rows(
    uni_model,
    agesex_model,
    agesex_ageatreg_model,
    agesex_ethnic_model,
    agesex_imd_model) %>%
    filter(!var %in% covar_names) 
  overall_model_df$var = str_remove(overall_model_df$var,pattern="model_data\\[\\[var_to_test\\]\\]")
  overall_model_df$variable = var_to_test
  return(overall_model_df)
}

# run for ethnicity
make_case_control_model(var_to_test = "ethnicity_pc_hes_hybrid")

# run for different variables 
bmi_cat = make_case_control_model(var_to_test = "bmi_category_simple")
smoking = make_case_control_model(var_to_test = "smoking_status")
im = make_case_control_model(var_to_test = "im_status")
vd_def = make_case_control_model(var_to_test = "vd_def_status")
alcohol = make_case_control_model(var_to_test = "alcohol_status")
hi = make_case_control_model(var_to_test = "head_injury_status")

table(hes_ms_cohort$MS_status)
bmi_cont = make_case_control_model(var_to_test = "bmi_z")
bmi_cat_who = make_case_control_model(var_to_test = "bmi_category")

# combine results and plot
model_overall_df = bind_rows(bmi_cat,smoking,im,vd_def,alcohol,hi)
model_overall_df$y = paste0(model_overall_df$variable,":",model_overall_df$var)

model_overall_df = model_overall_df %>%
  mutate(bonf_sig = ifelse(pval < 0.05/6,
                           "*",
                           "NS"))

model_overall_df$y = factor(model_overall_df$y)
levels(model_overall_df$y) = c("Non-drinker",
                               "Overweight/obese",
                               "Prior head injury",
                               "Prior IM",
                               "Ever smoker",
                               "Vitamin D deficiency")
p=ggplot(model_overall_df,
         aes(or,y,col=model,label=bonf_sig))+
  geom_point(position=ggstance::position_dodgev(height=0.5))+
  geom_errorbarh(mapping=aes(xmin = lower_ci,xmax=upper_ci,y=y),height=0.1,position=ggstance::position_dodgev(height=0.5))+
  theme_minimal()+
  geom_vline(xintercept=1,alpha=0.2)+
  labs(x="Odds Ratio (for MS)",y="Exposure")+
  scale_color_brewer(palette="Paired")+
  scale_x_log10(limits=c(0.5,5))

png("./figs/case_control_plots.png",res=300,units="in",width=6,height=6)
p
dev.off()

# repeat with just main analysis (age + sex)
p=ggplot(model_overall_df %>% filter(model=="age + sex"),
         aes(or,y,label=bonf_sig))+
  geom_point()+
  geom_errorbarh(mapping=aes(xmin = lower_ci,xmax=upper_ci,y=y),height=0.1)+
  theme_minimal()+
  geom_vline(xintercept=1,alpha=0.2)+
  labs(x="Odds Ratio (for MS)",y="Exposure")+
  scale_color_brewer(palette="Paired")+
  geom_text(mapping = aes(upper_ci,y,hjust=-1))+
  scale_x_log10(limits=c(0.5,5),breaks = c(0.5,1,2,3,4,5))+
  theme(legend.position="none")

png("./figs/case_control_plots_age_sex.png",res=300,units="in",width=4,height=4)
p
dev.off()

# save results to a table
tbl = model_overall_df %>% 
  dplyr::select(variable, var,eff_n,or,lower_ci,upper_ci,model, pval, bonf_sig)%>% 
  mutate(
    or_ci = paste0(
      round(or,2)," (",round(lower_ci,2)," - ",round(upper_ci,2),")"
    )
  ) %>% 
  dplyr::select(-or,-lower_ci,-upper_ci) %>%
  dplyr::select(1,2,3,7,4,5,6)

tbl = tbl %>% pivot_wider(id_cols = c(variable,var), names_from = model, values_from = c(eff_n,or_ci,pval, bonf_sig))
tbl = tbl %>% dplyr::select(1,2,3,8,13,4,9,14,5,10,15,6,11,16,7,12,17,8,13,18,9,14,19)

# recode P vals 
tbl = tbl %>% mutate_at(vars(contains("pval")),simplify_pval)
write_csv(tbl,"./tables/case_control_results.csv")


# do sensitivity analyses

sensitivity_plots = list()
do_sensitivity_analysis = function(dataset_for_cc_analysis_name = "post_2006_cohort"){
  
  # evaluate name of arg to get df 
  dataset_for_cc_analysis = get(dataset_for_cc_analysis_name)
  
  # recode MS status with control as ref 
  dataset_for_cc_analysis$MS_status = relevel(factor(dataset_for_cc_analysis$MS_status),ref="Control")
  
  # run for different variables 
  bmi_cat = make_case_control_model(var_to_test = "bmi_category_simple", data_for_model = dataset_for_cc_analysis)
  bmi_cont = make_case_control_model(var_to_test = "bmi_z", data_for_model = dataset_for_cc_analysis)
  smoking = make_case_control_model(var_to_test = "smoking_status", data_for_model = dataset_for_cc_analysis)
  im = make_case_control_model(var_to_test = "im_status", data_for_model = dataset_for_cc_analysis)
  vd_def = make_case_control_model(var_to_test = "vd_def_status", data_for_model = dataset_for_cc_analysis)
  alcohol = make_case_control_model(var_to_test = "alcohol_status", data_for_model = dataset_for_cc_analysis)
  hi = make_case_control_model(var_to_test = "head_injury_status", data_for_model = dataset_for_cc_analysis)
  
  # combine results and plot
  model_overall_df = bind_rows(bmi_cat,smoking,im,vd_def,alcohol,hi)
  model_overall_df$y = paste0(model_overall_df$variable,":",model_overall_df$var)
  
  model_overall_df = model_overall_df %>%
    mutate(bonf_sig = ifelse(pval < 0.05/6,
                             "*",
                             "NS"))
  
  model_overall_df$y = factor(model_overall_df$y)
  levels(model_overall_df$y) = c("Non-drinker",
                                 "Overweight/obese",
                                 "Prior head injury",
                                 "Prior IM",
                                 "Ever smoker",
                                 "Vitamin D deficiency")
  p=ggplot(model_overall_df,
           aes(or,y,col=model,label=bonf_sig))+
    geom_point(position=ggstance::position_dodgev(height=0.5))+
    geom_errorbarh(mapping=aes(xmin = lower_ci,xmax=upper_ci,y=y),height=0.1,position=ggstance::position_dodgev(height=0.5))+
    theme_minimal()+
    geom_vline(xintercept=1,alpha=0.2)+
    labs(x="Odds Ratio (for MS)",y="Exposure")+
    scale_color_brewer(palette="Paired")+
    scale_x_log10(limits=c(0.5,10))+
    ggtitle(dataset_for_cc_analysis_name)
  
  png(paste0("./figs/case_control_plots_sensitivity_analysis",dataset_for_cc_analysis_name,".png"),res=300,units="in",width=6,height=4)
  print(p)
  dev.off()
  
  sensitivity_plots[[length(sensitivity_plots)+1]] <<- p
  # save results to a table
  tbl = model_overall_df %>% 
    dplyr::select(variable, var,eff_n,or,lower_ci,upper_ci,model, pval, bonf_sig)%>% 
    mutate(
      or_ci = paste0(
        round(or,2)," (",round(lower_ci,2)," - ",round(upper_ci,2),")"
      )
    ) %>% 
    dplyr::select(-or,-lower_ci,-upper_ci) %>%
    dplyr::select(1,2,3,7,4,5,6)
  
  tbl = tbl %>% pivot_wider(id_cols = c(variable,var), names_from = model, values_from = c(eff_n,or_ci,pval, bonf_sig))
  tbl = tbl %>% dplyr::select(1,2,3,8,13,4,9,14,5,10,15,6,11,16,7,12,17,8,13,18,9,14,19)
  
  # recode P vals 
  tbl = tbl %>% mutate_at(vars(contains("pval")),simplify_pval)
  
  
  write_csv(tbl,paste0("./tables/case_control_plots_sensitivity_analysis",dataset_for_cc_analysis_name,".csv"))
}
do_sensitivity_analysis("post_1997_cohort")
do_sensitivity_analysis("post_2006_cohort")
do_sensitivity_analysis("hes_ms_cohort")

plot_fx(do.call("grid.arrange",sensitivity_plots),
        "./figs/sensitivity_cc_plots.png",
        width=12,height=6)

# combined model 
table(patient$bmi_category_simple)
multivariable_model  = glm(
  data = patient,
  MS_status ~ index_age + gender + 
    bmi_category_simple +
    smoking_status +
    im_status +
    vd_def_status +
    alcohol_status +
    head_injury_status,
  family=binomial(link="logit")) %>% 
  extract_coefs_from_model_obj()%>%
  mutate(model = "all")

multivariable_model = multivariable_model %>% 
  dplyr::select(var,or,lower_ci,upper_ci,pval) %>%
  mutate(
    or_ci = paste0(
      round(or,2)," (",round(lower_ci,2)," - ",round(upper_ci,2),")"
    )
  ) %>% 
  dplyr::select(-or,-lower_ci,-upper_ci) %>%
  dplyr::select(var,or_ci,pval)
write_csv(multivariable_model,"multi_model.csv")

################################################
# look in more detail at categorical vars 
################################################

# bmi categories
patient$bmi_category = factor(patient$bmi_category,ordered = F)
patient$bmi_category = relevel(factor(patient$bmi_category),ref="Healthy")
plot_dat = make_case_control_model(var_to_test = "bmi_category")
plot_dat$var = factor(plot_dat$var,
                      levels = c("Underweight","Healthy","Overweight","Obese","Morbidly obese"),
                      ordered = T)
table(patient$bmi_category,patient$MS_status)

p1=make_categorical_cc_plot(plot_dat = plot_dat)+scale_y_continuous(limits=c(0.1,7.5))+ggtitle("A")


# bmi z score
make_case_control_model(var_to_test = "bmi_z")

# smoking 
patient$smoking_status_detailed = relevel(factor(patient$smoking_status_detailed),ref="never smoker")

plot_dat = make_case_control_model(var_to_test = "smoking_status_detailed")
p2=make_categorical_cc_plot(plot_dat = plot_dat)+
  scale_y_continuous(limits=c(0.5,2.5))+ggtitle("C")


# smoking detailed categories
patient$smoking_status_extra_detailed = relevel(factor(patient$smoking_status_extra_detailed),ref="never smoker")
table(patient$MS_status,patient$smoking_status_extra_detailed)

table(patient$MS_status,patient$smoking_status_detailed)
plot_dat = make_case_control_model(var_to_test = "smoking_status_extra_detailed")

plot_dat$var = factor(plot_dat$var,
                      levels = c(
                        "current_very_light",
                        "current_light",
                        "current_moderate",
                        "current_heavy",
                        "ex_very_light",
                        "ex_light",
                        "ex_moderate",
                        "ex_heavy"
                      ),
                      ordered = T)
p3=make_categorical_cc_plot(plot_dat = plot_dat)+theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))+labs(x="Smoking status")+
  scale_y_continuous(limits=c(0.5,2.5))+ggtitle("D")

table(patient$smoking_status_detailed)



plot_fx(grid.arrange(p1,p2,p3),
        filename = "./figs/categorical_cc_plots.png",
        width=8,height=8)



###################################
# ethnic variation in risk factors
###################################

# first compare demographics between ethnicities
make_demographics(
  input_data = patient,
  response_var = "ethnicity_pc_hes_hybrid",
  cat_vars = c("gender","dead","e2019_imd_5","e2011_urban_rural","in_hes"),
  cont_vars = c("yob","numeric_year_of_reg","age_at_registration","fu_time","index_age","years_data_before_index","years_data_after_index"),
  outfile = "ethnicities"
)

# first compare demographics between ethnicities
make_demographics(
  input_data = patient %>% filter(ethnicity_pc_hes_hybrid=="Asian"),
  response_var = "MS_status",
  cat_vars = c("gender","dead","e2019_imd_5","e2011_urban_rural","in_hes"),
  cont_vars = c("yob","numeric_year_of_reg","age_at_registration","fu_time","index_age","years_data_before_index","years_data_after_index"),
  outfile = "asian_demographics"
)
make_demographics(
  input_data = patient %>% filter(ethnicity_pc_hes_hybrid=="Black"),
  response_var = "MS_status",
  cat_vars = c("gender","dead","e2019_imd_5","e2011_urban_rural","in_hes"),
  cont_vars = c("yob","numeric_year_of_reg","age_at_registration","fu_time","index_age","years_data_before_index","years_data_after_index"),
  outfile = "black_demographics"
)

make_demographics(
  input_data = patient %>% filter(ethnicity_pc_hes_hybrid=="White"),
  response_var = "MS_status",
  cat_vars = c("gender","dead","e2019_imd_5","e2011_urban_rural","in_hes"),
  cont_vars = c("yob","numeric_year_of_reg","age_at_registration","fu_time","index_age","years_data_before_index","years_data_after_index"),
  outfile = "white_demographics"
)

just_ms_ethnicities_dat = patient %>% filter(MS_status=="Case") 
make_demographics(
  input_data = just_ms_ethnicities_dat,
  response_var = "ethnicity_pc_hes_hybrid",
  cat_vars = c("gender","dead","e2019_imd_5","e2011_urban_rural","in_hes"),
  cont_vars = c("yob","numeric_year_of_reg","age_at_registration","fu_time","index_age","years_data_before_index","years_data_after_index"),
  outfile = "just_ms_by_ethnicity"
)

just_control_ethnicities_dat = patient %>% filter(MS_status=="Control") 
make_demographics(
  input_data = just_control_ethnicities_dat,
  response_var = "ethnicity_pc_hes_hybrid",
  cat_vars = c("gender","dead","e2019_imd_5","e2011_urban_rural","in_hes"),
  cont_vars = c("yob","numeric_year_of_reg","age_at_registration","fu_time","index_age","years_data_before_index","years_data_after_index"),
  outfile = "just_controls_by_ethnicity"
)

# look in detail at year of reg
yor = patient %>% group_by(ethnicity_pc_hes_hybrid) %>% summarise(med_yor = median(numeric_year_of_reg))
p=ggplot(patient %>% filter(!is.na(ethnicity_pc_hes_hybrid)),aes(numeric_year_of_reg,fill=ethnicity_pc_hes_hybrid))+
  geom_density(alpha=0.3)+
  theme_minimal()+
  labs(x="Year of registration",fill="Ethnicity")
p2=ggplot(patient %>% filter(!is.na(ethnicity_pc_hes_hybrid)),aes(numeric_year_of_reg,fill=ethnicity_pc_hes_hybrid))+
  geom_density(alpha=0.3)+
  theme_minimal()+
  labs(x="Year of registration",fill="Ethnicity")+
  facet_wrap(~ethnicity_pc_hes_hybrid,ncol=1)+ggtitle("A")
p3=ggplot(patient %>% filter(!is.na(ethnicity_pc_hes_hybrid)),
          aes(ethnicity_pc_hes_hybrid,numeric_year_of_reg,fill=MS_status))+
  geom_boxplot(alpha=0.8)+
  theme_minimal()+
  labs(x="Year of registration",fill="Ethnicity")+
  coord_flip()+
  scale_fill_brewer(palette="Set2")+
  labs(y="Year of registration",x="Ethnicity",fill="MS status")+ggtitle("B")

patient %>% group_by(ethnicity_pc_hes_hybrid,MS_status) %>% summarise(mean(numeric_year_of_reg))
plot_fx(grid.arrange(p2,p3),"./figs/year_of_reg_ethnicity.png",height=8)


patient %>% filter(is.na(ethnicity_pc_hes_hybrid)) %>%
  filter(years_data_after_index<0) %>% dplyr::select(1,3,4,regstartdate,regenddate,MS_status,index_code_date,ms_code_date) %>%
  filter(MS_status=="Case")


####################################
# Run models 
####################################


overall_res = data.frame()
get_counts = function(x){
  message(x)
  rf = levels(factor(patient[[x]]))[2]
  non_rf = levels(factor(patient[[x]]))[1]
  
  res = patient %>%
    filter(!is.na(.data[[x]])) %>%
    filter(!is.na(ethnicity_pc_hes_hybrid)) %>%
    dplyr::count(MS_status,ethnicity_pc_hes_hybrid,.data[[x]]) %>% 
    pivot_wider(id_cols = c(x,ethnicity_pc_hes_hybrid),names_from = c(MS_status),values_from=n) %>%
    replace_na(list(Case = 0,Control = 0)) %>%
    mutate(odds = Case/Control) %>%
    pivot_wider(id_cols = ethnicity_pc_hes_hybrid,names_from=x,values_from = c(Control,Case,odds)) %>%
    mutate(OR = .data[[paste0("odds_",rf)]] / .data[[paste0("odds_",non_rf)]]  )
  
  colnames(res) = c("Ethnicity","Control_unexposed","Control_exposed","Case_unexposed","Case_exposed","Odds_unexposed","Odds_exposed","OR")
  res$Exposure = x
  res$exposed_level = rf
  overall_res <<- bind_rows(overall_res,res)
}


vars = c("bmi_category_simple",
         "smoking_status",
         "im_status",
         "vd_def_status",
         "head_injury_status")
sapply(vars,get_counts)
write_csv(overall_res,"./tables/raw_exposure_counts.csv")


####################################################
# Analysis 1 - white vs other 
####################################################

# define new binary ethnicity variable 
patient = patient %>%
  mutate(ethnicity_binary = ifelse(
    ethnicity_pc_hes_hybrid=="White",
    "White",
    "Diverse"
  ))

hes_ms_cohort = hes_ms_cohort %>%
  mutate(ethnicity_binary = ifelse(
    ethnicity_pc_hes_hybrid=="White",
    "White",
    "Diverse"
  ))

# initialise master results df
master_results_df = data.frame()


# counts 

# repeat with binary ethnicity
overall_res = data.frame()
get_counts = function(x){
  message(x)
  rf = levels(factor(patient[[x]]))[2]
  non_rf = levels(factor(patient[[x]]))[1]
  
  res = patient %>%
    filter(!is.na(.data[[x]])) %>%
    filter(!is.na(ethnicity_binary)) %>%
    dplyr::count(MS_status,ethnicity_binary,.data[[x]]) %>% 
    pivot_wider(id_cols = c(x,ethnicity_binary),names_from = c(MS_status),values_from=n) %>%
    replace_na(list(Case = 0,Control = 0)) %>%
    mutate(odds = Case/Control) %>%
    pivot_wider(id_cols = ethnicity_binary,names_from=x,values_from = c(Control,Case,odds)) %>%
    mutate(OR = .data[[paste0("odds_",rf)]] / .data[[paste0("odds_",non_rf)]]  )
  
  colnames(res) = c("Ethnicity","Control_unexposed","Control_exposed","Case_unexposed","Case_exposed","Odds_unexposed","Odds_exposed","OR")
  res$Exposure = x
  res$exposed_level = rf
  overall_res <<- bind_rows(overall_res,res)
}


vars = c("bmi_category_simple",
         "smoking_status",
         "im_status",
         "vd_def_status",
         "head_injury_status")
sapply(vars,get_counts)
write_csv(overall_res,"./tables/raw_exposure_counts_binary_ethnicity.csv")

# create model function 
make_case_control_model_by_ethnicity = function(
  data_for_model = patient,
  var_to_test,
  covars){
  
  # print parameters
  message("variable: ",var_to_test)
  
  # subset to non-missing obs & non-missing ethnicity
  model_data = data_for_model %>% 
    filter(!is.na(.data[[var_to_test]])) %>% 
    filter(!is.na(ethnicity_binary))
  message("subsetted")
  
  # if var is ordered, unorder is 
  model_data[[var_to_test]] = if(is.factor(model_data[[var_to_test]])){
    factor(model_data[[var_to_test]],ordered=F)
  } else{
    model_data[[var_to_test]]
  }
  
  # make model
  message("making model")
  
  formula = paste0("MS_status ~ ",paste0(covars,collapse="+"),"+ ethnicity_binary * ",var_to_test)
  message(formula)
  ethnicity_interaction_model  = glm(
    data = model_data,
    eval(formula),
    family=binomial(link="logit")) 
  message("made model")
  
  # get n
  eff_n = ethnicity_interaction_model$fitted.values %>% length
  summary(ethnicity_interaction_model)
  # compare interaction and main effects model 
  aov = anova(ethnicity_interaction_model,test="Chisq")
  pval = aov[nrow(aov),"Pr(>Chi)"]
  message("P value for interaction: ",pval)
  
  
  # get summaries 
  ethnicity_interaction_model_summary = ethnicity_interaction_model %>% 
    extract_coefs_from_model_obj() %>%
    mutate(model = "ethnicity_interaction")
  
  # stratified models
  formula2 = paste0("MS_status ~ ",paste0(covars,collapse="+"),"+ ",var_to_test)
  
  white_data = model_data %>% filter(ethnicity_binary=="White")
  white_model  = glm(
    data = white_data,
    eval(formula2),
    family=binomial(link="logit")) %>% 
    extract_coefs_from_model_obj() %>%
    mutate(model = "white")
  
  non_white_data = model_data %>% filter(ethnicity_binary!="White")
  non_white_model  = glm(
    data = non_white_data %>% filter(ethnicity_binary!="White"),
    eval(formula2),
    family=binomial(link="logit")) %>% 
    extract_coefs_from_model_obj() %>%
    mutate(model = "diverse")
  

  overall_model_df = bind_rows(
    ethnicity_interaction_model_summary,
    white_model,
    non_white_model) 
  overall_model_df$var = str_remove(overall_model_df$var,pattern="model_data\\[\\[var_to_test\\]\\]")
  overall_model_df$variable = var_to_test
  overall_model_df$pval_interaction = pval
  overall_model_df$eff_n = eff_n
  return(overall_model_df)
}

make_case_control_model_by_ethnicity(var_to_test = "bmi_z",
                                     covars = c("index_age","gender"))
make_case_control_model_by_ethnicity(var_to_test = "bmi_category",
                                     covars = c("index_age","gender"))

run_interaction_analysis = function(dat = patient, 
                                    plot_title = "whole_cohort", 
                                    covars_to_test = c("index_age","gender")){
  # run for different variables 
  bmi_cont = make_case_control_model_by_ethnicity(var_to_test = "bmi_z",data_for_model = dat, covars = covars_to_test)
  bmi_cat = make_case_control_model_by_ethnicity(var_to_test = "bmi_category",data_for_model = dat,covars = covars_to_test)
  bmi_cat_simple = make_case_control_model_by_ethnicity(var_to_test = "bmi_category_simple",data_for_model = dat,covars = covars_to_test)
  smoking = make_case_control_model_by_ethnicity(var_to_test = "smoking_status",data_for_model = dat,covars = covars_to_test)
  smoking_detailed = make_case_control_model_by_ethnicity(var_to_test = "smoking_status_detailed",data_for_model = dat,covars = covars_to_test)
  im = make_case_control_model_by_ethnicity(var_to_test = "im_status",data_for_model = dat,covars = covars_to_test)
  vd_def = make_case_control_model_by_ethnicity(var_to_test = "vd_def_status",data_for_model = dat,covars = covars_to_test)
  hi = make_case_control_model_by_ethnicity(var_to_test = "head_injury_status",data_for_model = dat,covars = covars_to_test)
  
  # combine results and plot
  model_overall_df = bind_rows(bmi_cont,bmi_cat,bmi_cat_simple,smoking,smoking_detailed,im,vd_def,hi)
  
  # look at interaction models 
  interaction_model_df = model_overall_df %>% 
    filter(grepl("interaction",model)) %>%
    filter(grepl(":",var)) %>%
    mutate(p_bonf_sig = ifelse(pval_interaction<0.05/5,"sig","NS"))
  
  # make simple forest comparing ethnicities
  model_overall_df = model_overall_df %>% 
    filter(model %in% c("white","diverse")) %>%
    filter(!var %in% c("(Intercept)","age_at_registration","index_age","genderMale","numeric_imd","fu_time","e2011_urban_ruralUrban")) %>%
    mutate(bonf_sig = ifelse(pval < 0.05/5,
                             "*",
                             "NS")) %>%
    mutate(y = paste0(variable,":",var)) %>%
    mutate(y = factor(y))
  levels(model_overall_df$y) = c("Overweight/obese",
                                 "BMI>40",
                                 "BMI 30-40",
                                 "BMI 25-30",
                                 "BMI < 18.5",
                                 "1-SD increase in BMI",
                                 "Prior head injury",
                                 "Prior IM",
                                 "Current smoker",
                                 "Ex-smoker",
                                 "Ever smoker",
                                 "Vitamin D deficiency")
  model_overall_df$p_int_bonf = ifelse(model_overall_df$pval_interaction<0.05/5,"*","NS")
  p=ggplot(model_overall_df,
           aes(or,y,fill=model,label=p_int_bonf))+
    geom_errorbarh(mapping=aes(xmin = lower_ci,xmax=upper_ci,y=y),height=0.01,position=ggstance::position_dodgev(height=0.4),alpha=0.7)+
    geom_point(alpha=0.7,shape=21,color="black",position=ggstance::position_dodgev(height=0.4))+
    theme_minimal()+
    geom_vline(xintercept=1,alpha=0.2)+
    labs(x="Odds Ratio (for MS)",y="Exposure",fill="Ethnicity")+
    scale_color_brewer(palette="Set2")+
    scale_x_log10(limits=c(0.01,100),breaks = c(0.1,1,10))+
    geom_text(mapping = aes(x=90,y=y,label=p_int_bonf),show.legend = F)
  
  png(paste0("./figs/exposure_het_ethnicity",plot_title,".png"),res=300,units="in",width=4,height=3)
  print(p)
  dev.off()
  
  # simple plot 
  p2=ggplot(model_overall_df %>% filter(y %in%  c("Overweight/obese",
                                                  "Prior head injury",
                                                  "Prior IM",
                                                  "Ever smoker",
                                                  "Vitamin D deficiency")),
            aes(or,y,fill=model))+
    geom_errorbarh(mapping=aes(xmin = lower_ci,xmax=upper_ci,y=y),height=0.01,position=ggstance::position_dodgev(height=0.4),alpha=0.7)+
    geom_point(alpha=0.7,shape=21,color="black",position=ggstance::position_dodgev(height=0.4))+
    theme_minimal()+
    geom_vline(xintercept=1,alpha=0.2)+
    labs(x="Odds Ratio (for MS)",y="Exposure",fill="Ethnicity")+
    scale_color_brewer(palette="Set2")+
    scale_x_log10(limits=c(0.05,50))
  png(paste0("./figs/exposure_simple_plot_het_ethnicity",plot_title,".png"),res=300,units="in",width=4,height=3)
  print(p2)
  dev.off()
  
  master_results_df <<- bind_rows(master_results_df, model_overall_df %>% mutate("cohort" = plot_title))
  strat_models = model_overall_df %>% mutate(estimate = 
                                               paste0(round(or,2),
                                                      " (",round(lower_ci,2),
                                                      " - ",
                                                      round(upper_ci,2),
                                                      ")"
                                               )) %>%
    dplyr::select(y,eff_n,estimate,pval,pval_interaction,model) %>%
    pivot_wider(id_cols = c("y","eff_n"),names_from = model,values_from = c("estimate","pval","pval_interaction"))
  strat_models = strat_models %>% mutate_at(vars(contains("pval")),simplify_pval)
  
  
  write_csv(strat_models,paste0("./tables/stratified_exposure_het_ethnicity",plot_title,".csv"))
  write_csv(interaction_model_df,paste0("./tables/interactions_exposure_het_ethnicity",plot_title,".csv"))
}



run_interaction_analysis(dat=patient,plot_title="whole_cohort")
run_interaction_analysis(dat=patient,plot_title="whole_cohort_imd",
                         covars_to_test = c("index_age","gender","numeric_imd"))
run_interaction_analysis(dat=hes_ms_cohort,plot_title="hes_cohort",
                         covars_to_test = c("index_age","gender","numeric_imd"))


###################################################


# create model function 
make_case_control_model_by_ethnicity = function(
  data_for_model = patient,
  var_to_test,
  covars){
  
  # print parameters
  message("variable: ",var_to_test)
  
  # subset to non-missing obs & non-missing ethnicity
  model_data = data_for_model %>% 
    filter(!is.na(.data[[var_to_test]])) %>% 
    filter(!is.na(ethnicity_pc_hes_hybrid))
  message("subsetted")
  
  # if var is ordered, unorder is 
  model_data[[var_to_test]] = if(is.factor(model_data[[var_to_test]])){
    factor(model_data[[var_to_test]],ordered=F)
  } else{
    model_data[[var_to_test]]
  }
  
  # make model
  message("making model")
  
  formula = paste0("MS_status ~ ",paste0(covars,collapse="+"),"+ ethnicity_pc_hes_hybrid * ",var_to_test)
  message(formula)
  ethnicity_interaction_model  = glm(
    data = model_data,
    eval(formula),
    family=binomial(link="logit")) 
  message("made model")
  
  # get n
  eff_n = ethnicity_interaction_model$fitted.values %>% length
  summary(ethnicity_interaction_model)
  # compare interaction and main effects model 
  aov = anova(ethnicity_interaction_model,test="Chisq")
  pval = aov[nrow(aov),"Pr(>Chi)"]
  message("P value for interaction: ",pval)
  
 
  # get summaries 
  ethnicity_interaction_model_summary = ethnicity_interaction_model %>% 
    extract_coefs_from_model_obj() %>%
    mutate(model = "ethnicity_interaction")
  
  # stratified models
  formula2 = paste0("MS_status ~ ",paste0(covars,collapse="+"),"+ ",var_to_test)
  
  white_data = model_data %>% filter(ethnicity_pc_hes_hybrid=="White")
  white_model  = glm(
    data = white_data,
    eval(formula2),
    family=binomial(link="logit")) %>% 
    extract_coefs_from_model_obj() %>%
    mutate(model = "white")
  
  black_data = model_data %>% filter(ethnicity_pc_hes_hybrid=="Black")
  black_model  = glm(
    data = model_data %>% filter(ethnicity_pc_hes_hybrid=="Black"),
    eval(formula2),
    family=binomial(link="logit")) %>% 
    extract_coefs_from_model_obj() %>%
    mutate(model = "black")
  
  asian_data = model_data %>% filter(ethnicity_pc_hes_hybrid=="Asian")
  asian_model  = glm(
    data = model_data %>% filter(ethnicity_pc_hes_hybrid=="Asian"),
    eval(formula2),
    family=binomial(link="logit")) %>% 
    extract_coefs_from_model_obj() %>%
    mutate(model = "asian")
  
  mixed_data = model_data %>% filter(ethnicity_pc_hes_hybrid=="Mixed/Other")
  mixed_model  = glm(
    data = model_data %>% filter(ethnicity_pc_hes_hybrid=="Mixed/Other"),
    eval(formula2),
    family=binomial(link="logit")) %>% 
    extract_coefs_from_model_obj() %>%
    mutate(model = "mixed")
  
  overall_model_df = bind_rows(
    ethnicity_interaction_model_summary,
    white_model,
    asian_model,
    black_model,
    mixed_model) 
  overall_model_df$var = str_remove(overall_model_df$var,pattern="model_data\\[\\[var_to_test\\]\\]")
  overall_model_df$variable = var_to_test
  overall_model_df$pval_interaction = pval
  overall_model_df$eff_n = eff_n
  return(overall_model_df)
}

make_case_control_model_by_ethnicity(var_to_test = "bmi_z",
                                     covars = c("index_age","gender"))
make_case_control_model_by_ethnicity(var_to_test = "bmi_category",
                                     covars = c("index_age","gender"))

run_interaction_analysis = function(dat = patient, 
                                    plot_title = "whole_cohort", 
                                    covars_to_test = c("index_age","gender")){
  # run for different variables 
  bmi_cont = make_case_control_model_by_ethnicity(var_to_test = "bmi_z",data_for_model = dat, covars = covars_to_test)
  bmi_cat = make_case_control_model_by_ethnicity(var_to_test = "bmi_category",data_for_model = dat,covars = covars_to_test)
  bmi_cat_simple = make_case_control_model_by_ethnicity(var_to_test = "bmi_category_simple",data_for_model = dat,covars = covars_to_test)
  smoking = make_case_control_model_by_ethnicity(var_to_test = "smoking_status",data_for_model = dat,covars = covars_to_test)
  smoking_detailed = make_case_control_model_by_ethnicity(var_to_test = "smoking_status_detailed",data_for_model = dat,covars = covars_to_test)
  
  # combine results and plot
  model_overall_df = bind_rows(bmi_cont,bmi_cat,bmi_cat_simple,smoking,smoking_detailed)
  
  # look at interaction models 
  interaction_model_df = model_overall_df %>% 
    filter(grepl("interaction",model)) %>%
    filter(grepl(":",var)) %>%
    mutate(p_bonf_sig = ifelse(pval_interaction<0.05/5,"sig","NS"))
  
  # make simple forest comparing ethnicities
  model_overall_df = model_overall_df %>% 
    filter(model %in% c("white","asian","black","mixed")) %>%
    filter(!var %in% c("(Intercept)","age_at_registration","index_age","genderMale","numeric_imd","fu_time","e2011_urban_ruralUrban")) %>%
    mutate(bonf_sig = ifelse(pval < 0.05/5,
                             "*",
                             "NS")) %>%
    mutate(y = paste0(variable,":",var)) %>%
    mutate(y = factor(y))
  levels(model_overall_df$y) = c("Overweight/obese",
                                 "BMI>40",
                                 "BMI 30-40",
                                 "BMI 25-30",
                                 "BMI < 18.5",
                                 "1-SD increase in BMI",
                                 "Current smoker",
                                 "Ex-smoker",
                                 "Ever smoker")
  model_overall_df$p_int_bonf = ifelse(model_overall_df$pval_interaction<0.05/6,"*","NS")
  p=ggplot(model_overall_df,
           aes(or,y,fill=model,label=p_int_bonf))+
    geom_errorbarh(mapping=aes(xmin = lower_ci,xmax=upper_ci,y=y),height=0.01,position=ggstance::position_dodgev(height=0.4),alpha=0.7)+
    geom_point(alpha=0.7,shape=21,color="black",position=ggstance::position_dodgev(height=0.4))+
    theme_minimal()+
    geom_vline(xintercept=1,alpha=0.2)+
    labs(x="Odds Ratio (for MS)",y="Exposure",fill="Ethnicity")+
    scale_color_brewer(palette="Set2")+
    scale_x_log10(limits=c(0.01,100),breaks = c(0.1,1,10))
  png(paste0("./figs/exposure_het_ethnicity_full",plot_title,".png"),res=300,units="in",width=6,height=4)
  print(p)
  dev.off()
  
  # simple plot 
  p2=ggplot(model_overall_df %>% filter(y %in%  c("Overweight/obese",
                                                "Ever smoker"
                                               )),
           aes(or,y,fill=model,label=p_int_bonf))+
    geom_errorbarh(mapping=aes(xmin = lower_ci,xmax=upper_ci,y=y),height=0.01,position=ggstance::position_dodgev(height=0.4),alpha=0.7)+
    geom_point(alpha=0.7,shape=21,color="black",position=ggstance::position_dodgev(height=0.4))+
    theme_minimal()+
    geom_vline(xintercept=1,alpha=0.2)+
    labs(x="Odds Ratio (for MS)",y="Exposure",fill="Ethnicity")+
    scale_color_brewer(palette="Set2")+
    scale_x_log10(limits=c(0.01,100),breaks = c(0.1,1,10))
  
  png(paste0("./figs/exposure_simple_plot_het_ethnicity_full",plot_title,".png"),res=300,units="in",width=4,height=3)
  print(p2)
  dev.off()
  
  master_results_df <<- bind_rows(master_results_df, model_overall_df %>% mutate("cohort" = plot_title))
  
  strat_models = model_overall_df %>% mutate(estimate = 
                                paste0(round(or,2),
                                       " (",round(lower_ci,2),
                                       " - ",
                                       round(upper_ci,2),
                                       ")"
                                       )) %>%
    dplyr::select(y,eff_n,estimate,pval,pval_interaction,model) %>%
    pivot_wider(id_cols = c("y","eff_n"),names_from = model,values_from = c("estimate","pval","pval_interaction"))
  strat_models = strat_models %>% mutate_at(vars(contains("pval")),simplify_pval)
  
  
  write_csv(strat_models,paste0("./tables/stratified_exposure_het_ethnicity_full",plot_title,".csv"))
  write_csv(interaction_model_df,paste0("./tables/interactions_exposure_het_ethnicityfull",plot_title,".csv"))
}



run_interaction_analysis(dat=patient,plot_title="whole_cohort_indiv_ethnicity")
run_interaction_analysis(dat=patient,plot_title="whole_cohort_imd_indiv_ethnicity",
                         covars_to_test = c("index_age","gender","numeric_imd"))
run_interaction_analysis(dat=hes_ms_cohort,plot_title="hes_ms_indiv_ethnicity")


plot_dat = master_results_df %>% filter(y %in%  c("Overweight/obese",
                                                  "Prior head injury",
                                                  "Prior IM",
                                                  "Ever smoker",
                                                  "Vitamin D deficiency"))
plot_dat = plot_dat %>%
  mutate(cohort = 
           case_when(
             cohort == "hes_cohort" ~ "C: HES-confirmed MS \n(Binary ethnicity)",
             cohort == "hes_ms_indiv_ethnicity" ~ "F: HES-confirmed MS \n(Census ethnicity)",
             cohort == "whole_cohort" ~ "A: Whole cohort \n(Binary ethnicity)",
             cohort == "whole_cohort_imd" ~ "B: Whole cohort \n(Binary ethnicity,\nadjusted for IMD)",
             cohort == "whole_cohort_imd_indiv_ethnicity" ~  "E: Whole cohort \n(Census ethnicity,\nadjusted for IMD)",
             cohort == "whole_cohort_indiv_ethnicity" ~ "D: Whole cohort \n(Census ethnicity)"
           ))

plot_dat$cohort = factor(plot_dat$cohort, levels = c(
  "A: Whole cohort \n(Binary ethnicity)",
  "B: Whole cohort \n(Binary ethnicity,\nadjusted for IMD)",
  "C: HES-confirmed MS \n(Binary ethnicity)",
  "D: Whole cohort \n(Census ethnicity)",
  "E: Whole cohort \n(Census ethnicity,\nadjusted for IMD)",
  "F: HES-confirmed MS \n(Census ethnicity)"
),ordered = T)

          
p=ggplot(plot_dat,
       aes(or,y,fill=model))+
  geom_errorbarh(mapping=aes(xmin = lower_ci,xmax=upper_ci,y=y),height=0.01,position=ggstance::position_dodgev(height=0.4),alpha=0.7)+
  geom_point(alpha=0.7,shape=21,color="black",position=ggstance::position_dodgev(height=0.4))+
  theme_minimal()+
  facet_wrap(~cohort,nrow=2)+
  geom_vline(xintercept=1,alpha=0.2)+
  labs(x="Odds Ratio (for MS)",y="Exposure",fill="Ethnicity")+
  scale_color_brewer(palette="Set2")+
  scale_x_log10(limits=c(0.5,20))
png("./figs/fig3_plots.png",res=600,units="in",width=8,height=4)
p
dev.off()
