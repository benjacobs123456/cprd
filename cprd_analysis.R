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

get_year_from_date = function(x){
  as.numeric(format(x,format="%Y"))
}
datify = function(x){
  as.Date(x,format="%d/%m/%Y")
}
delta_dates = function(x){
  as.numeric(x)/365.25
}
make_hist = function(x){
  med = median(patient[[x]],na.rm=T) %>% round(2)
  ggplot(
    patient,
    aes(patient[[x]])
  ) +
    geom_density(alpha=0.5,fill="lightblue") +
    geom_vline(color="red",xintercept = med)+
    labs(x = x, y= "Density") +
    theme_minimal()+
    annotate("label",x=Inf,y=Inf,vjust = "inward", hjust = "inward",label = paste0("Median: ",med))
}
make_grouped_hist = function(x,group="gender"){
  ggplot(
    patient,
    aes(patient[[x]])
  ) +
    geom_density(alpha=0.5,fill="lightblue") +
    labs(x = x, y= "Density") +
    facet_wrap(~.data[[group]],nrow=2) +
    theme_minimal()
}

# plot function 
plot_fx = function(plot,filename,width=6,height=6){
  png(filename = filename, width = width, height = height, units="in",res=300)
  print(plot)
  dev.off()
}
# get tidy proportions
get_prop = function(x, dat){
  dat %>%
    count(.data[[x]]) %>%
    mutate(prop = n/sum(n)) %>%
    mutate(total = sum(n)) %>%
    mutate(pct = paste0(round(prop*100,2),"%")) %>%
    mutate(pastable = paste0(n,"/",total,"(",pct,")"))
}
# IPW functions 
# define function to format model outputs
extract_coefs_from_model_obj = function(model_obj){
  summ_df = summary(model_obj)$coefficients %>% data.frame()
  res_df = data.frame(
    var = rownames(summ_df),
    eff_n = model_obj$model %>% na.omit %>% nrow(),
    beta = summ_df$Estimate,
    se = summ_df$`Std..Error`,
    or = exp(summ_df$Estimate),
    lower_ci = exp(summ_df$Estimate - 1.96 * summ_df$`Std..Error`),
    upper_ci = exp(summ_df$Estimate + 1.96 * summ_df$`Std..Error`),
    pval = summ_df$`Pr...z..`
  )
  return(res_df)
}
get_weights = function(input_data, missing_col, miss_pred_cols){
  
  # arg checking 
  if(!is.data.frame(input_data)){
    stop(input_data, "must be a data frame")
  }
  if(is.null(input_data[[missing_col]])){
    stop(missing_col," must be present in the input data")
  }
  
  for(i in c(1:length(miss_pred_cols))){
    if(is.null(input_data[[miss_pred_cols[i]]])){
      stop(miss_pred_cols[i]," must be present in the input data")
    }
  }
  
  # copy data to new df
  data_for_ipw = input_data
  
  # find NAs
  message("finding NAs in ",missing_col)
  data_for_ipw = dplyr::mutate(data_for_ipw,
                               is_na = ifelse(is.na(.data[[missing_col]]),
                                              "missing",
                                              "nonmissing"))
  
  
  message("found ",nrow(filter(data_for_ipw,is_na == "missing"))," NAs")
  message("total n: ",nrow(data_for_ipw))
  
  # make is_na a factor (for regression) & set ref to missing
  data_for_ipw$is_na = relevel(factor(data_for_ipw$is_na),ref="missing")
  
  # ensure that no missingness among predictors
  nrow_data = nrow(data_for_ipw)
  for(i in c(1:length(miss_pred_cols))){
    data_for_ipw <<- data_for_ipw[!is.na(data_for_ipw[[miss_pred_cols[i]]]),]
  }
  if(nrow(data_for_ipw) != nrow_data){
    warning("Some predictors have missing data themselves. Excluding from the model.")
  }
  
  # regression of missing status using all predictors
  rhs = stringr::str_c(miss_pred_cols,collapse="+")
  formula = paste0("is_na ~ ",rhs)
  model = glm(data = data_for_ipw,
              formula,
              family=binomial(link="logit"))
  weights = model$fitted.values
  res_df = cbind(model$data,"p_nonmissing" = weights, "ipw" = 1/weights)
  return(list("missingness_model" = model,
              "ipw_results" = res_df))
}
check_missingness_predictors = function(input_data, missing_col, miss_pred_cols){
  
  # arg checking 
  if(!is.data.frame(input_data)){
    stop(input_data, "must be a data frame")
  }
  if(is.null(input_data[[missing_col]])){
    stop(missing_col," must be present in the input data")
  }
  for(i in c(1:length(miss_pred_cols))){
    if(is.null(input_data[[miss_pred_cols[i]]])){
      stop(miss_pred_cols[i]," must be present in the input data")
    }
  }
  
  # copy data to new df
  data_for_ipw = input_data
  
  # find NAs
  message("finding NAs in ",missing_col)
  data_for_ipw = dplyr::mutate(data_for_ipw,
                               is_na = ifelse(is.na(.data[[missing_col]]),
                                              "missing",
                                              "nonmissing"))
  
  
  message("found ",nrow(filter(data_for_ipw,is_na == "missing"))," NAs")
  message("total n: ",nrow(data_for_ipw))
  
  # make is_na a factor (for regression) & set ref to missing
  data_for_ipw$is_na = relevel(factor(data_for_ipw$is_na),ref="missing")
  
  # ensure that no missingness among predictors
  nrow_data = nrow(data_for_ipw)
  for(i in c(1:length(miss_pred_cols))){
    data_for_ipw <<- data_for_ipw[!is.na(data_for_ipw[[miss_pred_cols[i]]]),]
  }
  if(nrow(data_for_ipw) != nrow_data){
    warning("Some predictors have missing data themselves. Excluding from the model.")
  }
  
  # loop through predictors 
  overall_res = list()
  for(i in c(1:length(miss_pred_cols))){
    model = glm(data = data_for_ipw,
                is_na ~ data_for_ipw[[miss_pred_cols[i]]],
                family=binomial(link="logit"))
    
    res = dplyr::select(extract_coefs_from_model_obj(model),var,eff_n,or,lower_ci,upper_ci,pval)      
    res = filter(res,var != "(Intercept)")
    res$predictor = miss_pred_cols[i]
    overall_res[[i]] = res
  }
  overall_res = do.call("bind_rows",overall_res)
  overall_res$var = stringr::str_remove(overall_res$var,"data_for_ipw\\[\\[miss_pred_cols\\[i\\]\\]\\]")
  overall_res$predictor = paste0(overall_res$predictor,":",overall_res$var)
  overall_res$sig = ifelse(overall_res$pval < 0.05/nrow(overall_res),"*"," ")
  ggplot(overall_res,
         aes(or,predictor,label=sig))+
    geom_point()+
    geom_vline(xintercept = 1)+
    labs(x="OR for having complete (non-missing) data",y="Factor")+
    geom_errorbarh(mapping = aes(xmin = lower_ci,xmax = upper_ci,y=predictor),height=0.1)+
    annotate("label",x=Inf,y=Inf,hjust=1,vjust=1,label="More likely to\n have complete data")+
    theme_minimal()+
    geom_text(hjust=-2)+
    scale_x_log10(limits=c(0.05,10))+
    geom_label(mapping = aes(x=0.1,y=predictor,label=paste0("OR=",round(or,3))))
}
get_weights_stepwise = function(input_data, missing_col, miss_pred_cols){
  
  # arg checking 
  if(!is.data.frame(input_data)){
    stop(input_data, "must be a data frame")
  }
  if(is.null(input_data[[missing_col]])){
    stop(missing_col," must be present in the input data")
  }
  for(i in c(1:length(miss_pred_cols))){
    if(is.null(input_data[[miss_pred_cols[i]]])){
      stop(miss_pred_cols[i]," must be present in the input data")
    }
  }
  
  # copy data to new df
  data_for_ipw = input_data
  
  # find NAs
  message("finding NAs in ",missing_col)
  data_for_ipw = dplyr::mutate(data_for_ipw,
                               is_na = ifelse(is.na(.data[[missing_col]]),
                                              "missing",
                                              "nonmissing"))
  
  
  message("found ",nrow(filter(data_for_ipw,is_na == "missing"))," NAs")
  message("total n: ",nrow(data_for_ipw))
  
  # make is_na a factor (for regression) & set ref to missing
  data_for_ipw$is_na = relevel(factor(data_for_ipw$is_na),ref="missing")
  
  # ensure that no missingness among predictors
  nrow_data = nrow(data_for_ipw)
  for(i in c(1:length(miss_pred_cols))){
    data_for_ipw <<- data_for_ipw[!is.na(data_for_ipw[[miss_pred_cols[i]]]),]
  }
  if(nrow(data_for_ipw) != nrow_data){
    warning("Some predictors have missing data themselves. Excluding from the model.")
  }
  
  
  # regression of missing status using all predictors
  rhs = stringr::str_c(miss_pred_cols,collapse="+")
  formula = paste0("is_na ~ ",rhs)
  model = glm(data = data_for_ipw,
              formula,
              family=binomial(link="logit"))
  
  # now do stepwise regression
  stepwise_model = MASS::stepAIC(model, direction = "forward", 
                                 trace = F)
  weights = stepwise_model$fitted.values
  res_df = cbind(stepwise_model$data,"p_nonmissing" = weights, "ipw" = 1/weights)
  return(list("missingness_model" = stepwise_model,
              "ipw_results" = res_df))
}
plot_from_model_output = function(x){
  x = extract_coefs_from_model_obj(x)
  ggplot(filter(x,var != "(Intercept)"),
         aes(or,var))+
    geom_point()+
    geom_vline(xintercept = 1)+
    labs(x="OR for having non-missing data",y="Factor")+
    geom_errorbarh(mapping = aes(xmin = lower_ci,xmax = upper_ci,y=var),height=0.1)+
    annotate("label",x=Inf,y=Inf,hjust=1,vjust=1,label="More likely to\n have non-missing data")+
    theme_minimal()+
    scale_x_log10()
}
plot_missingness_prob = function(x){
  p1=ggplot2::ggplot(x,
                     aes(is_na,ipw))+
    geom_violin()+
    scale_y_log10()+
    labs(y="Weight",x="Missing data")+
    theme_minimal()+
    stat_boxplot(width=0.1)
  
  p2=ggplot2::ggplot(x,
                     aes(is_na,1-p_nonmissing))+
    geom_violin()+
    stat_boxplot(width=0.1)+
    labs(y="Fitted probability\nof missing data",x="Missing data")+
    theme_minimal()+
    scale_y_continuous(limits=c(0,1))
  
  gridExtra::grid.arrange(p1,p2)
}
hinkley_test = function(x){
  model_dat = cbind(x$data,fitted_prob_sq = x$fitted.values^2)
  # add fitted probs ^ 2 to original model call 
  original_formula = as.character(x$call)
  new_formula = paste0(original_formula[2]," + fitted_prob_sq") 
  
  # make new model   
  new_model = glm(data = model_dat,
                  new_formula,
                  family=binomial(link="logit"))
  pval = summary(new_model)$coefficients[rownames(summary(new_model)$coefficients)=="fitted_prob_sq",'Pr(>|z|)']  
  sig = ifelse(pval<0.05,"significant - missingness model poorly specified","not significant -  - missingness model well-specified")
  rounded_pval = round(pval,5)
  message("Hinkley's test pval: ",rounded_pval,"\nThis is ",sig)  
}
ipw_regression = function(input_data, outcome_pred_cols, outcome_col){
  
  res_df = input_data
  
  # make outcome into factor 
  res_df[[outcome_col]] = if(!is.factor(res_df[[outcome_col]])){
    warning(outcome_col," is not a factor. Coercing to one.")
    factor(res_df[[outcome_col]])
  } else {
    res_df[[outcome_col]]
  }
  
  # now do weighted regression
  design = survey::svydesign(id =~1,
                             weights =~ ipw,
                             data = res_df)
  
  # regression
  rhs = stringr::str_c(outcome_pred_cols,collapse="+")
  formula = paste0(outcome_col," ~ ",rhs)
  glm_res = survey::svyglm(formula, design = design,family = quasibinomial(link="logit"))
  summ_weighted_reg = summary(glm_res)
  
  # unweighted reg
  glm_unweighted = glm(data = res_df, formula,family = binomial(link="logit"))
  summ_unweighted_reg = summary(glm_unweighted)
  
  list("weights" = res_df,
       "model" = model,
       "unweighted_reg" = summ_unweighted_reg,
       "weighted_reg" = summ_weighted_reg)
}
ipw_match = function(input_data, outcome_col,id_col){
  
  res_df = input_data
  
  # make outcome into factor 
  res_df[[outcome_col]] = if(!is.factor(res_df[[outcome_col]])){
    warning(outcome_col," is not a factor. Coercing to one.")
    factor(res_df[[outcome_col]])
  } else {
    res_df[[outcome_col]]
  }
  
  
  # find cases 
  cases = dplyr::filter(res_df,.data[[outcome_col]]=="Case")
  controls = dplyr::filter(res_df,.data[[outcome_col]]=="Control")
  message("Cases: ",nrow(cases))
  message("Controls: ",nrow(controls))
  
  # thin both dfs to make things quickers
  cases = cases[,c(id_col,"ipw","yob")]
  controls = controls[,c(id_col,"ipw","yob")]
  
  
  # loop through cases and for each one find closest propensity-score control with same yob
  controls_to_keep = list()
  for(i in c(1:nrow(cases))){
    message("Matching MS case ",i," of ",nrow(cases))
    this_case = cases[i,]
    controls_for_this_case = controls[!controls$patid %in% unlist(controls_to_keep) & controls$yob == this_case$yob,]
    controls_for_this_case = dplyr::mutate(controls_for_this_case,ipw_delta = abs(ipw - this_case$ipw))
    controls_for_this_case = dplyr::slice_min(controls_for_this_case,ipw_delta) 
    controls_to_keep[[i]] = controls_for_this_case$patid
  }
  
  # filter original dataset to matched controls 
  output_df = res_df[res_df$patid %in% unlist(controls_to_keep) | res_df$patid %in% cases$patid,]
  return(output_df)
}


###########################
# read in data
###########################

# patient data 
patient = read_tsv("../../Data/Primary_care/Aurum/21_000677_Extract_Patient_001.txt",
                   col_types="cccccccccccc") %>%
  dplyr::select(1,2,4,5,8,10,12)

# read in codes which aren't of interest
codes_to_exclude = read_csv("../codelists/codes_to_exclude_aurum.csv", col_types = "ccccccccc")

# observation data 
clinical = purrr::map(c(1:41),function(x){
  number = str_pad(x,width=2,side="left",pad="0")
  # read in 
  # get rid of codes recorded prior to registration 
  # get rid of some commmon useless codes (eg 'text message sent')
  read_tsv(paste0("../../Data/Primary_care/Aurum/21_000677_Extract_Observation_0",number,".txt"),
           col_types = cols_only(
             patid = col_character(),
             obsdate = col_character(),
             medcodeid = col_character(),
             value = col_character()
           )
  ) %>% 
    filter(!medcodeid %in% codes_to_exclude$MedCodeId)
  
})
clinical = do.call("bind_rows",clinical)


# read in dictionary
dictionary = read_tsv("../../Dictionary Browsers/CPRD_CodeBrowser_202205_Aurum/CPRDAurumMedical.txt",
                      col_types = "ccccccccc") %>%
  rename("medcodeid" = MedCodeId) %>% 
  dplyr::select(medcodeid,CleansedReadCode,Term)

# add dictionary terms to clinical dataset
clinical = clinical %>% 
  left_join(dictionary,by="medcodeid")

# add in demographics
clinical = clinical %>%
  dplyr::select(-regstartdate) %>% 
  left_join(patient,by="patid")

# add in IMD (NB not available for everyone)
imd = read_tsv("../../Data/Linked_data/Aurum/patient_2019_imd_21_000677.txt",col_types="ccc")
patient = patient %>% left_join(imd,by=c("patid","pracid"))

# recode gender 
patient$gender = recode(patient$gender,
                        "2" = "Female",
                        "1" = "Male")

# recode year of birth as numeric
patient$yob = as.numeric(patient$yob)


# get LSOA info 
urban_rural = read_tsv("../../Data/Linked_data/Aurum/patient_urban_rural_21_000677.txt",col_types="ccc")
patient = patient %>% left_join(urban_rural,by=c("patid","pracid"))
patient$e2011_urban_rural = recode(patient$e2011_urban_rural,
                                   "1"="Urban",
                                   "2"="Rural")



rm(codes_to_exclude)
rm(urban_rural)
rm(imd)
rm(dictionary)

#################################
# look at age distros           #
#################################

# sort out dates in patient file
patient = patient %>%
  mutate(cprd_ddate = datify(cprd_ddate)) %>%
  mutate(dead = ifelse(!is.na(cprd_ddate),"dead","alive")) %>%
  mutate(regenddate = ifelse(
    is.na(regenddate),
    "01/05/2022",
    regenddate)) %>%
  mutate(regenddate = datify(regenddate)) %>%
  mutate(regstartdate = datify(regstartdate)) %>%
  mutate(fu_time = delta_dates(regenddate - regstartdate)) %>%
  mutate(rough_dob = datify(paste0("01/01/",yob))) %>%
  mutate(age_at_registration = delta_dates(regstartdate - rough_dob)) %>%
  mutate(age_at_deregistration = delta_dates(regenddate - rough_dob)) %>%
  mutate(age_at_death = delta_dates(cprd_ddate - rough_dob))

# get numeric year of reg
patient$numeric_year_of_reg = get_year_from_date(patient$regstartdate)

# explore ages
cols_to_plot = list("age_at_death","age_at_registration","age_at_deregistration","fu_time","numeric_year_of_reg")
plots = lapply(cols_to_plot,make_hist)

# save plots
plot_fx(do.call("grid.arrange",plots),"./figs/sf1.png")

#################################
# define MS outcome             #
#################################

# read in MS codelist
codelist = read_tsv("../codelists/aurum_codelist.txt",col_types = "cccccccc")

# derive first recorded MS code date
ms_first_report = clinical %>%
  filter(medcodeid %in% codelist$MedCodeId) %>%
  mutate(ms_code_date = datify(obsdate)) %>%
  group_by(patid) %>%
  slice_min(as.numeric(ms_code_date),n=1,with_ties = F) %>%
  dplyr::select(patid,ms_code_date) %>%
  ungroup()

# combine with main dataset
patient = patient %>%
  left_join(ms_first_report,by="patid") %>%
  mutate(MS_status = ifelse(is.na(ms_code_date),"Control","Case"))

# check that MS diagnosis after registration date
patient %>%
  mutate(ms_dx_pre_enrolment = ifelse(ms_code_date<regstartdate,"Yes","No")) %>%
  count(ms_dx_pre_enrolment)

## age at diagnosis
patient = patient %>%
  mutate(age_at_dx = delta_dates(ms_code_date - rough_dob))
patient$age_at_dx %>% summary
patient %>% group_by(gender) %>% summarise(median(age_at_dx,na.rm=T))

p=make_grouped_hist(x="age_at_dx",group="gender")+
  scale_x_continuous(limits=c(0,110))+
  labs(x="Age at first MS diagnostic\ncode report")
plot_fx(p,"./figs/age_ms_by_sex.png", width=4,height=4)

## verify diagnoses (see who has >1 code)
ms_diagnoses_verification = clinical %>%
  filter(medcodeid %in% codelist$MedCodeId) %>%
  count(patid) %>%
  rename("ms_diag_counts" = n)
patient = patient %>%
  left_join(ms_diagnoses_verification,by="patid")

patient %>% filter(MS_status=="Case") %>% count(ms_diag_counts>1) %>% mutate(n/sum(n)*100)
patient %>% filter(MS_status=="Case") %>% count(ms_diag_counts>2) %>% mutate(n/sum(n)*100)
patient %>% filter(MS_status=="Case") %>% count(ms_diag_counts>3) %>% mutate(n/sum(n)*100)

stats = patient %>% filter(MS_status=="Case") %>% summarise_at(vars(ms_diag_counts),.funs=c("mean","median","IQR","sd"))
p=ggplot(patient %>% filter(MS_status=="Case"),aes(ms_diag_counts))+geom_density(fill="lightblue",alpha=0.8)+scale_x_log10()+theme_minimal()+labs(x="MS diagnostic code counts \nper individual")
patient$n_ms_diag_codes = cut2(patient$ms_diag_counts,cuts=c(1,2,5,10))
levels(patient$n_ms_diag_codes) = c("1","2-5","6-10",">10")

p1=ggplot(patient %>% filter(MS_status=="Case"),aes(n_ms_diag_codes))+
  geom_bar(fill="lightblue",alpha=0.8,color="black")+
  theme_minimal()+
  labs(x="MS diagnostic code counts \nper individual")

plot_fx(p1,width=3,height=3,"./figs/sf3.png")

# HES validation of MS diagnosis 

# HES-OP
ms_hes_codes_overall = list()
for(year in seq(2003,2020,by=1)){
  hes_df = read_tsv(paste0("../../Data/Linked_data/Aurum/hesop_clinical_",year,"_21_000677.txt"),
                    col_types = cols(.default = "c"))
  ms_hes_codes = hes_df %>% filter_at(
    vars(contains("diag")),
    any_vars(.=="G35X")) %>% 
    dplyr::select(patid,attendkey,HES_yr)
  ms_hes_codes_overall[[length(ms_hes_codes_overall)+1]] = ms_hes_codes
}
# bind
ms_hes_codes_overall = do.call("bind_rows",ms_hes_codes_overall)

# keep earliest record
ms_hes_codes_overall = ms_hes_codes_overall %>% 
  group_by(patid) %>%
  mutate(HES_yr = as.numeric(HES_yr)) %>%
  slice_min(HES_yr,with_ties = F)

# HES-APC 
hes_apc = read_tsv("../../Data/Linked_data/Aurum/hes_diagnosis_epi_21_000677.txt",
                   col_types = cols(.default = "c"))
hes_op_patlist = read_tsv("../../Data/Linked_data/Aurum/hesop_patient_21_000677.txt",
                          col_types = cols(.default = "c"))

# collate list of all unique patids in HES 

overall_hes_ids = hes_apc %>% 
  distinct(patid) %>%
  bind_rows(hes_op_patlist %>% distinct(patid)) %>%
  distinct(patid)

# see what % of people are in HES
in_hes_n = overall_hes_ids %>% nrow
total_n = patient %>% distinct(patid) %>% nrow
in_hes_n / total_n

ms_hes_codes_apc = hes_apc %>% 
  filter(grepl("G35",(ICD))) %>% 
  dplyr::select(patid,epiend) %>%
  mutate(HES_yr = datify(epiend)) %>%
  mutate(HES_yr = get_year_from_date(HES_yr)) %>%
  dplyr::select(patid,HES_yr) %>%
  group_by(patid) %>%
  slice_min(HES_yr,with_ties = F)

# merge HES-OP and HES-APC
ms_hes_codes_overall = ms_hes_codes_overall %>% 
  bind_rows(ms_hes_codes_apc) %>%
  mutate(HES_yr = as.numeric(HES_yr)) %>%
  dplyr::select(-attendkey) %>%
  group_by(patid) %>%
  slice_min(HES_yr, with_ties = F) 

p=ggplot(ms_hes_codes_overall,
       aes(HES_yr))+
  geom_bar(fill="lightblue",color="black")+
  scale_x_continuous(breaks = seq(1996,2022,by=2))+
  theme_minimal()+
  labs(x="Year",y="No. of new recorded MS diagnoses")+
  ggtitle("HES MS diagnostic codes per year")

plot_fx(p,"./figs/hes_ms_dx_per_year.png",width=6,height=4)  

# check how many MS patients have a code in HES
patient = patient %>% 
  mutate(ms_code_source =
           ifelse(MS_status=="Case" & patid %in% ms_hes_codes_overall$patid,
                  "primary_care_and_hes",
                  ifelse(MS_status=="Case" & !patid %in% ms_hes_codes_overall$patid,
                         "primary_care_only",NA)))

# add var indicating if person is in HES at all 
patient = patient %>% 
  mutate(in_hes = ifelse(patid %in% overall_hes_ids$patid,"in_hes","not_in_hes"))

# count what % of MS patients also have a HES code
patient %>% 
  filter(MS_status=="Case") %>% 
  count(ms_code_source) %>% 
  mutate(prop = n/sum(n),
         total = sum(n))

patient %>% 
  filter(MS_status=="Case" & in_hes == "in_hes") %>% 
  count(ms_code_source) %>% 
  mutate(prop = n/sum(n),
         total = sum(n))

ms_hes_not_pc = patient %>% 
  filter(MS_status=="Control" & patid %in% ms_hes_codes_overall$patid) %>%
  left_join(ms_hes_codes_overall,by="patid") %>%
  mutate(hes_dx_within_cprd_reg = ifelse(HES_yr < format(regstartdate,"%Y") | HES_yr > format(regenddate,"%Y") ,"dx_outside_cprd","dx_during_cprd"))

ms_hes_not_pc %>% dplyr::select(HES_yr,regstartdate,hes_dx_within_cprd_reg)
ms_hes_not_pc$hes_dx_within_cprd_reg %>% table
ms_hes_not_pc$fu_time %>% median
patient$fu_time %>% median

# exclude these individuals from all downstream analyses
patient = patient %>% filter(!patid %in% ms_hes_not_pc$patid)

# assign pseudodiagnosis date for each matched control
matching_file = read_tsv("../../Data/Primary_care/Aurum/21_000677_CPRD_Aurum_matching_file.txt",col_types="cccccccc") %>%
  dplyr::select(control_patid,case_index) %>%
  rename("patid" = control_patid)

patient = patient %>%
  left_join(matching_file,by="patid")

## age at pseudodiagnosis
patient = patient %>%
  mutate(index_code_date = datify(case_index)) %>% 
  mutate(age_at_pseudodx = delta_dates(index_code_date - rough_dob))

# check that only cases have a dx date & vice-versa
patient %>% 
  filter(is.na(age_at_pseudodx)) %>%
  count(MS_status)

patient %>% 
  filter(is.na(age_at_dx)) %>%
  count(MS_status)

# define variable which encompasses age at dx/pseudodx
# also define FU time prior to and after this date
patient = patient %>%
  mutate(index_age = ifelse(MS_status=="Case",age_at_dx,age_at_pseudodx)) %>%
  mutate(years_data_before_index = index_age - age_at_registration) %>%
  mutate(years_data_after_index = age_at_deregistration - index_age) 

# look at MS diagnoses over time 
patient %>% glimpse

p=ggplot(patient %>% filter(!is.na(ms_code_date)) %>%
           mutate(ms_code_source = ifelse(ms_code_source=="primary_care_only","PC","PC+HES")),
         aes(as.numeric(format(ms_code_date,format="%Y")),fill=ms_code_source))+
  geom_bar(color="black")+
  theme_minimal()+
  labs(x="Year of MS report",y="N",fill="MS diagnostic\ncode source")+
  scale_fill_brewer(palette="Paired")

plot_fx(p,"./figs/ms_codes_over_time.png",width=5,height=3)

# cross-check against referrals
neuro_ref_codelist = read_csv("../codelists/neuro_ref_codelist.csv")
neuro_referrals = clinical %>% 
  filter(medcodeid %in% neuro_ref_codelist$MedCodeId)

# find earliest referral
earliest_neuro_ref = neuro_referrals %>%
  mutate(refdate = as.numeric(datify(obsdate))) %>%
  group_by(patid) %>%
  slice_min(refdate,n=1,with_ties = F)

patient = patient %>%
  mutate(any_ref_to_neurologist = ifelse(patid %in% earliest_neuro_ref$patid,"ref_to_neuro","not_referred"))
patient %>%
  group_by(MS_status) %>%
  count(any_ref_to_neurologist) %>%
  mutate(prop = n/sum(n), total = sum(n))

# check that referral precedes diagnosis 

patient = patient %>% left_join(
  earliest_neuro_ref %>% 
    dplyr::select(patid,obsdate) %>%
    rename("neuro_ref_date" = obsdate),by="patid") %>%
  mutate(neuro_ref_date = datify(neuro_ref_date))


# calculate time from neuro ref to dx
patient = patient %>% 
  mutate(time_from_neuro_ref_to_dx = delta_dates(ms_code_date - neuro_ref_date)) %>%
  mutate(ms_code_after_neuro_ref = ifelse(ms_code_date > neuro_ref_date,"yes","no")) 

get_prop(x = "ms_code_after_neuro_ref", dat = patient %>% filter(MS_status=="Case"))

#################################
# define ethnicity              #
#################################

# read in HES ethnicity data
hes_ethnicity =  read_tsv("../../Data/Linked_data/Aurum/hes_patient_21_000677.txt", col_types="cccccc")
hes_ae_ethnicity =  read_tsv("../../Data/Linked_data/Aurum/hesae_patient_21_000677.txt", col_types="cccccc")
hes_op_ethnicity =  read_tsv("../../Data/Linked_data/Aurum/hesop_patient_21_000677.txt", col_types="cccccc")
overall_hes_ethnicity = bind_rows(hes_ethnicity,hes_op_ethnicity,hes_ae_ethnicity)

# get rid of NAs and unknown
overall_hes_ethnicity = overall_hes_ethnicity %>%
  filter(!is.na(gen_ethnicity)) %>%
  filter(gen_ethnicity != "Unknown") %>%
  distinct(patid,.keep_all = T)

# join with main dataset
patient = patient %>%
  left_join(
    overall_hes_ethnicity %>% 
      dplyr::select(patid,gen_ethnicity),
    by="patid")

# see ethnicity missingness in HES
get_prop(dat = patient %>% mutate(eth_miss = is.na(gen_ethnicity)) %>% group_by(MS_status),
         "eth_miss")

# recode ethnicity as simple
patient$gen_ethnicity = recode(patient$gen_ethnicity,
                               "Bangladesi" = "Asian",
                               "Bl_Afric" = "Black",
                               "Bl_Carib" = "Black",
                               "Bl_Other" = "Black",
                               "Chinese" = "Mixed/Other",
                               "Indian" = "Asian",
                               "Mixed" = "Mixed/Other",
                               "Oth_Asian" = "Asian",
                               "Other" = "Mixed/Other",
                               "Pakistani" = "Asian",
                               "White" = "White")


# see if missing ethnicity codes can be derived from primary care
ethnicity_codes = read_csv(
  "../codelists/ethnicity_codelist.csv",
  col_types="cccccc")

ethnicity_df = clinical %>%
  filter(medcodeid %in% ethnicity_codes$MedCodeId)
ethnicity_df %>% distinct(patid)

# take the earliest recorded code for each person
ethnicity_df = ethnicity_df  %>%
  left_join(ethnicity_codes %>% rename("medcodeid" = MedCodeId),by="medcodeid") %>%
  mutate(date_ethnicity_recorded = datify(obsdate)) %>%
  group_by(patid) %>% 
  slice_min(as.numeric(date_ethnicity_recorded),with_ties = F) %>% 
  dplyr::select(patid,ethnicity)

# join with main dataset
patient = patient %>% left_join(ethnicity_df,by="patid")

# define source of ethnicity report
patient = patient %>% 
  mutate(source_of_ethnicity_report = 
           ifelse(!is.na(gen_ethnicity) & !is.na(ethnicity),
                  "HES + PC",
                  ifelse(!is.na(gen_ethnicity) & is.na(ethnicity),
                         "HES only",
                         ifelse(is.na(gen_ethnicity) & !is.na(ethnicity),
                                "PC only",
                                "Neither"))))


get_prop(dat=patient,x="source_of_ethnicity_report")
get_prop(dat=patient %>% filter(numeric_year_of_reg>=2000),x="source_of_ethnicity_report")
get_prop(dat=patient %>% filter(numeric_year_of_reg>=2006),x="source_of_ethnicity_report")
get_prop(dat=patient %>% filter(numeric_year_of_reg>=1990),x="source_of_ethnicity_report")

# calculate agreement 
patient %>% 
  filter(!is.na(gen_ethnicity) & !is.na(ethnicity)) %>%
  group_by(gen_ethnicity) %>%
  count(ethnicity) %>%
  mutate(n/sum(n)) %>%
  filter(gen_ethnicity == ethnicity)

p=ggplot(patient,
         aes(MS_status,fill=source_of_ethnicity_report))+
  geom_bar(position="fill",colour="black")+
  theme_minimal()+
  scale_fill_brewer(palette = "Paired")+
  scale_y_continuous(labels = scales::percent)+
  labs(x="MS status",y="% of patients with \nethnicity data from \neach source",
       fill="source of ethnicity\n data")+
  ggtitle("A")


p1=ggplot(patient %>%
            filter(source_of_ethnicity_report == "HES + PC"),
          aes(gen_ethnicity,fill=ethnicity))+
  geom_bar(position="fill",colour="black")+
  theme_minimal()+
  scale_fill_brewer(palette = "Set3")+
  scale_y_continuous(labels = scales::percent)+
  labs(x="HES ethnicity",y="% of patients",
       fill="Primary care ethnicity")+
  theme(axis.text.x = element_text(angle=90,vjust=0.3))+
  facet_wrap(~MS_status)+
  ggtitle("B")

# look at missingness of CPRD ethnicity data 
missing_ethnicity_data_df = patient
missing_ethnicity_data_df$decade_of_recruitment = cut2(
  as.numeric(format(missing_ethnicity_data_df$regstartdate,format="%Y")),
  cuts=seq(1880,2020,by=10))
levels(missing_ethnicity_data_df$decade_of_recruitment) = paste0(seq(1880,2020,by=10),"s")

plot_dat = get_prop(dat = missing_ethnicity_data_df %>% group_by(decade_of_recruitment, MS_status),
         x = "source_of_ethnicity_report") 
p2=ggplot(plot_dat %>% filter(source_of_ethnicity_report=="Neither"),
          aes(decade_of_recruitment,prop*100,fill=source_of_ethnicity_report))+
  geom_col(color="black",fill="orange",alpha=0.8)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90,vjust = 0.3))+
  labs(x="Decade of CPRD registration",y="% of individuals\nwith missing ethnicity data",fill="Source of ethnicity report")+
  facet_wrap(~MS_status)+
  ggtitle("C")

mat = rbind(c(1,2),c(3,3))
plot_fx(grid.arrange(p,p1,p2,layout_matrix = mat),
        "./figs/sf5.png",height=6,width=8)


# define ethnicity as HES ethnicity if present, and if not as CPRD ethnicity

patient = patient %>%
  mutate(ethnicity_pc_hes_hybrid = 
           ifelse(!is.na(gen_ethnicity),
                  gen_ethnicity,
                  ethnicity))
table(patient$ethnicity_pc_hes_hybrid,patient$MS_status)

#################################
# descriptive stats
#################################

# define function which makes demographics table 
make_demographics = function(
  input_data = patient,
  response_var = "MS_status",
  cat_vars = c("gender","ethnicity_pc_hes_hybrid","dead","e2019_imd_5","e2011_urban_rural","in_hes"),
  cont_vars = c("yob","age_at_registration","fu_time","index_age","years_data_before_index","years_data_after_index"),
  outfile = "basic")
{
  
  # initialise empty list 
  datasets_to_bind <<- list()
  
  # check levels in outcome 
  response_levels = input_data[[response_var]] %>% factor %>% levels
  
  # get n  
  demog_df = input_data %>% 
    group_by(.data[[response_var]]) %>% 
    count() %>% 
    t() 
  
  demog_df = cbind(rownames(demog_df),demog_df)
  demog_df = cbind(demog_df,c("P_value","_"))
  datasets_to_bind[[1]] <<- demog_df
  
  # process categorical vars 
  for(i in c(1:length(cat_vars))){
    this_col = cat_vars[i]
    message("processing ",this_col)
    
    res = input_data %>% 
      group_by(.data[[response_var]]) %>% 
      filter(!is.na(.data[[this_col]])) %>%
      count(.data[[this_col]]) %>%
      mutate(prop = n/sum(n)*100) %>%
      mutate(x = paste0(n," (",round(prop,1),"%)")) %>%
      pivot_wider(id_cols = .data[[response_var]],
                  names_from = .data[[this_col]],
                  values_from = x) %>% 
      t()
    
    res_2 = table(input_data[[response_var]],input_data[[this_col]])
    pval = round(chisq.test(res_2)$p.value,2)
    
    res = res[-1,]
    res = cbind(rownames(res),res)
    
    res = cbind(res,c(pval,rep("_",nrow(res)-1)))
    print(res)
    rownames(res)=NULL
    datasets_to_bind[[length(datasets_to_bind)+1]] <<- res
    
  }  
  # process continuous vars 
  for(i in c(1:length(cont_vars))){
    this_col = cont_vars[i]
    message("processing ",this_col)
    res = input_data %>% 
      filter(!is.na(.data[[this_col]])) %>%
      group_by(.data[[response_var]]) %>% 
      summarise_at(this_col,list(mean = mean,sd = sd)) %>%
      mutate(x = paste0(round(mean,1)," (",round(sd,1),")")) %>%
      dplyr::select(1,4) %>% 
      t()
    
    a = input_data %>% 
      filter(.data[[response_var]] == response_levels[1]) 
    b = input_data %>% 
      filter(.data[[response_var]] == response_levels[2]) 
    pval = round(t.test(a[[this_col]],b[[this_col]])$p.value,2)    
    res = res[-1,]
    res = c(this_col,res,pval)
    print(res)
    datasets_to_bind[[length(datasets_to_bind)+1]] <<- res
    
  }
  message("Binding")
  
  demog_df = do.call("rbind",datasets_to_bind)
  message(nrow(demog_df))
  print(demog_df)
  
  # get into nice format
  message("formatting")
  demog_df = data.frame(demog_df)
  colnames(demog_df) = demog_df[1,]
  demog_df = demog_df[-1,]
  rownames(demog_df) = NULL
  
  
  message("saving")
  write_csv(demog_df,paste0("./tables/",outfile,"_demographics.csv"))
  print(demog_df)
}

# now make demographics tables
make_demographics(input_data = patient, outfile = "whole_cohort")

make_demographics(input_data = patient %>% filter(
  as.numeric(format(regstartdate,format="%Y"))>=2006
), outfile = "post_2006")
make_demographics(input_data = patient %>% filter(years_data_before_index>=5), outfile = "5y_pre_index")

# calculate person-years & data per person
patient %>%
  summarise(sum(fu_time))
patient %>%
  group_by(MS_status) %>% 
  summarise(sum(fu_time))
patient %>%
  group_by(MS_status) %>% 
  count(years_data_before_index>=5) %>%
  mutate(prop = n/sum(n))

#################################
# define other exposures
#################################


#################################
# 1. BMI
#################################

# filter observation files 
weight_recordings = clinical %>% filter(medcodeid=="253677014")
height_recordings = clinical %>% filter(medcodeid=="253669010")

# see how many readings for each 
weight_recordings %>% nrow
weight_recordings %>% distinct(patid) %>% nrow
height_recordings %>% nrow
height_recordings %>% distinct(patid) %>% nrow

# qc on height & weight recordings 
weight_recordings = weight_recordings %>%
  mutate(value = as.numeric(value)) %>%
  filter(!is.na(value)) %>%
  filter(value > 20)
height_recordings = height_recordings %>%
  mutate(value = as.numeric(value)) %>%
  filter(!is.na(value)) %>%
  filter(value >= 121 & value <= 214)

# join height and weight
height_weight = weight_recordings %>% 
  bind_rows(height_recordings) %>%
  dplyr::select(patid,obsdate,Term,value) %>% 
  pivot_wider(id_cols = c(patid,obsdate),
              names_from = Term, 
              values_from = value,
              values_fn = mean)

# exclude NAs 
height_weight_nonmissing = height_weight %>% 
  filter(!is.na(`Body weight`) & !is.na(`Standing height`))
height_weight_nonmissing %>%  nrow
height_weight_nonmissing %>% distinct(patid) %>% nrow

# find missing data
ht_missing = patient %>% filter(!patid %in% height_recordings$patid)
wt_missing = patient %>% filter(!patid %in% weight_recordings$patid)
both_missing = patient %>% filter(!patid %in% weight_recordings$patid & !patid %in% height_recordings$patid)
either_missing = patient %>% filter(!patid %in% weight_recordings$patid | !patid %in% height_recordings$patid)

# for people with missing data, impute with nearest height 
height_weight_missing = height_weight %>% 
  filter(!patid %in% height_weight_nonmissing$patid) %>%
  filter(!is.na(`Body weight`) & is.na(`Standing height`))

# join with height recordings
height_weight_missing = height_weight_missing %>%
  filter(patid %in% height_recordings$patid) %>%
  left_join(height_recordings,by="patid") %>%
  mutate(wt_date = as.Date(obsdate.x,format="%d/%m/%Y")) %>%
  mutate(ht_date = as.Date(obsdate.y,format="%d/%m/%Y")) %>%
  mutate(delta_date = abs(as.numeric(ht_date - wt_date))) %>% 
  filter(!is.na(value)) %>%
  group_by(patid,wt_date) %>%
  slice_min(delta_date,n=1) %>%
  ungroup() %>%
  dplyr::select(patid,obsdate.x,`Body weight`,value,delta_date) %>% 
  rename("obsdate" = obsdate.x,
         `Standing height` = value)

height_weight_missing %>% nrow
height_weight_missing %>% distinct(patid) %>% nrow
ggplot(height_weight_missing,aes(delta_date/365))+
  geom_histogram()+
  geom_vline(xintercept=5)


# calculate BMI 
bmi_overall = height_weight_nonmissing %>% 
  bind_rows(height_weight_missing %>% dplyr::select(-delta_date)) %>%
  mutate(bmi = `Body weight`/((`Standing height`/100)^2))
summary(bmi_overall$bmi)

# fill in missing BMIs with CPRD-derived BMIs 
bmi_raw_values = clinical %>% filter(medcodeid=="100716012")
bmi_raw_values %>% distinct(patid)
combo_bmi = bmi_raw_values %>% 
  dplyr::select(patid,obsdate,value) %>%
  filter(!patid %in% bmi_overall$patid) %>%
  rename("bmi" = value) %>%
  mutate(bmi = as.numeric(bmi)) %>%
  filter(!is.na(bmi)) %>%
  bind_rows(
    bmi_overall %>% dplyr::select(1,2,5)
  )

# get rid of unrealistic readings
combo_bmi = combo_bmi %>% filter(bmi >=5 & bmi <=100)
hist(combo_bmi$bmi)

# examine validity of BMI 
combo_bmi %>% distinct(patid) %>% nrow

# get earliest bmi recording 
timepoint_bmi = combo_bmi %>% 
  mutate(obs_dt = datify(obsdate)) %>%
  filter(obs_dt < datify("01/05/2022")) %>%
  left_join(patient %>% dplyr::select(patid,rough_dob,regstartdate),by="patid") %>%
  filter(obs_dt >= regstartdate) %>%
  mutate(age_at_bmi = delta_dates(obs_dt - rough_dob)) %>%
  filter(age_at_bmi > 0)

timepoint_bmi$bmi_epoch = cut2(timepoint_bmi$age_at_bmi,cuts = seq(0,80,by=8))

timepoint_bmi = timepoint_bmi %>%
  group_by(patid,bmi_epoch) %>%
  summarise(bmi = mean(bmi,na.rm=T)) %>%
  pivot_wider(id_cols = patid, names_from = bmi_epoch, values_from = bmi, names_prefix = "bmi")

# get earliest BMI for each person after age 16
earliest_bmi = combo_bmi %>% 
  mutate(obs_dt = datify(obsdate)) %>%
  filter(obs_dt < datify("01/05/2022")) %>%
  left_join(patient %>% dplyr::select(patid,rough_dob,regstartdate),by="patid") %>%
  filter(obs_dt >= regstartdate) %>%
  mutate(age_at_bmi = delta_dates(obs_dt - rough_dob)) %>% 
  filter(age_at_bmi>=16) %>%
  group_by(patid) %>%
  slice_min(age_at_bmi,n=1)


get_prop(dat = earliest_bmi %>% ungroup %>%
           mutate("bmi_before_30" = ifelse(age_at_bmi<30,"yes","no")),
         x = "bmi_before_30") 

# combine with main dataset 
patient = patient %>%
  left_join(
    earliest_bmi %>% 
      dplyr::select(patid,obs_dt,bmi,age_at_bmi) %>%
      rename("bmi_date" = obs_dt),
    by="patid"
  ) %>%
  left_join(timepoint_bmi,by="patid")

# see whether BMI recorded prior to index age (1 year buffer)
patient = patient %>%
  mutate(bmi_to_index = index_age - age_at_bmi) %>%
  mutate(bmi_before_index = ifelse(
    !is.na(bmi) & bmi_to_index>=1,
    "yes",
    "no"))
get_prop(x="bmi_before_index",dat = patient)


# recode as NA any BMI recordings after index date 
patient = patient %>%
  mutate(bmi = 
           ifelse(bmi_before_index == "no",
                  NA,
                  bmi))
get_prop(dat = patient,x="bmi_before_index")
make_hist("bmi_to_index")

patient %>% filter(bmi_before_index=="yes") %>%
  summarise_at(vars(bmi_to_index),.funs = c("median","IQR"))
# make BMI categories
patient$bmi_category = cut2(patient$bmi,
                            cuts = c(18.5,25,30,40))
levels(patient$bmi_category) = c("Underweight","Healthy","Overweight","Obese","Morbidly obese")
patient$bmi_category = factor(patient$bmi_category,ordered=T)
patient$bmi_category_simple = recode(
  patient$bmi_category,
  "Underweight" = "Healthy/Underweight" ,
  "Healthy" = "Healthy/Underweight",
  "Overweight" = "Overweight/obese",
  "Obese" = "Overweight/obese",
  "Morbidly obese" = "Overweight/obese",
)

# make BMI z score 
patient = patient %>% mutate(
  bmi_z = (bmi - mean(bmi,na.rm=T)) / sd(bmi,na.rm = T))

# look at % in each category
get_prop(dat = patient %>% filter(!is.na(bmi_category)), x= "bmi_category")

# make basic plots
p=make_grouped_hist(x="bmi", group="gender")+ggtitle("A")

p1=ggplot(patient %>% filter(!is.na(bmi_category)),
       aes(MS_status,fill=bmi_category))+
  geom_bar(position="fill",col="black")+
  theme_minimal()+
  scale_fill_brewer(palette="Paired")+
  labs(x="MS status",y="Proportion of individuals in\neach BMI category",fill="BMI")+
  ggtitle("B")

make_missingness_plot_over_time = function(x, dat = patient){
  plot_dat = dat
  plot_dat$decade_of_recruitment = cut2(
    as.numeric(format(plot_dat$regstartdate,format="%Y")),
    cuts=seq(1880,2020,by=10))
  levels(plot_dat$decade_of_recruitment) = paste0(seq(1880,2020,by=10),"s")
  
  plot_dat = plot_dat %>% mutate(
    missing_x = ifelse(is.na(.data[[x]]),"missing","nonmissing")
  )
  plot_dat_summ = get_prop(dat = plot_dat %>% group_by(decade_of_recruitment, MS_status),
                      x = "missing_x") 
  p=ggplot(plot_dat_summ %>% filter(missing_x == "missing"),
            aes(decade_of_recruitment,prop*100))+
    geom_col(color="black",fill="orange",alpha=0.8)+
    theme_minimal()+
    theme(axis.text.x = element_text(angle=90,vjust = 0.3))+
    labs(x="Decade of CPRD registration",y="% of individuals\nwith missing data")+
    facet_wrap(~MS_status)
  return(p)
}

p2=make_missingness_plot_over_time("bmi")+ggtitle("C")
mat = rbind(c(1,2,2),c(3,3,3))
plot_fx(grid.arrange(p,p1,p2,layout_matrix =mat),"./figs/bmi_missing.png")

#################################
# 2. Smoking 
#################################

# smoking
## define start time
smoking_codes = read_csv("../codelists/smoking_status.csv", col_types = cols(.default = "c") )

# see how many people have any smoking recorded
any_smoking_data = clinical %>%
  filter(medcodeid %in% smoking_codes$medcodeid) %>%
  left_join(smoking_codes,by="medcodeid") %>%
  distinct(patid,.keep_all = T)

get_prop(dat = patient %>% 
           mutate(smok_dat = ifelse(patid %in% any_smoking_data$patid,"yes","no")),
         x="smok_dat")
get_prop(dat = any_smoking_data,x="smoking_status")


smoking_df = clinical %>%
  filter(medcodeid %in% smoking_codes$medcodeid) %>%
  left_join(smoking_codes,by="medcodeid") %>%
  left_join(patient %>% dplyr::select(patid,regstartdate,rough_dob,index_age),by="patid") %>%
  filter(datify(obsdate)>=regstartdate) %>%
  mutate(age_at_smok_code = delta_dates(datify(obsdate) - rough_dob)) %>%
  filter(age_at_smok_code <= (index_age - 1))

# derive smoking code date
smoking_df = smoking_df %>%
  group_by(patid) %>%
  slice_max(age_at_smok_code,n=1,with_ties = F) %>%
  dplyr::select(patid,age_at_smok_code, smoking_status) %>%
  ungroup()

# join with main dataset - define smoking status prior to index date
patient = patient %>%
  left_join(smoking_df,by="patid")

# crudely estimate no. of smoking years prior to index 
patient = patient %>%
  mutate(smoking_prior_to_index = ifelse(smoking_status == "smoker",
                                         index_age - age_at_smok_code,
                                         NA))

get_prop(dat = patient %>% 
           mutate(smok_miss = ifelse(is.na(smoking_status),"miss","nonmiss")), 
         x="smok_miss")
get_prop(dat = patient %>% filter(!is.na(smoking_status)), x="smoking_status")


hist(patient$smoking_prior_to_index)

# inspect discordant smoking data
discordant_smok_data = patient %>% 
  filter(patid %in% never_smokers$patid & patid %in% ever_smokers$patid) %>%
  dplyr::select(patid,smoking_status) %>%
  left_join(never_smokers,by="patid") %>%
  left_join(ever_smokers,by="patid") %>% 
  filter(!is.na(never_smok_date)) %>%
  filter(!is.na(ever_smok_date)) %>%
  filter(!is.na(smoking_status)) %>%
  mutate(never_dat = datify(never_smok_date))%>%
  mutate(ever_dat = datify(ever_smok_date)) %>%
  mutate(delta = delta_dates(ever_dat - never_dat)) %>%
  mutate(never_before_ever = ifelse(never_dat < ever_dat,"never_first","ever_first"))

get_prop(x = "never_before_ever", dat = discordant_smok_data)
hist(discordant_smok_data$delta)
summary(discordant_smok_data$delta)
discordant_smok_data %>% arrange(desc(delta))

#################################
# 3. IM 
#################################

# read in codelist 
im_codes = read_csv("../codelists/im_codelist.csv",
                    col_types = "ccccccccc")

# im data 
im_df = clinical %>%
  filter(medcodeid %in% im_codes$MedCodeId) %>%
  filter(medcodeid %in% c("406382010","406380019","288197011","288198018","404460017","820971000006110")) %>%
  mutate(im_dx_date = as.Date(obsdate,format="%d/%m/%Y")) %>%
  filter(im_dx_date < as.Date("01/05/2022",format="%d/%m/%Y")) %>%
  left_join(patient %>% select(patid,rough_dob),by="patid") %>%
  mutate(age_at_im = as.numeric(im_dx_date - rough_dob)/365) %>%
  filter(age_at_im>0) %>%
  group_by(patid) %>%
  slice_min(age_at_im,n=1,with_ties = F) %>%
  select(patid,im_dx_date,age_at_im)

# join with main dataset
patient = patient %>%
  left_join(im_df,by="patid") %>%
  mutate(im_status = ifelse(!is.na(age_at_im) & 
                              age_at_im < (index_age-1),
                            "IM",
                            "no_IM"))

patient %>%
  count(im_status) %>%
  mutate(prop = n/sum(n))

patient %>% group_by(MS_status) %>%
  count(im_status) %>%
  mutate(prop = n/sum(n))

summary(patient$age_at_im)
hist(patient$age_at_im)



#################################
# 4. Vit D 
#################################

# read in codelist 
vitd_def_codes = read_csv("../codelists/vitd_def_codelist.csv",
                          col_types = "ccccccccc")

# vitd  data 
vd_df = clinical %>%
  filter(medcodeid %in% vitd_def_codes$MedCodeId) %>%
  mutate(vd_def_dx_date = as.Date(obsdate,format="%d/%m/%Y")) %>%
  filter(vd_def_dx_date < as.Date("01/05/2022",format="%d/%m/%Y")) %>%
  left_join(patient %>% select(patid,rough_dob),by="patid") %>%
  mutate(age_at_vd_def = as.numeric(vd_def_dx_date - rough_dob)/365) %>%
  filter(age_at_vd_def>0) %>%
  group_by(patid) %>%
  slice_min(age_at_vd_def,n=1,with_ties = F) %>%
  select(patid,vd_def_dx_date,age_at_vd_def)

# join with main dataset
patient = patient %>%
  left_join(vd_df,by="patid") %>%
  mutate(vd_def_status = ifelse(!is.na(age_at_vd_def) & 
                                  age_at_vd_def < (index_age-1),
                                "vd_def",
                                "vd_replete"))

table(patient$vd_def_status,patient$MS_status)

# vitamin D tests
vd_test_codes = c("457822017",
                  "145751000006117",
                  "457821012")
vitd_tests = clinical %>% 
  filter(medcodeid %in% vd_test_codes) %>%
  filter(!is.na(value)) %>%
  select(patid,obsdate,value)

hist(vitd_tests$value %>% as.numeric)

vd_test_validation = patient %>%
  left_join(vitd_tests,by="patid") %>% 
  filter(!is.na(value)) %>%
  filter(!is.na(vd_def_dx_date)) %>%
  filter(!is.na(obsdate)) %>%
  mutate(vd_test_date = as.Date(obsdate,format="%d/%m/%Y")) %>%
  mutate(vd_test_to_def_dx = abs(as.numeric((vd_test_date - vd_def_dx_date )/365 ))) %>% 
  filter(vd_test_to_def_dx<1/12)

ggplot(
  vd_test_validation %>% filter(vd_def_status=="vd_def"),
  aes(as.numeric(value)))+
  geom_density(fill="lightblue",alpha=0.5)+
  geom_vline(xintercept = 50)+
  scale_x_log10()+
  theme_minimal()

vd_test_validation %>% 
  group_by(vd_def_status) %>%
  mutate(value = as.numeric(value)) %>% 
  summarise_at(.vars = vars(value),.funs = c(median,IQR))

#################################
# 5. Alcohol 
#################################

# read in codelist 
alcohol_codes = read_csv("../codelists/alcohol_codelist.csv",
                         col_types = "ccccccccc")

# restrict to just 'disease' codes 
alcohol_binary_codes = alcohol_codes %>% filter(!EmisCodeCategoryId %in% c(3,37))
alcohol_score_codes = alcohol_codes %>% filter(EmisCodeCategoryId %in% c(3,37))

# get most recent alcohol  data >1yr before index
alcohol_df = clinical %>%
  filter(medcodeid %in% alcohol_binary_codes$MedCodeId) %>%
  mutate(alc_date = as.Date(obsdate,format="%d/%m/%Y")) %>%
  filter(alc_date < as.Date("01/05/2022",format="%d/%m/%Y")) %>%
  left_join(patient %>% select(patid,rough_dob,index_age),by="patid") %>%
  mutate(age_at_alc = as.numeric(alc_date - rough_dob)/365) %>%
  filter(age_at_alc>0) %>%
  group_by(patid) %>%
  filter(index_age - age_at_alc >= 1) %>%
  slice_max(age_at_alc,n=1,with_ties = F) %>% 
  mutate(value = as.numeric(value)) %>%
  mutate(alcohol_status = ifelse(
    (Term == "Alcohol consumption" & value == 0) |
      (Term == "Alcohol units consumed per week" & value == 0) |
      (Term == "Alcohol consumption NOS" & value == 0) |
      (Term == "Alcohol units consumed per day" & value == 0),
    "never_drinker",
    "ever_drinker"
  )) %>%
  select(patid,alc_date,age_at_alc,alcohol_status)


# join with main dataset
patient = patient %>%
  left_join(alcohol_df,by="patid")

table(patient$alcohol_status)

# validate using scores 
audit_c_scores = clinical %>% 
  filter(medcodeid == "408548014") %>%
  filter(!is.na(value)) %>%
  select(patid,obsdate,value)

# score must be 0-12
audit_c_scores = audit_c_scores %>%
  mutate(value = as.numeric(value)) %>%
  filter(value >=0 & value <= 12)

hist(audit_c_scores$value %>% as.numeric)

alcohol_validation = patient %>%
  left_join(audit_c_scores,by="patid") %>% 
  filter(!is.na(value)) %>%
  filter(!is.na(obsdate)) %>%
  mutate(auditc_date = as.Date(obsdate,format="%d/%m/%Y")) %>%
  mutate(alc_audit_delta = abs(as.numeric((alc_date - auditc_date )/365 ))) %>% 
  filter(alc_audit_delta<1/12) %>%
  filter( !(alcohol_status == "never_drinker" & age_at_alc > (index_age-1) ) ) %>%
  filter(!is.na(alcohol_status))

ggplot(alcohol_validation,
       aes(alcohol_status,value))+
  geom_boxplot(fill="lightblue",alpha=0.5)+
  theme_minimal()


alcohol_validation %>% 
  group_by(alcohol_status) %>%
  mutate(value = as.numeric(value)) %>% 
  summarise_at(.vars = vars(value),.funs = c(mean,median,IQR,range))

table(alcohol_validation$alcohol_status,!is.na(alcohol_validation$value))

#####################################
# 6. Head injury
#####################################

# read in codelist 
hi_codes = read_csv("../codelists/headinjury_codelist.csv",
                    col_types = "ccccccccc")

# head injury data 
hi_df = clinical %>%
  filter(medcodeid %in% hi_codes$MedCodeId) %>%
  mutate(hi_dx_date = as.Date(obsdate,format="%d/%m/%Y")) %>%
  filter(hi_dx_date < as.Date("01/05/2022",format="%d/%m/%Y")) %>%
  left_join(patient %>% select(patid,rough_dob),by="patid") %>%
  mutate(age_at_hi = as.numeric(hi_dx_date - rough_dob)/365) %>%
  filter(age_at_hi>0) %>%
  group_by(patid) %>%
  slice_min(age_at_hi,n=1,with_ties = F) %>%
  select(patid,hi_dx_date,age_at_hi)

# join with main dataset
patient = patient %>%
  left_join(hi_df,by="patid") %>%
  mutate(head_injury_status = ifelse(!is.na(age_at_hi) & 
                                       age_at_hi < (index_age-1),
                                     "head_injury",
                                     "no_head_injury"))

patient %>%
  count(head_injury_status) %>%
  mutate(prop = n/sum(n))

patient %>% group_by(MS_status) %>%
  count(head_injury_status) %>%
  mutate(prop = n/sum(n))

summary(patient$age_at_hi)
hist(patient$age_at_hi)


#####################################
# 7. Solvents
#####################################

# read in codelist 
solvent_codes = read_csv("../codelists/solvent_codelist.csv",
                         col_types = "ccccccccc")

# head injury data 
solvent_df = clinical %>%
  filter(medcodeid %in% solvent_codes$MedCodeId) %>%
  mutate(solvent_date = as.Date(obsdate,format="%d/%m/%Y")) %>%
  filter(solvent_date < as.Date("01/05/2022",format="%d/%m/%Y")) %>%
  left_join(patient %>% select(patid,rough_dob),by="patid") %>%
  mutate(age_at_solvent_exp = as.numeric(solvent_date - rough_dob)/365) %>%
  filter(age_at_solvent_exp>0) %>%
  group_by(patid) %>%
  slice_min(age_at_solvent_exp,n=1,with_ties = F) %>%
  select(patid,solvent_date,age_at_solvent_exp)

# join with main dataset
patient = patient %>%
  left_join(solvent_df,by="patid") %>%
  mutate(solvent_exposure_status = ifelse(!is.na(age_at_solvent_exp) & 
                                            age_at_solvent_exp < (index_age-1),
                                          "solvents",
                                          "no_solvents"))

patient %>%
  count(solvent_exposure_status) %>%
  mutate(prop = n/sum(n))

patient %>% group_by(MS_status) %>%
  count(solvent_exposure_status) %>%
  mutate(prop = n/sum(n))

summary(patient$age_at_hi)
hist(patient$age_at_hi)

#####################################
# save progress & create sub-cohorts
#####################################

saveRDS(patient,"cprd_aurum_cleaned_exposures_patient.rds")
saveRDS(clinical,"cprd_aurum_cleaned_exposures_clinical.rds")

#################################
# case-control study
#################################

patient %>% glimpse

# recode MS status with control as ref 
patient$MS_status = relevel(factor(patient$MS_status),ref="Control")


# create model function 
make_case_control_model = function(
  data_for_model = patient,
  var_to_test = "bmi_category"){
  
  # print parameters
  message("variable: ",var_to_test)
  
  # subset to non-missing obs 
  model_data = data_for_model %>% filter(!is.na(.data[[var_to_test]]))
  
  # if var is ordered, undorder is 
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
    MS_status ~ age_at_registration + gender + model_data[[var_to_test]],
    family=binomial(link="logit")) %>% 
    extract_coefs_from_model_obj() %>%
    mutate(model = "age + sex")
  
  agesex_ethnic_model  = glm(
    data = model_data,
    MS_status ~ age_at_registration + gender +
      ethnicity_pc_hes_hybrid + model_data[[var_to_test]],
    family=binomial(link="logit")) %>% 
    extract_coefs_from_model_obj() %>%
    mutate(model = "age + sex + ethnicity")
  
  agesex_imd_model  = glm(
    data = model_data,
    MS_status ~ age_at_registration + gender +
      e2019_imd_5 + model_data[[var_to_test]],
    family=binomial(link="logit")) %>% 
    extract_coefs_from_model_obj()%>%
    mutate(model = "age + sex + IMD")
  
  covar_names = c("(Intercept)",
                  "age_at_registration",
                  "genderMale",
                  "ethnicity_pc_hes_hybridBlack",
                  "ethnicity_pc_hes_hybridMixed/Other",
                  "ethnicity_pc_hes_hybridWhite",
                  "e2019_imd_52",
                  "e2019_imd_53",
                  "e2019_imd_54",
                  "e2019_imd_55") 
  
  overall_model_df = bind_rows(
    uni_model,
    agesex_model,
    agesex_ethnic_model,
    agesex_imd_model) %>%
    filter(!var %in% covar_names) 
  overall_model_df$var = str_remove(overall_model_df$var,pattern="model_data\\[\\[var_to_test\\]\\]")
  overall_model_df$variable = var_to_test
  return(overall_model_df)
}

# relevel variables to orient towards risk factor
patient$vd_def_status = relevel(factor(patient$vd_def_status),ref="vd_replete")
patient$im_status = relevel(factor(patient$im_status),ref="no_IM")
patient$head_injury_status = relevel(factor(patient$head_injury_status),ref="no_head_injury")
patient$bmi_category = relevel(factor(patient$bmi_category,ordered = F),ref="Healthy")
patient$alcohol_status = relevel(factor(patient$alcohol_status),ref="ever_drinker")

# run for different variables 
bmi_cat = make_case_control_model(var_to_test = "bmi_category")
bmi_cont = make_case_control_model(var_to_test = "bmi_z")
smoking = make_case_control_model(var_to_test = "smoking_status")
im = make_case_control_model(var_to_test = "im_status")
vd_def = make_case_control_model(var_to_test = "vd_def_status")
alcohol = make_case_control_model(var_to_test = "alcohol_status")
hi = make_case_control_model(var_to_test = "head_injury_status")

# combine results and plot
model_overall_df = bind_rows(bmi_cat,smoking,im,vd_def,alcohol,hi)
model_overall_df$y = paste0(model_overall_df$variable,":",model_overall_df$var)

model_overall_df = model_overall_df %>%
  mutate(bonf_sig = ifelse(pval < 0.05/6,
                           "*",
                           "NS"))

model_overall_df$y = factor(model_overall_df$y)
levels(model_overall_df$y) = c("Never drinker",
                               "BMI > 40",
                               "BMI 30-40",
                               "BMI 25-30",
                               "BMI < 18.5",
                               "Prior head injury",
                               "Prior IM",
                               "Ever smoker",
                               "Vitamin D deficiency")
p=ggplot(model_overall_df,
         aes(or,y,col=model,label=bonf_sig))+
  facet_wrap(~model)+
  geom_point()+
  geom_errorbarh(mapping=aes(xmin = lower_ci,xmax=upper_ci,y=y),height=0.1)+
  theme_minimal()+
  geom_vline(xintercept=1,alpha=0.2)+
  labs(x="Odds Ratio (for MS)",y="Exposure")+
  scale_color_brewer(palette="Paired")+
  geom_text(mapping = aes(upper_ci,y,hjust=-1))+
  scale_x_continuous(limits=c(0.5,5))

png("./figs/case_control_plots.png",res=300,units="in",width=6,height=6)
p
dev.off()

# repeat with just main analysis (age + sex)
p=ggplot(model_overall_df %>% filter(model=="age + sex"),
         aes(or,y,col=variable,label=bonf_sig))+
  geom_point()+
  geom_errorbarh(mapping=aes(xmin = lower_ci,xmax=upper_ci,y=y),height=0.1)+
  theme_minimal()+
  geom_vline(xintercept=1,alpha=0.2)+
  labs(x="Odds Ratio (for MS)",y="Exposure")+
  scale_color_brewer(palette="Paired")+
  geom_text(mapping = aes(upper_ci,y,hjust=-1))+
  scale_x_continuous(limits=c(0.5,3.5))+
  theme(legend.position="none")

png("./figs/case_control_plots_age_sex.png",res=300,units="in",width=4,height=4)
p
dev.off()

# save results to a table
tbl = model_overall_df %>% 
  select(variable, var,eff_n,or,lower_ci,upper_ci,model, pval, bonf_sig)%>% 
  mutate(
    or_ci = paste0(
      round(or,2)," (",round(lower_ci,2)," - ",round(upper_ci,2),")"
    )
  ) %>% 
  select(-or,-lower_ci,-upper_ci) %>%
  select(1,2,3,7,4,5,6)

tbl = tbl %>% pivot_wider(id_cols = c(variable,var), names_from = model, values_from = c(eff_n,or_ci,pval, bonf_sig))
write_csv(tbl,"./tables/case_control_results.csv")


# repeat for sensitivity analysis 
do_sensitivity_analysis = function(dataset_for_cc_analysis_name = "post_2007_cohort"){
  
  # evaluate name of arg to get df 
  dataset_for_cc_analysis = get(dataset_for_cc_analysis_name)
  
  # recode MS status with control as ref 
  dataset_for_cc_analysis$MS_status = relevel(factor(dataset_for_cc_analysis$MS_status),ref="Control")
  
  # relevel variables to orient towards risk factor
  dataset_for_cc_analysis$vd_def_status = relevel(factor(dataset_for_cc_analysis$vd_def_status),ref="vd_replete")
  dataset_for_cc_analysis$im_status = relevel(factor(dataset_for_cc_analysis$im_status),ref="no_IM")
  dataset_for_cc_analysis$head_injury_status = relevel(factor(dataset_for_cc_analysis$head_injury_status),ref="no_head_injury")
  dataset_for_cc_analysis$bmi_category = relevel(factor(dataset_for_cc_analysis$bmi_category,ordered = F),ref="Healthy")
  dataset_for_cc_analysis$alcohol_status = relevel(factor(dataset_for_cc_analysis$alcohol_status),ref="ever_drinker")
  
  # run for different variables 
  bmi_cat = make_case_control_model(var_to_test = "bmi_category", data_for_model = dataset_for_cc_analysis)
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
  levels(model_overall_df$y) = c("Never drinker",
                                 "BMI > 40",
                                 "BMI 30-40",
                                 "BMI 25-30",
                                 "BMI < 18.5",
                                 "Prior head injury",
                                 "Prior IM",
                                 "Ever smoker",
                                 "Vitamin D deficiency")
  p=ggplot(model_overall_df,
           aes(or,y,col=model,label=bonf_sig))+
    facet_wrap(~model)+
    geom_point()+
    geom_errorbarh(mapping=aes(xmin = lower_ci,xmax=upper_ci,y=y),height=0.1)+
    theme_minimal()+
    geom_vline(xintercept=1,alpha=0.2)+
    labs(x="Odds Ratio (for MS)",y="Exposure")+
    scale_color_brewer(palette="Paired")+
    geom_text(mapping = aes(upper_ci,y,hjust=-1))
  png(paste0("./figs/case_control_plots_sensitivity_analysis",dataset_for_cc_analysis_name,".png"),res=300,units="in",width=6,height=6)
  print(p)
  dev.off()
  
  # save results to a table
  tbl = model_overall_df %>% 
    select(variable, var,eff_n,or,lower_ci,upper_ci,model, pval, bonf_sig)%>% 
    mutate(
      or_ci = paste0(
        round(or,2)," (",round(lower_ci,2)," - ",round(upper_ci,2),")"
      )
    ) %>% 
    select(-or,-lower_ci,-upper_ci) %>%
    select(1,2,3,7,4,5,6)
  
  tbl = tbl %>% pivot_wider(id_cols = c(variable,var), names_from = model, values_from = c(eff_n,or_ci,pval, bonf_sig))
  write_csv(tbl,paste0("./tables/case_control_plots_sensitivity_analysis",dataset_for_cc_analysis_name,".csv"))
}
cohorts_to_analyse = c("post_2007_cohort","post_1990_cohort","antecedent_5_cohort","ms_sensitivity_cohort")
sapply(cohorts_to_analyse,do_sensitivity_analysis)

# look at effect of BMI during different time windows on later MS risk
bmi_lower = seq(16,40,by=8)
bmi_upper = seq(24,48,by=8)
bmi_cat_levels = paste0("bmi[ ",
                        bmi_lower,
                        ", ",
                        bmi_upper,
                        ")")
overall_res = list()
for(i in c(1:length(bmi_cat_levels))){
  bmi_level = bmi_cat_levels[i]
  bmi_ceiling = bmi_upper[i]
  data_for_model = patient %>% filter(index_age>bmi_ceiling)
  model = glm(data = data_for_model,
              factor(MS_status)=="Case" ~ age_at_registration + gender + z_score(data_for_model[[bmi_level]]),
              family=binomial(link="logit")) %>%
    extract_coefs_from_model_obj() %>%
    mutate(bmi_window = paste0("Age < ",bmi_ceiling)) %>%
    dplyr::select(var, or,lower_ci,upper_ci,pval,bmi_window,eff_n) %>%
    filter(var == "z_score(data_for_model[[bmi_level]])")
  overall_res[[i]] = model   
}
overall_res = do.call("bind_rows",overall_res)
ggplot(overall_res,aes(bmi_window,or))+
  geom_point()+
  theme_minimal()+
  geom_errorbar(mapping = aes(x = bmi_window,ymin = lower_ci,ymax = upper_ci),width=0.3)+
  geom_hline(yintercept = 1,alpha=0.2)

#################################
# colliding & ethnicity 
#################################

# define some useful new / transformed predictors for missingness
patient$numeric_imd = as.numeric(patient$e2019_imd_5)
patient$registered_pre_2007 = ifelse(patient$numeric_year_of_reg<2007,"pre_2007","post_2007")
patient$registered_pre_2007 = relevel(factor(patient$registered_pre_2007),ref="pre_2007")
patient$yob_sq = patient$yob^2
patient$fu_time_sqrt = sqrt(patient$fu_time)
patient$fu_time_sqrt %>% hist

output_df %>%
  group_by(MS_status) %>%
  count(is_na) %>%
  mutate(prop = n/sum(n))

res_df %>%
  group_by(MS_status) %>%
  count(is_na) %>%
  mutate(prop = n/sum(n))

check_missingness_predictors(input_data = patient,
                             missing_col = "ethnicity_pc_hes_hybrid",
                             miss_pred_cols = c("registered_pre_2007",
                                                "gender",
                                                "fu_time",
                                                "fu_time_sqrt",
                                                "yob",
                                                "yob_sq",
                                                "numeric_imd",
                                                "e2011_urban_rural",
                                                "MS_status"))

check_missingness_predictors(input_data = output_df,
                             missing_col = "ethnicity_pc_hes_hybrid",
                             miss_pred_cols = c("registered_pre_2007",
                                                "gender",
                                                "fu_time",
                                                "fu_time_sqrt",
                                                "yob",
                                                "yob_sq",
                                                "numeric_imd",
                                                "e2011_urban_rural",
                                                "MS_status"))

check_missingness_predictors(input_data = post_2007_cohort,
                             missing_col = "ethnicity_pc_hes_hybrid",
                             miss_pred_cols = c("gender",
                                                "fu_time",
                                                "yob",
                                                "e2011_urban_rural",
                                                "MS_status"))

ipw_res_step = get_weights_stepwise(input_data = patient,
                                    missing_col = "ethnicity_pc_hes_hybrid",
                                    miss_pred_cols = c("registered_pre_2007",
                                                       "gender",
                                                       "fu_time_sqrt",
                                                       "yob"))

plot_from_model_output(ipw_res_step$missingness_model)
plot_missingness_prob(ipw_res_step$ipw_results)
hinkley_test(ipw_res_step$missingness_model)


weighted_reg_res = ipw_regression(input_data = ipw_res_step$ipw_results,
                                  outcome_col = "MS_status",
                                  outcome_pred_cols = c("ethnicity_pc_hes_hybrid"))




weighted_reg_res$weighted_reg

z_score = function(x){
  (x - mean(x,na.rm=T)) / sd(x, na.rm = T)
}


make_hist("numeric_year_of_reg")
make_hist("numeric_year_of_reg_sq")
make_hist("numeric_year_of_reg_log10_z")

#################################
# timing of risk factors
#################################


###################################
# ethnic variation in risk factors
###################################

# create model function 
make_case_control_model_by_ethnicity = function(
  data_for_model = patient,
  var_to_test = "bmi_category"){
  
  # print parameters
  message("variable: ",var_to_test)
  
  # subset to non-missing obs 
  model_data = data_for_model %>% filter(!is.na(.data[[var_to_test]]))
  
  # filter out non-missing ethnicity
  model_data = data_for_model %>% filter(!is.na(ethnicity_pc_hes_hybrid))
  
  # if var is ordered, undorder is 
  model_data[[var_to_test]] = if(is.factor(model_data[[var_to_test]])){
    factor(model_data[[var_to_test]],ordered=F)
  } else{
    model_data[[var_to_test]]
  }
  
  ethnicity_interaction_model  = glm(
    data = model_data,
    MS_status ~ age_at_registration + gender +
      ethnicity_pc_hes_hybrid * model_data[[var_to_test]],
    family=binomial(link="logit")) %>% 
    extract_coefs_from_model_obj() %>%
    mutate(model = "ethnicity_interaction")
  
  ethnicity_interaction_model_imd  = glm(
    data = model_data,
    MS_status ~ age_at_registration + gender + e2019_imd_5 +
      ethnicity_pc_hes_hybrid * model_data[[var_to_test]],
    family=binomial(link="logit")) %>% 
    extract_coefs_from_model_obj() %>%
    mutate(model = "ethnicity_imd_interaction")
  
  white_data = model_data %>% filter(ethnicity_pc_hes_hybrid=="White")
  white_model  = glm(
    data = white_data,
    MS_status ~ age_at_registration + gender + white_data[[var_to_test]],
    family=binomial(link="logit")) %>% 
    extract_coefs_from_model_obj() %>%
    mutate(model = "white")
  
  black_data = model_data %>% filter(ethnicity_pc_hes_hybrid=="Black")
  black_model  = glm(
    data = model_data %>% filter(ethnicity_pc_hes_hybrid=="Black"),
    MS_status ~ age_at_registration + gender + black_data[[var_to_test]],
    family=binomial(link="logit")) %>% 
    extract_coefs_from_model_obj() %>%
    mutate(model = "black")
  
  asian_data = model_data %>% filter(ethnicity_pc_hes_hybrid=="Asian")
  asian_model  = glm(
    data = model_data %>% filter(ethnicity_pc_hes_hybrid=="Asian"),
    MS_status ~ age_at_registration + gender + asian_data[[var_to_test]],
    family=binomial(link="logit")) %>% 
    extract_coefs_from_model_obj() %>%
    mutate(model = "asian")
  
  white_model_imd  = glm(
    data = white_data,
    MS_status ~ age_at_registration + gender + e2019_imd_5 + white_data[[var_to_test]],
    family=binomial(link="logit")) %>% 
    extract_coefs_from_model_obj() %>%
    mutate(model = "white_imd")
  
  black_model_imd  = glm(
    data = model_data %>% filter(ethnicity_pc_hes_hybrid=="Black"),
    MS_status ~ age_at_registration + gender + e2019_imd_5 + black_data[[var_to_test]],
    family=binomial(link="logit")) %>% 
    extract_coefs_from_model_obj() %>%
    mutate(model = "black_imd")
  
  asian_model_imd  = glm(
    data = model_data %>% filter(ethnicity_pc_hes_hybrid=="Asian"),
    MS_status ~ age_at_registration + gender + e2019_imd_5 +asian_data[[var_to_test]],
    family=binomial(link="logit")) %>% 
    extract_coefs_from_model_obj() %>%
    mutate(model = "asian_imd")
  
  
  overall_model_df = bind_rows(
    ethnicity_interaction_model,
    ethnicity_interaction_model_imd,
    white_model,
    asian_model,
    black_model,
    white_model_imd,
    asian_model_imd,
    black_model_imd) 
  overall_model_df$var = str_remove(overall_model_df$var,pattern="model_data\\[\\[var_to_test\\]\\]")
  overall_model_df$variable = var_to_test
  return(overall_model_df)
}

# run for different variables 
bmi_cat = make_case_control_model_by_ethnicity(var_to_test = "bmi_category")
bmi_cont = make_case_control_model_by_ethnicity(var_to_test = "bmi_z")
smoking = make_case_control_model_by_ethnicity(var_to_test = "smoking_status")
im = make_case_control_model_by_ethnicity(var_to_test = "im_status")
vd_def = make_case_control_model_by_ethnicity(var_to_test = "vd_def_status")
alcohol = make_case_control_model_by_ethnicity(var_to_test = "alcohol_status")
hi = make_case_control_model_by_ethnicity(var_to_test = "head_injury_status")

# combine results and plot
model_overall_df = bind_rows(bmi_cat,smoking,im,vd_def,alcohol,hi)
model_overall_df$model

# look at interaction models 
interaction_model_df = model_overall_df %>% 
  filter(grepl("interaction",model)) %>%
  filter(grepl(":",var))
interaction_model_df %>% filter(pval<0.005)

# make simple forest comparing ethnicities
model_overall_df = model_overall_df %>% 
  filter(model %in% c("white_imd","asian_imd","black_imd")) %>%
  filter(!var %in% c("(Intercept)","age_at_registration","genderMale")) %>%
  filter(!grepl("e2019_imd",var)) %>%
  tidyr::separate(var, sep = "\\]]",into=c("eth","variable_level")) %>%
  select(-eth) %>%
  mutate(bonf_sig = ifelse(pval < 0.05/6,
                           "*",
                           "NS")) %>%
  mutate(y = paste0(variable,":",variable_level)) %>%
  mutate(y = factor(y))

levels(model_overall_df$y) = c("Never drinker",
                               "BMI > 40",
                               "BMI 30-40",
                               "BMI 25-30",
                               "BMI < 18.5",
                               "Prior head injury",
                               "Prior IM",
                               "Ever smoker",
                               "Vitamin D deficiency")
p=ggplot(model_overall_df,
         aes(or,y,col=model,label=bonf_sig))+
  geom_point(position=ggstance::position_dodgev(height=1))+
  geom_errorbarh(mapping=aes(xmin = lower_ci,xmax=upper_ci,y=y),height=0.1,position=ggstance::position_dodgev(height=1))+
  theme_minimal()+
  facet_wrap(~model)+
  geom_vline(xintercept=1,alpha=0.2)+
  labs(x="Odds Ratio (for MS)",y="Exposure")+
  scale_color_brewer(palette="Set2")+
  scale_x_log10(limits=c(0.1,20))+
  geom_text(mapping = aes(upper_ci,y,hjust=-1))+
  theme(legend.position="none")

png("./figs/exposure_het_ethnicity.png",res=300,units="in",width=6,height=6)
p
dev.off()

# meta-analysis 

pvals = c()
for(exposure in levels(model_overall_df$y)){
  te = model_overall_df[model_overall_df$y == exposure,][["beta"]]
  labs = model_overall_df[model_overall_df$y == exposure,][["model"]]
  se = model_overall_df[model_overall_df$y == exposure,][["se"]]
  meta_res = meta::metagen(TE = te, seTE = se, studlab = labs)
  pvals <<- c(pvals,meta_res$pval.Q)
}
data.frame(exposure = levels(model_overall_df$y), p = pvals) %>%
  mutate(p_adj = p.adjust(p,method="bonf")) %>%
  mutate(bonf_sig = ifelse(p_adj<0.05,"*","NS"))

# repeat with just main analysis (age + sex)
p=ggplot(model_overall_df %>% filter(model=="age + sex"),
         aes(or,y,col=variable,label=bonf_sig))+
  geom_point()+
  geom_errorbarh(mapping=aes(xmin = lower_ci,xmax=upper_ci,y=y),height=0.1)+
  theme_minimal()+
  geom_vline(xintercept=1,alpha=0.2)+
  labs(x="Odds Ratio (for MS)",y="Exposure")+
  scale_color_brewer(palette="Paired")+
  geom_text(mapping = aes(upper_ci,y,hjust=-1))+
  scale_x_continuous(limits=c(0.5,3.5))+
  theme(legend.position="none")

png("./figs/case_control_plots_age_sex.png",res=300,units="in",width=4,height=4)
p
dev.off()

# save results to a table
tbl = model_overall_df %>% 
  select(variable, var,eff_n,or,lower_ci,upper_ci,model, pval, bonf_sig)%>% 
  mutate(
    or_ci = paste0(
      round(or,2)," (",round(lower_ci,2)," - ",round(upper_ci,2),")"
    )
  ) %>% 
  select(-or,-lower_ci,-upper_ci) %>%
  select(1,2,3,7,4,5,6)

tbl = tbl %>% pivot_wider(id_cols = c(variable,var), names_from = model, values_from = c(eff_n,or_ci,pval, bonf_sig))
write_csv(tbl,"./tables/case_control_results.csv")


# repeat for sensitivity analysis 
do_sensitivity_analysis = function(dataset_for_cc_analysis_name = "post_2007_cohort"){
  
  # evaluate name of arg to get df 
  dataset_for_cc_analysis = get(dataset_for_cc_analysis_name)
  
  # recode MS status with control as ref 
  dataset_for_cc_analysis$MS_status = relevel(factor(dataset_for_cc_analysis$MS_status),ref="Control")
  
  # relevel variables to orient towards risk factor
  dataset_for_cc_analysis$vd_def_status = relevel(factor(dataset_for_cc_analysis$vd_def_status),ref="vd_replete")
  dataset_for_cc_analysis$im_status = relevel(factor(dataset_for_cc_analysis$im_status),ref="no_IM")
  dataset_for_cc_analysis$head_injury_status = relevel(factor(dataset_for_cc_analysis$head_injury_status),ref="no_head_injury")
  dataset_for_cc_analysis$bmi_category = relevel(factor(dataset_for_cc_analysis$bmi_category,ordered = F),ref="Healthy")
  dataset_for_cc_analysis$alcohol_status = relevel(factor(dataset_for_cc_analysis$alcohol_status),ref="ever_drinker")
  
  # run for different variables 
  bmi_cat = make_case_control_model(var_to_test = "bmi_category", data_for_model = dataset_for_cc_analysis)
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
  levels(model_overall_df$y) = c("Never drinker",
                                 "BMI > 40",
                                 "BMI 30-40",
                                 "BMI 25-30",
                                 "BMI < 18.5",
                                 "Prior head injury",
                                 "Prior IM",
                                 "Ever smoker",
                                 "Vitamin D deficiency")
  p=ggplot(model_overall_df,
           aes(or,y,col=model,label=bonf_sig))+
    facet_wrap(~model)+
    geom_point()+
    geom_errorbarh(mapping=aes(xmin = lower_ci,xmax=upper_ci,y=y),height=0.1)+
    theme_minimal()+
    geom_vline(xintercept=1,alpha=0.2)+
    labs(x="Odds Ratio (for MS)",y="Exposure")+
    scale_color_brewer(palette="Paired")+
    geom_text(mapping = aes(upper_ci,y,hjust=-1))
  png(paste0("./figs/case_control_plots_sensitivity_analysis",dataset_for_cc_analysis_name,".png"),res=300,units="in",width=6,height=6)
  print(p)
  dev.off()
  
  # save results to a table
  tbl = model_overall_df %>% 
    select(variable, var,eff_n,or,lower_ci,upper_ci,model, pval, bonf_sig)%>% 
    mutate(
      or_ci = paste0(
        round(or,2)," (",round(lower_ci,2)," - ",round(upper_ci,2),")"
      )
    ) %>% 
    select(-or,-lower_ci,-upper_ci) %>%
    select(1,2,3,7,4,5,6)
  
  tbl = tbl %>% pivot_wider(id_cols = c(variable,var), names_from = model, values_from = c(eff_n,or_ci,pval, bonf_sig))
  write_csv(tbl,paste0("./tables/case_control_plots_sensitivity_analysis",dataset_for_cc_analysis_name,".csv"))
}
cohorts_to_analyse = c("post_2007_cohort","post_1990_cohort","antecedent_5_cohort","ms_sensitivity_cohort")
sapply(cohorts_to_analyse,do_sensitivity_analysis)



#################################
# ethnic variation in basic demog
#################################

# filter to just ms
ms = patient %>% filter(MS_status=="Case")

# age of onset
p=ggplot(ms %>%
           filter(!is.na(ethnicity_pc_hes_hybrid)),
         aes(ethnicity_pc_hes_hybrid, 
             age_at_dx,
             fill=ethnicity_pc_hes_hybrid))+
  geom_violin(alpha=0.8)+
  stat_boxplot(alpha=0.2,width=0.1)+
  theme_minimal()+
  scale_fill_brewer(palette="Set3")+
  labs(x="Ethnic background",y="Age at MS diagnosis")
p+facet_wrap(~source_of_ethnicity_report)
p+facet_wrap(~gender)
p+facet_wrap(~e2019_imd_5)

# mortality 

table(ms$dead, ms$ethnicity_pc_hes_hybrid)
ms = ms %>% filter(years_data_after_index>1) %>%
  filter(regstartdate>"1990-01-01")
ms = ms %>%
  mutate(
    survival_time = ifelse(!is.na(age_at_death),
                           age_at_death - age_at_dx,
                           age_at_deregistration  - age_at_dx)) %>%
  mutate(status = ifelse(dead=="dead",
                         2,1))

ms$ethnicity_pc_hes_hybrid = relevel(factor(ms$ethnicity_pc_hes_hybrid),ref="White")
fit <- survfit(Surv(survival_time, status) ~ ethnicity_pc_hes_hybrid, data = ms)
cox_fit_raw <- coxph(Surv(survival_time, status) ~ ethnicity_pc_hes_hybrid, data = ms)
cox_fit_full <- coxph(Surv(survival_time, status) ~ gender + e2019_imd_5 + ethnicity_pc_hes_hybrid, data = ms)
cox_fit_imd <- coxph(Surv(survival_time, status) ~ gender + e2019_imd_5, data = ms)

ggsurvplot(fit, 
           data = ms, 
           censor.shape="|", 
           censor.size = 2,
           risk.table = T,
           legend.labs = c("Asian","Black","Mixed/Other","White"),
           fun="cumhaz",
           conf.int = F)

# IMD cumhaz 
fit <- survfit(Surv(survival_time, status) ~ e2019_imd_5, data = ms)

patient_plot = patient %>%
  filter(regstartdate>"1990-01-01") %>%
  filter(years_data_after_index>1) %>%
  filter(index_code_date > "2000-01-01" | ms_code_date > "2000-01-01") %>%
  mutate(
    survival_time = ifelse(!is.na(age_at_death),
                           age_at_death - index_age,
                           age_at_deregistration  - index_age)) %>%
  mutate(status = ifelse(dead=="dead",2,1))

Surv(patient_plot$survival_time, patient_plot$status)
fit <- survfit(Surv(survival_time, status) ~ factor(MS_status), data = patient_plot )

patient_plot$ethnicity_pc_hes_hybrid = relevel(factor(patient_plot$ethnicity_pc_hes_hybrid),ref="White")
cox_ethnic_nodep <- coxph(Surv(survival_time, status) ~ age_at_registration + gender + ethnicity_pc_hes_hybrid, data = patient_plot%>% filter(MS_status=="Case"))
cox_ethnic_dep <- coxph(Surv(survival_time, status) ~ age_at_registration + gender + ethnicity_pc_hes_hybrid + e2019_imd_5, data = patient_plot%>% filter(MS_status=="Case"))

p2=ggsurvplot(fit, 
              data = patient_plot, 
              censor.shape="|", 
              censor.size = 2,
              fun="cumhaz",
              legend.labs = c("Cont","MS"),
              conf.int = F)+ggtitle("post-1990")

p2
gridExtra::grid.arrange(p1,p2)
