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
# read in data
###########################

# set master index age filter - this will be used to define exposures 
master_index_age_window = 5

# patient data 
patient = read_tsv("../../Data/Primary_care/Aurum/21_000677_Extract_Patient_001.txt",
                   col_types="cccccccccccc") %>%
  dplyr::select(1,2,4,5,8,10,12)

# read in codes which aren't of interest
codes_to_exclude = read_csv("../codelists/codes_to_exclude_aurum.csv", col_types = "ccccc")


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
cols_to_plot = list("age_at_registration","numeric_year_of_reg")
plots = lapply(cols_to_plot,make_hist)

# save plots
plot_fx(do.call("grid.arrange",plots),"./figs/sf1.png")

#################################
# define MS outcome             #
#################################

# read in MS codelist
codelist = read_csv("/data/Wolfson-UKBB-Dobson/CPRD/Lookups/Aurum/ms_aurum_codelist_1608.csv",col_types = "cccccccc")

# separate into definite & probable/possible
definite_ms_codes = codelist %>% filter(Confidence=="Definite")
other_ms_codes = codelist %>% filter(Confidence!="Definite")

# remove cases with a non-definite MS code 
ms_code_obs = clinical %>%
  mutate(definite_ms_code = ifelse(medcodeid %in% definite_ms_codes$MedCodeId,"definite_MS",NA)) %>%
  mutate(non_definite_ms_code = ifelse(medcodeid %in% other_ms_codes$MedCodeId,"non_definite_MS",NA)) %>% 
  filter(!is.na(definite_ms_code) | !is.na(non_definite_ms_code))

definite_code_counts = ms_code_obs %>%
  group_by(patid) %>%
  dplyr::count(definite_ms_code) %>%
  ungroup %>%
  na.omit
nrow(definite_code_counts)

non_definite_code_counts = ms_code_obs %>%
  group_by(patid) %>%
  dplyr::count(non_definite_ms_code) %>%
  ungroup %>%
  na.omit

# 17405 with >=1 definite code 


# exclude 
patient = patient %>%
  mutate(MS_status = ifelse(
    patid %in% definite_code_counts$patid,"Case",
    ifelse(
      !(patid %in% non_definite_code_counts$patid) & !(patid %in% definite_code_counts$patid),
      "Control",NA
    ))) 

patient %>% 
  filter(patid %in% non_definite_code_counts$patid & !patid %in% definite_code_counts$patid) %>%
  dplyr::count(MS_status)
get_prop(patient,x = "MS_status")

# filter 
patient = patient %>%
  filter(!is.na(MS_status))

get_prop(patient,x = "MS_status")

# filter clinical table to included participants
clinical = clinical %>%
  filter(patid %in% patient$patid)

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
  left_join(ms_first_report,by="patid")

table(patient$MS_status)

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

missing_date_dx = patient %>% filter(MS_status=="Case" & is.na(ms_code_date))
clinical %>%  filter(medcodeid %in% codelist$MedCodeId) %>%
  filter(patid %in% missing_date_dx$patid)

## there are five people with MS codes but missing dates - exclude 
patient = patient %>% 
  filter(!(MS_status=="Case" & is.na(ms_code_date)))

## age at diagnosis
patient = patient %>%
  mutate(age_at_dx = delta_dates(ms_code_date - rough_dob))
patient$age_at_dx %>% summary
patient %>% group_by(gender) %>% summarise(median(age_at_dx,na.rm=T))

# check that only cases have a dx date & vice-versa
patient %>% 
  filter(is.na(age_at_pseudodx)) %>%
  count(MS_status)

patient %>% 
  filter(is.na(age_at_dx)) %>%
  count(MS_status)

# check that this derived index date matched supplied index date 
study_pop = read_table("../../Data/Primary_care/Aurum/CPRD_Aurum_study_pop_file_21_000677.txt", col_types = cols(.default = "c"))

patient %>% 
  filter(MS_status=="Case") %>%
  dplyr::select(patid,MS_status,ms_code_date) %>%
  left_join(study_pop,by="patid") %>%
  mutate(cprd_index = datify(indexdate)) %>%
  mutate(delta = delta_dates(cprd_index - ms_code_date)) %>%
  arrange(desc(delta))

# define variable which encompasses age at dx/pseudodx
# also define FU time prior to and after this date
patient = patient %>%
  mutate(index_age = ifelse(MS_status=="Case",age_at_dx,age_at_pseudodx)) %>%
  mutate(years_data_before_index = index_age - age_at_registration) %>%
  mutate(years_data_after_index = age_at_deregistration - index_age) 

# remove controls without 5y antecedent data 
hist(patient$years_data_before_index)
patient %>% filter(MS_status == "Case"  & years_data_before_index < 5) %>% dplyr::select(years_data_before_index)
summary(patient$years_data_before_index)
get_prop(patient %>% group_by(MS_status) %>% mutate(keep = ifelse(years_data_before_index>=4.99,"Y","N")),x="keep")
get_prop(patient %>% group_by(MS_status) %>% mutate(keep = ifelse(years_data_before_index>=9.99,"Y","N")),x="keep")

patient = patient %>%
  filter(years_data_before_index>=4.99)

# count 
get_prop(patient,x="MS_status")

# get rid of controls whose case has been removed - run without this 3/11/22
#matching_file = read_tsv("../../Data/Primary_care/Aurum/21_000677_CPRD_Aurum_matching_file.txt",col_types="cccccccc") %>% rename("patid" = control_patid)
#just_ms_cases = patient %>% filter(MS_status=="Case")
#matching_file = matching_file %>% filter(case_patid %in% just_ms_cases$patid)
#patient = patient %>% filter(patid %in% matching_file$patid | patid %in% just_ms_cases$patid)

# check that MS diagnosis after registration date
patient %>%
  mutate(ms_dx_pre_enrolment = ifelse(ms_code_date<regstartdate,"Yes","No")) %>%
  count(ms_dx_pre_enrolment)

patient %>% group_by(gender) %>% summarise(median(age_at_dx,na.rm=T))

p=make_grouped_hist(x="age_at_dx",group="gender")+
  scale_x_continuous(limits=c(0,110))+
  labs(x="Age at first MS diagnostic\ncode report")
plot_fx(p,"./figs/age_ms_by_sex.png", width=4,height=4)

# filter observation dataset
clinical =clinical %>% filter(patid %in% patient$patid)

## verify diagnoses (see who has >1 code)
ms_diagnoses_verification = clinical %>%
  filter(medcodeid %in% definite_ms_codes$MedCodeId) %>%
  count(patid) %>%
  rename("ms_diag_counts" = n)
patient = patient %>%
  left_join(ms_diagnoses_verification,by="patid")

patient %>% filter(MS_status=="Case") %>% count(ms_diag_counts>1) %>% mutate(n/sum(n)*100)

## repeat with other codes 
ms_diagnoses_verification_probable = clinical %>%
  filter(medcodeid %in% codelist$MedCodeId) %>%
  count(patid) %>%
  rename("ms_diag_counts_any_ms_code" = n)
patient = patient %>%
  left_join(ms_diagnoses_verification_probable,by="patid")

patient %>% filter(MS_status=="Case") %>% count(ms_diag_counts_any_ms_code>1) %>% mutate(n/sum(n)*100)

stats = patient %>% filter(MS_status=="Case") %>% summarise_at(vars(ms_diag_counts),.funs=c("mean","median","IQR","sd"))
p=ggplot(patient %>% filter(MS_status=="Case"),aes(ms_diag_counts))+geom_density(fill="lightblue",alpha=0.8)+scale_x_log10()+theme_minimal()+labs(x="MS diagnostic code counts \nper individual")
patient$n_ms_diag_codes = cut2(patient$ms_diag_counts,cuts=c(1,2,5,10))
levels(patient$n_ms_diag_codes) = c("1","2-5","6-10",">10")

p1=ggplot(patient %>% filter(MS_status=="Case"),aes(n_ms_diag_codes))+
  geom_bar(fill="lightblue",alpha=0.8,color="black")+
  theme_minimal()+
  labs(x="MS diagnostic code counts \nper individual")

plot_fx(p1,width=3,height=3,"./figs/ms_diag_counts.png")

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

get_prop(x = "ms_code_after_neuro_ref", dat = patient %>% filter(MS_status=="Case" & !is.na(ms_code_after_neuro_ref)))
get_prop(x = "MS_status", dat = patient)

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
get_prop(dat=patient %>% filter(numeric_year_of_reg>=1997),x="source_of_ethnicity_report")

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
# define other exposures
#################################


#################################
# 1. BMI
#################################

NumUnit = read.delim("/data/Wolfson-UKBB-Dobson/CPRD/Lookups/Aurum/NumUnit.txt")
clinical$unit    = ""
clinical$unit[clinical$numunitid != ""] = NumUnit$Description[match(clinical$numunitid[clinical$numunitid != ""], NumUnit$numunitid)]

# filter observation files 
weight_recording = clinical %>% filter(medcodeid %in% c("253677014","253688015","402435016","411922013","7879461000006116"))
height_recording = clinical %>% filter(medcodeid %in% c("253669010","11905141000006112","253676017","3315721000006119"))

weight_recording = weight_recording %>%
  mutate(value = as.numeric(value)) %>%
  filter(!is.na(value)) 

#convert stone to kg, if units in stone or if weight is between 10 to 20
weight_recording$value[which(weight_recording$unit== "decimal stones")]  = 
  weight_recording$value[which(weight_recording$unit== "decimal stones")]*6.35029318

weight_recording$value[which(weight_recording$value>=10 & weight_recording$value< 20)]=
  weight_recording$value[which(weight_recording$value>=10 & weight_recording$value< 20)]*6.35029318

weight_recording = weight_recording[weight_recording$value>20,]


# convert m to cm, if units in m and if height is between 1.21 to 2.14
height_recordings = height_recording %>%
  mutate(value = as.numeric(value)) %>%
  filter(!is.na(value)) 

height_recordings$value[height_recordings$unit %in%c("m","metres") & height_recordings$value<=2.14 & height_recordings$value>=1.21] =
  height_recordings$value[height_recordings$unit %in%c("m","metres") & height_recordings$value<=2.14 & height_recordings$value>=1.21]*100

height_recordings$value[height_recordings$value<=2.14 & height_recordings$value>=1.21] = 
  height_recordings$value[height_recordings$value<=2.14 & height_recordings$value>=1.21]*100

height_recordings$value[height_recordings$value<=21.4 & height_recordings$value>=12.1] = 
  height_recordings$value[height_recordings$value<=21.4 & height_recordings$value>=12.1]*10

height_recordings$value[height_recordings$value<=2140 & height_recordings$value>=1210] = 
  height_recordings$value[height_recordings$value<=2140 & height_recordings$value>=1210]/10

height_recordings$value[height_recordings$value<=21400 & height_recordings$value>=12100] = 
  height_recordings$value[height_recordings$value<=21400 & height_recordings$value>=12100]/100

height_recordings = height_recordings[height_recordings$value<=214 & height_recordings$value>=121,]


# join height and weight
height_weight = weight_recording %>% 
  bind_rows(height_recordings) %>%
  dplyr::select(patid,obsdate,Term,value) %>% 
  pivot_wider(id_cols = c(patid,obsdate),
              names_from = Term, 
              values_from = value,
              values_fn = mean) # takes average if repeated in 1 day

height_weight$weight = ifelse(!is.na(height_weight$`Body weight`),height_weight$`Body weight`,
                              ifelse(!is.na(height_weight$`O/E - weight NOS`), height_weight$`O/E - weight NOS`,
                                     ifelse(!is.na(height_weight$Weight), height_weight$Weight,
                                            ifelse(!is.na(height_weight$`O/E - weight 10-20% over ideal`),height_weight$`O/E - weight 10-20% over ideal`,
                                                   ifelse(!is.na(height_weight$`O/E - overweight`), height_weight$`O/E - overweight`, NA)))))

height_weight$height = ifelse(!is.na(height_weight$`Standing height`),height_weight$`Standing height`,
                              ifelse(!is.na(height_weight$`O/E - height`), height_weight$`O/E - height`,
                                     ifelse(!is.na(height_weight$`O/E - height NOS`), height_weight$`O/E - height NOS`,
                                            ifelse(!is.na(height_weight$`Body height`),height_weight$`Body height`, NA))))

height_weight = height_weight[order(height_weight$patid,height_weight$obsdate), c("patid","obsdate","height","weight")]
height_weight = as.data.frame(height_weight)

#------ last value carried forward for height, if no last value available use most recent value
height_weight <- data.frame(height_weight %>%
                              group_by(patid) %>%
                              fill(height, .direction = "downup") %>%
                              fill(height, .direction = "updown"))

head(height_weight)

# exclude NAs 
height_weight_nonmissing = height_weight %>% 
  filter(!is.na(`weight`) & !is.na(`height`))
height_weight_nonmissing %>%  nrow
height_weight_nonmissing %>% distinct(patid) %>% nrow

# find missing data
ht_missing = patient %>% filter(!patid %in% height_recordings$patid)
wt_missing = patient %>% filter(!patid %in% weight_recording$patid)
both_missing = patient %>% filter(!patid %in% weight_recording$patid & !patid %in% height_recordings$patid)
either_missing = patient %>% filter(!patid %in% weight_recording$patid | !patid %in% height_recordings$patid)

# for people with missing data, impute with nearest height 
height_weight_missing = height_weight %>% 
  filter(!patid %in% height_weight_nonmissing$patid) %>%
  filter(!is.na(`weight`) & is.na(`height`))

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
  dplyr::select(patid,obsdate.x,`weight`,value,delta_date) %>% 
  rename("obsdate" = obsdate.x,
         `height` = value)

# calculate BMI 
bmi_overall = height_weight_nonmissing %>% 
  bind_rows(height_weight_missing %>% dplyr::select(-delta_date)) %>%
  mutate(bmi = `weight`/((`height`/100)^2))
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
      distinct(patid,.keep_all=T) %>%
      dplyr::select(patid,obs_dt,bmi,age_at_bmi) %>%
      rename("bmi_date" = obs_dt),
    by="patid"
  ) %>%
  left_join(timepoint_bmi,by="patid")

# see whether BMI recorded prior to index age
patient = patient %>%
  mutate(bmi_to_index = index_age - age_at_bmi) %>%
  mutate(bmi_before_index = ifelse(
    !is.na(bmi) & bmi_to_index>=master_index_age_window,
    "yes",
    "no"))


# recode as NA any BMI recordings after index date 
patient = patient %>%
  mutate(bmi = 
           ifelse(bmi_before_index == "no",
                  NA,
                  bmi))
get_prop(dat = patient,x="bmi_before_index")
make_hist("bmi_to_index")
summary(patient$age_at_bmi)
patient %>% filter(bmi_before_index=="yes") %>%
  summarise_at(vars(bmi_to_index, age_at_bmi),.funs = c("median","IQR"))

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


p2=make_missingness_plot_over_time("bmi")+ggtitle("C")
mat = rbind(c(1,2,2),c(3,3,3))
plot_fx(grid.arrange(p,p1,p2,layout_matrix =mat),"./figs/bmi_missing.png")

get_prop_miss("bmi_z",patient)

#################################
# 2. Smoking 
#################################

# smoking
## define start time
smoking_codelist = read_csv("../codelists/smoking_status.csv", col_types = cols(.default = "c") )
smoking_codes = smoking_codelist %>% filter(smoking_status=="smoker")
non_smoking_codes = smoking_codelist %>% filter(smoking_status=="never_smoker")

# now extract smoking status & age code recorded
smoking_df = clinical %>%
  filter(medcodeid %in% smoking_codelist$medcodeid) %>%
  left_join(smoking_codelist,by="medcodeid") %>%
  left_join(patient %>% dplyr::select(patid,rough_dob,index_age),by="patid") %>%
  filter(datify(obsdate)>=regstartdate) %>%
  mutate(age_at_smok_code = delta_dates(datify(obsdate) - rough_dob))

# restrict to codes before index date or never-smoking codes
smoking_df = smoking_df %>%
  filter(age_at_smok_code <= (index_age - master_index_age_window) | smoking_status=="never_smoker")

ever_smok = smoking_df %>%
  filter(smoking_status=="smoker")
never_smok = smoking_df %>%
  filter(smoking_status=="never_smoker")
ever_smok_latest = ever_smok %>% 
  group_by(patid) %>%
  slice_max(age_at_smok_code, with_ties = F)
ever_smok_latest_ex = ever_smok_latest %>% filter(smoking_status_detailed=="ex_smoker")
ever_smok_latest_current = ever_smok_latest %>% filter(smoking_status_detailed=="current_smoker")

# define smoking status
patient = patient %>%
  mutate(smoking_status = ifelse(patid %in% ever_smok$patid,"smoker",
                                 ifelse(patid %in% never_smok$patid,"never_smoker",NA)))



# add in extra detailed smoking data 
smokers_main_dataset = patient %>% filter(smoking_status=="smoker")
ever_smok_latest_detailed = ever_smok_latest %>% 
  filter(patid %in% smokers_main_dataset$patid) %>%
  filter(!is.na(smoking_status_extra_detailed))
ever_smok_latest_detailed = ever_smok_latest_detailed %>% ungroup %>% dplyr::select(patid,smoking_status_extra_detailed)


# define detailed smoking status
patient = patient %>%
  mutate(smoking_status_detailed = ifelse(patid %in% ever_smok_latest_current$patid,"current smoker",
                                 ifelse(patid %in% ever_smok_latest_ex$patid,"ex-smoker",
                                        ifelse(patid %in% never_smok$patid,"never smoker",NA))))


get_prop(dat = patient %>% filter(!is.na(smoking_status)), x="smoking_status")
get_prop(dat = patient %>% filter(!is.na(smoking_status)), x="smoking_status_detailed")
get_prop_miss(dat = patient , x="smoking_status")

# add in very detailed smoking status 
patient = patient %>% left_join(ever_smok_latest_detailed,by="patid")

# make plots
p=make_missingness_plot_over_time("smoking_status")+  ggtitle("C")

p1=ggplot(patient %>% filter(!is.na(smoking_status)),
       aes(MS_status,fill=smoking_status_detailed))+
  geom_bar(position="fill",color="black")+
  scale_y_continuous(labels = scales::percent)+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()+
  labs(x="MS status", y = "Smoking status (%)",fill="Smoking status")+
  ggtitle("B")


# look at pack years 
packyears = clinical %>% filter(medcodeid %in% c("1780396011","460127014"))

packyear_df = patient %>%
  left_join(packyears,by="patid") %>%
  filter(!is.na(value)) %>%
  mutate(packyears = as.numeric(value)) %>%
  filter(packyears<1000) %>%
  mutate(py_date = datify(obsdate)) %>%
  mutate(py_age = delta_dates(py_date - rough_dob)) %>%
  filter(py_age <= (index_age - master_index_age_window)) %>%
  group_by(patid) %>%
  slice_max(packyears,with_ties = F)

p2=ggplot(packyear_df %>% filter(!is.na(smoking_status)),
       aes(smoking_status_detailed,packyears))+
  geom_violin(fill="lightblue")+
  stat_boxplot(width=0.1,alpha=0.5)+
  theme_minimal()+
  labs(x="Smoking status",y="Estimated pack-years \nprior to index date")+
  ggtitle("A")

ggplot(packyear_df %>% filter(!is.na(smoking_status_extra_detailed)),
          aes(smoking_status_extra_detailed,packyears))+
  geom_violin(fill="lightblue")+
  stat_boxplot(width=0.1,alpha=0.5)+
  theme_minimal()+
  labs(x="Smoking status",y="Estimated pack-years \nprior to index date")+
  ggtitle("A")+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust = 1))

# fill in never smokers 
patient = patient %>%
  mutate(smoking_status_extra_detailed = ifelse(
    smoking_status == "never_smoker",
    "never smoker",
    ifelse(smoking_status == "smoker" & is.na(smoking_status_extra_detailed),
           NA,
           smoking_status_extra_detailed)
  ))

get_prop(packyear_df %>% ungroup,x="smoking_status")
packyear_df %>% ungroup %>% group_by(smoking_status_detailed) %>%
  summarise(median(packyears))
packyear_df %>% filter(smoking_status_detailed=="never smoker") %>% filter(packyears>0)

get_prop_miss(patient %>% filter(numeric_year_of_reg>=2004),x="smoking_status")
get_prop_miss(patient %>% filter(numeric_year_of_reg<=2003),x="smoking_status")

mat = rbind(c(1,1,2,2),c(3,3,3,3))
plot_fx(grid.arrange(p2,p1,p,layout_matrix =mat),"./figs/smoking_validation.png",width=10)

# crudely estimate no. of smoking years prior to index 
earliest_smok_code = ever_smok %>%
  group_by(patid) %>%
  slice_min(age_at_smok_code,with_ties = F) %>%
  dplyr::select(patid,age_at_smok_code)

patient = patient %>%
  left_join(earliest_smok_code %>% ungroup,by="patid") %>%
  mutate(smoking_prior_to_index = ifelse(smoking_status == "smoker",
                                         index_age - age_at_smok_code,
                                         NA))

 


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
  mutate(im_dx_date = datify(obsdate)) %>%
  filter(im_dx_date < datify("01/05/2022")) %>%
  left_join(patient %>% dplyr::select(patid,rough_dob,age_at_registration),by="patid") %>%
  mutate(age_at_im = delta_dates(im_dx_date - rough_dob)) %>%
  filter(age_at_im>0) %>%
  group_by(patid) %>%
  slice_min(age_at_im,n=1,with_ties = F) %>%
  dplyr::select(patid,im_dx_date,age_at_im)

# join with main dataset
patient = patient %>%
  left_join(im_df,by="patid") %>%
  mutate(im_status = ifelse(!is.na(age_at_im) & 
                              age_at_im < (index_age-master_index_age_window),
                            "IM",
                            "no_IM"))

get_prop(dat=patient,x="im_status")
summary(patient$age_at_im)
p= make_hist(x="age_at_im")+labs(x="Age at IM")
plot_fx(p,filename = "./figs/im_age.png",width=3,height=3)

#################################
# 4. Vit D 
#################################

# read in codelist 
vitd_def_codes = read_csv("../codelists/vitd_def_codelist.csv",
                          col_types = "ccccccccc")

# vitd  data 
vd_df = clinical %>%
  filter(medcodeid %in% vitd_def_codes$MedCodeId) %>%
  mutate(vd_def_dx_date = datify(obsdate)) %>%
  filter(vd_def_dx_date < datify("01/05/2022")) %>%
  left_join(patient %>% dplyr::select(patid,rough_dob),by="patid") %>%
  mutate(age_at_vd_def = delta_dates(vd_def_dx_date - rough_dob)) %>%
  filter(age_at_vd_def>0) %>%
  group_by(patid) %>%
  slice_min(age_at_vd_def,n=1,with_ties = F) %>%
  dplyr::select(patid,vd_def_dx_date,age_at_vd_def)

# join with main dataset
patient = patient %>%
  left_join(vd_df,by="patid") %>%
  mutate(vd_def_status = ifelse(!is.na(age_at_vd_def) & 
                                  age_at_vd_def < (index_age-master_index_age_window),
                                "vd_def",
                                "vd_replete"))

get_prop("vd_def_status",patient)

# vitamin D tests
vd_test_codes = c("457822017",
                  "145751000006117",
                  "457821012")
vitd_tests = clinical %>% 
  filter(medcodeid %in% vd_test_codes) %>%
  filter(!is.na(value)) %>%
  dplyr::select(patid,obsdate,value)


vd_test_validation = patient %>%
  left_join(vitd_tests,by="patid") %>% 
  filter(!is.na(value)) %>%
  filter(!is.na(vd_def_dx_date)) %>%
  filter(!is.na(obsdate)) %>%
  mutate(vd_test_date = datify(obsdate)) %>%
  mutate(vd_test_to_def_dx = abs(delta_dates((vd_test_date - vd_def_dx_date )))) %>% 
  filter(vd_test_to_def_dx<1/12)

p=ggplot(
  vd_test_validation %>% filter(vd_def_status=="vd_def"),
  aes(as.numeric(value)))+
  geom_density(fill="lightblue",alpha=0.5)+
  geom_vline(xintercept = 50)+
  scale_x_log10()+
  theme_minimal()+
  labs(x="Measured vitamin D (nM)")

plot_fx(p,filename="./figs/vd_validation.png",width=3,height=3)

vd_test_validation %>% 
  group_by(vd_def_status) %>%
  mutate(value = as.numeric(value)) %>% 
  summarise_at(vars(value),.funs = c("median","IQR")) 

#################################
# 5. Alcohol 
#################################

# read in codelist 
alcohol_codes = read_csv("../codelists/alcohol_codelist_updated.csv",
                         col_types = "ccccccccc")

# restrict to just 'disease' codes 
alcohol_binary_codes = alcohol_codes %>% filter(!EmisCodeCategoryId %in% c(3,37))
alcohol_score_codes = alcohol_codes %>% filter(EmisCodeCategoryId %in% c(3,37))

# get most recent alcohol  data >2yr before index
alcohol_df = clinical %>%
  filter(medcodeid %in% alcohol_binary_codes$MedCodeId) %>%
  mutate(alc_date = datify(obsdate)) %>%
  filter(alc_date < datify("01/05/2022")) %>%
  left_join(patient %>% dplyr::select(patid,rough_dob,index_age),by="patid") %>%
  mutate(age_at_alc = delta_dates(alc_date - rough_dob)) %>%
  filter(age_at_alc>0) %>%
  group_by(patid) %>%
  filter(index_age - age_at_alc >= master_index_age_window) %>%
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
  dplyr::select(patid,alc_date,age_at_alc,alcohol_status)


# join with main dataset
patient = patient %>%
  left_join(alcohol_df,by="patid")


get_prop_miss(dat = patient, x="alcohol_status")
get_prop(dat = patient %>% filter(!is.na(alcohol_status)), x="alcohol_status")
p=make_missingness_plot_over_time("alcohol_status")+ggtitle("A")

# validate using scores 
audit_c_scores = clinical %>% 
  filter(medcodeid == "408548014") %>%
  filter(!is.na(value)) %>%
  dplyr::select(patid,obsdate,value)

# score must be 0-12
audit_c_scores = audit_c_scores %>%
  mutate(value = as.numeric(value)) %>%
  filter(value >=0 & value <= 12)


alcohol_validation = patient %>%
  left_join(audit_c_scores,by="patid") %>% 
  filter(!is.na(value)) %>%
  filter(!is.na(obsdate)) %>%
  mutate(auditc_date = datify(obsdate)) %>%
  mutate(alc_audit_delta = abs(delta_dates((alc_date - auditc_date )))) %>% 
  filter(alc_audit_delta<1/12) %>%
  filter( !(alcohol_status == "never_drinker" & age_at_alc > (index_age-master_index_age_window) ) ) %>%
  filter(!is.na(alcohol_status))

p1=ggplot(alcohol_validation %>% mutate(alcohol_status = ifelse(alcohol_status=="ever_drinker","Drinker","Non-drinker")),
       aes(alcohol_status,value))+
  stat_boxplot(fill="lightblue",alpha=0.5)+
  theme_minimal()+
  labs(y="AUDIT-C score (0-12)",x="Alcohol drinker status")+ggtitle("B")

plot_fx(grid.arrange(p,p1),"./figs/alcohol_validation.png")

alcohol_validation %>% 
  group_by(alcohol_status) %>%
  mutate(value = as.numeric(value)) %>% 
  summarise_at(.vars = vars(value),.funs = c(mean,median,IQR,range))


#####################################
# 6. Head injury
#####################################

# read in codelist 
hi_codes = read_csv("../codelists/headinjury_codelist_updated.csv",
                    col_types = "ccccccccc")

# head injury data 
hi_df = clinical %>%
  filter(medcodeid %in% hi_codes$MedCodeId) %>%
  mutate(hi_dx_date = datify(obsdate)) %>%
  filter(hi_dx_date < datify("01/05/2022")) %>%
  left_join(patient %>% dplyr::select(patid,rough_dob),by="patid") %>%
  mutate(age_at_hi = delta_dates(hi_dx_date - rough_dob)) %>%
  filter(age_at_hi>0) %>%
  group_by(patid) %>%
  slice_min(age_at_hi,n=1,with_ties = F) %>%
  dplyr::select(patid,hi_dx_date,age_at_hi)

# join with main dataset
patient = patient %>%
  left_join(hi_df,by="patid") %>%
  mutate(head_injury_status = ifelse(!is.na(age_at_hi) & 
                                       age_at_hi < (index_age-master_index_age_window),
                                     "head_injury",
                                     "no_head_injury"))

get_prop(dat=patient,"head_injury_status")
summary(patient$age_at_hi)
make_hist("age_at_hi")


#####################################
# 7. Solvents
#####################################

# read in codelist 
solvent_codes = read_csv("../codelists/solvent_codelist.csv",
                         col_types = "ccccccccc")

# solvent  data 
solvent_df = clinical %>%
  filter(medcodeid %in% solvent_codes$MedCodeId) %>%
  mutate(solvent_date = datify(obsdate)) %>%
  filter(solvent_date < datify("01/05/2022")) %>%
  left_join(patient %>% dplyr::select(patid,rough_dob),by="patid") %>%
  mutate(age_at_solvent_exp = delta_dates(solvent_date - rough_dob)) %>%
  filter(age_at_solvent_exp>0) %>%
  group_by(patid) %>%
  slice_min(age_at_solvent_exp,n=1,with_ties = F) %>%
  dplyr::select(patid,solvent_date,age_at_solvent_exp)

# join with main dataset
patient = patient %>%
  left_join(solvent_df,by="patid") %>%
  mutate(solvent_exposure_status = ifelse(!is.na(age_at_solvent_exp) & 
                                            age_at_solvent_exp < (index_age-master_index_age_window),
                                          "solvents",
                                          "no_solvents"))

patient %>%
  count(solvent_exposure_status) %>%
  mutate(prop = n/sum(n))

patient %>% group_by(MS_status) %>%
  count(solvent_exposure_status) %>%
  mutate(prop = n/sum(n))


#####################################
# make sensitivity cohorts
#####################################

# recode MS status with control as ref 
patient$MS_status = relevel(factor(patient$MS_status),ref="Control")
# relevel variables to orient towards risk factor
patient$vd_def_status = relevel(factor(patient$vd_def_status),ref="vd_replete")
patient$ethnicity_pc_hes_hybrid = relevel(factor(patient$ethnicity_pc_hes_hybrid),ref="White")
patient$im_status = relevel(factor(patient$im_status),ref="no_IM")
patient$head_injury_status = relevel(factor(patient$head_injury_status),ref="no_head_injury")
patient$bmi_category = relevel(factor(patient$bmi_category,ordered = F),ref="Healthy")
patient$alcohol_status = relevel(factor(patient$alcohol_status),ref="ever_drinker")
patient$smoking_status_detailed = relevel(factor(patient$smoking_status_detailed),ref="never smoker")
patient$numeric_imd = as.numeric(patient$e2019_imd_5)

# make sensitivity cohorts 
# repeat for sensitivity analysis 
post_2006_cohort = patient %>% filter(get_year_from_date(regstartdate)>=2006)
post_1997_cohort = patient %>% filter(get_year_from_date(regstartdate)>=1997)

# define HES-confirmed MS cohort 
matching_file = read_tsv("../../Data/Primary_care/Aurum/21_000677_CPRD_Aurum_matching_file.txt",col_types="cccccccc") %>%
  rename("patid" = control_patid)
just_hes_ms_cases = patient %>% filter(ms_code_source=="primary_care_and_hes")
matching_file = matching_file %>% filter(case_patid %in% just_hes_ms_cases$patid)
hes_ms_cohort = patient %>% filter(patid %in% matching_file$patid | patid %in% just_hes_ms_cases$patid)


# 10y antecedent
ten_year_antecedent_cohort = patient %>% filter(years_data_before_index >= 9.99)

cohort_list = list("primary" = patient,
                   "hes_ms" = hes_ms_cohort,
                   "post_97" = post_1997_cohort,
                   "post_06" = post_2006_cohort,
                   "ten_y_antecedent" = ten_year_antecedent_cohort)

#####################################
# save progress 
#####################################
saveRDS(cohort_list,"./datasets/cprd_aurum_cleaned_exposures_patient.rds")

# matched cohort - not run
#ms_cases = patient %>% filter(MS_status=="Case") %>% dplyr::select(patid,yob,gender,ethnicity_pc_hes_hybrid)
#controls = patient %>% filter(MS_status!="Case") %>% dplyr::select(patid,yob,gender,ethnicity_pc_hes_hybrid)

#overall_control_df = data.frame()
#overall_case_df = data.frame()
#for(i in c(1:nrow(ms_cases))){
  message(i," of ",nrow(ms_cases))
  this_case = ms_cases[i,]
  
  matched_controls = controls %>%
      filter(!patid %in% overall_control_df$patid & yob == this_case$yob & gender == this_case$gender & ethnicity_pc_hes_hybrid == this_case$ethnicity_pc_hes_hybrid)
  
  if(nrow(matched_controls) >=4){
    matched_controls = matched_controls %>%
      sample_n(size=4,replace=F)
    overall_control_df <<- bind_rows(overall_control_df,matched_controls)
    overall_case_df <<- bind_rows(overall_case_df,this_case)
    controls <<- controls %>% filter(!patid %in% overall_control_df$patid)
    ms_cases <<- ms_cases %>% filter(!patid %in% overall_case_df$patid)
  }
}
