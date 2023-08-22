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

# transform variables
ms_post_1997 = cohorts$post_97 %>% filter(MS_status=="Case")
ms_post_1997$year_of_ms_date = get_year_from_date(ms_post_1997$ms_code_date)
ms_post_1997 = ms_post_1997 %>%
  mutate(age_of_dx_z_score = z_score(age_at_dx))
ms_post_1997$ethnicity_pc_hes_hybrid = relevel(factor(ms_post_1997$ethnicity_pc_hes_hybrid),ref="White")
ms_all = cohorts$primary %>% filter(MS_status=="Case")
ms_all$year_of_ms_date = get_year_from_date(ms_all$ms_code_date)
ms_all = ms_post_1997 %>%
  mutate(age_of_dx_z_score = z_score(age_at_dx))
ms_all$ethnicity_pc_hes_hybrid = relevel(factor(ms_all$ethnicity_pc_hes_hybrid),ref="White")


#################################
# age of onset
#################################

overall_res_df = list()

# plot
age_at_dx_plot =ggplot(ms_post_1997 %>%
           filter(!is.na(ethnicity_pc_hes_hybrid)),
         aes(ethnicity_pc_hes_hybrid, 
             age_at_dx,
             fill=gender))+
  geom_violin(alpha=0.8)+
  stat_boxplot(alpha=0.2,width=0.1,show.legend = F,position = position_dodge(width=0.9))+
  theme_minimal()+
  scale_fill_brewer(palette="Set3")+
  labs(x="Ethnic background",y="Age at MS diagnosis")+
  ggtitle("A")
plot_fx(p,"./figs/age_at_dx_ethnic.png",width=5,height=4)


# print basics 
ms_post_1997 %>%
  filter(!is.na(ethnicity_pc_hes_hybrid)) %>%
  group_by(ethnicity_pc_hes_hybrid) %>%
  summarise_at(vars(age_at_dx),.funs = c("mean","sd"))
# make function to extract info from model

get_info_from_lm = function(model_to_test, dat = ms_post_1997, outcome = "age_at_dx"){
  sd_1_unit = sd(dat[[outcome]],na.rm=T)
  model_summ = model_to_test %>% summary()
  dat = data.frame(model_summ$coefficients)
  dat$est = dat$Estimate * sd_1_unit
  dat$lower = (dat$Estimate - 1.96* dat$Std..Error) * sd_1_unit
  dat$upper = (dat$Estimate + 1.96* dat$Std..Error) * sd_1_unit
  dat = dat %>% dplyr::select(est,lower,upper,Pr...t..)
  dat = dat %>% mutate(estimate = paste0(round(est,2)," (95% CI ",round(lower,2)," to ",round(upper,2),", p",simplify_pval(Pr...t..),")"))
  dat$var = rownames(dat)
  rownames(dat) = NULL
  dat$var = str_remove(dat$var,"ethnicity_pc_hes_hybrid")
  return(dat)
}

# univariable model
mod1 = glm(data = ms_post_1997,
            age_of_dx_z_score ~ ethnicity_pc_hes_hybrid) %>%
  get_info_from_lm() %>%
  mutate(model = "univariable", cohort = "Post-1997")


# with sex model
mod2 = glm(data = ms_post_1997,
    age_of_dx_z_score ~ gender + ethnicity_pc_hes_hybrid) %>%
  get_info_from_lm() %>%
  mutate(model = "sex", cohort = "Post-1997")


# with sex & IMD model
mod3 = glm(data = ms_post_1997,
    age_of_dx_z_score ~ e2019_imd_5 + gender + ethnicity_pc_hes_hybrid) %>%
  get_info_from_lm() %>%
  mutate(model = "sex + imd", cohort = "Post-1997")


# with sex & IMD model
mod4 = glm(data = ms_post_1997,
    age_of_dx_z_score ~ e2011_urban_rural + e2019_imd_5 + gender + ethnicity_pc_hes_hybrid) %>%
  get_info_from_lm() %>%
  mutate(model = "sex + imd + urban/rural", cohort = "Post-1997")


# repeat with whole cohort
mod5 = glm(data = ms_all,
    age_of_dx_z_score ~ e2011_urban_rural + e2019_imd_5 + gender + ethnicity_pc_hes_hybrid) %>%
  get_info_from_lm(dat = ms_all) %>%
  mutate(model = "sex + imd + urban/rural", cohort = "Whole Cohort")


# repeat with just HES-confirmed dx
post_1997_pc_hes = ms_post_1997 %>% filter(ms_code_source == "primary_care_and_hes")
mod6 = glm(data = post_1997_pc_hes,
    age_of_dx_z_score ~ e2011_urban_rural + e2019_imd_5 + gender + ethnicity_pc_hes_hybrid) %>%
  get_info_from_lm(dat = post_1997_pc_hes) %>%
  mutate(model = "sex + imd + urban/rural", cohort = "Post-1997 (HES-confirmed)")



# look at age of neuro ref
ms_post_1997$year_neuro_ref = get_year_from_date(ms_post_1997$neuro_ref_date)
ms_post_1997 = ms_post_1997 %>% mutate(age_neuro_ref = year_neuro_ref - yob)

# plot
p=ggplot(ms_post_1997 %>%
           filter(!is.na(ethnicity_pc_hes_hybrid) & age_neuro_ref>0) ,
         aes(ethnicity_pc_hes_hybrid, 
             age_neuro_ref,
             fill=gender))+
  geom_violin(alpha=0.8)+
  stat_boxplot(alpha=0.2,width=0.1,show.legend = F,position = position_dodge(width=0.9))+
  theme_minimal()+
  scale_fill_brewer(palette="Set3")+
  labs(x="Ethnic background",y="Age at first neurological referral")
plot_fx(p,"./figs/age_at_neuro_ref_ethnic.png",width=5,height=4)

ms_post_97_age_neuro = ms_post_1997 %>%
  filter(!is.na(ethnicity_pc_hes_hybrid) & age_neuro_ref>0) %>%
  mutate(age_neuro_z = z_score(age_neuro_ref))

get_prop(dat = ms_post_97_age_neuro %>% mutate(neuro_pre_ms = 
                                           ifelse(age_neuro_ref<age_at_dx,
                                                  "neuro_pre_ms",
                                                  "neuro_post_ms")),
                                         x="neuro_pre_ms")

  
mod7 = glm(data = ms_post_97_age_neuro,
    age_neuro_z ~ e2011_urban_rural + e2019_imd_5 + gender + ethnicity_pc_hes_hybrid) %>%
  get_info_from_lm(dat = ms_post_97_age_neuro, outcome = "age_neuro_ref") %>%
  mutate(model = "sex + imd + urban/rural", cohort = "Post-1997", proxy = "Neurology referral")

ggplot(ms_post_97_age_neuro,
       aes(age_neuro_z,age_at_dx))+
  geom_point()

# look in HES data
hes_attendances =  read_tsv("../../Data/Linked_data/Aurum/hesae_attendance_21_000677.txt",
                   col_types = cols(.default = "c")) %>% dplyr::select(1:3)
hes_ae =  read_tsv("../../Data/Linked_data/Aurum/hesae_diagnosis_21_000677.txt",
                   col_types = cols(.default = "c")) %>%
  filter(diagscheme==1) %>% 
  filter(diag2=="24") %>%
  filter(patid %in% ms_all$patid) %>%
  left_join(hes_attendances,by=c("aekey","patid")) %>%
  mutate(hes_ae_date = datify(arrivaldate)) %>%
  group_by(patid) %>%
  slice_min(hes_ae_date,n=1) %>%
  dplyr::select(patid,hes_ae_date)


ms_post_1997 = ms_post_1997 %>% left_join(hes_ae,by="patid") %>%
  mutate(has_hes_encounter_pre_ms = ifelse(hes_ae_date < ms_code_date,"yes","no"))

get_prop(ms_post_1997,x="has_hes_encounter_pre_ms")

ms_post_1997_with_hes = ms_post_1997 %>% filter(!is.na(hes_ae_date)) %>% 
  mutate(age_at_hes_neuro_pres = delta_dates ( hes_ae_date - rough_dob ) ) %>% 
  mutate(hes_neuro_pres_to_ms_date = delta_dates ( hes_ae_date - ms_code_date ) ) %>%
  mutate(age_hes_neuro_pres_z = z_score (age_at_hes_neuro_pres))

hist(ms_post_1997_with_hes$age_at_hes_neuro_pres)
hist(ms_post_1997_with_hes$hes_neuro_pres_to_ms_date)

mod8 = glm(data = ms_post_1997_with_hes,
    age_hes_neuro_pres_z ~ e2011_urban_rural + e2019_imd_5 + gender + ethnicity_pc_hes_hybrid) %>%
  get_info_from_lm(dat = ms_post_1997_with_hes, outcome = "age_at_hes_neuro_pres") %>%
  mutate(model = "sex + imd + urban/rural", cohort = "Post-1997", proxy = "A&E attendance")

# combine results and plot
res_df = bind_rows(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8)
res_df$proxy = replace_na(res_df$proxy,"Age of first MS code")


res_df = res_df %>% filter(var %in% c("Asian","Black","Mixed/Other"))
                  
p=ggplot(res_df,
       aes(est,var,shape=model,col = cohort))+
  geom_point(position=ggstance::position_dodgev(height=0.5))+
  facet_wrap(~proxy,nrow=3)+
  geom_errorbarh(mapping = aes(xmin = lower, xmax = upper, y = var),height=0.1,
                 position=ggstance::position_dodgev(height=0.5))+
  theme_minimal()+
  scale_color_brewer(palette = "Set2")+
  labs(x="Age at MS onset relative to\nindividuals of White ethnicity (years)",
       y="Ethnicity")+
  ggtitle("B")


mat = matrix(c(1,1,2,2,2,2),ncol=1)
plot_fx(
  grid.arrange(age_at_dx_plot,p,layout_matrix = mat),
  height=8,
  width=8,
  "./figs/age_at_dx_plots.png")



#################################
# sex ratio 
#################################


plot_dat = ms_post_1997 %>% 
  filter(!is.na(ethnicity_pc_hes_hybrid)) %>%
  group_by(ethnicity_pc_hes_hybrid) %>%
  get_prop(x="gender")
p=ggplot(plot_dat,
       aes(ethnicity_pc_hes_hybrid,prop,fill=gender,label=pct))+
  geom_col(color="black",position=position_dodge())+
  scale_fill_brewer(palette="YlOrRd")+
  scale_y_continuous(labels = scales::percent,limits=c(0,1))+
  labs(x="Ethnicity",y="Percentage")+
  theme_minimal()+
  coord_flip()+
  geom_text(position = position_dodge(width=1),hjust=-0.2)


ms_post_1997 %>% 
  group_by(ethnicity_pc_hes_hybrid) %>% 
  get_prop(x="gender") %>%
  pivot_wider(values_from = prop, names_from = gender,id_cols = ethnicity_pc_hes_hybrid) %>%
  mutate(f_m_ratio = Female / Male)

ms_all %>% 
  group_by(ethnicity_pc_hes_hybrid) %>% 
  get_prop(x="gender") %>%
  pivot_wider(values_from = prop, names_from = gender,id_cols = ethnicity_pc_hes_hybrid) %>%
  mutate(f_m_ratio = Female / Male)


ms_all$ms_dx_year = cut2(get_year_from_date(ms_all$ms_code_date),cuts = seq(1990,2020,by=5))
sex_ratio_by_year = ms_all %>% 
  group_by(ethnicity_pc_hes_hybrid,ms_dx_year) %>% 
  get_prop(x="gender") %>%
  pivot_wider(values_from = prop, names_from = gender,id_cols = c(ethnicity_pc_hes_hybrid,ms_dx_year,total)) %>%
  mutate(f_m_ratio = Female / Male) %>%
  ungroup() 
sex_ratio_by_year$ms_dx_year = str_remove(sex_ratio_by_year$ms_dx_year,"\\[") 
sex_ratio_by_year = sex_ratio_by_year %>% separate(ms_dx_year,sep=",",into=c("start","end")) %>%
  mutate(start = as.numeric(start))


# function to get CIs of prop 
get_se_prop = function(prop,n){
se_prop = sqrt ( (prop * (1-prop) )/ n )
return(se_prop)
}

sex_ratio_by_year = sex_ratio_by_year %>% mutate(se = get_se_prop(prop = Female,n = total)) %>%
  mutate(lower = Female - 1.96 * se) %>%
  mutate(upper = Female + 1.96 * se)
  
sex_ratio_by_year %>% arrange(total)
p2=ggplot(sex_ratio_by_year %>% filter(!is.na(ethnicity_pc_hes_hybrid)) %>%
         filter(ethnicity_pc_hes_hybrid %in% c("White","Black","Asian")),
       aes(start,Female))+
  geom_point()+
  geom_line()+
  geom_errorbar(mapping = aes(x = start,ymin = lower,ymax = upper),width=0.1)+
  facet_wrap(~ethnicity_pc_hes_hybrid,nrow=1)+
  theme_minimal()+
  labs(y="Proportion of females",x="Year of MS diagnostic code report")

plot_fx(grid.arrange(p,p2,nrow=2),
        "./figs/sex_ratio.png",
        width=8)

# test gender ratio with regression
glm(data = ms_post_1997,
    factor(gender)=="Male" ~ ethnicity_pc_hes_hybrid,
    family=binomial(link="logit")) %>% extract_coefs_from_model_obj() %>%
  mutate(p = simplify_pval(pval))

glm(data = ms_post_1997,
    factor(gender)=="Male" ~ yob + ethnicity_pc_hes_hybrid,
    family=binomial(link="logit")) %>% extract_coefs_from_model_obj() %>%
  mutate(p = simplify_pval(pval))

glm(data = ms_post_1997,
    factor(gender)=="Male" ~ yob + index_age + ethnicity_pc_hes_hybrid,
    family=binomial(link="logit")) %>% extract_coefs_from_model_obj() %>%
  mutate(p = simplify_pval(pval))

glm(data = ms_post_1997,
    factor(gender)=="Male" ~ e2019_imd_5 + ethnicity_pc_hes_hybrid,
    family=binomial(link="logit")) %>% extract_coefs_from_model_obj() %>%
  mutate(p = simplify_pval(pval))

glm(data = ms_post_1997,
    factor(gender)=="Male" ~  yob + e2019_imd_5 + ethnicity_pc_hes_hybrid,
    family=binomial(link="logit")) %>% extract_coefs_from_model_obj() %>%
  mutate(p = simplify_pval(pval))


glm(data = ms_post_1997,
    factor(gender)=="Male" ~  e2011_urban_rural +  yob + e2019_imd_5 + ethnicity_pc_hes_hybrid,
    family=binomial(link="logit")) %>% extract_coefs_from_model_obj() %>%
  mutate(p = simplify_pval(pval))

glm(data = ms_all,
    factor(gender)=="Male" ~ ethnicity_pc_hes_hybrid,
    family=binomial(link="logit")) %>% extract_coefs_from_model_obj() %>%
  mutate(p = simplify_pval(pval))
# look at changes in sex ratio over time

plot_dat = ms_all %>% 
  filter(!is.na(ethnicity_pc_hes_hybrid)) %>%
  group_by(ethnicity_pc_hes_hybrid) %>%
  get_prop(x="gender")
p=ggplot(plot_dat,
         aes(ethnicity_pc_hes_hybrid,prop,fill=gender,label=pct))+
  geom_col(color="black",position=position_dodge())+
  scale_fill_brewer(palette="YlOrRd")+
  scale_y_continuous(labels = scales::percent,limits=c(0,1))+
  labs(x="Ethnicity",y="Percentage")+
  theme_minimal()+
  coord_flip()+
  geom_text(position = position_dodge(width=1),hjust=-0.2)

#################################
# mortality 
#################################

get_prop(ms_all %>% group_by(ethnicity_pc_hes_hybrid),x="dead")

ms_all$diagnosed_pre_2000 = ifelse(ms_all$year_of_ms_date<=2000,"pre_2000","post_2000")
ggplot(ms_all,
       aes(ethnicity_pc_hes_hybrid,years_data_after_index))+
         geom_violin()


x=20
glm(data = ms_all %>% filter(years_data_after_index<x),
    z_score(years_data_after_index) ~
      ethnicity_pc_hes_hybrid) %>% get_info_from_lm()


get_prop(ms_all %>% 
           mutate(data_post_index_10 = ifelse(years_data_after_index>=10 | (dead=="dead" & age_at_death <= (age_at_dx+10)),
                                              "yes","no")) %>%
           group_by(ethnicity_pc_hes_hybrid),
         x="data_post_index_10")

get_prop(ms_post_1997 %>% 
           mutate(data_post_index_10 = ifelse(years_data_after_index>=10 | (dead=="dead" & age_at_death <= (age_at_dx+10)),
                                              "yes","no")) %>%
           group_by(ethnicity_pc_hes_hybrid),
         x="data_post_index_10")

# define survival time (10 years after index)
ms_all = ms_all %>%
  mutate(
    survival_time = ifelse(!is.na(age_at_death) & age_at_death<=(age_at_dx+10),
                           age_at_death - age_at_dx,
                           ifelse(
                             age_at_deregistration - age_at_dx >=10,
                             10,
                             age_at_deregistration - age_at_dx
                           ))) %>%
  mutate(status = ifelse(dead=="dead" & age_at_death <= (age_at_dx+10),
                         1,0)) %>%
  filter(survival_time>0)


# basic km - case vs control 
patient = patient %>%
  mutate(
    survival_time = ifelse(!is.na(age_at_death) & age_at_death<=(index_age+10),
                           age_at_death - index_age,
                           ifelse(
                             age_at_deregistration - index_age >=10,
                             10,
                             age_at_deregistration - index_age
                           ))) %>%
  mutate(status = ifelse(dead=="dead" & age_at_death <= (index_age+10),
                         1,0)) %>%
  filter(survival_time>0)


# get basic descriptives 
get_prop(patient %>% group_by(MS_status),x="status")
get_prop(ms_all %>% group_by(ethnicity_pc_hes_hybrid),x="status")


fit <- survfit(Surv(survival_time, status) ~ MS_status, 
               data = patient)
p00=ggsurvplot(fit, 
              data = patient, 
              censor.shape="|", 
              censor.size = 1,
              conf.int = T,
              fun = "pct",
              legend = c(0.3,0.3),
              ylim = c(50,100),
              palette = RColorBrewer::brewer.pal(name = "Set1",n=5),
              ylab = "Survival probability (%)",
              xlab = "Years from diagnosis")

p00
patient$MS_status = as.factor(patient$MS_status)
cox_fit <- coxph(Surv(survival_time, status) ~ yob + index_age + gender + MS_status, 
               data = patient %>% filter(!is.na(MS_status)))

summary(cox_fit)
fixed_ggadjustedcurves(fit = cox_fit,variable = "MS_status")

# basic km plots
# imd 
fit <- survfit(Surv(survival_time, status) ~ e2019_imd_5, 
               data = ms_all)
p0=ggsurvplot(fit, 
              data = ms_all, 
              censor.shape="|", 
              censor.size = 1,
              conf.int = F,
              fun = "pct",
              legend = c(0.3,0.3),
              legend.title = "IMD",
              legend.labs = c("Least deprived","2nd quintile","3rd quintile","4th quintile","Most deprived"),
              ylim = c(80,100),
              palette = RColorBrewer::brewer.pal(name = "Set1",n=5),
              ylab = "Survival probability (%)",
              xlab = "Years from diagnosis")+
  ggtitle("A")

fit <- survfit(Surv(survival_time, status) ~ gender, 
               data = ms_all)
p1=ggsurvplot(fit, 
              data = ms_all, 
              censor.shape="|", 
              censor.size = 1,
              conf.int = F,
              fun = "pct",
              legend = c(0.3,0.3),
              legend.title = "Sex",
              legend.labs = c("Female","Male"),
              ylim = c(80,100),
              palette = RColorBrewer::brewer.pal(name = "Set1",n=5),
              ylab = "Survival probability (%)",
              xlab = "Years from diagnosis")+
  ggtitle("B")

fit <- survfit(Surv(survival_time, status) ~ diagnosed_pre_2000, 
               data = ms_all)
p2=ggsurvplot(fit, 
              data = ms_all, 
              censor.shape="|", 
              censor.size = 1,
              conf.int = F,
              fun = "pct",
              legend = c(0.3,0.3),
              legend.title = "Diagnosis year",
              legend.labs = c("Post-2000","Pre-2000"),
              ylim = c(80,100),
              palette = RColorBrewer::brewer.pal(name = "Set1",n=5),
              ylab = "Survival probability (%)",
              xlab = "Years from diagnosis")+
  ggtitle("C")

ms_all$diagnosed_pre_age_40 = ifelse(ms_all$age_at_dx <= 40,"Pre-40","Post-40")

fit <- survfit(Surv(survival_time, status) ~ diagnosed_pre_age_40, 
               data = ms_all)
p3=ggsurvplot(fit, 
           data = ms_all, 
           censor.shape="|", 
           censor.size = 1,
           conf.int = F,
           fun = "pct",
           legend = c(0.3,0.3),
           legend.title = "Age at diagnosis",
           ylim = c(80,100),
           legend.labs = c(">40 at diagnosis","<40 at diagnosis"),
           palette = RColorBrewer::brewer.pal(name = "Set1",n=4),
           ylab = "Survival probability (%)",
           xlab = "Years from diagnosis")+
  ggtitle("D")

fit <- survfit(Surv(survival_time, status) ~ ethnicity_pc_hes_hybrid, 
               data = ms_all)
p4=ggsurvplot(fit, 
           data = ms_all, 
           censor.shape="|", 
           censor.size = 1,
           conf.int = F,
           fun = "pct",
           legend = c(0.3,0.3),
           legend.title = "Ethnicity",
           ylim = c(80,100),
           legend.labs = c("White","Asian","Black","Mixed/Other"),
           palette = RColorBrewer::brewer.pal(name = "Set1",n=4),
           ylab = "Survival probability (%)",
           xlab = "Years from diagnosis") +
  ggtitle("Ethnicity")

plot_fx(grid.arrange(p0$plot,p1$plot,p2$plot,p3$plot),
        "./figs/predictors_mortality.png",
        width=8,height=8)

p4$plot


# univariable cox models

ms_all$numeric_imd = as.numeric(ms_all$e2019_imd_5)

overall_res = list()
do_univariable_cox_reg = function(var){
model = coxph(Surv(survival_time, status) ~ 
                ms_all[[var]], data =  ms_all)
summ = summary(model)
phtest = cox.zph(model)
pval = phtest$table[2,3]
df = data.frame(summ$coefficients, var = var, ph_pval = pval)
df$hr_rounded = round(df$`exp.coef.`,2)
df$lower_ci = exp(df$coef - 1.96 * df$`se.coef.`)
df$upper_ci = exp(df$coef + 1.96 * df$`se.coef.`)
df$pastable = paste0("HR: ",df$hr_rounded," (95% CI ",round(df$lower_ci,2)," - ",round(df$upper_ci,2),")")
overall_res[[length(overall_res)+1]] <<- df
}

vars = c("year_of_ms_date","numeric_imd","gender","index_age","ethnicity_pc_hes_hybrid")
sapply(vars,do_univariable_cox_reg)
overall_res = do.call("bind_rows",overall_res)
overall_res = overall_res %>% dplyr::select(var,pastable,Pr...z..,ph_pval)
rownames(overall_res) = NULL
overall_res$pval = simplify_pval(overall_res$Pr...z..)

# look at linearity 
plots = list()
make_linearity_plots = function(var){
  model_dat <<- ms_all %>% filter(!is.na(.data[[var]]))
  model = coxph(Surv(survival_time, status) ~ 
                  model_dat[[var]], data =  model_dat)
  plots[[length(plots)+1]] <<- ggcoxfunctional(model)
}
vars = c("year_of_ms_date","numeric_imd","index_age")
sapply(vars,make_linearity_plots)
do.call("grid.arrange",plots)


# make multivariable models
model0 = coxph(Surv(survival_time, status) ~ 
                 year_of_ms_date + gender + numeric_imd + index_age + tt(index_age) + ethnicity_pc_hes_hybrid, data =  ms_all)

summary(model0)

# adjust for risk factors 
model001 = coxph(Surv(survival_time, status) ~ 
                 year_of_ms_date + gender + smoking_status + numeric_imd + index_age + tt(index_age) + ethnicity_pc_hes_hybrid, data =  ms_all)

summary(model001)


summary(model01)
table(ms_all$dead,ms_all$ethnicity_pc_hes_hybrid)
model = coxph(Surv(survival_time, status) ~ 
                year_of_ms_date + gender + numeric_imd + tt(index_age) + ethnicity_pc_hes_hybrid, data =  ms_all)

model2 = coxph(Surv(survival_time, status) ~ 
                 gender + numeric_imd + tt(index_age) + ethnicity_pc_hes_hybrid, data =  ms_all)

model3 = coxph(Surv(survival_time, status) ~ 
                 diagnosed_pre_2000 + gender + e2019_imd_5 + tt(index_age) + ethnicity_pc_hes_hybrid, data =  ms_all)

model4 = coxph(Surv(survival_time, status) ~ 
                 tt(year_of_ms_date) + gender + numeric_imd + tt(index_age) + ethnicity_pc_hes_hybrid, data =  ms_all)


summary(model4)

# look at whether censoring is equal 
ms_all %>% group_by(ethnicity_pc_hes_hybrid) %>%
  summarise(median(years_data_after_index))

x=20
ms_all %>% group_by(ethnicity_pc_hes_hybrid) %>%
  mutate(fu_more_than_x = years_data_after_index>x) %>%
  get_prop(x="fu_more_than_x")


# calculate crude mortality rates
deaths = get_prop(ms_all %>%
                    group_by(ethnicity_pc_hes_hybrid),
                  x="dead") %>% filter(dead=="dead")
ms_all %>%
  group_by(ethnicity_pc_hes_hybrid) %>%
  summarise(fu_years = sum(years_data_after_index)) %>%
  left_join(deaths,by="ethnicity_pc_hes_hybrid") %>%
  mutate(mortality_rate = n/fu_years*100000)



