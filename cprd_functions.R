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

# get tidy proportions
get_prop_miss = function(x, dat){
  dat %>%
    count(is.na(.data[[x]])) %>%
    mutate(prop = n/sum(n)) %>%
    mutate(total = sum(n)) %>%
    mutate(pct = paste0(round(prop*100,2),"%")) %>%
    mutate(pastable = paste0(n,"/",total,"(",pct,")"))
}

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
    scale_x_log10(limits=c(0.001,1000))+
    geom_label(mapping = aes(x=0.01,y=predictor,label=paste0("OR=",round(or,3))))
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
  stepwise_model = MASS::stepAIC(model, direction = "both", 
                                 trace = T)
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
  message(pval)
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


# define function which makes demographics table 
make_demographics = function(
  input_data = patient,
  response_var = "MS_status",
  cat_vars = c("gender","ethnicity_pc_hes_hybrid","e2019_imd_5","e2011_urban_rural"),
  cont_vars = c("yob","age_at_registration","index_age","years_data_before_index"),
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
    input_data[[response_var]] = factor(input_data[[response_var]])
    res_2 = table(input_data[[response_var]],input_data[[this_col]])
    pval = round(chisq.test(res_2)$p.value,5)
    
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
    pval = round(t.test(a[[this_col]],b[[this_col]])$p.value,5)    
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

z_score = function(x){
  (x - mean(x,na.rm=T) ) / sd(x,na.rm = T)
}

make_categorical_cc_plot  = function(plot_dat){
  p=ggplot(plot_dat %>% filter(model=="age + sex"),
           aes(var,or))+
    geom_point()+
    geom_errorbar(mapping = aes(ymin = lower_ci, ymax= upper_ci,x=var),width=0.1)+
    theme_minimal()+
    geom_hline(yintercept = 1,alpha=0.2)+
    labs(x="Exposure",y="Odds ratio for MS")
  return(p)
}

fixed_ggadjustedcurves = function(fit,
                                  variable = NULL,
                                  data = NULL,
                                  reference = NULL,
                                  method = "conditional",
                                  fun = NULL,
                                  palette = "hue",
                                  ylab = "Survival rate", size = 1,
                                  ggtheme = theme_survminer(), ...) {
  stopifnot(method %in% c("marginal", "average", "conditional", "single"))
  ylim <- NULL
  if (is.null(fun)) ylim <- c(0, 1)
  
  curve <- surv_adjustedcurves(fit = fit,
                               variable = variable,
                               data = data,
                               reference = reference,
                               method = method,
                               size = size,
                               ...)
  time <- surv <- NULL
  pl <- ggplot(curve, aes(x = time, y = surv, color = variable)) +
    geom_step(size = size) + ggtheme +
    scale_y_continuous(limits = ylim) +
    ylab(ylab)
  
  if (method == "single")
    pl <- pl + theme(legend.position = "none")
  
  ggpubr::ggpar(pl,  palette = palette, ...)
  
}


surv_adjustedcurves <- function(fit,
                                variable = NULL,
                                data = NULL,
                                reference = NULL,
                                method = "conditional",
                                size = 1,
                                ...) {
  stopifnot(method %in% c("marginal", "average", "conditional", "single"))
  data <- .get_data(fit, data)
  # deal with default arguments
  # reference = NULL
  if (is.null(reference))
    reference <- data
  
  # variable = NULL
  if (is.null(variable)) {
    # is there a 'strata' component?
    term.labels <- attr(terms(fit$formula), "term.labels")
    strata.term.labels <- grep(term.labels, pattern = "strata(", fixed = TRUE, value = TRUE)
    if (length(strata.term.labels) > 0) {
      variable <- gsub(
        gsub(
          strata.term.labels,
          pattern = "strata(", replacement = "", fixed = TRUE)[1],
        pattern = "[\\) ]", replacement = "")
      cat("The variable argument is missing. Using", variable, "as extracted from strata\n")
    } else {
      # if not then leave variable = NULL
      method = "single"
    }
  }
  
  curve <- switch(method,
                  single = ggadjustedcurves.single(data, fit, size = size),
                  average =  ggadjustedcurves.average(data, fit, variable, size = size),
                  conditional = ggadjustedcurves.conditional(data, fit, variable, size = size),
                  marginal = ggadjustedcurves.marginal(data, fit, variable, reference, size = size))
  
  curve
}


ggadjustedcurves.single <- function(data, fit, size = 1) {
  time <- surv <- variable <- NULL
  
  pred <- survexp(~1, data = data, ratetable = fit)
  
  curve <- data.frame(time = c(0,pred$time),
                      variable = "total",
                      surv = c(1, pred$surv))
  
  # plot moved to ggadjustedcurves
  # ggplot(curve, aes(x = time, y = surv, color = variable)) +
  #   geom_step(size = size) + theme(legend.position = "none")
  curve
}

ggadjustedcurves.average <- function(data, fit, variable, size = 1) {
  time <- surv <- NULL
  
  lev <- sort(unique(pull(data, variable)))
  
  pred <- survexp(as.formula(paste("~", variable)), data = data,
                  ratetable = fit)
  
  curve <- data.frame(time = rep(c(0,pred$time), length(lev)),
                      variable = factor(rep(lev, each=1+length(pred$time))),
                      surv = c(rbind(1, pred$surv)))
  
  # plot moved to ggadjustedcurves
  # ggplot(curve, aes(x = time, y = surv, color=variable)) +
  #   geom_step(size = size)
  curve
}

ggadjustedcurves.marginal <- function(data, fit, variable, reference, size = 1) {
  time <- surv <- NULL
  
  lev <- sort(unique(pull(data, variable)))
  reference[,variable] = "_reference_"
  df0 <- reference
  form <- paste(variable, "~", gsub(as.character(formula(fit))[3], pattern="\\+ *strata.*[^\\)].", replacement=""))
  
  allRes <- list()
  rwt <- numeric(nrow(data))
  for (level in lev) {
    indexes <- which(data[,variable] == level)
    if (length(indexes) > 0) {
      df1 <- data[indexes, ]
      ndf <- rbind(df0, df1)
      ndf[,variable] <- factor(ndf[,variable])
      model <- glm(as.formula(form), ndf, family="binomial")
      allRes[[level]] <- predict(model, newdata = data, type = "response")
      rwt[indexes] <- 1/allRes[[level]][indexes]
    }
  }
  
  nform <- paste(as.character(formula(fit))[2], "~", variable)
  nfit <- coxph(as.formula(nform), data = data, weights = rwt)
  
  pred <- survexp(as.formula(paste("~", variable)), data = data, ratetable = nfit)
  
  # remove leading zeros
  # while survexp returns non monotonic results
  if (length(dim(pred$surv))==2) {
    for (i in 1:ncol(pred$surv))
      for (j in nrow(pred$surv):2)
        if (pred$surv[j,i] > pred$surv[j - 1,i])
          pred$surv[j - 1,i] <- 1
  }
  
  curve <- data.frame(time = rep(c(0,pred$time), length(lev)),
                      variable = factor(rep(lev, each=1+length(pred$time))),
                      surv = c(rbind(1, pred$surv)))
  
  # plot moved to ggadjustedcurves
  # ggplot(curve, aes(x = time, y = surv, color=variable)) +
  #   geom_step(size = size)
  curve
}

ggadjustedcurves.conditional <- function(data, fit, variable, size = 1) {
  time <- surv <- NULL
  
  lev <- sort(unique(pull(data, variable)))
  
  ndata <- data[rep(1:nrow(data), each=length(lev)),
                setdiff(colnames(data), variable)]
  ndata[,variable] = rep(lev, nrow(data))
  
  pred <- survexp(as.formula(paste("~", variable)), data = ndata,
                  ratetable = fit)
  # remove leading zeros
  # while survexp returns non monotonic results
  if (length(dim(pred$surv)) == 2) {
    for (i in 1:ncol(pred$surv))
      for (j in nrow(pred$surv):2)
        if (pred$surv[j,i] > pred$surv[j - 1,i])
          pred$surv[j - 1,i] <- 1
  }
  
  curve <- data.frame(time = rep(c(0,pred$time), length(lev)),
                      variable = factor(rep(lev, each=1+length(pred$time))),
                      surv = c(rbind(1, pred$surv)))
  # plot moved to ggadjustedcurves
  # ggplot(curve, aes(x = time, y = surv, color=variable)) +
  #   geom_step(size = size)
  curve
}

.get_data <- function(fit, data = NULL, complain = TRUE) {
  if(is.null(data)){
    if (complain)
      warning ("The `data` argument is not provided. Data will be extracted from model fit.")
    data <- eval(fit$call$data)
    if (is.null(data))
      stop("The `data` argument should be provided either to ggsurvfit or survfit.")
  }
  data
}

simplify_pval = function(x){
  y = ifelse(x<0.0001,"<0.0001",round(x,3))
  return(y)
}

