#Support functions


#Corrected version of the Log-rank test to avoid integer overflow

LogrankTest_num<-function (sample.list, parameter) 
{
  call = (parameter[[1]] == "Description")
  if (call == FALSE | is.na(call)) {
    if (is.na(parameter[[2]])) {
      larger = TRUE
    }
    else {
      if (!all(names(parameter[[2]]) %in% c("larger"))) 
        stop("Analysis model: LogRankTest test: this function accepts only one argument (larger)")
      if (!is.logical(parameter[[2]]$larger)) 
        stop("Analysis model: LogRankTest test: the larger argument must be logical (TRUE or FALSE).")
      larger = parameter[[2]]$larger
    }
    outcome1 = sample.list[[1]][, "outcome"]
    outcome1.complete = outcome1[stats::complete.cases(outcome1)]
    event1 = !sample.list[[1]][, "patient.censor.indicator"]
    event1.complete = event1[stats::complete.cases(outcome1)]
    n1 = length(outcome1.complete)
    outcome2 = sample.list[[2]][, "outcome"]
    outcome2.complete = outcome2[stats::complete.cases(outcome2)]
    event2 = !sample.list[[2]][, "patient.censor.indicator"]
    event2.complete = event2[stats::complete.cases(outcome2)]
    n2 = length(outcome2.complete)
    outcome = c(outcome1.complete, outcome2.complete)
    event = c(event1.complete, event2.complete)
    treatment = c(rep(0, n1), rep(1, n2))
    data = data.frame(time = outcome, event = event, treatment = treatment)
    data = data[order(data$time), ]
    data$event1 = data$event * (data$treatment == 0)
    data$event2 = data$event * (data$treatment == 1)
    data$eventtot = data$event1 + data$event2
    data$n.risk1.prior = length(outcome1) - cumsum(data$treatment == 
                                                     0) + (data$treatment == 0)
    data$n.risk2.prior = length(outcome2) - cumsum(data$treatment == 
                                                     1) + (data$treatment == 1)
    data$n.risk.prior = data$n.risk1.prior + data$n.risk2.prior
    data$e1 = data$n.risk1.prior * data$eventtot/data$n.risk.prior
    data$u1 = data$event1 - data$e1
    data$v1 = ifelse(data$n.risk.prior > 1, 
                     (as.numeric(data$n.risk1.prior) * data$n.risk2.prior * data$eventtot * (data$n.risk.prior - data$eventtot))/(data$n.risk.prior^2 * (data$n.risk.prior - 1)),
                     0)
    stat.test = sum(data$u1)/sqrt(sum(data$v1))
    result = stats::pnorm(stat.test, lower.tail = !larger)
  }
  else if (call == TRUE) {
    result = list("Log-rank test")
  }
  return(result)
}

#Random sample following an Asymptotic exponential distribution
rexp_asymp<-function(n, asymptote, hr) {
  y<-runif(n, min=0, max=1)
  z<-rep(NA, n)
  
  zero_df<-function(t ,y, asymptote, hr) {
    val<-asymptote*(1-exp(-hr*t))-y
    return(val)
  }
  
  my.uniroot<-function(y, asymptote, hr) uniroot(zero_df, c(0,1000),asymptote=asymptote, hr=hr, y=y)$root
  z[y<asymptote]<-sapply(y[y<asymptote],FUN=my.uniroot, asymptote=asymptote, hr=hr)
  z[y>=asymptote]<-100
  return(z)
}

#Custom function for survival (using a asymptotic regression curve)
asymp_exp=function(parameter) {
  call = (parameter[[1]] == "description")
  
  # Generate random variables
  if (call == FALSE) {
    # The number of observations to generate
    n = parameter[[1]]
    
    # Distribution-specific parameters
    asymptote = parameter[[2]]$asymptote #Asymptote
    hr = parameter[[2]]$hr #Hazard rate
    
    if (n%%1 != 0)
      stop("Data model: TemplateDist distribution: Number of observations must be an integer.")
    if (n <= 0)
      stop("Data model: TemplateDist distribution: Number of observations must be positive.")
    
    # Observations are generated using the "rexp_asymp" function and assigned to the "result" object
    result = rexp_asymp(n = n, asymptote = asymptote, hr = hr)
    ##############################################################
    
  } else {
    # Provide information about the distribution function
    if (call == TRUE) {
      
      result = list(list(asymptote = "asymptote", hr = "hr"),
                    list("Template"))
    }
  }
  return(result)
  
}

#Custom function to extract standard errors

HR_se<-function (sample.list, parameter) 
{
  call = (parameter[[1]] == "Description")
  if (call == FALSE | is.na(call)) {
    if (length(sample.list) != 2) 
      stop("Analysis model: Two samples must be specified in the HazardRatioCoxStat statistic.")
    outcome1 = sample.list[[1]][, "outcome"]
    outcome1.complete = outcome1[stats::complete.cases(outcome1)]
    event1 = !sample.list[[1]][, "patient.censor.indicator"]
    event1.complete = event1[stats::complete.cases(outcome1)]
    n1 = length(outcome1.complete)
    outcome2 = sample.list[[2]][, "outcome"]
    outcome2.complete = outcome2[stats::complete.cases(outcome2)]
    event2 = !sample.list[[2]][, "patient.censor.indicator"]
    event2.complete = event2[stats::complete.cases(outcome2)]
    n2 = length(outcome2.complete)
    outcome = c(outcome1.complete, outcome2.complete)
    event = c(event1.complete, event2.complete)
    treatment = c(rep(0, n1), rep(1, n2))
    result = summary(survival::coxph(survival::Surv(outcome, 
                                                    event) ~ treatment))$coef[, "se(coef)"]
  }
  else if (call == TRUE) {
    result = list("Hazard Ratio (Cox)")
  }
  return(result)
}






