library(Mediana)

#Load required functions to make this thing work:
#LogrankTest_num: Corrected version of the log-rank test to allow large sample sizes
#rexp_asymp: Random samples following an asymptotic logistic distribution
#asymp_exp: Cumulative Asymptotic logistic distribution
#HR_se: Custom function to extract standard errors from a cox model in the simulation

source("IFIH1_support_functions.R")

#Model parameters

a_Cx<-0.283 #Asymptotic mortality in Cx variants without steroids (25.7% at day 28 * 1.1)
RR_dexa_Cx<-0.7 #RR of dexamethasone in patients with a Cx variant
RR_TTCx_no_dexa<-0.821 #RR of a TT variant compared to Cx, with no steroids
RR_dexa_TT<-1.487 ##RR of dexamethasone in patients with a TT variant
a_Cx_dexa<-a_Cx*RR_dexa_Cx #Mortality in Cx variants with steroids
a_TT<-a_Cx*RR_TTCx_no_dexa #Mortality in TT variants without steroids
a_TT_dexa<-a_TT*RR_dexa_TT #Mortality in TT variants with steroids

c<-0.08564 #Curvature parameter, assuming asymptotic mortality is 1.1 time higher than the 28-day mortality

# Outcome parameter set 1
death.placebo.Cx = parameters(asymptote=a_Cx, hr=c)
death.placebo.TT = parameters(asymptote=a_TT, hr=c)
death.dexa.Cx = parameters(asymptote=a_Cx_dexa, hr=c)
death.dexa.TT = parameters(asymptote=a_TT_dexa, hr=c)

# Sample size parameters
sample.size.total = c(6000)
TT_prop<-seq(0.13, 0.61, 0.16)^2
sample.size.placebo.Cx = as.list(as.integer((1-TT_prop)/2 * sample.size.total))
sample.size.placebo.TT = as.list(as.integer(TT_prop/2 * sample.size.total))
sample.size.dexa.Cx = as.list(as.integer((1-TT_prop)/2 * sample.size.total))
sample.size.dexa.TT = as.list(as.integer(TT_prop/2 * sample.size.total))

#Data model
IFIH1.data.model = DataModel() +
  OutcomeDist(outcome.dist = "asymp_exp", outcome.type = "event") +
  Design(enroll.period = 365,            
         followup.period = 28,
         enroll.dist = "UniformDist",
         dropout.dist = "UniformDist",
         dropout.dist.par = parameters(prop=0.00001) )+
  Sample(id = "Placebo Cx",
         sample.size = sample.size.placebo.Cx,
         outcome.par = parameters(death.placebo.Cx)) +
  Sample(id = "Placebo TT",
         sample.size = sample.size.placebo.TT,
         outcome.par = parameters(death.placebo.TT)) +
  Sample(id = "Dexamethasone Cx",
         sample.size = sample.size.dexa.Cx,
         outcome.par = parameters(death.dexa.Cx)) +
  Sample(id = "Dexamethasone TT",
         sample.size = sample.size.dexa.TT,
         outcome.par = parameters(death.dexa.TT))

# Analysis model
IFIH1.analysis.model = AnalysisModel() +
  Test(id = "Overall population",
       samples = samples(c("Placebo Cx", "Placebo TT"),
                         c("Dexamethasone Cx", "Dexamethasone TT")),
       method = "LogrankTest_num") +
  Statistic(id="HR Overall population", 
            samples=samples(c("Placebo Cx", "Placebo TT"),
                            c("Dexamethasone Cx", "Dexamethasone TT")),
            method = "HazardRatioStat",
            par = parameters(method = "Cox")) +
  Statistic(id="HR Overall population_SE", 
            samples=samples(c("Placebo Cx", "Placebo TT"),
                            c("Dexamethasone Cx", "Dexamethasone TT")),
            method = "HR_se") +
  Test(id = "Cx population",
       samples = samples("Placebo Cx", "Dexamethasone Cx"),
       method = "LogrankTest_num") +
  Statistic(id="HR Cx population", 
            samples=samples("Placebo Cx", "Dexamethasone Cx"),
            method = "HazardRatioStat",
            par = parameters(method = "Cox")) +
  Statistic(id="HR Cx population_SE", 
            samples=samples("Placebo Cx", "Dexamethasone Cx"),
            method = "HR_se") +
  Test(id = "TT population",
       samples = samples("Placebo TT", "Dexamethasone TT"),
       method = "LogrankTest_num") +
  Statistic(id="HR TT population", 
            samples=samples("Placebo TT", "Dexamethasone TT"),
            method = "HazardRatioStat",
            par = parameters(method = "Cox")) +
  Statistic(id="HR TT population_SE", 
            samples=samples("Placebo TT", "Dexamethasone TT"),
            method = "HR_se") +
  Test(id = "No steroids population",
       samples = samples("Placebo Cx", "Placebo TT"),
       method = "LogrankTest_num")+
  Statistic(id="HR No steroids population", 
            samples=samples("Placebo Cx", "Placebo TT"),
            method = "HazardRatioStat",
            par = parameters(method = "Cox"))+
    Statistic(id="HR No steroids population_SE", 
              samples=samples("Placebo Cx", "Placebo TT"),
              method = "HR_se") +
  Statistic(id="Log-rank overall",
            samples=samples(c("Placebo Cx", "Placebo TT"),
                            c("Dexamethasone Cx", "Dexamethasone TT")),
            method="LogrankTest_num",
            par=parameters(larger=TRUE))+
  Statistic(id="Log-rank Cx",
            samples = samples("Placebo Cx", "Dexamethasone Cx"),
            method = "LogrankTest_num")+
  Statistic(id = "Log-rank TT",
       samples = samples("Placebo TT", "Dexamethasone TT"),
       method = "LogrankTest_num",
       par=parameters(larger=FALSE)) +
  Statistic(id = "Log-rank No steroids",
       samples = samples("Placebo Cx", "Placebo TT"),
       method = "LogrankTest_num")

# Evaluation model
IFIH1.evaluation.model = EvaluationModel() +
  Criterion(id = "Marginal power",
            method = "MarginalPower",
            tests = tests("Overall population",
                          "Cx population",
                          "TT population",
                          "No steroids population"),
            labels = c("Overall population test",
                       "Cx population test", 
                       "TT population test", 
                       "No steroids population test"),
            par = parameters(alpha = 0.025)) +
  Criterion(id = "Hazard Ratio",
            method = "MeanSumm",
            statistics = tests ("HR Overall population",
                                "HR Overall population_SE",
                                "HR Cx population",
                                "HR Cx population_SE",
                                "HR TT population",
                                "HR TT population_SE",
                                "HR No steroids population",
                                "HR No steroids population_SE"),
            labels=c("HR Overall population",
                     "HR Overall population_SE",
                     "HR Cx population",
                     "HR Cx population_SE",
                     "HR TT population",
                     "HR TT population_SE",
                     "HR No steroids population",
                     "HR No steroids population_SE"))+
  Criterion(id =" Log-rank test", 
            method="MeanSumm",
            statistics=tests ("Log-rank overall",
                              "Log-rank Cx", 
                              "Log-rank TT",
                              "Log-rank No steroids"),
            labels=c("Log-rank p Overall population",
                     "Log-rank p Cx population", 
                     "Log-rank p TT population",
                     "Log-rank p No steroids population"))

# Simulation Parameters
IFIH1.sim.parameters =  SimParameters(n.sims = 1000,
                                            proc.load = "full",
                                            seed = 1990760)

# Perform clinical scenario evaluation
IFIH1.results = CSE(IFIH1.data.model,
                          IFIH1.analysis.model,
                          IFIH1.evaluation.model,
                          IFIH1.sim.parameters)

summary(IFIH1.results)

#Plotting hazard ratios
library(ggplot2)
library(dplyr)
library(tidyr)

trial_results<-summary(IFIH1.results)
trial_results$MAF<-rep(sqrt(TT_prop), each=16)

HR_table<-trial_results[trial_results$criterion=="Hazard Ratio",]
HR_table$SE<-rep(HR_table[grepl("_SE", HR_table$test.statistic), "result"], each=2)
HR_table<-HR_table[!grepl("_SE", HR_table$test.statistic),]
HR_table$lci<-HR_table$result/exp(2*HR_table$SE)
HR_table$uci<-HR_table$result*exp(2*HR_table$SE)

ggplot(HR_table[HR_table$test.statistic=="HR No steroids population",], aes(x=MAF, y=result, ymin=lci, ymax=uci, col=test.statistic))+
  geom_hline(yintercept=1, linetype=2)+
  geom_pointrange(position=position_dodge2(width=0.05), size=1, col="black")+
  scale_x_continuous(breaks=sqrt(TT_prop))+
  labs(y="Hazard Ratio (TT vs CC/CT)", x="Minor allele frequency")+
  theme_light(base_size = 24)+
  theme(legend.position="none", aspect.ratio=1.618)
  
HR_table$test.statistic<-factor(HR_table$test.statistic, levels=c("HR No steroids population", "HR Cx population", "HR TT population", "HR Overall population"))
ggplot(HR_table[HR_table$test.statistic!="HR No steroids population",], aes(x=MAF, y=result, ymin=lci, ymax=uci, col=test.statistic))+
  geom_hline(yintercept=1, linetype=2)+
  geom_pointrange(position=position_dodge2(width=0.05), size=1)+
  scale_x_continuous(breaks=sqrt(TT_prop))+
  labs(y="Hazard Ratio (Dexamethasone vs std care)", x="Minor allele frequency")+
  theme_light(base_size = 24)+
  theme(legend.position="none", aspect.ratio=1/1.618)+
  scale_color_manual(values=c("#BC3C29FF","#0072B5FF", "black"))+
  scale_y_log10()

#Plotting survival curves (equal sample sizes, all groups)
set.seed(1990760)

n<-1000

Cx<-rexp_asymp(n, a_Cx, c)
Cx_st<-rexp_asymp(n, a_Cx_dexa, c)
TT<-rexp_asymp(n, a_TT, c)
TT_st<-rexp_asymp(n, a_TT_dexa, c)
trial<-data.frame(group=rep(c("Cx", "Cx_st", "TT", "TT_st"), each=n), fu_time=c(Cx, Cx_st, TT, TT_st))
trial$status<-1

library(survival)
library(ggfortify)
surv.obj<-Surv(trial$fu_time, trial$status)
icu.model<-survfit(surv.obj~group, data=trial)
autoplot(icu.model, fun="event", conf.int = FALSE, surv.size = 1, censor=FALSE)+
  theme_light(base_size = 24)+
  theme(aspect.ratio = 1/1.618, legend.box.background = element_rect(), legend.position = "none")+
  scale_color_manual(values=c("#BC3C29FF", "#BC3C29FF","#0072B5FF", "#0072B5FF"))+
  labs(x="Time (days)", y="Hospital deaths (%)")+
  coord_cartesian(xlim=c(0,28))
survdiff(surv.obj~group, data=trial)

#Plotting survival curves in different MAF scenarios

# Sample size parameters
sample.size.total = c(6000)
TT_prop<-c(0.13^2, 0.4^2, 0.5^2, 0.6^2)

sample.size.placebo.Cx = as.list(as.integer((1-TT_prop)/2 * sample.size.total))
sample.size.placebo.TT = as.list(as.integer(TT_prop/2 * sample.size.total))
sample.size.dexa.Cx = as.list(as.integer((1-TT_prop)/2 * sample.size.total))
sample.size.dexa.TT = as.list(as.integer(TT_prop/2 * sample.size.total))

Cx<-lapply(sample.size.placebo.Cx, FUN=function(x) cbind(rexp_asymp(x, a_Cx, c), "Placebo", "Cx"))
Cx_st<-lapply(sample.size.dexa.Cx, FUN=function(x) cbind(rexp_asymp(x, a_Cx_dexa, c), "Dexamethasone", "Cx"))
TT<-lapply(sample.size.placebo.TT, FUN=function(x) cbind(rexp_asymp(x, a_TT, c), "Placebo", "TT"))
TT_st<-lapply(sample.size.dexa.TT, FUN=function(x) cbind(rexp_asymp(x, a_TT_dexa, c), "Dexamethasone", "TT"))

trial_steroids_1<-as.data.frame(rbind(Cx[[1]], Cx_st[[1]], TT[[1]], TT_st[[1]]))
trial_steroids_2<-as.data.frame(rbind(Cx[[2]], Cx_st[[2]], TT[[2]], TT_st[[2]]))
trial_steroids_3<-as.data.frame(rbind(Cx[[3]], Cx_st[[3]], TT[[3]], TT_st[[3]]))
trial_steroids_4<-as.data.frame(rbind(Cx[[4]], Cx_st[[4]], TT[[4]], TT_st[[4]]))
trial_steroids_1$MAF<-TT_prop[[1]]
trial_steroids_2$MAF<-TT_prop[[2]]
trial_steroids_3$MAF<-TT_prop[[3]]
trial_steroids_4$MAF<-TT_prop[[4]]
trial_steroids<-rbind(trial_steroids_1, trial_steroids_2, trial_steroids_3, trial_steroids_4)
trial_steroids$status<-1
colnames(trial_steroids)<-c("fu_time", "group", "variant", "MAF", "status")
trial_steroids$MAF<-as.factor(trial_steroids$MAF)

chunk_surv<-function(x, y) {
  surv.obj<-Surv(as.numeric(x$fu_time), x$status)
  icu.model<-survfit(surv.obj~x[[y]], data=x)
  sdf<-survdiff(surv.obj~x[[y]], data=x)
  p.val <- round(1 - pchisq(sdf$chisq, length(sdf$n) - 1), 3)
  p.val<-ifelse(p.val<0.001, "<0.001", paste("=",p.val, sep=""))
  dibujo<-autoplot(icu.model, fun="event", conf.int = FALSE, surv.size = 1, surv.linetype = y)+
    annotate("text", x=0, y=0.9, hjust=0,size=6, label=paste("Log-rank p", p.val, sep=""))+
    coord_cartesian(xlim=c(0,28))+
    theme_light(base_size = 24)+
    theme(aspect.ratio = 1/1.618, legend.box.background = element_rect(), legend.position = "none")+
    scale_linetype_manual(values=c(2,1))+
    scale_color_manual(values=c("#BC3C29FF", "#BC3C29FF","#0072B5FF", "#0072B5FF"))+
    labs(x="Time (days)", y="Hospital deaths")
  return(dibujo)
}

effects_steroids<-by(trial_steroids, trial_steroids$MAF, FUN=function(x) chunk_surv(x, "group"))
effects_steroids[[4]]
ggsave("sim_surv_21_4.pdf")

# Get the data generated in the CSE (beware!, large object)
IFIH1.data.stack = DataStack(data.model = IFIH1.data.model,
                                   sim.parameters = IFIH1.sim.parameters)
