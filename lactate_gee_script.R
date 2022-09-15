
#### ANALYSIS SCRIPT FOR:
## PAPER: Serum lactate is associated with increased illness severity in immunocompromised pediatric hematology oncology patients presenting to the emergency department with fever
## JOURNAL: Frontiers in Oncology
## AUTHOR(S): Leonora Rose Slatnick, Kristen Miller, Halden F Scott, Michele Loi, Adam J Esbenshade, Anna Franklin and Alisa B Lee-Sherick
## BIOSTATISTICIAN: Kristen Miller, University of Colorado 

library(gee)
library(MuMIn)
# dataset: not publically available, but format is as follows:
# long format (1 row per ED visit - multiple rows per patient if repeat visits)
# study id variable: id_n
# outcomes: CDE.within.48.hrs_01 (CDE within 48 hours, 0=no, 1=yes), IBI_yn_48hr_01 (IBI within 48 hours, 0=no, 1=yes)
# covariaets: listed below in "covariates_all" object
dat<-read.csv('K:/CCBD/Slatnick/Lactate/Data/lactatemanuscript_dataset_6.9.2022.csv')

#list of covariates for univariate analyses:
covariates_all<-c("Lactate....initial.within.2.hours","lactate_2","lactate_2_4",
                  "Age..years.","Sex","leukemia","lymphoma","cns","solid",
                  "BMT.last.6.months","Type.of.line","maxtemp",
                  "Hypotension_ED.2hrs","tachycardic_goldstein",
                  "tachypneic_goldstein","Chills.or.rigors","URI.symptoms",
                  "WBC...initial","ANC...initial",
                  "Neutropenic","Initial.AMC","Initial.ALC",
                  "Plt...initial_10","Hemoglobin...initial",
                  "intensity_grp")

#GEE univariate model function: fits the GEE for a given covariate and outcome, formats output into table
logit_gee<- function(var1,outcome) {

  mod<-gee(outcome ~ var1,data = dat, id = id_n, family = binomial,
           corstr = "exchangeable")
  mod_sum<-data.frame(summary(mod)$coefficients)
  mod_sum$p<-2 * pnorm(abs(mod_sum$Robust.z), lower.tail = FALSE)
  mod_sum$Var<-row.names(mod_sum)
  mod_sum<-mod_sum[-1,-c(2,3,5)]
  p<-round(mod_sum$p,4)
  p[p<0.001]<-"<0.001"
    #table
  if (is.factor(var1)){
    tab.1<-data.frame(Var=paste0(label(var1),": ",levels(as.factor(var1))[-1]),
                      ref=levels(as.factor(var1))[1],
                      OR=paste0(round(exp(mod_sum[,1]),2)," 95% CI: (",round(exp(mod_sum[,1]-1.96*mod_sum[,2]),2),
                                ", ",round(exp(mod_sum[,1]+1.96*mod_sum[,2]),2),")"),
                      pval=p,
                      nobs=mod$nobs,
                      row.names=NULL)
  }
  if (!is.factor(var1)){
    tab.1<-data.frame(Var=label(var1),
                      ref=NA,
                      OR=paste0(round(exp(mod_sum[,1]),2)," 95% CI: (",round(exp(mod_sum[,1]-1.96*mod_sum[,2]),2),
                                ", ",round(exp(mod_sum[,1]+1.96*mod_sum[,2]),2),")"),
                      pval=p,
                      nobs=mod$nobs,
                      row.names=NULL)
  }
  return(tab.1)
}

#call GEE function to produce univariate model results for both outcomes
CDE_uni_gee<-data.frame(do.call(rbind,lapply(dat[covariates_all],logit_gee,outcome=dat$CDE.within.48.hrs_01)),row.names=NULL)
BACT_uni_gee<-data.frame(do.call(rbind,lapply(dat[covariates_all],logit_gee,outcome=dat$IBI_yn_48hr_01)),row.names=NULL)

# MULTIVARIABLE MODEL BUILDING: for CDE
#1. fit full model with all potential predictors
CDE_full<-gee(CDE.within.48.hrs_01 ~
                Lactate....initial.within.2.hours+
                Neutropenic+
                Hypotension_ED.2hrs+
                tachycardic_goldstein+
                Chills.or.rigors+
                Age..years.+
                intensity_grp,
              data = dat, id = id_n, family = binomial,corstr = "exchangeable")

#2. select subset of covariates that result in model with lowest QIC
CDE_select<-dredge(CDE_full,beta="none",rank="QIC",subset=~Lactate....initial.within.2.hours)
#3. fit model with best subset of covariates for reporting
CDE_final<-gee(CDE.within.48.hrs_01 ~
                 Lactate....initial.within.2.hours+
                 Age..years.+
                 Hypotension_ED.2hrs+
                 intensity_grp_new+
                 tachycardic_goldstein
               
               ,
               data = dat, id = id_n, family = binomial,corstr = "exchangeable")
QIC(CDE_final)
#4. compare this model with the same model but without lactate
CDE_final_NOLACTATE<-gee(CDE.within.48.hrs_01 ~
                           Age..years.+
                           Hypotension_ED.2hrs+
                           intensity_grp+
                           tachycardic_goldstein,
                         data = dat, id = id_n, family = binomial,corstr = "exchangeable")
QIC(CDE_final_NOLACTATE)

# MULTIVARIABLE MODEL BUILDING: for IBI
#1. fit full model with all potential predictors
INF_full<-gee(IBI_yn_48hr_01 ~
                Lactate....initial.within.2.hours+
                Neutropenic+
                Chills.or.rigors+
                intensity_grp+
                Hypotension_ED.2hrs+
                tachycardic_goldstein+
                Type.of.line,
              data = dat, id = id_n, family = binomial,corstr = "exchangeable")
#2. select subset of covariates that result in model with lowest QIC
INF_select<-dredge(INF_full,beta="none",rank="QIC",subset=~Lactate....initial.within.2.hours)
#3. fit model with best subset of covariates for reporting
INF_final<-gee(IBI_yn_48hr_01 ~
                 Lactate....initial.within.2.hours+
                 Chills.or.rigors+
                 Neutropenic+
                 Type.of.line,
               data = dat, id = id_n, family = binomial,corstr = "exchangeable")
QIC(INF_final)
#4. compare this model with the same model but without lactate
INF_final_NOLACTATE<-gee(IBI_yn_48hr_01 ~
                           Chills.or.rigors+
                           Neutropenic+
                           Type.of.line,
                         data = dat, id = id_n, family = binomial,corstr = "exchangeable")
QIC(INF_final_NOLACTATE)