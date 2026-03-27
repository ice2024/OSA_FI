library(haven)
library(plyr)
library(dplyr) 
library(arsenal) 
library(survey)
library(gtsummary)
library(survey)
library(haven)
library(tableone)
library(plyr)
library(dplyr) 
library(tidyverse)
library(arsenal) 
demo.a <- read_xpt("DEMO/demo.xpt")
demo.b <- read_xpt("DEMO/demo_b.xpt")
demo.c <- read_xpt("DEMO/demo_c.xpt")
demo.d <- read_xpt("DEMO/demo_d.xpt")
demo.e <- read_xpt("DEMO/demo_e.xpt")
demo.f <- read_xpt("DEMO/demo_f.xpt")
demo.g <- read_xpt("DEMO/demo_g.xpt")
demo.h <- read_xpt("DEMO/demo_h.xpt")
demo.i <- read_xpt("DEMO/demo_i.xpt")
demo.j <- read_xpt("DEMO/DEMO_J.XPT")
demo.data.file <- dplyr::bind_rows(list(demo.a,demo.b,demo.c,demo.d,demo.e,demo.f,demo.g, demo.h,demo.i,demo.j))
demo.data <- demo.data.file[,c('SEQN', 'RIDAGEYR', 'RIAGENDR', 'RIDRETH1', 'DMDEDUC2','DMDCITZN','DMDMARTL','INDFMPIR')]
tab1 <- tableby( ~ RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) + factor(DMDEDUC2) + factor(DMDCITZN)+factor(DMDMARTL)+INDFMPIR, 
smq.a <- read_xpt("SMQ/smq.xpt")
smq.b <- read_xpt("SMQ/smq_b.xpt")
smq.c <- read_xpt("SMQ/smq_c.xpt")
smq.d <- read_xpt("SMQ/smq_d.xpt")
smq.e <- read_xpt("SMQ/smq_e.xpt")
smq.f <- read_xpt("SMQ/smq_f.xpt")
smq.g <- read_xpt("SMQ/smq_g.xpt")
smq.h <- read_xpt("SMQ/smq_h.xpt")
smq.i <- read_xpt("SMQ/smq_i.xpt")
smq.j <- read_xpt("SMQ/SMQ_J.XPT")

smq.data.file <- dplyr::bind_rows(list(smq.a,smq.b,smq.c,smq.d,smq.e,smq.f,smq.g, smq.h,smq.i,smq.j))
colnames(smq.data.file)
smq.data <- smq.data.file[,c('SEQN', 'SMQ020', 'SMQ040', 'SMQ050Q', 'SMQ050U')]
alq.a <- read_xpt("ALQ/alq.xpt")
colnames(alq.a)
colnames(alq.a)[2]="ALQ101"
alq.b <- read_xpt("ALQ/alq_b.xpt")
colnames(alq.b)
colnames(alq.b)[2]="ALQ101"

alq.c <- read_xpt("ALQ/alq_c.xpt")
alq.d <- read_xpt("ALQ/alq_d.xpt")
alq.e <- read_xpt("ALQ/alq_e.xpt")
alq.f <- read_xpt("ALQ/alq_f.xpt")
alq.g <- read_xpt("ALQ/alq_g.xpt")
alq.h <- read_xpt("ALQ/alq_h.xpt")
alq.i <- read_xpt("ALQ/alq_i.xpt")
alq.j <- read_xpt("ALQ/ALQ_J.XPT")
colnames(alq.j)[2]="ALQ101"

alq.data.file <- dplyr::bind_rows(list(alq.a,alq.b,alq.c,alq.d,alq.e,alq.f,alq.g, alq.h,alq.i,alq.j))
colnames(alq.data.file)
alq.data <- alq.data.file[,c('SEQN', 'ALQ101', 'ALQ110', 'ALQ120Q', 'ALQ120U')]
tab1 <- tableby(~factor(ALQ101) + factor(ALQ110) + ALQ120Q + factor(ALQ120U), 
                data = alq.data)
bmx.a <- read_xpt("BMX/bmx.xpt")
bmx.b <- read_xpt("BMX/bmx_b.xpt")
bmx.c <- read_xpt("BMX/bmx_c.xpt")
bmx.d <- read_xpt("BMX/bmx_d.xpt")
bmx.e <- read_xpt("BMX/bmx_e.xpt")
bmx.f <- read_xpt("BMX/bmx_f.xpt")
bmx.g <- read_xpt("BMX/bmx_g.xpt")
bmx.h <- read_xpt("BMX/bmx_h.xpt")
bmx.i <- read_xpt("BMX/bmx_i.xpt")
bmx.j <- read_xpt("BMX/BMX_J.XPT")
bmx.data.file <- dplyr::bind_rows(list(bmx.a,bmx.b,bmx.c,bmx.d,bmx.e,bmx.f,bmx.g, bmx.h,bmx.i,bmx.j))
bmx.data <- bmx.data.file[,c('SEQN', 'BMXBMI')]
dr1tot.a <- read_xpt("DR1TOT/DRXTOT_1999_2000.XPT")
dr1tot.a=dr1tot.a[,c("SEQN","DRXTALCO","WTDRD1")]
colnames(dr1tot.a)[2]="DR1TALCO"
dr1tot.b <- read_xpt('DR1TOT/DRXTOT_B_2001_2002.XPT')
dr1tot.b=dr1tot.b[,c("SEQN","DRXTALCO","WTDRD1")]
colnames(dr1tot.b)[2]="DR1TALCO"
dr1tot.c <- read_xpt('DR1TOT/DR1TOT_C.XPT')
dr1tot.d <- read_xpt('DR1TOT/DR1TOT_D.XPT')
dr1tot.e <- read_xpt('DR1TOT/DR1TOT_E.XPT')
dr1tot.f <- read_xpt('DR1TOT/DR1TOT_F.XPT')
dr1tot.g <- read_xpt('DR1TOT/DR1TOT_G.XPT')
dr1tot.h <- read_xpt('DR1TOT/DR1TOT_H.XPT')
dr1tot.i <- read_xpt('DR1TOT/DR1TOT_I.XPT')
dr1tot.j <- read_xpt('DR1TOT/DR1TOT_J.XPT')
dr1tot.data.file <- dplyr::bind_rows(list(dr1tot.a,dr1tot.b,dr1tot.c,dr1tot.d,dr1tot.e,dr1tot.f,dr1tot.g, dr1tot.h,dr1tot.i,dr1tot.j))
dr1tot.data <- dr1tot.data.file[,c('SEQN', 'DR1TALCO')]

# 第2天
dr2tot.a <- read_xpt("DR2TOT/DRXTOT_1999_2000.XPT")
dr2tot.a=dr2tot.a[,c("SEQN","DRXTALCO","WTDRD1")]
colnames(dr2tot.a)[2]="DR2TALCO"
dr2tot.b <- read_xpt('DR2TOT/DRXTOT_B_2001_2002.XPT')
dr2tot.b=dr2tot.b[,c("SEQN","DRXTALCO","WTDRD1")]
colnames(dr2tot.b)[2]="DR2TALCO"

dr2tot.c <- read_xpt('DR2TOT/DR2TOT_C.XPT')
dr2tot.d <- read_xpt('DR2TOT/DR2TOT_D.XPT')
dr2tot.e <- read_xpt('DR2TOT/DR2TOT_E.XPT')
dr2tot.f <- read_xpt('DR2TOT/DR2TOT_F.XPT')
dr2tot.g <- read_xpt('DR2TOT/DR2TOT_G.XPT')
dr2tot.h <- read_xpt('DR2TOT/DR2TOT_H.XPT')
dr2tot.i <- read_xpt('DR2TOT/DR2TOT_I.XPT')
dr2tot.j <- read_xpt('DR2TOT/DR2TOT_J.XPT')
dr2tot.data.file <- dplyr::bind_rows(list(dr2tot.a,dr2tot.b,dr2tot.c,dr2tot.d,dr2tot.e,dr2tot.f,dr2tot.g, dr2tot.h,dr2tot.i, dr2tot.j))
dr2tot.data <- dr2tot.data.file[,c('SEQN', 'DR2TALCO')]

dr.data <- merge(dr2tot.data, dr1tot.data)
View(dr.data)
mcq.a <- read_xpt('MCQ/mcq.xpt')
mcq.a$tmp=paste(mcq.a$MCQ160A,mcq.a$MCQ190,sep = "_")
s=which(mcq.a$tmp=="1_1"|mcq.a$tmp=="2_NA")
mcq.a=mcq.a[s,]
mcq.a=mcq.a[,c('SEQN', 'MCQ160A')]##1
table(mcq.a$MCQ160A)

mcq.b <- read_xpt('MCQ/mcq_b.xpt')
mcq.b$tmp=paste(mcq.b$MCQ160A,mcq.b$MCQ190,sep = "_")
s=which(mcq.b$tmp=="1_1"|mcq.b$tmp=="2_NA")
mcq.b=mcq.b[s,]
mcq.b=mcq.b[,c('SEQN', 'MCQ160A')]

mcq.c <- read_xpt('MCQ/mcq_c.xpt')
mcq.c$tmp=paste(mcq.c$MCQ160A,mcq.c$MCQ190,sep = "_")
table(mcq.c$tmp)
s=which(mcq.c$tmp=="1_1"|mcq.c$tmp=="2_NA")
mcq.c=mcq.c[s,]
mcq.c=mcq.c[,c('SEQN', 'MCQ160A')]

mcq.d <- read_xpt('MCQ/mcq_d.xpt')
mcq.d$tmp=paste(mcq.d$MCQ160A,mcq.d$MCQ190,sep = "_")
s=which(mcq.d$tmp=="1_1"|mcq.d$tmp=="2_NA")
mcq.d=mcq.d[s,]
mcq.d=mcq.d[,c('SEQN', 'MCQ160A')]

mcq.e <- read_xpt('MCQ/mcq_e.xpt')
mcq.e$tmp=paste(mcq.e$MCQ160A,mcq.e$MCQ190,sep = "_")
s=which(mcq.e$tmp=="1_1"|mcq.e$tmp=="2_NA")
mcq.e=mcq.e[s,]
mcq.e=mcq.e[,c('SEQN', 'MCQ160A')]

mcq.f <- read_xpt('MCQ/mcq_f.xpt')
mcq.f$tmp=paste(mcq.f$MCQ160A,mcq.f$MCQ191,sep = "_")
table(mcq.f$tmp)
s=which(mcq.f$tmp=="1_1"|mcq.f$tmp=="2_NA")
mcq.f=mcq.f[s,]
mcq.f=mcq.f[,c('SEQN', 'MCQ160A')]
mcq.g <- read_xpt('MCQ/mcq_g.xpt')
mcq.g$tmp=paste(mcq.g$MCQ160A,mcq.g$MCQ195,sep = "_")
table(mcq.g$tmp)
s=which(mcq.g$tmp=="1_2"|mcq.g$tmp=="2_NA")
mcq.g=mcq.g[s,]
mcq.g=mcq.g[,c('SEQN', 'MCQ160A')]

mcq.h <- read_xpt('MCQ/mcq_h.xpt')
mcq.h$tmp=paste(mcq.h$MCQ160A,mcq.h$MCQ195,sep = "_")
table(mcq.h$tmp)
s=which(mcq.h$tmp=="1_2"|mcq.h$tmp=="2_NA")
mcq.h=mcq.h[s,]
mcq.h=mcq.h[,c('SEQN', 'MCQ160A')]

mcq.i <- read_xpt('MCQ/mcq_i.xpt')
mcq.i$tmp=paste(mcq.i$MCQ160A,mcq.i$MCQ195,sep = "_")
table(mcq.i$tmp)
s=which(mcq.i$tmp=="1_2"|mcq.i$tmp=="2_NA")
mcq.i=mcq.i[s,]
mcq.i=mcq.i[,c('SEQN', 'MCQ160A')]

mcq.j <- read_xpt('MCQ/MCQ_J.XPT')
mcq.j$tmp=paste(mcq.j$MCQ160A,mcq.j$MCQ195,sep = "_")
table(mcq.j$tmp)
s=which(mcq.j$tmp=="1_2"|mcq.j$tmp=="2_NA")
mcq.j=mcq.j[s,]
mcq.j=mcq.j[,c('SEQN', 'MCQ160A')]

mcq.data.file <- dplyr::bind_rows(list(mcq.a,mcq.b,mcq.c,mcq.d,mcq.e,mcq.f,mcq.g, mcq.h,mcq.i,mcq.j))
mcq.data <- mcq.data.file[,c('SEQN', 'MCQ160A')]
table(mcq.data$MCQ160A)

weight.data <- dr1tot.data.file[,c('SEQN', 'WTDRD1')]
weight.data$WTDRD1 <- weight.data$WTDRD1/10
survey.design.data <- demo.data.file[,c('SEQN', 'SDMVPSU', 'SDMVSTRA')]
slq.d <- read_xpt("SLQ/SLQ_D.xpt")
slq.e <- read_xpt("SLQ/SLQ_E.xpt")
slq.f <- read_xpt("SLQ/SLQ_F.xpt")
slq.g <- read_xpt("SLQ/SLQ_G.xpt")
slq.h <- read_xpt("SLQ/SLQ_H.xpt")
slq.i <- read_xpt("SLQ/SLQ_I.xpt")
slq.j <- read_xpt("SLQ/P_SLQ.xpt")
slq.data.file <- dplyr::bind_rows(list(slq.d,slq.e,slq.f,slq.g, slq.h,slq.i,slq.j))
colnames(slq.data.file)
slq.data <- slq.data.file[,c('SEQN', 'SLQ030', 'SLQ040','SLQ120')]

output <- plyr::join_all(list(demo.data, smq.data, alq.data, bmx.data,slq.data,
                              dr.data, mcq.data,weight.data,survey.design.data),
                         by='SEQN', type='left')
output1=output[-which(output$SEQN%in%demo.l$SEQN),]

dim(output1)
paper.data <- subset.data.frame(output1, RIDAGEYR >= 20 )
dim(paper.data)
table(paper.data$SLQ030)
table(paper.data$SLQ040)
table(paper.data$SLQ120)
colnames(paper.data)
s1=which(paper.data$SLQ030=="7"|paper.data$SLQ030=="9")
s2=which(paper.data$SLQ040=="7"|paper.data$SLQ040=="9")
s3=which(paper.data$SLQ120=="7"|paper.data$SLQ120=="9")
paper.data=paper.data[-c(s1,s2,s3),]
write.csv(paper.data,"paper.data.csv")
paper.data$Sex <- ifelse(paper.data$RIAGENDR == 1, 'male', 'female')
paper.data$Age=paper.data$RIDAGEYR
paper.data$Age.group <- ifelse(paper.data$RIDAGEYR >= 20 & paper.data$RIDAGEYR < 65, '20-65 years',
                               '65+ years')#)

table(paper.data$DMDMARTL)
Marital <- recode_factor(paper.data$DMDMARTL, 
                         `1` = 'Married',
                         `2` = 'Never married',
                         `3` = 'Never married',
                         `4` = 'Never married',
                         `5` = 'Never married',
                         `6` = 'Married'
                         
)
paper.data$Marital <- Marital
table(paper.data$Marital)

table(paper.data$ALQ101)
alq <- recode_factor(paper.data$ALQ101, 
                     `2` = 'drinker',
                     `1` = 'Non-drinker'                
)
paper.data$Alq.group <- alq

paper.data$BMI=paper.data$BMXBMI
paper.data$BMI.group <- ifelse(paper.data$BMXBMI <18.5, 'Underweight(<18.5)',
                               ifelse(paper.data$BMXBMI >=18.5 & paper.data$BMXBMI < 25, 'Normal(18.5 to <25)',
                                      ifelse(paper.data$BMXBMI >=25 & paper.data$BMXBMI < 30, 'Overweight(25 to <30)',
                                             'Obese(30 or greater)')))

paper.data <- mutate(paper.data, Smoke.group = case_when(
  SMQ020 == 2 ~ 'Never smoker',
  SMQ020 == 1 & SMQ040 == 3 ~ 'Former smoker',
  SMQ020 == 1 & SMQ040 <= 2 ~ 'Current smoker'
))

table(paper.data$DMDEDUC2)
education.attainment <- recode_factor(paper.data$DMDEDUC2, 
                                      `1` = 'Below high school',
                                      `2` = 'Below high school',
                                      `3`= 'High School',
                                      `4`= 'Above High school',
                                      `5`= 'Above High school'
)

paper.data$Education.attainment <- education.attainment
table(paper.data$SLQ030)
SLQ030.fenzu <- recode_factor(paper.data$SLQ030, 
                                      `0` = '0',
                                      `1` = '0',
                                      `2`= '0',                                      
                                      `3`= '1'
)

paper.data$SLQ030.group <- SLQ030.fenzu
table(paper.data$SLQ040)
SLQ040.fenzu <- recode_factor(paper.data$SLQ040, 
                              `0` = '0',
                              `1` = '0',
                              `2`= '0',
                              `3`= '1'
)
paper.data$SLQ040.group <- SLQ040.fenzu
table(paper.data$SLQ120)
SLQ120.fenzu <- recode_factor(paper.data$SLQ120, 
                              `0` = '0',
                              `1` = '0',
                              `2`= '0',
                              `3`= '0',
                              `4`= '1'
)
paper.data$SLQ120.group <- SLQ120.fenzu

paper.data <- mutate(paper.data, Type = case_when(
  SLQ030.group == 0 & SLQ040.group == 0 & SLQ120.group == 0 ~ 'Control',
  SLQ030.group == 1 | SLQ040.group == 1 | SLQ120.group == 1 ~ 'OSA'
))

colnames(paper.data)
table(paper.data$Type)
fi=read.csv("Frailty_index.csv",header = T,stringsAsFactors = F,check.names = F)
fi=fi[,c(2,53)]
paper.data1=merge(paper.data,fi,by="SEQN")
colnames(paper.data1)
colnames(paper.data1)[38]="Frailty.index"
paper.data2 <- subset.data.frame(paper.data1, 
                                (!is.na(Type)) & 
                                  (!is.na(Frailty.index)) &
                                  (!is.na(WTINT2YR))&
                                  (!is.na(Smoke.group))& 
                                  (!is.na(Marital))&
                                  (!is.na(Sex)) & 
                                  (!is.na(Education.attainment))&
                                  (!is.na(Age.group)) & 
                                  (!is.na(Alq.group)) &
                                  (!is.na(BMI.group))) 
dim(paper.data2)
table(paper.data2$Type)
colnames(paper.data2)
analyze.variable <- c("SEQN", "WTINT2YR", "SDMVPSU", "SDMVSTRA",
                      "Sex", "Age", "Age.group", 
                      "Alq.group","BMI", "BMI.group","Smoke.group","Education.attainment","Marital","Frailty.index","Type"
)
paper.data2 <- paper.data2[, analyze.variable]
write.csv(paper.data2,"paper.data_result.csv")

all_variable <- c("Sex", "Age", "Age.group", 
                  "Alq.group","BMI", "BMI.group","Smoke.group","Education.attainment","Marital","Frailty.index")
factor_variable <- c('Sex', 'Age.group',  
                     'Alq.group', 'BMI.group', 'Smoke.group','Education.attainment','Marital')

table_statistics <- CreateTableOne(vars = all_variable, 
                                   strata = 'Type', 
                                   data = paper.data2, 
                                   factorVars = factor_variable)
table_statistics <- print(table_statistics, 
                          showAllLevels = TRUE 
)
write.csv(table_statistics, '01.table_statistics.csv')
paper.data=read.csv("paper.data_result.csv",header = T,stringsAsFactors = FALSE,check.names = FALSE,row.names = 1)
paper.data1=paper.data
paper.data1$Frailty.index.var <- cut(paper.data1$Frailty.index,
                           breaks = quantile(paper.data1$Frailty.index),
                           labels = c('Q1', 'Q2', 'Q3', 'Q4'))
paper.data1$Frailty.index.var[which(is.na(paper.data1$Frailty.index.var))] <- 'Q1'

table(paper.data1$Frailty.index.var)

write.csv(paper.data1,"paper.data_result.Frailty.index.var.csv")
paper.data1=paper.data1[which(paper.data1$Type=="OSA"),]

all_variable <- c("Sex", "Age", "Age.group", 
                  "Alq.group","BMI", "BMI.group","Smoke.group",
                  "Education.attainment",
                  "Marital"
                  )
factor_variable <- c('Sex', 'Age.group',
                     'Alq.group', 'BMI.group', 'Smoke.group',
                     'Education.attainment',"Marital")
table_statistics <- CreateTableOne(vars = all_variable, 
                                   strata = 'Frailty.index.var', 
                                   data = paper.data1, 
                                   factorVars = factor_variable)
table_statistics <- print(table_statistics, 
                          showAllLevels = TRUE 
)
write.csv(table_statistics, '02.table_statistics.csv')

nhance_data <- read.csv('paper.data_result.Frailty.index.var.csv', row.names = 1, check.names = F)
nhance_data$Type=ifelse(nhance_data$Type == "OSA", '1', '0')
nhance_data$Type <- as.numeric(nhance_data$Type)
colnames(nhance_data)
exposure <- 'Frailty.index'
model_first <- c('Age.group','Sex')
model_second <- c('Alq.group','Smoke.group','BMI.group','Education.attainment','Marital')
library(survey)
stamoch_lead <- svydesign(ids = ~ SDMVPSU, 
                          strata = ~ SDMVSTRA,
                          weights = ~WTINT2YR, 
                          nest = TRUE, 
                          data = nhance_data)

formula_first <- as.formula(paste0 ("Type~",
                                    paste0(c(exposure),
                                           collapse = "+")))
options(survey.lonely.psu="adjust")
model_svyglm_first <- svyglm(formula_first, 
                             design = stamoch_lead,family = 'binomial')
first_model_df <- data.frame(ID = exposure, tidy(model_svyglm_first))

B_value <- coef(model_svyglm_first)[2] ## B值
OR_value <- exp(coef(model_svyglm_first))[2] ## OR值
OR_value
CI <- exp(confint(model_svyglm_first))[2, ] ## 95%CI
CI
pvalue <- tidy(model_svyglm_first)[2, ]$p.value ## p值

model1 <- data.frame(exposure = exposure, 
                     OR_model1 = paste0(format(sprintf("%.3f",OR_value), scientific = F, digits = 3, width = 5, trim = TRUE),
                                        '(',
                                        format(sprintf("%.3f",CI[1]), scientific = F, digits = 3, width = 5, trim = TRUE),
                                        '-',
                                        format(sprintf("%.3f",CI[2]),scientific = F, digits = 3, width = 5, trim = TRUE),
                                        ')'))
colnames(model1) <- c('exposure', 'model1_OR(95%_CI)')
p_model1 <- data.frame(exposure = 'p_value', 
                       OR_model1 =format(pvalue, scientific = TRUE, digits = 3, width = 5, trim = TRUE))
colnames(p_model1) <- c('exposure', 'model1_OR(95%_CI)')
finally_model1 <- rbind(model1, p_model1)
formula_second <- as.formula(paste0 ("Type~",
                                     paste0(c(exposure, model_first),
                                            collapse = "+")))
options(survey.lonely.psu="adjust")
model_svyglm_second <- svyglm(formula_second, 
                              design = stamoch_lead, family = 'binomial')

B_value <- coef(model_svyglm_second)[2] ## B值
OR_value <- exp(coef(model_svyglm_second))[2] ## OR值
CI <- exp(confint(model_svyglm_second))[2, ] ## 95%CI
pvalue <- tidy(model_svyglm_second)[2, ]$p.value ## p值

model2 <- data.frame(exposure = exposure, 
                     OR_model2 = paste0(format(sprintf("%.3f",OR_value), scientific = F, digits = 3, width = 5, trim = TRUE),
                                        '(',
                                        format(sprintf("%.3f",CI[1]), scientific = F, digits = 3, width = 5, trim = TRUE),
                                        '-',
                                        format(sprintf("%.3f",CI[2]), scientific = F, digits = 3, width = 5, trim = TRUE),
                                        ')'))
colnames(model2) <- c('exposure', 'model2_OR(95%_CI)')
p_model2 <- data.frame(exposure = 'p_value', 
                       OR_model2 =format(pvalue, scientific = TRUE, digits = 3, width = 5, trim = TRUE))
colnames(p_model2) <- c('exposure', 'model2_OR(95%_CI)')

finally_model2 <- rbind(model2, p_model2)
formula_third <- as.formula(paste0 ("Type~",
                                    paste0(c(exposure, model_first, model_second),
                                           collapse = "+")))
options(survey.lonely.psu="adjust")
model_svyglm_third <- svyglm(formula_third, 
                             design = stamoch_lead, family = 'binomial')

B_value <- coef(model_svyglm_third)[2] ## B值
OR_value <- exp(coef(model_svyglm_third))[2] ## OR值
CI <- exp(confint(model_svyglm_third))[2, ] ## 95%CI
pvalue <- tidy(model_svyglm_third)[2, ]$p.value ## p值

model3 <- data.frame(exposure = exposure, 
                     OR_model3 = paste0(format(sprintf("%.3f",OR_value), scientific = F, digits = 3, width = 5, trim = TRUE),
                                        '(',
                                        format(sprintf("%.3f",CI[1]), scientific = F, digits = 3, width = 5, trim = TRUE),
                                        '-',
                                        format(sprintf("%.3f",CI[2]), scientific = F, digits = 3, width = 5, trim = TRUE),
                                        ')'))
colnames(model3) <- c('exposure', 'model3_OR(95%_CI)')
p_model3 <- data.frame(exposure = 'p_value', 
                       OR_model3 =format(pvalue, scientific = TRUE, digits = 3, width = 5, trim = TRUE))
colnames(p_model3) <- c('exposure', 'model3_OR(95%_CI)')

finally_model3 <- rbind(model3, p_model3)

table_data <- full_join(finally_model1, finally_model2, by="exposure")
table_data <- full_join(table_data, finally_model3, by="exposure")

write.csv(table_data, '01.association_res.csv')

library(questionr)
model1_OR <- odds.ratio(model_svyglm_first)
model2_OR <- odds.ratio(model_svyglm_second)
model3_OR <- odds.ratio(model_svyglm_third)
model1_OR$`OR (95% CI)` <- ifelse(is.na(model1_OR$`2.5 %`), "",
                                  sprintf("%.3f (%.3f-%.3f)",
                                          model1_OR$OR, model1_OR$`2.5 %`, model1_OR$`97.5 %`))
model1_OR$`p.value` <- ifelse(is.na(model1_OR$`2.5 %`), "",ifelse(model1_OR$p<0.001, "<0.001",sprintf("%.3f", model1_OR$p)))

model2_OR$`OR (95% CI)` <- ifelse(is.na(model2_OR$`2.5 %`), "",
                                  sprintf("%.3f (%.3f-%.3f)",
                                          model2_OR$OR, model2_OR$`2.5 %`, model2_OR$`97.5 %`))
model2_OR$`p.value` <- ifelse(is.na(model2_OR$`2.5 %`), "",ifelse(model2_OR$p<0.001, "<0.001",sprintf("%.3f", model2_OR$p)))

model3_OR$`OR (95% CI)` <- ifelse(is.na(model3_OR$`2.5 %`), "",
                                  sprintf("%.3f (%.3f-%.3f)",
                                          model3_OR$OR, model3_OR$`2.5 %`, model3_OR$`97.5 %`))
model3_OR$`p.value` <- ifelse(is.na(model3_OR$`2.5 %`), "",ifelse(model3_OR$p<0.001, "<0.001",sprintf("%.3f", model3_OR$p)))
write.csv(model1_OR, '01.model1_res.csv')
write.csv(model2_OR, '02.model2_res.csv')
write.csv(model3_OR, '03.model3_res.csv')

colnames(nhance_data)

exposure <- 'Frailty.index.var'
model_first <- c('Age.group', 'Sex')
model_second <- c('Alq.group','Smoke.group','BMI.group','Education.attainment','Marital')

library(survey)
stamoch_lead <- svydesign(ids = ~ SDMVPSU, 
                          strata = ~ SDMVSTRA,
                          weights = ~WTINT2YR, 
                          nest = TRUE, 
                          # survey.lonely.psu="adjust",  #抽样单元为1时不报错
                          data = nhance_data)

## 1.1 模型1-----------------
formula_first <- as.formula(paste0 ("Type~",
                                    paste0(c(exposure),
                                           collapse = "+")))
options(survey.lonely.psu="adjust")
model_svyglm_first <- svyglm(formula_first, 
                             design = stamoch_lead,family = 'binomial')
first_model_df <- data.frame(ID = exposure, tidy(model_svyglm_first))

#write.csv(finally_model1,"model1.csv")
## 1.2 模型2-----------------
formula_second <- as.formula(paste0 ("Type~",
                                     paste0(c(exposure, model_first),
                                            collapse = "+")))
options(survey.lonely.psu="adjust")
model_svyglm_second <- svyglm(formula_second, 
                              design = stamoch_lead, family = 'binomial')


## 1.3 模型3-----------------
formula_third <- as.formula(paste0 ("Type~",
                                    paste0(c(exposure, model_first, model_second),
                                           collapse = "+")))
options(survey.lonely.psu="adjust")
model_svyglm_third <- svyglm(formula_third, 
                             design = stamoch_lead, family = 'binomial')


library(questionr)
model1_OR <- odds.ratio(model_svyglm_first)
model2_OR <- odds.ratio(model_svyglm_second)
model3_OR <- odds.ratio(model_svyglm_third)

model1_OR$`OR (95% CI)` <- ifelse(is.na(model1_OR$`2.5 %`), "",
                                  sprintf("%.3f (%.3f-%.3f)",
                                          model1_OR$OR, model1_OR$`2.5 %`, model1_OR$`97.5 %`))
model1_OR$`p.value` <- ifelse(is.na(model1_OR$`2.5 %`), "",ifelse(model1_OR$p<0.001, "<0.001",sprintf("%.3f", model1_OR$p)))

model2_OR$`OR (95% CI)` <- ifelse(is.na(model2_OR$`2.5 %`), "",
                                  sprintf("%.3f (%.3f-%.3f)",
                                          model2_OR$OR, model2_OR$`2.5 %`, model2_OR$`97.5 %`))
model2_OR$`p.value` <- ifelse(is.na(model2_OR$`2.5 %`), "",ifelse(model2_OR$p<0.001, "<0.001",sprintf("%.3f", model2_OR$p)))

model3_OR$`OR (95% CI)` <- ifelse(is.na(model3_OR$`2.5 %`), "",
                                  sprintf("%.3f (%.3f-%.3f)",
                                          model3_OR$OR, model3_OR$`2.5 %`, model3_OR$`97.5 %`))
model3_OR$`p.value` <- ifelse(is.na(model3_OR$`2.5 %`), "",ifelse(model3_OR$p<0.001, "<0.001",sprintf("%.3f", model3_OR$p)))

write.csv(model1_OR, '01.model1_res_var.csv')
write.csv(model2_OR, '02.model2_res_var.csv')
write.csv(model3_OR, '03.model3_res_var.csv')




