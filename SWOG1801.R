## How to reconstruct data from published KM curves

### everything explained here: https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-021-01308-8

### install package
install.packages("IPDfromKM")
library(IPDfromKM)

## First step, go on digitilzest and for each arm (important, use 100 as the maximum in y axis), export as excel files on your Desktop (EXP for the experimental arm, and CON for the control arm)
## use the Import Dataset function in R Studio to import excel files. 
## then put your EXP data in E : 
E <- EXP
## attribute new column names to E
colnames(E) <-c("time", "survival probability")

## the same with the control arm (called C)
C <- CON
colnames(C) <-c("time", "survival probability")

##time risk, separate by comma (no space allowed)
trisk <- c(0,6,12,18,24,30,36)

## number at risk experimental arm (no space allowed)
nrisk.E <- c(154,96,69,46,25,17,1)

## number at risk control arm (no space allowed)
nrisk.C <- c(159,98,67,40,22,10,2)

## preprocess the data
pre_E <- preprocess(dat=E, trisk=trisk, nrisk=nrisk.E, maxy=100)
pre_C <- preprocess(dat=C, trisk=trisk, nrisk=nrisk.C, maxy=100)

## then individual patient data = IPD, (experimental treatment is 1, control group treatment is 0)
est_E_OS <- getIPD(prep=pre_E, armID=1, tot.events = NULL)
est_C_OS <- getIPD(prep=pre_C, armID=0, tot.events = NULL)

## you can isolate the IPD part (time, status, treat) from est_E_OS and est_C_OS
est_E_OS$IPD
est_C_OS$IPD

### calculate percentage censored at each time interval in both arms
est_E_OS$riskmat$censor.perc <- round((est_E_OS$riskmat$censor.hat / est_E_OS$riskmat$nrisk.hat) * 100, 1)
est_C_OS$riskmat$censor.perc <- round((est_C_OS$riskmat$censor.hat / est_C_OS$riskmat$nrisk.hat) * 100, 1)

## the summary function allow to have the data of events, censoring, percentage censoring, etc. 
summary(est_E_OS)
summary (est_C_OS)

## percentages to be compared easily : 
print(est_E_OS$riskmat$censor.perc)
print(est_C_OS$riskmat$censor.perc)

#### SURVIVAL ANALYIS
install.packages("rms")
library(rms)

## here you can create a KM curve based on the novel IPD data generated based on the data that have been digitilzed
## you can also recapitulate the Cox analysis. 
surv_data<- rbind(est_E_OS$IPD, est_C_OS$IPD)
names(surv_data) <- c("time", "status", "arm")
x <- surv_data$time
y <- surv_data$status
z <- surv_data$arm

## Kaplan-Meier Curve
par(mar = c(4.1, 4.4, 4.1, 1.9), xaxs = "i", yaxs = "i")
plot(survfit(Surv(x,y)~z), xlim = c(0, 30), ylim = c(0, 1), 
     col=c("#66CC66", "#6699CC"), lwd=3, xaxt='n', bty = "l", 
     las = 1, cex.axis = 1.5, tcl  = -0.5)
axis(side=1, at=seq(0, 30, by=3), cex.axis = 1.5)

## Cox, HR and p-value
summary(coxph(Surv(x,y)~z))

#------------------------------------------------------------------------------

## simulation of a specific "censoring scenario"
## modify some values of the reconstructed IPD data

newdata <- rbind(est_E_OS$IPD, est_C_OS$IPD)

# Set the seed for reproducibility
set.seed(123)

#---------------------------------EXPERIMENTAL ARM

# Parameters for the experimental arm (censored for additional toxicity, more likely to present the event)
# here we will take the censored patients from 0 to 6 months
# then randomly pick a percentage of those (to pick the 8 additionnal censored patients during the first time point, we will set the percentage at 29.6%)
# and change their status between being censored to having the event
EstartTime <- 0 # start time
EendTime <- 6 # end time
Eperc <- 29.6 # Percentage of "0" status to be changed to "1"

# Subset rows with treat = 1 and time within [EstartTime, EendTime]
subset_treat_1 <- newdata[newdata$treat == 1 & newdata$time >= EstartTime & newdata$time <= EendTime,]

# Further subset those with status = 0
subset_treat_1_status_0 <- subset_treat_1[subset_treat_1$status == 0,]

# Calculate number of rows to change based on percentage
e_change <- round(nrow(subset_treat_1_status_0) * (Eperc / 100))

# Randomly select rows to change
rows_to_change <- sample(nrow(subset_treat_1_status_0), e_change)

# Change status to 1 for these rows
subset_treat_1_status_0$status[rows_to_change] <- 1

# Replace the original rows in newdata
newdata[newdata$treat == 1 & newdata$time >= EstartTime & newdata$time <= EendTime & newdata$status == 0,] <- subset_treat_1_status_0

## the number of patients with data being change in the experimental arm : 
e_change

## -------------------------"Modeled" SURVIVAL ANALYIS
##preparation
x <- newdata$time
y <- newdata$status
z <- newdata$treat

## Kaplan-Meier Curve
par(mar = c(4.1, 4.4, 4.1, 1.9), xaxs = "i", yaxs = "i")
plot(survfit(Surv(x,y)~z), xlim = c(0, 30), ylim = c(0, 1), 
     col=c("#66CC66", "#6699CC"), lwd=3, xaxt='n', bty = "l", 
     las = 1, cex.axis = 1.5, tcl  = -0.5)
axis(side=1, at=seq(0, 30, by=3), cex.axis = 1.5)

## Cox, HR and p-value
summary(coxph(Surv(x,y)~z))