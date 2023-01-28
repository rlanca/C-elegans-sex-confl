

library(gdata)
library(coxme)
library(multcomp)
library(lme4)
library(survival)
library(dplyr)
library(RColorBrewer)

# Add replicate column in sheets 

# Read .csv data into R

setwd("~/Desktop/scripts_and_data/csv files")
original <- read.csv(file="RL_031719.csv", header=TRUE, sep=",")

lifespan1 <- subset(original, select = c("Replicate", "Condition.1", "Condition.2", "Age.at.Death", "Censored", "Plate.Column", "Plate.Row", "Strain"))
names(lifespan1) <- c("replicate", "condition.1", "condition.2", "death.age", "censored", "column", "row", "strain")
lifespan1$treatment <- with(lifespan1, interaction(condition.1, condition.2, strain))
lifespan1$plate <- with(lifespan1, interaction(column, row))
lifespan1$dead <- 1 

lifespan2 <- read.csv("RL_010419CC.csv", header = T)
lifespan2 <- subset(lifespan2, select = c("Replicate", "Condition1", "Condition2", "Ageofdeath", "Censored", "Plate.Column", "Plate.Row", "Strain", "Dead"))
names(lifespan2) <- c("replicate", "condition.1", "condition.2", "death.age", "censored", "column", "row", "strain", "dead")
lifespan2$treatment <- with(lifespan2, interaction(condition.1, condition.2, strain))
lifespan2$plate <- with(lifespan2, interaction(replicate, column, row))
head(lifespan2)

#individuals assay

JK574.equal <- subset(lifespan2, treatment == "mated.equal.JK574", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
JK574.skew <- subset(lifespan2, treatment == "mated.skew.JK574", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
JK574.virgin <- subset(lifespan2, treatment == "virgin.null.JK574", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))

combined.JK574 <- rbind(JK574.skew, JK574.virgin, JK574.equal)
head(combined.JK574)

combreml <- lmer(death.age ~ treatment + (1|plate), data = combined.JK574)
combprofile <- profile(combreml)
summary(combreml)
confint(combprofile)
summary(glht(combreml, linfct = mcp(treatment = c("mated.skew.JK574 - mated.equal.JK574 =0",
                                                  "mated.skew.JK574 - virgin.null.JK574 =0",
                                                  "virgin.null.JK574 - mated.equal.JK574=0"))))


combined.JK574$SurvObj <- with(combined.JK574, Surv(death.age, censored == 0))
model.plot2 <- survfit(SurvObj ~ treatment, data = combined.JK574, conf.type = "log-log")
summary(model.plot2)

mycols2 <- c("blue", "red", "green")

plot(model.plot2, las = 1, col = mycols2, lty = c(1,1,1), lwd = 2, xlab = "Age (Days)", xlim = c(0, 40),
     ylab = "Survival", bty = "l")
legend("bottomleft", c("equal", "virgin", "skew)"),
       lty = c(1,1,1), lwd = 2, cex = 1, col = mycols2, bty = "n")


#combined JK574 replicates 

lifespan3 <- read.csv("RL_031719.csv", header = T)
lifespan3 <- subset(lifespan3, select = c("Replicate", "Condition.1", "Condition.2", "Age.at.Death", "Censored", "Plate.Column", "Plate.Row", "Strain"))
names(lifespan3) <- c("replicate", "condition.1", "condition.2", "death.age", "censored", "column", "row", "strain")
lifespan3$treatment <- with(lifespan3, interaction(condition.1, condition.2, strain, replicate))
lifespan3$plate <- with(lifespan3, interaction(replicate, column, row))
lifespan3$dead <- 1 
head(lifespan3)

lifespan2 <- read.csv("RL_010419CC.csv", header = T)
lifespan2 <- subset(lifespan2, select = c("Replicate", "Condition1", "Condition2", "Ageofdeath", "Censored", "Plate.Column", "Plate.Row", "Strain", "Dead"))
names(lifespan2) <- c("replicate", "condition.1", "condition.2", "death.age", "censored", "column", "row", "strain", "dead")
lifespan2$treatment <- with(lifespan2, interaction(condition.1, condition.2, strain, replicate))
lifespan2$plate <- with(lifespan2, interaction(replicate, column, row))
head(lifespan2)

JK574.mated.individual <- subset(lifespan3, treatment == "mated.skew.JK574.3", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
JK574.virgin.individual <- subset(lifespan3, treatment == "virgin.null.JK574.3", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
JK574.mated.group <- subset(lifespan2, treatment == "mated.skew.JK574.2", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
JK574.virgin.group <- subset(lifespan2, treatment == "virgin.null.JK574.2", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))

combined.JK574 <- rbind(JK574.mated.individual, JK574.virgin.individual, JK574.mated.group, JK574.virgin.group)
head(combined.JK574)  

combined.JK574$SurvObj <- with(combined.JK574, Surv(death.age, censored == 0))
model.plot2 <- survfit(SurvObj ~ treatment, data = combined.JK574, conf.type = "log-log")
summary(model.plot2)

mycols2 <- c("#501800", "#501800", "#ff802b", "#ff802b")

plot(model.plot2, las = 1, col = mycols2, lty = c(1, 3, 1, 3), lwd = 2, xlab = "Age (Days)", xlim = c(0, 40),
     ylab = "Survival", bty = "l")
legend("bottomleft", c("Individual virgins", "Individual mated 3:1", "Group virgins", "Group mated 3:1"),
       lty = c(1, 3, 1, 3), lwd = 2, cex = 1, col = mycols2, bty = "n")


#estimate random effects variance 

combreml <- lmer(death.age ~ treatment + (1|replicate / plate), data = combined.JK574)
combprofile <- profile(combreml)
summary(combreml)
confint(combprofile)
summary(glht(combreml, linfct = mcp(treatment = c("mated.skew.JK574.3 - virgin.null.JK574.3 = 0", 
                                                  "virgin.null.JK574.3 - mated.skew.JK574.2 = 0",
                                                  "mated.skew.JK574.2 - virgin.null.JK574.2 = 0",
                                                  "virgin.null.JK574.2 - mated.skew.JK574.3 = 0",
                                                  "virgin.null.JK574.3 - virgin.null.JK574.2 = 0", 
                                                  "mated.skew.JK574.3 - JK574.mated.skew.2 = 0"))))

summary(glht(combreml, linfct = mcp(treatment = c(
                                                  "virgin.null.JK574.3 - virgin.null.JK574.2 = 0", 
                                                  "mated.skew.JK574.3 - mated.skew.JK574.2 = 0"))))

#virgin vs mated JK574 is significant for both runs (group p<0.001, individual p=0.0146)

# cox hazard 

combined.JK574$dead <- 1
combfit <- coxme(Surv(death.age, dead) ~ treatment + (1|replicate / plate), data = combined.JK574)
print(combfit)
summary(glht(combfit, linfct = mcp(strain.treatment = c("mated.skew.JK574.3 - virgin.null.JK574.3", 
                                                        "virgin.null.JK574.3 - mated.skew.JK574.2",
                                                        "mated.skew.JK574.2 - virgin.null.JK574.2",
                                                        "virgin.null.JK574.2 - mated.skew.JK574.3"))))

#no starting estimate successful 
  