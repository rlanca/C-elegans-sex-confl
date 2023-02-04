setwd("C:/Users/Laura/Desktop/C-elegans-sex-confl-2023")
library(gdata)
library(multcomp)
library(lme4)
library(survival)
library(coxme)
library(dplyr)
library(RColorBrewer)

#Figure A:  Canonical strain analysis. 

lifespan3 <- read.csv("RL_031719 PX632 vs JK574 all mating ratios.csv", header = T)
lifespan3 <- subset(lifespan3, select = c("Replicate", "Condition.1", "Condition.2", "Age.at.Death", "Censored", "Plate.Column", "Plate.Row", "Strain"))
names(lifespan3) <- c("replicate", "condition.1", "condition.2", "death.age", "censored", "column", "row", "strain")
lifespan3$treatment <- with(lifespan3, interaction(condition.1, condition.2, strain, replicate))
lifespan3$plate <- with(lifespan3, interaction(replicate, column, row))
lifespan3$dead <- 1 
head(lifespan3)

lifespan2 <- read.csv("RL_010419CC PX632 vs JK574 all mating ratios.csv", header = T)
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

#Figure B: Wild isolate analyses. 

#clear previous plots
try(dev.off(dev.list()["RStudioGD"]), silent=TRUE)

lifespan1 <- read.csv("RL_110918COMB PX632 Auxin vs NGM.csv", header = T)
lifespan1 <- subset(lifespan1, select = c("Replicate", "Strain",
                                          "Age.at.Death..d..Raw", "Censored", "Plate.Column", "Plate.Row", "Genotype"))
names(lifespan1) <- c("replicate", "strain", "death.age", "censored", "column", "row", "treatment")
lifespan1$plate <- with(lifespan1, interaction(replicate, column, row))
lifespan1$dead <- 1
head(lifespan1)


lifespan2 <- read.csv("RL_010419CC PX632 vs JK574 all mating ratios.csv", header = T)
lifespan2 <- subset(lifespan2, select = c("Replicate", "Condition1", "Condition2", "Ageofdeath", "Censored", "Plate.Column", "Plate.Row", "Strain", "Dead"))
names(lifespan2) <- c("replicate", "condition.1", "condition.2", "death.age", "censored", "column", "row", "strain", "dead")
lifespan2$treatment <- with(lifespan2, interaction(condition.1, condition.2, strain))
lifespan2$plate <- with(lifespan2, interaction(replicate, column, row))
head(lifespan2)

lifespan3 <- read.csv("RL_031719 PX632 vs JK574 all mating ratios.csv", header = T)
lifespan3 <- subset(lifespan3, select = c("Replicate", "Condition.1", "Condition.2", "Age.at.Death", "Censored", "Plate.Column", "Plate.Row", "Strain"))
names(lifespan3) <- c("replicate", "condition.1", "condition.2", "death.age", "censored", "column", "row", "strain")
lifespan3$treatment <- with(lifespan3, interaction(condition.1, condition.2, strain))
lifespan3$plate <- with(lifespan3, interaction(replicate, column, row))
lifespan3$dead <- 1 
head(lifespan3)

lifespan1.mated <- subset(lifespan1, treatment == "mated.equal.PX632", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
lifespan1.virgin <- subset(lifespan1, treatment == "virgin.null.PX632", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
lifespan2.mated <- subset(lifespan2, treatment == "mated.equal.PX632", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
lifespan2.virgin <- subset(lifespan2, treatment == "virgin.null.PX632", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
lifespan2.skew <- subset(lifespan2, treatment == "mated.skew.PX632", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
lifespan3.mated <- subset(lifespan3, treatment == "mated.equal.PX632", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
lifespan3.virgin <- subset(lifespan3, treatment == "virgin.null.PX632", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
lifespan3.skew <- subset(lifespan3, treatment == "mated.skew.PX632", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))

combined.PX632 <- rbind(lifespan1.mated, lifespan1.virgin, lifespan2.mated, lifespan2.virgin, lifespan2.skew, lifespan3.mated, lifespan3.virgin, lifespan3.skew) 
head(lifespan2.mated)
#Plot
combined.PX632$SurvObj <- with(combined.PX632, Surv(death.age, censored == 0))
model.plot2 <- survfit(SurvObj ~ treatment, data = combined.PX632, conf.type = "log-log")
summary(model.plot2)
mycols2 <- c("#4286f4", "#6ba3ff","#01317f")
plot(model.plot2, las = 1, col = mycols2, lty = c(1,3,1), lwd = 2, xlab = "Age (Days)", xlim = c(0, 40),
     ylab = "Survival", bty = "l")
legend("bottomleft", c("mated.equal.PX632","virgin.null.PX632 ","mated.skew.PX632"),
       lty = c(1,3,1), lwd = 2, cex = 1, col = mycols2, bty = "n")

#Figure C

#Clear plots
try(dev.off(dev.list()["RStudioGD"]), silent=TRUE)

lifespan1 <- read.csv("RL_110918COMB PX632 Auxin vs NGM.csv", header = T)
lifespan1 <- subset(lifespan1, select = c("Replicate", "Strain", "Age.at.Death..d..Raw", "Censored", "Plate.Column", "Plate.Row", "Genotype"))
names(lifespan1) <- c("replicate", "strain", "death.age", "censored", "column", "row", "treatment")
lifespan1$plate <- with(lifespan1, interaction(replicate, column, row))
lifespan1$dead <- 1
head(lifespan1)

lifespan2 <- read.csv("RL_010419CC PX632 vs JK574 all mating ratios.csv", header = T)
lifespan2 <- subset(lifespan2, select = c("Replicate", "Condition1", "Condition2", "Ageofdeath", "Censored", "Plate.Column", "Plate.Row", "Strain", "Dead"))
names(lifespan2) <- c("replicate", "condition.1", "condition.2", "death.age", "censored", "column", "row", "strain", "dead")
lifespan2$treatment <- with(lifespan2, interaction(condition.1, condition.2, strain))
lifespan2$plate <- with(lifespan2, interaction(replicate, column, row))
head(lifespan2)

lifespan3 <- read.csv("RL_031719 PX632 vs JK574 all mating ratios.csv", header = T)
lifespan3 <- subset(lifespan3, select = c("Replicate", "Condition.1", "Condition.2", "Age.at.Death", "Censored", "Plate.Column", "Plate.Row", "Strain"))
names(lifespan3) <- c("replicate", "condition.1", "condition.2", "death.age", "censored", "column", "row", "strain")
lifespan3$treatment <- with(lifespan3, interaction(condition.1, condition.2, strain))
lifespan3$plate <- with(lifespan3, interaction(replicate, column, row))
lifespan3$dead <- 1 
head(lifespan3)

lifespan1.mated <- subset(lifespan1, treatment == "mated.equal.PX632", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
lifespan1.virgin <- subset(lifespan1, treatment == "virgin.null.PX632", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
lifespan2.mated <- subset(lifespan2, treatment == "mated.equal.PX632", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
lifespan2.virgin <- subset(lifespan2, treatment == "virgin.null.PX632", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
lifespan3.mated <- subset(lifespan3, treatment == "mated.equal.PX632", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
lifespan3.virgin <- subset(lifespan3, treatment == "virgin.null.PX632", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
combined.PX632 <- rbind(lifespan1.mated, lifespan1.virgin,  lifespan2.mated, lifespan2.virgin, lifespan3.mated, lifespan3.virgin) 


lifespan2.mated.J <- subset(lifespan2, treatment == "mated.equal.JK574", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
lifespan2.virgin.J <- subset(lifespan2, treatment == "virgin.null.JK574", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
combined.JK574 <- rbind(lifespan2.mated.J, lifespan2.virgin.J)

worm <-rbind(combined.JK574, combined.PX632)

#Create a survival object for plotting survival curves
worm$SurvObj <- with(worm, Surv(death.age, censored == 0))
model.plot2 <- survfit(SurvObj ~ treatment, data = worm, conf.type = "log-log")
summary(model.plot2)

#Add plot details
mycols <- c("#ff9216","#ffb159","#5495ff","#8cb8ff","#003a99")
plot(model.plot2, las = 1, col = mycols, lwd = 2, lty = c(1,3,1,3,1), xlab = "Age (Days)", ylab = "Survival", bty = "l")
legend("bottomleft", c("mated.equal.JK574", "virgin.null.JK574", "mated.equal.PX632", "virgin.null.PX632"), lty = c(1,3,1,3,1,3), lwd = 2, cex = 0.95, col = mycols, bty = "n")

#Supplementary Figure 1. JK574 group mating

#Clear plots 
try(dev.off(dev.list()["RStudioGD"]), silent=TRUE)

lifespan2.equal.J <- subset(lifespan2, treatment == "mated.equal.JK574", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
lifespan2.skew.J <- subset(lifespan2, treatment == "mated.skew.JK574", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
lifespan2.virgin.J <- subset(lifespan2, treatment == "virgin.null.JK574", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
worm <- rbind(lifespan2.skew.J, lifespan2.virgin.J, lifespan2.equal.J)
worm$SurvObj <- with(worm, Surv(death.age, censored == 0))
model.plot2 <- survfit(SurvObj ~ treatment, data = worm, conf.type = "log-log")
summary(model.plot2)
mycols <- c("#ff9216","#ffb159","#703d03")
plot(model.plot2, las = 1, col = mycols, lwd = 2, lty = c(1,3,1,3,1), xlab = "Age (Days)",
     ylab = "Survival", bty = "l")
legend("bottomleft", c("mated.equal.JK574", "virgin.null.JK574", "mated.skew.JK574"), lty = c(1,3,1,3,1,3), lwd = 2, cex = 0.95,
       col = mycols, bty = "n")
