library(gdata)
library(coxme)
library(multcomp)
library(lme4)
library(survival)
library(dplyr)
library(RColorBrewer)

setwd("~/Desktop/scripts_and_data/csv files")

lifespanN2 <- read.csv("N2 - Sheet1.csv", header = T)
lifespanN2 <- subset(lifespanN2, select = c("Age.at.Death", "Censored", "Plate.Column", "Plate.Row", "Strain"))
names(lifespanN2) <- c("death.age", "censored", "column", "row", "treatment")
lifespanN2$plate <- with(lifespanN2, interaction(column, row))


lifespanN2 <- subset(lifespanN2, treatment == "N2", select = c("treatment", "death.age", "censored", "plate"))
head(lifespanN2)

lifespan2 <- read.csv("RL_010419CC.csv", header = T)
lifespan2 <- subset(lifespan2, select = c("Replicate", "Condition1", "Condition2", "Ageofdeath", "Censored", "Plate.Column", "Plate.Row", "Strain", "Dead"))
names(lifespan2) <- c("replicate", "condition.1", "condition.2", "death.age", "censored", "column", "row", "strain", "dead")
lifespan2$treatment <- with(lifespan2, interaction(condition.1, condition.2, strain))
lifespan2$plate <- with(lifespan2, interaction(replicate, column, row))


#lifespan2.mated.J <- subset(lifespan2, treatment == "mated.equal.JK574", select = c("death.age", "censored", "plate", "treatment"))
lifespan2.virgin.J <- subset(lifespan2, treatment == "virgin.null.JK574", select = c("death.age", "censored", "plate", "treatment"))
combined.J <- rbind(lifespan2.mated.J, lifespan2.virgin.J)
head(combined.J)

worm <-rbind(lifespanN2,lifespan2.virgin.J)


worm$SurvObj <- with(worm, Surv(death.age, censored == 0))
model.plot2 <- survfit(SurvObj ~ treatment, data = worm, conf.type = "log-log")
summary(model.plot2)

mycols <- c("red","green")




plot(model.plot2, las = 1, col = mycols, lwd = 2, lty = c(3,3), xlab = "Age (Days)",
     ylab = "Survival", bty = "l")

legend("bottomleft", c("N2", "JK574 virgin"), lty = c(3,3), lwd = 2, cex = 0.95,
       col = mycols, bty = "n")




combreml <- lmer(death.age ~ treatment + (1|plate), data = worm)
combprofile <- profile(combreml)
summary(combreml)
confint(combprofile)
summary(glht(combreml, linfct = mcp(treatment = c("N2 - mated.equal.JK574 = 0",
                                                  "N2 - virgin.null.JK574 = 0",
                                                  "mated.equal.JK574 - virgin.null.JK574 = 0"))))

N2 - mated.equal.JK574 == 0                  0.1248    
N2 - virgin.null.JK574 == 0                  0.0266 *  
mated.equal.JK574 - virgin.null.JK574 == 0   <0.001 ***
  
  
