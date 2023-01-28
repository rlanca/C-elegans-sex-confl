library(gdata)
library(coxme)
library(multcomp)
library(lme4)
library(survival)
library(dplyr)
library(RColorBrewer)

## In Excel: add "Replicate"

setwd("~/Desktop/scripts_and_data/csv files")

lifespan1 <- read.csv("RL_110918UP.csv", header = T)
lifespan1 <- subset(lifespan1, select = c("Replicate", "Strain",
                                          "Age.at.Death..d..Raw", "Censored", "Plate.Column", "Plate.Row", "Genotype"))
names(lifespan1) <- c("replicate", "strain", "death.age", "censored", "column", "row", "treatment")
lifespan1$plate <- with(lifespan1, interaction(replicate, column, row))
lifespan1$dead <- 1
head(lifespan1)


lifespan2 <- read.csv("RL_010419CC.csv", header = T)
lifespan2 <- subset(lifespan2, select = c("Replicate", "Condition1", "Condition2", "Ageofdeath", "Censored", "Plate.Column", "Plate.Row", "Strain", "Dead"))
names(lifespan2) <- c("replicate", "condition.1", "condition.2", "death.age", "censored", "column", "row", "strain", "dead")
lifespan2$treatment <- with(lifespan2, interaction(condition.1, condition.2, strain, replicate))
lifespan2$plate <- with(lifespan2, interaction(replicate, column, row))
head(lifespan2)

lifespan3 <- read.csv("RL_031719.csv", header = T)
lifespan3 <- subset(lifespan3, select = c("Replicate", "Condition.1", "Condition.2", "Age.at.Death", "Censored", "Plate.Column", "Plate.Row", "Strain"))
names(lifespan3) <- c("replicate", "condition.1", "condition.2", "death.age", "censored", "column", "row", "strain")
lifespan3$treatment <- with(lifespan3, interaction(condition.1, condition.2, strain, replicate))
lifespan3$plate <- with(lifespan3, interaction(replicate, column, row))
lifespan3$dead <- 1 
head(lifespan3)

lifespan1.mated <- subset(lifespan1, treatment == "NGM_mated", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
lifespan1.virgin <- subset(lifespan1, treatment == "NGM_virgin", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
lifespan2.mated <- subset(lifespan2, treatment == "mated.equal.PX632.2", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
lifespan2.virgin <- subset(lifespan2, treatment == "virgin.null.PX632.2", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
lifespan3.mated <- subset(lifespan3, treatment == "mated.equal.PX632.3", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
lifespan3.virgin <- subset(lifespan3, treatment == "virgin.null.PX632.3", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))

head(lifespan1.mated)
head(lifespan2.virgin)
head(lifespan3.mated)

onetwo <- rbind(lifespan1.mated, lifespan1.virgin,  lifespan2.mated, lifespan2.virgin)
three <- rbind(lifespan3.mated, lifespan3.virgin)
combined.PX632 <- rbind(onetwo, three)
head(combined.PX632)
tail(combined.PX632)

## Plot
combined.PX632$SurvObj <- with(combined.PX632, Surv(death.age, censored == 0))
model.plot2 <- survfit(SurvObj ~ treatment, data = combined.PX632, conf.type = "log-log")
summary(model.plot2)

# "NGM_mated", "NGM_virgin", "PX632.mated.equal", "PX632.virgin.null"
mycols2 <- c("blue", "blue", "red", "red", "green", "green")

plot(model.plot2, las = 1, col = mycols2, lty = c(1, 3, 1, 3, 1, 3), lwd = 2, xlab = "Age (Days)", xlim = c(0, 40),
     ylab = "Survival", bty = "l")
legend("bottomleft", c("Rep1 mated (n = 183)", "Rep1 virgin (n = 229)", "Rep2 mated (n = 70)", "Rep2 'virgin' (n = 70)", "Rep3 'mated' (n = 180)", "Rep3 'virgin' (n = 177)"),
       lty = c(1, 3, 1, 3), lwd = 2, cex = 1, col = mycols2, bty = "n")



## Estimation of random effects variance
combreml <- lmer(death.age ~ treatment + (1|replicate / plate), data = combined.PX632)
combprofile <- profile(combreml)
summary(combreml)
confint(combprofile)
summary(glht(combreml, linfct = mcp(treatment = c("NGM_mated - NGM_virgin = 0",
                                                         "NGM_virgin - mated.equal.PX632.2 = 0",
                                                         "mated.equal.PX632.2 - virgin.null.PX632.2 = 0",
                                                         "virgin.null.PX632.2 - mated.equal.PX632.3 = 0",
                                                         "mated.equal.PX632.3 - virgin.null.PX632.3 = 0",
                                                         "virgin.null.PX632.3 - NGM_mated = 0"))))

control=lmerControl(check.nlev.gtr.1="ignore")

#no significant differences

## Cox Proportional Hazards treatment effects using random effects model within each strain
combined.PX632$dead <- 1
combfit <- coxme(Surv(death.age, dead) ~ treatment + (1|replicate / plate), data = combined.PX632)
print(combfit)
summary(glht(combfit, linfct = mcp(strain.treatment = c("NGM_mated - NGM_virgin = 0",
                                                        "NGM_virgin - mated.equal.PX632.2 = 0",
                                                        "mated.equal.PX632.2 - virgin.null.PX632.2 = 0",
                                                        "virgin.null.PX632.2 - mated.equal.PX632.3 = 0",
                                                        "mated.equal.PX632.3 - virgin.null.PX632.3 = 0",
                                                        "virgin.null.PX632.3 - NGM_mated = 0"))))
#model not converging

