library(gdata)
library(coxme)
library(multcomp)
library(lme4)
library(survival)
library(dplyr)
library(RColorBrewer)

setwd("~/Desktop/scripts_and_data/csv files")


lifespan1 <- read.csv("RL_110918COMB.csv", header = T)
lifespan1 <- subset(lifespan1, select = c("Replicate", "Strain",
                                          "Age.at.Death..d..Raw", "Censored", "Plate.Column", "Plate.Row", "Genotype"))
names(lifespan1) <- c("replicate", "strain", "death.age", "censored", "column", "row", "treatment")
lifespan1$plate <- with(lifespan1, interaction(replicate, column, row))
lifespan1$dead <- 1
head(lifespan1)


lifespan2 <- read.csv("RL_010419CC.csv", header = T)
lifespan2 <- subset(lifespan2, select = c("Replicate", "Condition1", "Condition2", "Ageofdeath", "Censored", "Plate.Column", "Plate.Row", "Strain", "Dead"))
names(lifespan2) <- c("replicate", "condition.1", "condition.2", "death.age", "censored", "column", "row", "strain", "dead")
lifespan2$treatment <- with(lifespan2, interaction(condition.1, condition.2, strain))
lifespan2$plate <- with(lifespan2, interaction(replicate, column, row))
head(lifespan2)

lifespan3 <- read.csv("RL_031719.csv", header = T)
lifespan3 <- subset(lifespan3, select = c("Replicate", "Condition.1", "Condition.2", "Age.at.Death", "Censored", "Plate.Column", "Plate.Row", "Strain"))
names(lifespan3) <- c("replicate", "condition.1", "condition.2", "death.age", "censored", "column", "row", "strain")
lifespan3$treatment <- with(lifespan3, interaction(condition.1, condition.2, strain))
lifespan3$plate <- with(lifespan3, interaction(replicate, column, row))
lifespan3$dead <- 1 
head(lifespan3)

lifespan1.mated <- subset(lifespan1, treatment == "mated.equal.PX632", select = c( "strain", "death.age", "censored", "plate", "treatment", "dead"))
lifespan1.virgin <- subset(lifespan1, treatment == "virgin.null.PX632", select = c("strain", "death.age", "censored", "plate", "treatment", "dead"))
lifespan2.mated <- subset(lifespan2, treatment == "mated.equal.PX632", select = c( "strain", "death.age", "censored", "plate", "treatment", "dead"))
lifespan2.virgin <- subset(lifespan2, treatment == "virgin.null.PX632", select = c( "strain", "death.age", "censored", "plate", "treatment", "dead"))
lifespan3.mated <- subset(lifespan3, treatment == "mated.equal.PX632", select = c("strain", "death.age", "censored", "plate", "treatment","dead"))
lifespan3.virgin <- subset(lifespan3, treatment == "virgin.null.PX632", select = c("strain", "death.age", "censored", "plate", "treatment", "dead"))

head(lifespan1.mated)
head(lifespan2.virgin)
head(lifespan3.mated)

combined.PX632 <- rbind(lifespan1.virgin, lifespan2.virgin, lifespan3.virgin)
head(combined.PX632)
tail(combined.PX632)

lifespan4 <-  read.csv(file="JU2526_lifespan.csv", header=T)
lifespan4 <- subset(lifespan4, select = c("Age.at.Death..d..Raw", "Censored", "Plate.Column", "Plate.Row", "Strain"))
names(lifespan4) <- c("death.age", "censored", "column", "row", "strain")
lifespan4$treatment <- with(lifespan4, interaction(strain))
lifespan4$plate <- with(lifespan4, interaction(column, row))
lifespan4$dead <- 1 
head(lifespan4)

lifespan4 <- subset(lifespan4, treatment == "JU2526", select = c( "strain", "death.age", "censored", "plate", "treatment", "dead"))

combined.2 <- rbind(combined.PX632, lifespan4)
head(combined.2)
tail(combined.2)




combined.2$SurvObj <- with(combined.2, Surv(death.age, censored == 0))
model.plot2 <- survfit(SurvObj ~ strain, data = combined.2, conf.type = "log-log")
summary(model.plot2)

# "NGM_mated", "NGM_virgin", "PX632.mated.equal", "PX632.virgin.null"
mycols2 <- c("green", "red")

plot(model.plot2, las = 1, col = mycols2, lty = c(3, 3), lwd = 2, xlab = "Age (Days)", xlim = c(0, 40),
     ylab = "Survival", bty = "l")
legend("bottomleft", c("PX632", "JU2526"),
       lty = c(3, 3), lwd = 2, cex = 1, col = mycols2, bty = "n")






combreml <- lmer(death.age ~ strain + (1| plate), data = combined.2)
combprofile <- profile(combreml)
summary(combreml)
confint(combprofile)
summary(glht(combreml, linfct = mcp(strain = c("PX632 - JU2526 = 0"))))


control=lmerControl(check.nlev.gtr.1="ignore")




combined.2$dead <- 1
combfit <- coxme(Surv(death.age, dead) ~ strain + (1| plate), data = combined.2)
print(combfit)
summary(glht(combfit, linfct = mcp(strain.treatment = c("PX632 - JU2526 = 0"))))

#model not converging