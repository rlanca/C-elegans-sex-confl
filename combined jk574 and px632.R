


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

#separate 1:1, 3:1

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


worm$SurvObj <- with(combined.JK574, Surv(death.age, censored == 0))
model.plot2 <- survfit(SurvObj ~ treatment, data = worm, conf.type = "log-log")
summary(model.plot2)

mycols <- c("#ff9216","#ffb159","#5495ff","#8cb8ff","#003a99")




plot(model.plot2, las = 1, col = mycols, lwd = 2, lty = c(1,3,1,3,1), xlab = "Age (Days)",
     ylab = "Survival", bty = "l")

legend("bottomleft", c("mated.equal.JK574", "virgin.null.JK574", "mated.equal.PX632", "virgin.null.PX632"), lty = c(1,3,1,3,1,3), lwd = 2, cex = 0.95,
       col = mycols, bty = "n")






combreml <- lmer(death.age ~ treatment + (1|replicate / plate), data = worm)
combprofile <- profile(combreml)
summary(combreml)
confint(combprofile)
summary(glht(combreml, linfct = mcp(treatment = c("mated.equal.PX632 - virgin.null.PX632 = 0",
                                                  "virgin.null.PX632 - virgin.null.JK574 = 0",
                                                  "virgin.null.JK574 - mated.equal.JK574 = 0",
                                                  "mated.equal.PX632 - mated.equal.JK574 =0"))))


Linear Hypotheses:
  Estimate Std. Error z value Pr(>|z|)    
mated.equal.PX632 - virgin.null.PX632 == 0 -0.14058    0.54121  -0.260    0.992    
virgin.null.PX632 - virgin.null.JK574 == 0  3.91015    0.70935   5.512   <1e-04 ***
  virgin.null.JK574 - mated.equal.JK574 == 0 -3.69062    0.82164  -4.492   <1e-04 ***
  mated.equal.PX632 - mated.equal.JK574 == 0  0.07894    0.68178   0.116    0.999  




surv_object <- Surv(time = worm$death.age, event = worm$dead)
surv_object

fit1 <- survfit(surv_object ~ treatment, data = worm)
summary(fit1)


fit <- coxme(Surv(death.age, dead) ~ treatment + (1|replicate / plate), data = worm)
print(combfit)
summary(glht(combfit, linfct = mcp(treatment = c(""))))



#JK574 group mating

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

combreml <- lmer(death.age ~ treatment + (1|replicate / plate), data = worm)
combprofile <- profile(combreml)
summary(combreml)
confint(combprofile)
summary(glht(combreml, linfct = mcp(treatment = c("mated.equal.PX632 - virgin.null.PX632 = 0",
                                                  "virgin.null.PX632 - virgin.null.JK574 = 0",
                                                  "virgin.null.JK574 - mated.equal.JK574 = 0",
                                                  "mated.equal.PX632 - mated.equal.JK574 =0"))))

