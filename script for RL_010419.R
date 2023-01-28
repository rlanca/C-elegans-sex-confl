

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
original <- read.csv(file="RL_010419CC.csv", header=TRUE, sep=",")

lifespan1 <- subset(original, select = c("Replicate", "Condition1", "Condition2", "Ageofdeath", "Censored", "Plate.Column", "Plate.Row", "Strain", "Dead"))
names(lifespan1) <- c("replicate", "condition.1", "condition.2", "death.age", "censored", "column", "row", "strain", "dead")
lifespan1$treatment <- with(lifespan1, interaction(condition.1, condition.2, strain))
lifespan1$plate <- with(lifespan1, interaction(column, row))
#lifespan1$dead <- 1 added to this specific .csv file already




##Median Lifespan

jk_equal <- subset(lifespan1, treatment == "mated.equal.JK574", select = "death.age")
jk_skew <- subset(lifespan1, treatment == "mated.skew.JK574", select = "death.age")
jk_null <- subset(lifespan1, treatment == "virgin.null.JK574", select = "death.age")
px_equal <- subset(lifespan1, treatment == "mated.equal.PX632", select = "death.age")
px_skew <- subset(lifespan1, treatment == "mated.skew.PX632", select = "death.age")
px_null <- subset(lifespan1, treatment == "virgin.null.PX632", select = "death.age")

median (jk_equal$death.age)
median (jk_skew$death.age)
median (jk_null$death.age)
median (px_equal$death.age)
median (px_skew$death.age)
median (px_null$death.age)

##> median (jk_equal$death.age)
##[1] 21.64205
##> median (jk_skew$death.age)
##[1] 21.71495
##> median (jk_null$death.age)
##[1] 16.8295
##> median (px_equal$death.age)
##[1] 23.90245
##> median (px_skew$death.age)
##[1] 22.89725
##> median (px_null$death.age)
##[1] 20.78785


#04.16.2019: why are these different?????????? 
#> median (jk_equal$death.age)
#[1] 21.6629
#> median (jk_skew$death.age)
#[1] 21.7566
# median (jk_null$death.age)
#[1] 17.6004
#> median (px_equal$death.age)
#[1] 23.8816
#> median (px_skew$death.age)
#[1] 23.4233
#> median (px_null$death.age)
#[1] 20.4649


## Estimation of random effects variance
genreml <- lmer(death.age ~ treatment + (1|plate), data = lifespan1)
genprofile <- profile(genreml)
# genprofile generates warnings?
summary(genreml)
confint(genprofile)
#bad fit. 
summary(glht(genreml, linfct = mcp(treatment = c("mated.equal.JK574 - mated.skew.JK574 = 0","mated.equal.JK574 - virgin.null.JK574 = 0",
                                                 "virgin.null.JK574 - mated.equal.PX632 = 0" ,
                                                 "mated.equal.PX632 - mated.skew.PX632 = 0",
                                                 "mated.skew.PX632 - virgin.null.PX632 = 0",
                                                 "virgin.null.PX632 - mated.skew.JK574 = 0",
                                                 "mated.skew.JK574 - virgin.null.JK574 = 0"))))



#mated.equal.JK574 - mated.skew.JK574 = 0
#mated.equal.JK574 - virgin.null.JK574 = 0
#virgin.null.JK574 - mated.skew.JK574 = 0
#mated.equal.PX632 - mated.skew.PX632 = 0
#mated.equal.PX632 - virgin.null.PX632 = 0
#virgin.null.PX632 - mated.skew.PX632 = 0
#mated.skew.JK574 - mated.skew.PX632 = 0


## Cox Proportional Hazards treatment effects using random effects model within each strain
fit <- coxme(Surv(death.age, dead) ~ condition.1 + (1|plate), data = lifespan1)
# no starting estimate successful - see A below
print(fit)
summary(glht(fit, linfct = mcp(treatment = c("mated.equal.JK574 - mated.skew.JK574 = 0","mated.equal.JK574 - virgin.null.JK574 = 0",
                                          "virgin.null.JK574 - mated.equal.PX632 = 0" ,
                                          "mated.equal.PX632 - mated.skew.PX632 = 0",
                                          "mated.skew.PX632 - virgin.null.PX632 = 0",
                                          "virgin.null.PX632 - mated.skew.JK574 = 0"))))


# [A] change treatment to strain and it works - it's something about "treatment" that isn't working
fit <- coxme(Surv(death.age, dead) ~ strain + (1|plate), data = lifespan1)

# [B] example from coxme package online resource
fit2 <- coxme(Surv(time, status) ~ age + sex + (1|ph.ecog), lung)


## Plot Survival Curve

lifespan1$SurvObj <- with(lifespan1, Surv(death.age, censored == 0))


model.plot1 <- survfit(SurvObj ~ treatment, data = lifespan1, conf.type = "log-log")
summary(model.plot1)

mycols <- c("#ff9216","#ffb159","#a35803", "#5495ff","#8cb8ff","#003a99")


## each genotype is the same color lighter/darker 

plot(model.plot1, las = 1, col = mycols, lwd = 2, lty = c(1,3,1,1,3,1), xlab = "Age (Days)",
     ylab = "Survival", bty = "l")
legend("bottomleft", c("Canonical mated 1:1 (n=147)", "Canonical virgin (n=117)", "Canonical mated 3:1 (n=146)", "Wild isolate mated 1:1 (n=70)", "Wild isolate virgin (n=70)", "Wild isolate mated 3:1 (n=86)" ), lty = c(1,3,1,1,3,1), lwd = 2, cex = 0.95,
       col = mycols, bty = "n")


# Pick a color palette 
display.brewer.all(colorblindFriendly = T)
brewer.pal(#number of colors, "name of palette") 
brewer.pal(3,"BrBG")
  
  ##Save plot to pdf or png => export => save as pdf or png 
  ## running code will save to Finder => all files
  
  dev.copy(pdf,"curve.pdf", width=4, height=4) 

