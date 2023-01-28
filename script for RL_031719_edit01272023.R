#Install R packages if necessary
install.packages("Rtools")
install.packages("gdata")
install.packages("coxme")
install.packages("multcomp")
install.packages("lme4")
install.packages("survival")
install.packages("dplyr")
install.packages("RColorBrewer")

#Load into current session
library(gdata)
library(coxme)
library(multcomp)
library(lme4)
library(survival)
library(dplyr)
library(RColorBrewer)

#Read data into R
setwd("~/Documents/scripts_and_data/csv files")
original <- read.csv(file="RL_031719.csv", header=TRUE, sep=",")

#
lifespan1 <- subset(original, select = c("Replicate", "Condition.1", "Condition.2", "Age.at.Death", "Censored", "Plate.Column", "Plate.Row", "Strain"))
names(lifespan1) <- c("replicate", "condition.1", "condition.2", "death.age", "censored", "column", "row", "strain")
lifespan1$treatment <- with(lifespan1, interaction(condition.1, condition.2, strain))
lifespan1$plate <- with(lifespan1, interaction(column, row))
lifespan1$dead <- 1 
head(lifespan1)


# only JK574 individuals



lifespan2.skew.J <- subset(lifespan1, treatment == "mated.skew.JK574", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
lifespan2.virgin.J <- subset(lifespan1, treatment == "virgin.null.JK574", select = c("replicate", "strain", "death.age", "censored", "plate", "treatment"))
worm <- rbind(lifespan2.skew.J, lifespan2.virgin.J)
worm$SurvObj <- with(worm, Surv(death.age, censored == 0))
model.plot2 <- survfit(SurvObj ~ treatment, data = worm, conf.type = "log-log")
summary(model.plot2)
mycols <- c("#ffb159","#ff9216")
plot(model.plot2, las = 1, col = mycols, lwd = 2, lty = c(3,1), xlab = "Age (Days)",
     ylab = "Survival", bty = "l")
legend("bottomleft", c("virgin.null.JK574.ind", "mated.skew.JK574.ind"), lty = c(3,1), lwd = 2, cex = 0.95,
       col = mycols, bty = "n")
##Median Lifespan


jk_skew <- subset(lifespan1, treatment == "mated.skew.JK574", select = "death.age")
jk_null <- subset(lifespan1, treatment == "virgin.null.JK574", select = "death.age")
px_equal <- subset(lifespan1, treatment == "mated.equal.PX632", select = "death.age")
px_skew <- subset(lifespan1, treatment == "mated.skew.PX632", select = "death.age")
px_null <- subset(lifespan1, treatment == "virgin.null.PX632", select = "death.age")


median (jk_skew$death.age)
median (jk_null$death.age)
median (px_equal$death.age)
median (px_skew$death.age)
median (px_null$death.age)


#results
#> median (jk_skew$death.age)
#[1] 20.2519
#> median (jk_null$death.age)
#[1] 16.2415
#> median (px_equal$death.age)
#[1] 19.6998
#> median (px_skew$death.age)
#[1] 20.26235
#> median (px_null$death.age)
#[1] 20.82485


## Estimation of random effects variance
genreml <- lmer(death.age ~ treatment + (1|plate), data = lifespan1)
genprofile <- profile(genreml)

# genprofile generates warnings?
summary(genreml)
confint(genprofile)
summary(glht(genreml, linfct = mcp(treatment = c("virgin.null.JK574 - mated.skew.JK574 = 0",
                                                 "mated.skew.JK574 - mated.equal.PX632 = 0" ,
                                                 "mated.skew.JK574 - mated.skew.PX632 = 0",
                                                 "mated.equal.PX632 - mated.skew.PX632 = 0",
                                                 "mated.skew.PX632 - virgin.null.PX632 = 0",
                                                 "virgin.null.PX632 - virgin.null.JK574 = 0"))))



#Linear Hypotheses:
#  Estimate Std. Error z value Pr(>|z|)    
#virgin.null.JK574 - mated.skew.JK574 == 0  -0.92905    0.92268  -1.007    0.795    
#mated.skew.JK574 - mated.equal.PX632 == 0  -1.68856    0.83450  -2.023    0.190    
#mated.skew.JK574 - mated.skew.PX632 == 0   -1.67263    0.85837  -1.949    0.221    
#mated.equal.PX632 - mated.skew.PX632 == 0   0.01593    0.72959   0.022    1.000    
#mated.skew.PX632 - virgin.null.PX632 == 0  -0.98910    0.75411  -1.312    0.598    
#virgin.null.PX632 - virgin.null.JK574 == 0  3.59079    0.82657   4.344   <0.001 ***


## Cox Proportional Hazards treatment effects using random effects model within each strain
fit <- coxme(Surv(death.age, dead) ~ treatment + (1|plate), data = lifespan1)
print(fit)
summary(glht(fit, linfct = mcp(treatment = c("virgin.null.JK574 - mated.skew.JK574 = 0",
                                             "mated.skew.JK574 - mated.equal.PX632 = 0" ,
                                             "mated.skew.JK574 - mated.skew.PX632 = 0",
                                             "mated.equal.PX632 - mated.skew.PX632 = 0",
                                             "mated.skew.PX632 - virgin.null.PX632 = 0",
                                             "virgin.null.PX632 - virgin.null.JK574 = 0"))))


## Plot Survival Curve

lifespan1$SurvObj <- with(lifespan1, Surv(death.age, censored == 0))


model.plot1 <- survfit(SurvObj ~ treatment, data = lifespan1, conf.type = "log-log")
summary(model.plot1)

mycols <- c("#ff9216","#ffb159","#5495ff","#8cb8ff","#003a99")


## each genotype is the same color lighter/darker 

plot(model.plot1, las = 1, col = mycols, lwd = 2, lty = c(1,3,1,1,3,1), xlab = "Age (Days)",
     ylab = "Survival", bty = "l")
legend("bottomleft", c("Canonical strain virgins (n = )", "Canonical strain mated 1:1 (n = )", "Wild isolate mated 1:1 (n = )", "Wild isolate virgins (n = )", "Wild isolate mated 3:1 (n=)"), lty = c(1,3,1,1,3,1), lwd = 2, cex = 0.95,
       col = mycols, bty = "n")


# Pick a color palette 
display.brewer.all(colorblindFriendly = T)
brewer.pal(#number of colors, "name of palette") 
  brewer.pal(3,"BrBG")
  
  ##Save plot to pdf or png => export => save as pdf or png 
  ## running code will save to Finder => all files
  
  dev.copy(pdf,"curve.pdf", width=4, height=4) 
  
  
  
 
  