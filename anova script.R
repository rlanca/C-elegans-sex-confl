
library(dplyr)
library(ggpubr)

original <- read.csv(file="RL_010419.csv", header=TRUE, sep=",")

lifespan1 <- subset(original, select = c("Condition1", "Condition2", "Ageofdeath", "Strain"))
names(lifespan1) <- c("condition.1", "condition.2", "death.age", "strain")
lifespan1$conditions <- with(lifespan1, interaction(condition.1, condition.2))
lifespan1$treatment <- with (lifespan1, interaction(conditions, strain))
lifespan1$dead <- 1

print(lifespan1)

lifespan2 <- subset(lifespan1, select = c(death.age, treatment))
print(lifespan2)

dplyr::sample_n(lifespan2, 10)
levels(lifespan2$treatment)

group_by(lifespan2, treatment) %>%
  summarise(
    count = n(),
    mean = mean(death.age, na.rm = TRUE),
    sd = sd(death.age, na.rm = TRUE)
  )

ggboxplot(lifespan2, x = "treatment", y = "death.age", 
          color = "treatment", palette = c("#FC4E07", "#FC4E07", "#FC4E07", "#00AFBB", "#00AFBB", "#00AFBB"),
          order = c("mated.equal.JK574", "virgin.null.JK574", "mated.skew.JK574", "mated.equal.PX632", "virgin.null.PX632", "mated.skew.PX632"), 
          ylab = "Age of Death", xlab = "Treatment")


ggline(lifespan2, x = "treatment", y = "death.age", 
       add = c("mean_se", "jitter"), 
       order = c("mated.equal.JK574", "virgin.null.JK574", "mated.skew.JK574", "mated.equal.PX632", "virgin.null.PX632", "mated.skew.PX632"),
       ylab = "Age of Death", xlab = "Treatment")

res.aov <- aov(death.age ~ treatment, data = lifespan2)
summary(res.aov)
TukeyHSD(res.aov)


#virgin.null.JK574-mated.equal.JK574 0.0000137
#mated.skew.JK574-virgin.null.JK574 0.0017520
#mated.equal.PX632-virgin.null.JK574 0.0000004
#mated.skew.PX632-virgin.null.JK574  0.0001180

