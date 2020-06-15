#R Student competition, Mikolaj Kacki, June 2019.

install.packages("tidyverse")
install.packages("Metrics")

library(rgl)
library(corrplot)
library(plot3D)#packages for 3d plotting
library(MASS)#selection 
library(standardize)#standarization for PCA

library(mlbench)#caret 
library(caret)

library(Boruta)
library(randomForest)

library(ggplot2)
library(factoextra)
library(RColorBrewer)

library(Metrics)

#missing values in data
crystal <- read.table("crystal.txt", header=TRUE)
head(crystal)
sum(is.na(crystal)) 
#no missing values

summary(crystal$y)
#melt.point ranges from 68.17 to 220.22 with a mean of 142.05 and median of 142.07

#distribution plot
ggplot(crystal, aes(y))+
  geom_histogram(binwidth = 1, color="red")+
  ggtitle("Distribution of melting point")+
  xlab("Melting point temperature")+
  ylab("")

unique(crystal$y)

#Correlation plot -----------------------------------------------------------------------------------------------------------------

C <- cor(crystal)
corrplot(C, type = "upper",col = brewer.pal(n = 8, name = "RdBu"),method = "number", number.cex=0.7)

#y is strongly correlated with x3 (0.66), x4 (0.82), x7 (0.55) 

#REGRESSION ANALYSIS ---------------------------------------------------------------------------------------------------------------

#fitting the complete model with all explanatory variables
crystal.lm.all <- lm(y ~ ., crystal)
summary(crystal.lm.all)

#intercept only model
crystal.lm.itcpt <- lm(y ~ 1, crystal)
summary(crystal.lm.itcpt)

#Model Selection ----------------------------

#stepwise algorithm chooses by AIC criterion
fit1 <- step(crystal.lm.itcpt, direction= 'forward', scope=formula(crystal.lm.all))
summary(fit1)

#GOODNESS of fit ----------------------------

#from the output R^2=0.82 adjR^2=0.8

rmse(crystal$y, predict(fit1))#15.15497 

rmse(crystal$y, predict(crystal.lm.all))#11.84387
rmse(crystal$y, predict(crystal.lm.itcpt))#35.76147

#Analysis------------------------------------

plot(crystal$x4, crystal$y)
plot(crystal$x3, crystal$y)
#plots indicate clear linear dependence between these two variables and y

var(crystal$x3) #1023.528
var(crystal$x4) #98.57
#x4 has 10 times smaller variance then x3

cor(crystal$x3, crystal$x4)#0.9159431

#huge corelation also between these two what makes sense
#as they represent similar variable 

ggplot(crystal, aes(x7, y))+
  geom_point(color="blue")+
  xlab("Melting point")+
  ylab("Polar surface area")

unique(crystal$x7) #x7 has 8 distinct values, but plot shows
#linear dependence excluding one outlier of x7=136.817

enthalpy.2 <- crystal$x4
enthalpy.1 <- crystal$x3
melting.point <- crystal$y

scatter3D(enthalpy.2, melting.point , enthalpy.1, grid=TRUE,
          surface = TRUE, col = "blue", pch = 21, xlab="Enthalpy 2", ylab="Enthalpy 2",
          zlab="Melting Point") 

crystal$enthalpy <- enthalpy.1+enthalpy.2

ggplot(crystal, aes(melting.point, enthalpy))+
  geom_point(color="blue", size=2)+ 
  geom_smooth(se=FALSE, color="red")+
  ggtitle("Dependence between melting point and summed enthalpy")+
  xlab("Melting point")+
  ylab("Summed enthalpies")

#BORUTA--------------------------------------------------------------------------------- 

crystal$enthalpy <- NULL

set.seed(123) #generating random numbers
boruta1 <- Boruta(y ~., data= crystal, doTrace=3, maxRuns = 500)
plot(boruta1, cex.axis =0.5, main="Boruta output")
#confirmed big impact of x3,4 and 7

#PCA ------------------------------------------------------------------------------------

crystal$enthalpy <- NULL
crystal_scaled <- scale(crystal)
#I will scale it to make sure differences between variances don't change pca

pc <- princomp(crystal_scaled, cor=TRUE, score=TRUE)
summary(pc)
plot(pc)

pr_var <- pc$sdev^2
prop_varex <- pr_var/sum(pr_var) #proportion of variance explained by each component
plot(cumsum(prop_varex), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained", type="b")

summary(pc)

#running PCA on crystal indicates that PC1-PC11 explain 96.3% of variability.