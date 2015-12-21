# This file is part of the Actuarial Statistics Project.
#
# Actuarial Statistics project is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Actuarial Statistics project is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Acturial Statistics Project  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2015 Jellen Vermeir 
# jellenvermeir@gmail.com

# Set the working directory
# setwd("..")

basefreq = read.table("baseFREQ.csv",header=TRUE,sep=";")
basefreq[1:10,]
# summary(basefreq)

n <- basefreq$nbsin

# Calculate the exposures, expressed in years
expo <- as.numeric(difftime(strptime(basefreq$fin_pol,
                                    format="%d/%m/%Y"),
        strptime(basefreq$debut_pol,format="%d/%m/%Y"),
                                  units="days")/365)

# copy categorical covariates to separate variables
sex <- basefreq$sexe
paymentFrequency <- basefreq$freq_paiement
profession <- basefreq$type_prof
nutrition <- basefreq$alimentation
terrain <- basefreq$type_territoire
usage <- basefreq$utilisation
alarm <- basefreq$presence_alarme

#####################################################
## Quantile binning:                              ###
## change continuous into categorical covariates  ###
#####################################################

# install.packages("Hmisc")
library("Hmisc")  
# Split a vector of continuous x-values in distinct categories
# 'split' is a vector that defines the category borders
# 'categoryNames' is a vector that defines the categorynames
categorize <- function(x,split,categoryNames)
{
  categories <- rep(categoryNames[1],length(x))
  for(i in seq(2,length(split)-2)){
    indices <- which(x >= split[i] & x < split[i+1])
    categories[indices] <- categoryNames[i]
  }
  # tail..
  indices <- which(x >= split[length(split)-1])
  categories[indices] <- tail(categoryNames,1)
  return(categories)
}


# Quantile bins
bins <- c(0,0.025,0.10,0.25,0.5,0.75,0.9,0.975,1)
split.Age <- wtd.quantile(basefreq$age,probs=bins)
split.Duration <- wtd.quantile(basefreq$duree_permis,
                                          probs=bins)
split.VehicleYear <- wtd.quantile(basefreq$annee_vehicule,
                                          probs=bins)

split.Age
# define category names for age category
ageCat.names <- c("18-25","26-30","31-35","36-43",
                  "44-54","54-63","63-70","71-87")
ageCat <- as.factor(categorize(basefreq$age,
              split=split.Age, categoryNames = ageCat.names))

split.Duration
durationCat.names <- c("0-4","05-09","10-14","15-21",
                       "22-29","30-32","33-37","38-77")
durationCat <- as.factor(categorize(basefreq$duree_permis,
        split=split.Duration, categoryNames = durationCat.names))

split.VehicleYear
vehicleYearCat.names <- c("1940-1884","1985-1987",
                          "1988-1990","1991-1994",
                          "1995-1997","1998-1999",
                          "2000","2001-2003")
vehicleYearCat <- as.factor(categorize(basefreq$annee_vehicule,
  split=split.VehicleYear, categoryNames = vehicleYearCat.names))

vehicleCat <- rep("Regular",nrow(basefreq))
luxuryidx <- which(basefreq$marque_voiture %in% 
                     c("BMW","MERCEDES-BENZ","ROVER","PORSCHE"))
vehicleCat[luxuryidx] <- "Luxury"
vehicleCat <- as.factor(vehicleCat)


###########################################################
############# DESCRIPTIVE STATISTICS ######################
###########################################################

freq <- table(n)
yearlyExposure <- xtabs(expo ~ n)

weights <- as.numeric(names(freq))
meanFreq <- as.numeric((freq %*% weights)/sum(expo))
# Annual claim frequency
meanFreq


# Risk classification for different age classes
age.claims <- xtabs(n ~ ageCat)
age.expo <- xtabs(expo ~ ageCat)
# Annualized claimfrequency per age class
claimFrequency.age <- age.claims/age.expo

# Risk classification for gender category
sex.claims <- xtabs(n ~ sex)
sex.expo <- xtabs(expo ~ sex)
claimFrequency.sex <- sex.claims/sex.expo

# Risk classification for different year classes
year.claims <- xtabs(n ~ vehicleYearCat)
year.expo <- xtabs(expo ~ vehicleYearCat)
claimFrequency.year <- year.claims/year.expo

# Risk classification for different terrain classes
terrain.claims <- xtabs(n ~ terrain)
terrain.expo <- xtabs(expo ~ terrain)
claimFrequency.terrain <- terrain.claims/terrain.expo

# Risk classification for different duration classes
duration.claims <- xtabs(n ~ durationCat)
duration.expo <- xtabs(expo ~ durationCat)
claimFrequency.duration <- duration.claims/duration.expo

# Risk classification for different profession classes
profession.claims <- xtabs(n ~ profession)
profession.expo <- xtabs(expo ~ profession)
claimFrequency.profession <- profession.claims/profession.expo

# Risk classification for different vehicle classes 
# (luxary vs regular)
vehicleCat.claims <- xtabs(n ~ vehicleCat)
vehicleCat.expo <- xtabs(expo ~ vehicleCat)
claimFrequency.vehicleCat <- vehicleCat.claims/vehicleCat.expo

# Risk classification for different usage classes
usage.claims <- xtabs(n ~ usage)
usage.expo <- xtabs(expo ~ usage)
claimFrequency.usage <- usage.claims/usage.expo

# Risk classification for different alarm classes
alarm.claims <- xtabs(n ~ alarm)
alarm.expo <- xtabs(expo ~ alarm)
claimFrequency.alarm <- alarm.claims/alarm.expo

# Risk classification for different paymentfrequency classes
paymentFrequency.claims <- xtabs(n ~ paymentFrequency)
paymentFrequency.expo <- xtabs(expo ~ paymentFrequency)
claimFrequency.paymentFrequency <- paymentFrequency.claims/
                                      paymentFrequency.expo


coll <- c(1,2,3,4,5)
layout(matrix(c(1,1,1,2,2,2,3,3,4,4,5,5), 
                      2, 6, byrow = TRUE))
barplot(age.claims,xlab="Age",
        main="Policyholders (Age)",
        ylab="Number of policyholders",
        beside=TRUE,col=coll[1])
barplot(year.claims,xlab="vehicle year",
        main="Policyholders (vehicle year)",
        ylab="Number of policyholders",
        beside=TRUE,col=coll[2])

barplot(sex.claims,xlab="Gender",
        main="Policyholders (Gender)",
        ylab="Number of policyholders",
        beside=TRUE,col=coll[3])
barplot(terrain.claims,xlab="Terrain",
        main="Policyholders (Terrain)",
        ylab="Number of policyholders",
        beside=TRUE,col=coll[4])
barplot(paymentFrequency.claims,xlab="Frequency",
        main="Policyholders (Frequency)",
        ylab="Number of policyholders",
        beside=TRUE,col=coll[5])

coll <- c(6,7,8,9,10)
layout(matrix(c(1,1,1,2,2,2,3,3,4,4,5,5), 
              2, 6, byrow = TRUE))
barplot(duration.claims,xlab="Duration",
        main="Policyholders (Duration)",
        ylab="Number of policyholders",
        beside=TRUE,col=coll[1])
barplot(profession.claims,xlab="Profession",
        main="Policyholders (Profession)",
        ylab="Number of policyholders",
        beside=TRUE,col=coll[2])

barplot(vehicleCat.claims,xlab="Category",
        main="Policyholders (Category)",
        ylab="Number of policyholders",
        beside=TRUE,col=coll[3])
barplot(usage.claims,xlab="Usage",
        main="Policyholders (Usage)",
        ylab="Number of policyholders",
        beside=TRUE,col=coll[4])
barplot(alarm.claims,xlab="Alarm",
        main="Policyholders (Alarm)",
        ylab="Number of policyholders",
        beside=TRUE,col=coll[5])

coll <- c(1,2,3,4,5)
layout(matrix(c(1,1,1,2,2,2,3,3,4,4,5,5), 
                        2, 6, byrow = TRUE))
barplot(claimFrequency.age,xlab="Age",
        main="Claim Frequency (Age)",
        ylab="Annual claim frequency",
        beside=TRUE,col=coll[1])
barplot(claimFrequency.year,xlab="vehicle year",
        main="Claim Frequency (vehicle year)",
        ylab="Annual claim frequency",
        beside=TRUE,col=coll[2])

barplot(claimFrequency.sex,xlab="Gender",
        main="Claim Frequency (Gender)",
        ylab="Annual claim frequency",
        beside=TRUE,col=coll[3])
barplot(claimFrequency.terrain,xlab="Terrain",
        main="Claim Frequency (Terrain)",
        ylab="Annual claim frequency",
        beside=TRUE,col=coll[4])
barplot(claimFrequency.paymentFrequency,xlab="paymentFrequency",
        main="Claim Frequency (Frequency)",
        ylab="Annual claim frequency",
        beside=TRUE,col=coll[5])

coll <- c(6,7,8,9,10)
layout(matrix(c(1,1,1,2,2,2,3,3,4,4,5,5), 
                        2, 6, byrow = TRUE))
barplot(claimFrequency.duration,xlab="Duration",
        main="Claim Frequency (Duration)",
        ylab="Annual claim frequency",
        beside=TRUE,col=coll[1])
barplot(claimFrequency.profession,xlab="Profession",
        main="Claim Frequency (Profession)",
        ylab="Annual claim frequency",
        beside=TRUE,col=coll[2])

barplot(claimFrequency.vehicleCat,xlab="Vehicle Category",
        main="Claim Frequency (Category)",
        ylab="Annual claim frequency",
        beside=TRUE,col=coll[3])
barplot(claimFrequency.usage,xlab="Usage",
        main="Claim Frequency (Usage)",
        ylab="Annual claim frequency",
        beside=TRUE,col=coll[4])
barplot(claimFrequency.alarm,xlab="Alarm",
        main="Claim Frequency (alarm)",
        ylab="Annual claim frequency",
        beside=TRUE,col=coll[5])


par(mfrow=c(2,1))
barplot(age.claims,xlab="Age Category",
        main="Descriptive statistics (Age)",
        ylab="Number of policyholders",
        beside=TRUE)
barplot(claimFrequency.age,xlab="Age Category",
        main="Descriptive statistics (Age)",
        ylab="Annual claim frequency",
        beside=TRUE)

# Risk classification for each sex
sex.claims <- xtabs(n ~ sex)
sex.expo <- xtabs(expo ~ sex)
# Annualized claimfrequency per gender
claimFrequency.sex <- sex.claims/sex.expo

par(mfrow=c(1,2))
barplot(sex.claims,xlab="Gender",
        main="Descriptive statistics (Gender)",
        ylab="Number of policyholders",
        beside=TRUE)
barplot(claimFrequency.sex,xlab="Gender",
        main="Descriptive statistics (Gender)",
        ylab="Annual claim frequency",
        beside=TRUE)


# Risk classification: Possible interaction between age and sex?
combo.claims <- xtabs(n ~ sex + ageCat)
combo.expo <- xtabs(expo ~ sex + ageCat)
claimFrequency.combo <- combo.claims/combo.expo

par(mfrow=c(2,1))
barplot(combo.claims,xlab="Age Category",
        main="Policyholders - Interaction age and gender",
        ylab="Number of policyholders",
        beside=TRUE,legend = rownames(claimFrequency.combo))
barplot(claimFrequency.combo,xlab="Age Category",
        main="Claim Frequency - Interaction age and gender",
        ylab="Annual claim frequency",
        beside=TRUE,legend = rownames(claimFrequency.combo))


###############################################################
###### CALIBRATING A simple GLM ###############################
## Verify significance of interaction between age and gender ##
###############################################################

# Trivial model
g1 <- glm(n ~ 1, fam=poisson(link=log) , offset=log(expo))
anova(g1,"Chisq")

# Only sex
g2 <- glm(n ~ 1 + sex, 
          fam=poisson(link=log), offset=log(expo))
summary(g2)
# Q-square test on drop in deviance
anova(g1,g2,"Chisq")
# H0 rejected in favor of more complexe model H1
1-pchisq(28.917,1) 

# Add age
g3 <- glm(n ~ 1 + sex + ageCat, 
          fam=poisson(link=log), offset=log(expo))
summary(g3)
anova(g2,g3)
# the sex-only model rejected in favor of age+sex model
1-pchisq(187.69,7)


# Add interactions between age and sex
g4 <- glm(n ~ 1 + sex*ageCat, 
          fam=poisson(link=log), offset=log(expo))
summary(g4)
anova(g3,g4)
1-pchisq(15.745,7) # 0.02, rejected
# !! H_0 is rejected in favor of the interaction model !!

anova(g4)


# Create a dummy variable that holds the interaction
# between age and gender
ageSexCat <- rep(0,nrow(basefreq))
ageSexCat[sex=="F"] <- categorize(basefreq$age[sex=="F"],
                    split=split.Age,categoryNames=seq(1,8))
ageSexCat[sex=="M"] <- categorize(basefreq$age[sex=="M"],
                    split=split.Age,categoryNames=seq(9,16))
ageSexCat <- as.factor(ageSexCat)


#############################################################
################# GLM - Include ALL covariates ##############
#############################################################

# sex <- basefreq$sexe
# paymentFrequency <- basefreq$freq_paiement
# profession <- basefreq$type_prof
# nutrition <- basefreq$alimentation
# terrain <- basefreq$type_territoire
# usage <- basefreq$utilisation
# alarm <- basefreq$presence_alarme


g1 <- glm(n ~ 1 + ageSexCat + paymentFrequency + profession + 
  terrain + usage + alarm + durationCat +
  vehicleYearCat + vehicleCat, fam=poisson(link=log), 
                                          offset=log(expo))
summary(g1)

##############################################################
############# Check AgeSexCat ################################
##############################################################
library(gmodels)
fit.contrast(g1,"ageSexCat",c(rep(0,8),1,-1,rep(0,6))) # 0.056

# put variable 9 and 10 together
ageSexCat5 <- rep(0,nrow(basefreq))
ageSexCat5[sex=="F"] <- categorize(basefreq$age[sex=="F"],
                        split=split.Age,categoryNames=seq(1,8))
ageSexCat5[sex=="M"] <- categorize(basefreq$age[sex=="M"],
          split=split.Age,categoryNames=c(9,9,10,11,12,13,14,15))
ageSexCat5 <- as.factor(ageSexCat5)

g2<- glm(n ~ 1 + ageSexCat5 + paymentFrequency + profession + 
  terrain + usage + alarm + durationCat +
  vehicleYearCat + vehicleCat, fam=poisson(link=log), 
                                        offset=log(expo))
summary(g2)
fit.contrast(g2,"ageSexCat5",c(rep(0,8),1,-1,rep(0,5))) # 0.34

# put variable 9 and 10 together
ageSexCat6 <- rep(0,nrow(basefreq))
ageSexCat6[sex=="F"] <- categorize(basefreq$age[sex=="F"],
                      split=split.Age,categoryNames=seq(1,8))
ageSexCat6[sex=="M"] <- categorize(basefreq$age[sex=="M"],
        split=split.Age,categoryNames=c(9,9,9,10,11,12,13,14))
ageSexCat6 <- as.factor(ageSexCat6)

g3 <- glm(n ~ 1 + ageSexCat6 + paymentFrequency + profession + 
  terrain + usage + alarm + durationCat +
  vehicleYearCat + vehicleCat, fam=poisson(link=log), 
                                        offset=log(expo))
summary(g3)
fit.contrast(g3,"ageSexCat6",c(rep(0,8),1,-1,rep(0,4))) # 0.22
fit.contrast(g3,"ageSexCat6",c(rep(0,8),1,0,-1,rep(0,3))) # 0.97

# Add var 9,10 and 11 together
ageSexCat7 <- rep(0,nrow(basefreq))
ageSexCat7[sex=="F"] <- categorize(basefreq$age[sex=="F"],
                    split=split.Age,categoryNames=seq(1,8))
ageSexCat7[sex=="M"] <- categorize(basefreq$age[sex=="M"],
            split=split.Age,categoryNames=c(9,9,9,9,9,10,11,12))
ageSexCat7 <- as.factor(ageSexCat7)

g4 <- glm(n ~ 1 + ageSexCat7 + paymentFrequency + profession + 
  terrain + usage + alarm + durationCat +
  vehicleYearCat + vehicleCat, fam=poisson(link=log), 
                                          offset=log(expo))
summary(g4)
fit.contrast(g4,"ageSexCat7",c(0,1,-1,rep(0,9))) # 0.83
fit.contrast(g4,"ageSexCat7",c(0,1,0,-1,rep(0,8))) # 0.58
fit.contrast(g4,"ageSexCat7",c(0,1,0,0,-1,rep(0,7))) # 0.99
fit.contrast(g4,"ageSexCat7",c(0,1,0,0,0,-1,rep(0,6))) # 0.74

# Add var 2,3,4,5,6 together

# Add var 9,10 and 11 together
ageSexCat8 <- rep(0,nrow(basefreq))
ageSexCat8[sex=="F"] <- categorize(basefreq$age[sex=="F"],
            split=split.Age,categoryNames=c(1,2,2,2,2,2,3,4))
ageSexCat8[sex=="M"] <- categorize(basefreq$age[sex=="M"],
            split=split.Age,categoryNames=c(5,5,5,5,5,6,7,8))
ageSexCat8 <- as.factor(ageSexCat8)

g5 <- glm(n ~ 1 + ageSexCat8 + paymentFrequency + profession + 
    terrain + usage + alarm + durationCat +
    vehicleYearCat + vehicleCat, fam=poisson(link=log), 
                                          offset=log(expo))
summary(g5)

# add all the insignificant vars togther 
ageSexCat9 <- rep(0,nrow(basefreq))
ageSexCat9[sex=="F"] <- categorize(basefreq$age[sex=="F"],
            split=split.Age,categoryNames=c(1,1,1,1,1,1,1,1))
ageSexCat9[sex=="M"] <- categorize(basefreq$age[sex=="M"],
            split=split.Age,categoryNames=c(2,2,2,2,2,1,1,1))
ageSexCat9 <- as.factor(ageSexCat9)

g5 <- glm(n ~ 1 + ageSexCat9 + paymentFrequency + profession + 
    terrain + usage + alarm + durationCat +
    vehicleYearCat + vehicleCat, fam=poisson(link=log), 
                                            offset=log(expo))
summary(g5)
# ageSexCat_1 --> females all ages AND males > 54 years old
# ageSexCat_2 --> males aged 18-54 ..
# .. positive value --> higher risk than the reference class

##############################################################
############# Check VehicleYearCat ###########################
##############################################################
claimFrequency.year

fit.contrast(g5,"vehicleYearCat",c(0,0,0,0,1,-1,0,0)) # 0.78
fit.contrast(g5,"vehicleYearCat",c(0,0,0,0,1,0,-1,0)) # 0.82
fit.contrast(g5,"vehicleYearCat",c(0,0,0,1,0,0,0,-1))

# Add var 5,6 and 7 together
vehicleYearCat2 <- as.factor(categorize(basefreq$annee_vehicule,
                          split=split.VehicleYear, 
                          categoryNames = c(1,2,3,4,5,5,5,6)))

g6 <- glm(n ~ 1 + ageSexCat9 + paymentFrequency + profession + 
  terrain + usage + alarm + durationCat +
  vehicleYearCat2 + vehicleCat, fam=poisson(link=log), 
                                          offset=log(expo))
summary(g6)
fit.contrast(g6,"vehicleYearCat2",c(0,1,-1,0,0,0)) # 0.21

# Add var 2 and 3 together
vehicleYearCat3 <- as.factor(categorize(basefreq$annee_vehicule,
                            split=split.VehicleYear, 
                            categoryNames = c(1,2,2,3,4,4,4,5)))

g7 <- glm(n ~ 1 + ageSexCat9 + paymentFrequency + profession + 
  terrain + usage + alarm + durationCat +
  vehicleYearCat3 + vehicleCat, fam=poisson(link=log), 
                                          offset=log(expo))
summary(g7)
fit.contrast(g7,"vehicleYearCat3",c(0,1,0,0,-1)) # 0.54

# Add var 2 and 5 together
vehicleYearCat4 <- as.factor(categorize(basefreq$annee_vehicule,
                            split=split.VehicleYear, 
                            categoryNames = c(2,1,1,3,4,4,4,1)))

g8 <- glm(n ~ 1 + ageSexCat9 + paymentFrequency + profession + 
  terrain + usage + alarm + durationCat +
  vehicleYearCat4 + vehicleCat, fam=poisson(link=log), 
                                          offset=log(expo))
summary(g8)
fit.contrast(g8,"vehicleYearCat4",c(0,1,0,0)) # 0.54

# param 2 insignificant (p > 0.05)
vehicleYearCat5 <- as.factor(categorize(basefreq$annee_vehicule,
                          split=split.VehicleYear, 
                          categoryNames = c(1,1,1,2,3,3,3,1)))

g9 <- glm(n ~ 1 + ageSexCat9 + paymentFrequency + profession + 
  terrain + usage + alarm + durationCat +
  vehicleYearCat5 + vehicleCat, fam=poisson(link=log), 
                                          offset=log(expo))
summary(g9)
fit.contrast(g9,"vehicleYearCat5",c(0,1,-1)) # <0.0001


split.VehicleYear
claimFrequency.year
# B_1: year 1945-1991(1-3) + 2001-2003 (8)
# B_2: year 1991-1995 (4)
# B_3: year 1995-2001 (5-7)


##############################################################
############# Check Terrain        ###########################
##############################################################
# Risk classification for different year classes
terrain.claims <- xtabs(n ~ terrain)
terrain.expo <- xtabs(expo ~ terrain)
# Annualized claimfrequency per age class
claimFrequency.terrain <- terrain.claims/terrain.expo

fit.contrast(g9,"terrain",c(0,-1,1)) 

claimFrequency.terrain
# B_1 = rural --> default
# B_2 = semi_urban --> lower claim frequency
# B_3 = urbal --> lowest claim frequency

##############################################################
############# Check durationCat        #######################
##############################################################
claimFrequency.duration

fit.contrast(g9,"durationCat",c(0,0,1,-1,0,0,0,0)) # 0.23
fit.contrast(g9,"durationCat",c(0,0,1,0,-1,0,0,0)) # 0.44
fit.contrast(g9,"durationCat",c(0,0,1,0,0,-1,0,0)) # 0.89

durationCat2 <- as.factor(categorize(basefreq$duree_permis,
    split=split.Duration, categoryNames = c(1,2,3,3,3,3,4,5)))
g10 <- glm(n ~ 1 + ageSexCat9 + paymentFrequency + profession + 
  terrain + usage + alarm + durationCat2 +
  vehicleYearCat5 + vehicleCat, fam=poisson(link=log), 
                                          offset=log(expo))
summary(g10)

fit.contrast(g10,"durationCat2",c(0,0,1,0,-1)) # 0.65

durationCat3 <- as.factor(categorize(basefreq$duree_permis,
  split=split.Duration, categoryNames = c(1,2,3,3,3,3,4,3)))
g11 <- glm(n ~ 1 + ageSexCat9 + paymentFrequency + profession + 
  terrain + usage + alarm + durationCat3 +
  vehicleYearCat5 + vehicleCat, fam=poisson(link=log), 
                                        offset=log(expo))
summary(g11)
fit.contrast(g11,"durationCat3",c(1,-1,0,0)) # 0.017 (rejected)

split.Duration
claimFrequency.duration
# B_1 -> 0-5
# B_2 -> 5-10
# B_3 -> 10-33 + 38-70
# B_4 -> 33-38


##############################################################
############# Check profession and vehiclecategory ###########
##############################################################
claimFrequency.profession
claimFrequency.vehicleCat

# not significant.. remove

g12 <- glm(n ~ 1 + ageSexCat9 + paymentFrequency + 
             terrain + usage + alarm + durationCat3 +
             vehicleYearCat5, fam=poisson(link=log), 
                                      offset=log(expo))
summary(g12)

##############################################################
############# Check usage              #######################
##############################################################
claimFrequency.usage

fit.contrast(g12,"usage",c(0,-1,1)) # 0.006 (rejected)

claimFrequency.usage
# B_1 - leasure (default)
# B_2 - occasional work --> lower claim frequency
# B_3 - daily work --> lowest claim frequency

##############################################################
############# Check alarm              #######################
##############################################################
claimFrequency.alarm
# B_1 no (default)
# B_2: yes --> negative param, lower claim frequency

##############################################################
############# Check paymentFrequency   #######################
##############################################################
claimFrequency.paymentFrequency
# B_1 yearly (default)
# B_2 monthly --> positve.. higher claim frequency

# Coefficients:
#  Estimate Std. Error z value Pr(>|z|)    
# (Intercept)              -0.97268    0.12649  -7.690 1.47e-14 ***
#  ageSexCat92               0.19718    0.02632   7.490 6.86e-14 ***
#  paymentFrequencymensuel   0.26366    0.02938   8.974  < 2e-16 ***
#  terrainSemi-urbain       -0.20518    0.06925  -2.963 0.003050 ** 
#  terrainUrbain            -0.39113    0.07094  -5.514 3.51e-08 ***
#  usageTravail-occasionnel -0.17976    0.04883  -3.682 0.000232 ***
#  usageTravail-quotidien   -0.26577    0.05401  -4.921 8.62e-07 ***
#  alarmoui                 -0.28418    0.02866  -9.914  < 2e-16 ***
#  durationCat32            -0.22900    0.09761  -2.346 0.018973 *  
#  durationCat33            -0.46655    0.09025  -5.170 2.35e-07 ***
#  durationCat34            -0.73206    0.10158  -7.207 5.72e-13 ***
#  vehicleYearCat52          0.20694    0.03895   5.313 1.08e-07 ***
#  vehicleYearCat53          0.36438    0.03636  10.021  < 2e-16 ***

claimFrequency.combo 
barplot(claimFrequency.combo,xlab="Age Category",
        main="Interaction age and gender",
        ylab="Annual claim frequency",
        beside=TRUE,legend = rownames(claimFrequency.combo))
# ageSexCat_1 --> females all ages AND males > 54 years old
# ageSexCat_2 --> males aged 18-54 ..
# .. positive value --> higher risk than 

split.VehicleYear
claimFrequency.year
# B_1: year 1945-1991(1-3) + 2001-2003 (8)
# B_2: year 1991-1995 (4)
# B_3: year 1995-2001 (5-7)

claimFrequency.terrain
# B_1 = rural --> default
# B_2 = semi_urban --> lower claim frequency
# B_3 = urbal --> lowest claim frequency

split.Duration
claimFrequency.duration
# B_1 -> 0-5 (default. highest claim frequency)
# B_2 -> 5-10 (negative param, lower frequency)
# B_3 -> 10-33 + 38-70 (even lower frequency)
# B_4 -> 33-38 (lowest claim frequency)

claimFrequency.usage
# B_1 - leasure (default)
# B_2 - occasional work --> negative param, lower claim frequency
# B_3 - daily work --> lowest claim frequency

claimFrequency.paymentFrequency
# B_1 yearly (default)
# B_2 monthly --> positve.. higher claim frequency

claimFrequency.alarm
# B_1 no (default)
# B_2: yes --> negative param, lower claim frequency