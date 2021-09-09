#====================================================================.
# Advances in completely automated vowel analysis for sociophonetics:
# Using end-to-end speech recognition systems with DARLA                                                    
# Rolando Coto-Solano (rolando.a.coto.solano@dartmouth.edu)                                                   
# James Stanford (james.n.stanford@dartmouth.edu)
# Last updated: September 9th, 2021
#====================================================================.

library(tidyverse)
library(sp)
library(adehabitatHR)
library(vowels)
library(lme4)
library(optimx)
library(effects)
library(sjPlot)
library(lmerTest)
library(MuMIn)
library(emmeans)

#====================================================================.
#== Read and prepare raw data                                    ====
#====================================================================.


folder = ""
idea <- read.csv(paste(folder, "idea-dialect-data-2020.csv", sep = ""), header=T, sep=",")
lev  <- read.csv(paste(folder, "idea-lev-dist-2020.csv", sep = ""), header=T, sep=",")


# Remove all tokens that precede a liquid, which is standard practice for NCS, etc.
idea <- subset(idea, idea$fol_seg!="R")
idea <- subset(idea, idea$fol_seg!="L")


# Make a new column anae_type for managing regional generalizations 

idea$anae_type <- "other"
idea$anae_type[idea$anae=="border blue/pink"] <- "border blue/pink"
idea$anae_type[idea$anae=="border red/orange"] <- "border red/orange"
idea$anae_type[idea$anae=="GeneralSouth_red"] <- "GeneralSouth_red"
idea$anae_type[idea$anae=="InlandNorth_dark_blue"] <- "InlandNorth_dark_blue"
idea$anae_type[idea$anae=="North_light_blue_not_InlandNorth"] <- "North_light_blue_not_InlandNorth"


# Make a new column binaryNorthSouth for analyzing SVS as north versus south:

idea$binaryNorthSouth <- "other"
idea$binaryNorthSouth[idea$anae_type == "GeneralSouth_red" & idea$state != "FLORIDA"] <- "generalSouth"
idea$binaryNorthSouth[idea$anae_type != "GeneralSouth_red"
                      & idea$anae_type != "border red/orange" & idea$state != "FLORIDA"] <- "generalNorth"


# Make a new column compareNorths for analyzing NCS as inland north versus other north
idea$compareNorths <- "other"
idea$compareNorths[idea$anae_type == "InlandNorth_dark_blue"] <- "inlandNorth"
idea$compareNorths [idea$anae_type != "GeneralSouth_red" & idea$anae_type != "border red/orange" & idea$state != "FLORIDA" & idea$anae_type != "InlandNorth_dark_blue"] <- "northWithoutInlandNorth"


# Double-check the counts of these regions, in case of overlapping subsets:

nrow(idea[idea$binaryNorthSouth == "generalSouth",])          # 28225
nrow(idea[idea$binaryNorthSouth == "generalNorth",])          # 55607
nrow(idea[idea$binaryNorthSouth == "other",])                 # 4680

28225 + 55607 + 4680
nrow(idea)

nrow(idea[idea$compareNorths == "inlandNorth",])              # 14427
nrow(idea[idea$compareNorths == "northWithoutInlandNorth",])  # 41180
nrow(idea[idea$compareNorths == "other",])                    # 32905

14427 + 41180 + 32905
nrow(idea)


# Divide fol_seg into a subcategories so that we can test
# following environment for: nasal, voicing, and coda/no coda

idea$fol_type <- "voicedObstruent"
idea$fol_type[idea$fol_seg =="N" | idea$fol_seg =="NG" |  idea$fol_seg =="M"] <- "nasal"
idea$fol_type[idea$fol_seg=="T" | 
                idea$fol_seg=="P" | 
                idea$fol_seg=="K" | 
                idea$fol_seg=="S" |
                idea$fol_seg=="F" |
                idea$fol_seg=="CH" |
                idea$fol_seg=="SH" | 
                idea$fol_seg=="TH" |
                idea$fol_seg=="HH"  	
] <- "voicelessObstruent"

summary(as.factor(idea$fol_type))
mydata = idea


# prepare dataset with just the fields that we need:
forRbrul <- c()
forRbrul$speaker <- idea$nameAudioFile
forRbrul$F1 <- idea$F1_LobanovNormed_unscaled
forRbrul$F2 <- idea$F2_LobanovNormed_unscaled
forRbrul$vowel <- idea$vowel
forRbrul$word <- idea$word
forRbrul$fol_type <- idea$fol_type
forRbrul$typeTranscription <- idea$typeTranscription
forRbrul$gender <- as.character(idea$gender)
forRbrul$estimatedYOB <- idea$estimatedYOB
forRbrul$compareNorths <- idea$compareNorths
forRbrul$binaryNorthSouth <- idea$binaryNorthSouth
forRbrul <- data.frame(forRbrul)
forRbrul$zEstimatedYOB = (forRbrul$estimatedYOB - mean(forRbrul$estimatedYOB))/sd(forRbrul$estimatedYOB)
t = forRbrul


#====================================================================.
#== Calculate Bhattacharyya distance                         ====
#====================================================================.

i = 1
batDist = data.frame(speaker=character(), anae=character(), vowel=character(), typeTranscription=character(), bat=double(), gender=character(), estimatedYOB=integer())
for (inSpk in unique(t$speaker)) {
  print(inSpk)
  print(i)
  inAnae = as.character(unique(subset(t,speaker==inSpk)$binaryNorthSouth))
  inGender = as.character(unique(subset(t,speaker==inSpk)$gender))
  inEstimatedYOB = as.character(unique(subset(t,speaker==inSpk)$estimatedYOB))
  for (inVowel in unique(t$vowel)) {
    tempVowelsDS = subset(t,speaker==inSpk&vowel==inVowel&(typeTranscription == "deepspeech"|typeTranscription == "groundtruth"))
    tempVowelsIH = subset(t,speaker==inSpk&vowel==inVowel&(typeTranscription == "inhouse"|typeTranscription == "groundtruth"))
    rowsDS = nrow(subset(tempVowelsDS,typeTranscription=="deepspeech"))
    rowsDSGT = nrow(subset(tempVowelsDS,typeTranscription=="groundtruth"))
    rowsIH = nrow(subset(tempVowelsIH,typeTranscription=="inhouse"))
    rowsIHGT = nrow(subset(tempVowelsIH,typeTranscription=="inhouse"))
    if ( rowsDS > 4 & rowsIH > 4 & rowsDSGT > 4 & rowsIHGT > 4) {
      dfDS = SpatialPointsDataFrame(cbind(tempVowelsDS$F1, tempVowelsDS$F2), data.frame(factor(tempVowelsDS$typeTranscription)))
      tempDS = kerneloverlap(dfDS, method = 'BA')[1,2]
      dfIH = SpatialPointsDataFrame(cbind(tempVowelsIH$F1, tempVowelsIH$F2), data.frame(factor(tempVowelsIH$typeTranscription)))
      tempIH = kerneloverlap(dfIH, method = 'BA')[1,2]
      tempRow = c(speaker=inSpk, anae=inAnae, vowel=inVowel, typeTranscription="deepspeech", bat=tempDS, gender=inGender, estimatedYOB=inEstimatedYOB)
      batDist = rbind(batDist, tempRow)
      tempRow = c(speaker=inSpk, anae=inAnae, vowel=inVowel, typeTranscription="inhouse", bat=tempIH, gender=inGender, estimatedYOB=inEstimatedYOB)
      batDist = rbind(batDist, tempRow)
      i = i+1
    }
  }
}
colnames(batDist) = c("speaker", "anae", "vowel", "typeTranscription", "bat", "gender", "estimatedYOB")
batDist$vowel = as.factor(batDist$vowel)
batDist$typeTranscription = as.factor(batDist$typeTranscription)
batDist$bat = as.numeric(batDist$bat)


#====================================================================.
#== Sec3.1 CER LMER                                              ====
#====================================================================.

#== Prepare data

names = subset(idea, typeTranscription=="groundtruth")$name
names = unique(str_replace(names, "-gt",""))
lev = lev[lev$id %in% names,]
unique(lev$id)

lev$percImprov = lev$diffLev / lev$levIH

lev$anaeBigGroups = "other"
lev$anaeBigGroups[lev$anae=="GeneralSouth_red"] = "south"
lev$anaeBigGroups[lev$anae=="border red/orange"] = "south"
lev$anaeBigGroups[lev$anae=="InlandNorth_dark_blue"] = "north"
lev$anaeBigGroups[lev$anae=="border blue/pink"] = "other"
lev$anaeBigGroups[lev$anae=="border North_light_blue_not_InlandNorth/pink"] = "other"

unique(lev$anae)
unique(lev$anaeBigGroups)

lev$gender = ""
lev$estimatedAge = -1
lev$estimatedYOB = -1

for (i in c(1:length(lev$id))){
  id = lev$id[[i]]
  id = paste(id,"-gt",sep="")
  print(id)
  g = idea$gender[idea$name==id][1]
  print(g)
  ag = idea$estimatedAge[idea$name==id][1]
  print(ag)
  yob = idea$estimatedYOB[idea$name==id][1]
  print(yob)
  lev$gender[[i]] = g
  lev$estimatedAge[[i]] = ag
  lev$estimatedYOB[[i]] = yob
  print(lev$gender[[i]])
  print(lev$estimatedAge[[i]])
  print(lev$estimatedYOB[[i]])
  
}

t1 = lev
t1$cer = t1$cerDS
t1$typeTranscription = "DeepSpeech"
t2 = lev
t2$cer = t2$cerIH
t2$typeTranscription = "Sphinx"
t = rbind(t1,t2)
t$anaeBigGroups = factor(t$anaeBigGroups, levels = c("other","north","south"))
t = subset(t, cer<3)

t = subset(t, !is.na(estimatedYOB))

# Centering age

t$zEstimatedYOB = (t$estimatedYOB - mean(t$estimatedYOB))/sd(t$estimatedYOB)
t$estimatedYOB
t$zEstimatedYOB

#== LMER

hist(t$cer)
hist(log(t$cer))

# 
#mLev.maximal <- lmer(log(cer) ~ typeTranscription*gender + anaeBigGroups +  estimatedYOB + (1 | id), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
#summary(mLev)

mLev.minimal <- lmer(log(cer) ~ typeTranscription*gender + anaeBigGroups +  zEstimatedYOB + (1 | id), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(mLev.minimal)


# Assumptions:
plot(fitted(mLev.minimal),resid(mLev.minimal))
abline(h=0, col="blue")
# Normality of residuals: the distribution is unimodal with similar tails on both sides.
# Moreover, there are 5730 datapoints, so, using the Central Limit Theorem, we can 
# assume that the residuals behave in a manner similar to the Normal distribution.
hist(resid(mLev.minimal))
qqnorm(resid(mLev.minimal))
qqline(resid(mLev.minimal), col = "steelblue", lwd = 2)
length(resid(mLev.minimal))

plot_model(mLev.minimal, type='diag') # you can ajust type (see package info: ?plot_model)

# Descriptive statistics

wT = with(t, tapply(cer, list(typeTranscription, gender), median))
wT

wT = with(t, tapply(cer, list(typeTranscription, gender), mean))
wT

wT = with(t, tapply(cer, list(typeTranscription, gender), sd))
wT

wt = with(t, tapply(cer, typeTranscription, median))
wT

wt = with(t, tapply(cer, typeTranscription, mean))
wt

wt = with(t, tapply(cer, gender, mean))
wt

wt = with(t, tapply(cer, anaeBigGroups, mean))
wt


#====================================================================.
#== Fig6 South Vowel Triangles (GT)                              ====
#====================================================================.

x <- subset(idea, mydata$binaryNorthSouth == "generalSouth" | mydata$binaryNorthSouth == "generalNorth"  )
x <- subset(x, x$typeTranscription == "groundtruth")

x$binaryNorthSouth = as.character(x$binaryNorthSouth)
x$binaryNorthSouth[x$binaryNorthSouth == "generalSouth"] = "Southern"
x$binaryNorthSouth[x$binaryNorthSouth == "generalNorth"] = "General North"
threeColorChoice = c("#FF0000", "#0000FF", "#00FF00", "#FF0000", "#FF3D00", "#FF7A00", "#FF0000", "#FF3D00", "#FF7A00")

to.plot <- c()
to.plot$speaker <- as.character(x$binaryNorthSouth)
to.plot$vowel <- as.character(x$vowel)
to.plot$context <- "NA"
to.plot$F1 <- as.numeric(x$F1_LobanovNormed_unscaled)
to.plot$F2 <- as.numeric(x$F2_LobanovNormed_unscaled )
to.plot$F3 <- "NA"
to.plot$gl.F1 <- "NA"
to.plot$gl.F2 <- "NA"
to.plot$gl.F3 <- "NA"
to.plot <- data.frame(to.plot)

# Export size: 600 x 600
vowelplot(compute.means(to.plot, separate=T), color="speakers",color.choice=threeColorChoice,  size=0, label="vowels", title="Southern versus General North\nGround truth only")


#====================================================================.
#== Fig7 South Vowel Triangles (All three types)                 ====
#====================================================================.

# SVS all three

x <- subset(idea, mydata$binaryNorthSouth == "generalSouth")

x$typeTranscription = as.character(x$typeTranscription)
x$typeTranscription[x$typeTranscription=="groundtruth"] = "Ground truth"
x$typeTranscription[x$typeTranscription=="inhouse"] = "Sphinx"
x$typeTranscription[x$typeTranscription=="deepspeech"] = "DeepSpeech"
threeColorChoice = c("#0000FF", "#FF0000", "#00FF00", "#FF0000", "#FF3D00", "#FF7A00", "#FF0000", "#FF3D00", "#FF7A00")

to.plot <- c()
to.plot$speaker <- as.character(x$typeTranscription)
to.plot$vowel <- as.character(x$vowel)
to.plot$context <- "NA"
to.plot$F1 <- as.numeric(x$F1_LobanovNormed_unscaled)
to.plot$F2 <- as.numeric(x$F2_LobanovNormed_unscaled )
to.plot$F3 <- "NA"
to.plot$gl.F1 <- "NA"
to.plot$gl.F2 <- "NA"
to.plot$gl.F3 <- "NA"
to.plot <- data.frame(to.plot)

# Export size: 600x600
vowelplot(compute.means(to.plot, separate=T), color="speakers", color.choice=threeColorChoice, size=0, label="vowels", title="Southern vowels in\nall three transcription types")


#====================================================================.
#== Fig8  BhatDist South (Bhattacharyya distance)                ====
#====================================================================.

batSouth = subset(batDist, anae=="generalSouth")

tBatSouth = batSouth
tBatSouth$typeTranscription = as.character(tBatSouth$typeTranscription)
tBatSouth$typeTranscription[tBatSouth$typeTranscription=="deepspeech"] = "DeepSpeech"
tBatSouth$typeTranscription[tBatSouth$typeTranscription=="inhouse"] = "Sphinx"
tBatSouth$typeTranscription = as.factor(tBatSouth$typeTranscription)

# export size: 1200 x 350
ggplot(tBatSouth, aes(x=vowel, y=bat, fill=typeTranscription)) +
  geom_boxplot()+
  labs(title="Bhattacharyya's Affinity for Southern vowels by type of transcription",x="Vowel", y = "Bhattacharyya's Affinity", fill="Transcription")+
  ylim(0.6, 1)

hist(batSouth$bat)
hist(batSouth$logBat)


#====================================================================.
#== Table1 BhatDist South Vowels                                 ====
#====================================================================.

medBhatDistSouth = with(batSouth, tapply(bat, list(typeTranscription, vowel), median))
medBhatDistSouth = t(medBhatDistSouth)
medBhatDistSouth = data.frame(medBhatDistSouth)
medBhatDistSouth$diff = medBhatDistSouth$deepspeech - medBhatDistSouth$inhouse

medBhatDistSouth <- medBhatDistSouth[order(medBhatDistSouth$diff),]
medBhatDistSouth


#====================================================================.
#== Fig9 South Vowel Triangle (DS, Sphinx)                       ====
#====================================================================.

# Plot8a: Ground truth versus DS

x <- subset(idea, mydata$binaryNorthSouth == "generalSouth" | mydata$binaryNorthSouth == "generalNorth"  )
x <- subset(x, x$typeTranscription == "deepspeech")

x$binaryNorthSouth = as.character(x$binaryNorthSouth)
x$binaryNorthSouth[x$binaryNorthSouth == "generalSouth"] = "Southern"
x$binaryNorthSouth[x$binaryNorthSouth == "generalNorth"] = "General North"
threeColorChoice = c("#FF0000", "#0000FF", "#00FF00", "#FF0000", "#FF3D00", "#FF7A00", "#FF0000", "#FF3D00", "#FF7A00")

to.plot <- c()
to.plot$speaker <- as.character(x$binaryNorthSouth)
to.plot$vowel <- as.character(x$vowel)
to.plot$context <- "NA"
to.plot$F1 <- as.numeric(x$F1_LobanovNormed_unscaled)
to.plot$F2 <- as.numeric(x$F2_LobanovNormed_unscaled )
to.plot$F3 <- "NA"
to.plot$gl.F1 <- "NA"
to.plot$gl.F2 <- "NA"
to.plot$gl.F3 <- "NA"
to.plot <- data.frame(to.plot)

# Export size: 600 x 600
vowelplot(compute.means(to.plot, separate=T), color="speakers",color.choice=threeColorChoice,  size=0, label="vowels", title="Southern versus General North\nDeepSpeech only")

# Plot8b: Ground truth versus Sphinx

x <- subset(idea, mydata$binaryNorthSouth == "generalSouth" | mydata$binaryNorthSouth == "generalNorth"  )
x <- subset(x, x$typeTranscription == "inhouse")

x$binaryNorthSouth = as.character(x$binaryNorthSouth)
x$binaryNorthSouth[x$binaryNorthSouth == "generalSouth"] = "Southern"
x$binaryNorthSouth[x$binaryNorthSouth == "generalNorth"] = "General North"
x$binaryNorthSouth = as.factor(x$binaryNorthSouth)
threeColorChoice = c("#0000FF", "#FF0000", "#00FF00", "#FF0000", "#FF3D00", "#FF7A00", "#FF0000", "#FF3D00", "#FF7A00")

to.plot <- c()
to.plot$speaker <- as.character(x$binaryNorthSouth)
to.plot$vowel <- as.character(x$vowel)
to.plot$context <- "NA"
to.plot$F1 <- as.numeric(x$F1_LobanovNormed_unscaled)
to.plot$F2 <- as.numeric(x$F2_LobanovNormed_unscaled )
to.plot$F3 <- "NA"
to.plot$gl.F1 <- "NA"
to.plot$gl.F2 <- "NA"
to.plot$gl.F3 <- "NA"
to.plot <- data.frame(to.plot)

# Export size: 600 x 600
vowelplot(compute.means(to.plot, separate=T), color="speakers",color.choice=threeColorChoice,  size=0, label="vowels", title="Southern versus General North\nSphinx only")


#====================================================================.
#== Fig10 Table2 Data for South vowels                              ====
#== SVS: Southern Vowel Shift
#====================================================================.

#== EY_SVS

EY_SVS <- subset(forRbrul, binaryNorthSouth=="generalNorth" | binaryNorthSouth=="generalSouth")
EY_SVS <- subset(EY_SVS, vowel=="EY")
EY_SVS <- subset(EY_SVS, typeTranscription != "inhouse")
EY_SVS$typeTranscription <- factor(EY_SVS$typeTranscription) 
# calculate EYshifting using Labov's formula F2- 2 x F1
EY_SVS$EYshifting <- EY_SVS$F2 - 2*EY_SVS$F1
EY_SVS = na.omit(EY_SVS)

nrow(EY_SVS)
mean(EY_SVS$EYshifting)
hist(EY_SVS$EYshifting)
qqnorm(EY_SVS$EYshifting, pch = 1, frame = FALSE)
qqline(EY_SVS$EYshifting, col = "steelblue", lwd = 2)

EY_SVS_south_ds = subset(EY_SVS, binaryNorthSouth == "generalSouth" & typeTranscription == "deepspeech")
mean_EY_SVS_south_ds = mean(EY_SVS_south_ds$EYshifting)
EY_SVS_south_gt = subset(EY_SVS, binaryNorthSouth == "generalSouth" & typeTranscription == "groundtruth")
mean_EY_SVS_south_gt = mean(EY_SVS_south_gt$EYshifting)
EY_SVS_north_ds = subset(EY_SVS, binaryNorthSouth == "generalNorth" & typeTranscription == "deepspeech")
mean_EY_SVS_north_ds = mean(EY_SVS_north_ds$EYshifting)
EY_SVS_north_gt = subset(EY_SVS, binaryNorthSouth == "generalNorth" & typeTranscription == "groundtruth")
mean_EY_SVS_north_gt = mean(EY_SVS_north_gt$EYshifting)
delta_EY_SVS_ds = abs(mean_EY_SVS_south_ds - mean_EY_SVS_north_ds)
delta_EY_SVS_gt = abs(mean_EY_SVS_south_gt - mean_EY_SVS_north_gt)
delta_EY_SVS_ds
delta_EY_SVS_gt

plot(EY_SVS$estimatedYOB, EY_SVS$EYshifting)
abline(lm(EY_SVS$EYshifting ~ EY_SVS$estimatedYOB))

#== EH_SVS

EH_SVS <- subset(forRbrul, binaryNorthSouth=="generalNorth" | binaryNorthSouth=="generalSouth")
EH_SVS <- subset(EH_SVS, vowel=="EH")
EH_SVS <- subset(EH_SVS, typeTranscription != "inhouse")
EH_SVS$typeTranscription <- factor(EH_SVS$typeTranscription) 
# calculate EHshifting using Labov's formula F2- 2 x F1
EH_SVS$EHshifting <- EH_SVS$F2 - 2*EH_SVS$F1
EH_SVS = na.omit(EH_SVS)

nrow(EH_SVS)
mean(EH_SVS$EHshifting)
hist(EH_SVS$EHshifting)
qqnorm(EH_SVS$EHshifting, pch = 1, frame = FALSE)
qqline(EH_SVS$EHshifting, col = "steelblue", lwd = 2)

ehSVSgender = with(EH_SVS, tapply(EHshifting, list(gender), mean))
ehSVSgender

plot(EH_SVS$estimatedYOB, EH_SVS$EHshifting)
abline(lm(EH_SVS$EHshifting ~ EH_SVS$estimatedYOB))
lm(EH_SVS$EHshifting ~ EH_SVS$estimatedYOB)

EH_SVS_south_ds = subset(EH_SVS, binaryNorthSouth == "generalSouth" & typeTranscription == "deepspeech")
mean_EH_SVS_south_ds = mean(EH_SVS_south_ds$EHshifting)
EH_SVS_south_gt = subset(EH_SVS, binaryNorthSouth == "generalSouth" & typeTranscription == "groundtruth")
mean_EH_SVS_south_gt = mean(EH_SVS_south_gt$EHshifting)
EH_SVS_north_ds = subset(EH_SVS, binaryNorthSouth == "generalNorth" & typeTranscription == "deepspeech")
mean_EH_SVS_north_ds = mean(EH_SVS_north_ds$EHshifting)
EH_SVS_north_gt = subset(EH_SVS, binaryNorthSouth == "generalNorth" & typeTranscription == "groundtruth")
mean_EH_SVS_north_gt = mean(EH_SVS_north_gt$EHshifting)
delta_EH_SVS_ds = abs(mean_EH_SVS_south_ds - mean_EH_SVS_north_ds)
delta_EH_SVS_gt = abs(mean_EH_SVS_south_gt - mean_EH_SVS_north_gt)
delta_EH_SVS_ds
delta_EH_SVS_gt

#== IY_SVS

IY_SVS <- subset(forRbrul, binaryNorthSouth=="generalNorth" | binaryNorthSouth=="generalSouth")
IY_SVS <- subset(IY_SVS, vowel=="IY")
IY_SVS <- subset(IY_SVS, typeTranscription != "inhouse")
IY_SVS$typeTranscription <- factor(IY_SVS$typeTranscription) 
# calculate IYshifting using Labov's formula F2- 2 x F1
IY_SVS$IYshifting <- IY_SVS$F2 - 2*IY_SVS$F1
IY_SVS = na.omit(IY_SVS)

nrow(IY_SVS)
mean(IY_SVS$IYshifting)
hist(IY_SVS$IYshifting)
qqnorm(IY_SVS$IYshifting, pch = 1, frame = FALSE)
qqline(IY_SVS$IYshifting, col = "steelblue", lwd = 2)

plot(IY_SVS$estimatedYOB, IY_SVS$IYshifting)
abline(lm(IY_SVS$IYshifting ~ IY_SVS$estimatedYOB))
lm(IY_SVS$IYshifting ~ IY_SVS$estimatedYOB)

IY_SVS_south_ds = subset(IY_SVS, binaryNorthSouth == "generalSouth" & typeTranscription == "deepspeech")
mean_IY_SVS_south_ds = mean(IY_SVS_south_ds$IYshifting)
IY_SVS_south_gt = subset(IY_SVS, binaryNorthSouth == "generalSouth" & typeTranscription == "groundtruth")
mean_IY_SVS_south_gt = mean(IY_SVS_south_gt$IYshifting)
IY_SVS_north_ds = subset(IY_SVS, binaryNorthSouth == "generalNorth" & typeTranscription == "deepspeech")
mean_IY_SVS_north_ds = mean(IY_SVS_north_ds$IYshifting)
IY_SVS_north_gt = subset(IY_SVS, binaryNorthSouth == "generalNorth" & typeTranscription == "groundtruth")
mean_IY_SVS_north_gt = mean(IY_SVS_north_gt$IYshifting)
delta_IY_SVS_ds = abs(mean_IY_SVS_south_ds - mean_IY_SVS_north_ds)
delta_IY_SVS_gt = abs(mean_IY_SVS_south_gt - mean_IY_SVS_north_gt)
delta_IY_SVS_ds
delta_IY_SVS_gt

#== IH_SVS

IH_SVS <- subset(forRbrul, binaryNorthSouth=="generalNorth" | binaryNorthSouth=="generalSouth")
IH_SVS <- subset(IH_SVS, vowel=="IH")
IH_SVS <- subset(IH_SVS, typeTranscription != "inhouse")
IH_SVS$typeTranscription <- factor(IH_SVS$typeTranscription) 
# calculate IHshifting using Labov's formula F2- 2 x F1
IH_SVS$IHshifting <- IH_SVS$F2 - 2*IH_SVS$F1
IH_SVS = na.omit(IH_SVS)

nrow(IH_SVS)
mean(IH_SVS$IHshifting)
hist(IH_SVS$IHshifting)
qqnorm(IH_SVS$IHshifting, pch = 1, frame = FALSE)
qqline(IH_SVS$IHshifting, col = "steelblue", lwd = 2)

plot(IH_SVS$estimatedYOB, IH_SVS$IHshifting)
abline(lm(IH_SVS$IHshifting ~ IH_SVS$estimatedYOB))
lm(IH_SVS$IHshifting ~ IH_SVS$estimatedYOB)

IH_SVS_south_ds = subset(IH_SVS, binaryNorthSouth == "generalSouth" & typeTranscription == "deepspeech")
mean_IH_SVS_south_ds = mean(IH_SVS_south_ds$IHshifting)
IH_SVS_south_gt = subset(IH_SVS, binaryNorthSouth == "generalSouth" & typeTranscription == "groundtruth")
mean_IH_SVS_south_gt = mean(IH_SVS_south_gt$IHshifting)
IH_SVS_north_ds = subset(IH_SVS, binaryNorthSouth == "generalNorth" & typeTranscription == "deepspeech")
mean_IH_SVS_north_ds = mean(IH_SVS_north_ds$IHshifting)
IH_SVS_north_gt = subset(IH_SVS, binaryNorthSouth == "generalNorth" & typeTranscription == "groundtruth")
mean_IH_SVS_north_gt = mean(IH_SVS_north_gt$IHshifting)
delta_IH_SVS_ds = abs(mean_IH_SVS_south_ds - mean_IH_SVS_north_ds)
delta_IH_SVS_gt = abs(mean_IH_SVS_south_gt - mean_IH_SVS_north_gt)
delta_IH_SVS_ds
delta_IH_SVS_gt

#== AE_SVS

AE_SVS <- subset(forRbrul, binaryNorthSouth=="generalNorth" | binaryNorthSouth=="generalSouth")
AE_SVS <- subset(AE_SVS, vowel=="AE")
AE_SVS <- subset(AE_SVS, typeTranscription != "inhouse")
AE_SVS$typeTranscription <- factor(AE_SVS$typeTranscription) 
# calculate AEshifting using Labov's formula F2- 2 x F1
AE_SVS$AEshifting <- AE_SVS$F2 - 2*AE_SVS$F1
AE_SVS = na.omit(AE_SVS)

nrow(AE_SVS)
mean(AE_SVS$AEshifting)
hist(AE_SVS$AEshifting)
qqnorm(AE_SVS$AEshifting, pch = 1, frame = FALSE)
qqline(AE_SVS$AEshifting, col = "steelblue", lwd = 2)

AE_SVS_south_ds = subset(AE_SVS, binaryNorthSouth == "generalSouth" & typeTranscription == "deepspeech")
mean_AE_SVS_south_ds = mean(AE_SVS_south_ds$AEshifting)
AE_SVS_south_gt = subset(AE_SVS, binaryNorthSouth == "generalSouth" & typeTranscription == "groundtruth")
mean_AE_SVS_south_gt = mean(AE_SVS_south_gt$AEshifting)
AE_SVS_north_ds = subset(AE_SVS, binaryNorthSouth == "generalNorth" & typeTranscription == "deepspeech")
mean_AE_SVS_north_ds = mean(AE_SVS_north_ds$AEshifting)
AE_SVS_north_gt = subset(AE_SVS, binaryNorthSouth == "generalNorth" & typeTranscription == "groundtruth")
mean_AE_SVS_north_gt = mean(AE_SVS_north_gt$AEshifting)
delta_AE_SVS_ds = abs(mean_AE_SVS_south_ds - mean_AE_SVS_north_ds)
delta_AE_SVS_gt = abs(mean_AE_SVS_south_gt - mean_AE_SVS_north_gt)
delta_AE_SVS_ds
delta_AE_SVS_gt

#== AW_SVS

AW_SVS <- subset(forRbrul, binaryNorthSouth=="generalNorth" | binaryNorthSouth=="generalSouth")
AW_SVS <- subset(AW_SVS, vowel=="AW")
AW_SVS <- subset(AW_SVS, typeTranscription != "inhouse")
AW_SVS$typeTranscription <- factor(AW_SVS$typeTranscription) 
AW_SVS = na.omit(AW_SVS)

nrow(AW_SVS)
mean(AW_SVS$F2)
hist(AW_SVS$F2)
qqnorm(AW_SVS$F2, pch = 1, frame = FALSE)
qqline(AW_SVS$F2, col = "steelblue", lwd = 2)

plot(AW_SVS$estimatedYOB, AW_SVS$F2)
abline(lm(AW_SVS$F2 ~ AW_SVS$estimatedYOB))
lm(AW_SVS$F2 ~ AW_SVS$estimatedYOB)

AW_SVS_south_ds = subset(AW_SVS, binaryNorthSouth == "generalSouth" & typeTranscription == "deepspeech")
mean_AW_SVS_south_ds = mean(AW_SVS_south_ds$F2)
AW_SVS_south_gt = subset(AW_SVS, binaryNorthSouth == "generalSouth" & typeTranscription == "groundtruth")
mean_AW_SVS_south_gt = mean(AW_SVS_south_gt$F2)
AW_SVS_north_ds = subset(AW_SVS, binaryNorthSouth == "generalNorth" & typeTranscription == "deepspeech")
mean_AW_SVS_north_ds = mean(AW_SVS_north_ds$F2)
AW_SVS_north_gt = subset(AW_SVS, binaryNorthSouth == "generalNorth" & typeTranscription == "groundtruth")
mean_AW_SVS_north_gt = mean(AW_SVS_north_gt$F2)
delta_AW_SVS_ds = abs(mean_AW_SVS_south_ds - mean_AW_SVS_north_ds)
delta_AW_SVS_gt = abs(mean_AW_SVS_south_gt - mean_AW_SVS_north_gt)
delta_AW_SVS_ds
delta_AW_SVS_gt

#== OW_SVS

OW_SVS <- subset(forRbrul, binaryNorthSouth=="generalNorth" | binaryNorthSouth=="generalSouth")
OW_SVS <- subset(OW_SVS, vowel=="OW")
OW_SVS <- subset(OW_SVS, typeTranscription != "inhouse")
OW_SVS$typeTranscription <- factor(OW_SVS$typeTranscription) 
OW_SVS = na.omit(OW_SVS)

nrow(OW_SVS)
mean(OW_SVS$F2)
hist(OW_SVS$F2)
qqnorm(OW_SVS$F2, pch = 1, frame = FALSE)
qqline(OW_SVS$F2, col = "steelblue", lwd = 2)

OW_SVS_south_ds = subset(OW_SVS, binaryNorthSouth == "generalSouth" & typeTranscription == "deepspeech")
mean_OW_SVS_south_ds = mean(OW_SVS_south_ds$F2)
OW_SVS_south_gt = subset(OW_SVS, binaryNorthSouth == "generalSouth" & typeTranscription == "groundtruth")
mean_OW_SVS_south_gt = mean(OW_SVS_south_gt$F2)
OW_SVS_north_ds = subset(OW_SVS, binaryNorthSouth == "generalNorth" & typeTranscription == "deepspeech")
mean_OW_SVS_north_ds = mean(OW_SVS_north_ds$F2)
OW_SVS_north_gt = subset(OW_SVS, binaryNorthSouth == "generalNorth" & typeTranscription == "groundtruth")
mean_OW_SVS_north_gt = mean(OW_SVS_north_gt$F2)
delta_OW_SVS_ds = abs(mean_OW_SVS_south_ds - mean_OW_SVS_north_ds)
delta_OW_SVS_gt = abs(mean_OW_SVS_south_gt - mean_OW_SVS_north_gt)
delta_OW_SVS_ds
delta_OW_SVS_gt

#== UW_SVS

UW_SVS <- subset(forRbrul, binaryNorthSouth=="generalNorth" | binaryNorthSouth=="generalSouth")
UW_SVS <- subset(UW_SVS, vowel=="UW")
UW_SVS <- subset(UW_SVS, typeTranscription != "inhouse")
UW_SVS$typeTranscription <- factor(UW_SVS$typeTranscription) 
UW_SVS = na.omit(UW_SVS)

nrow(UW_SVS)
mean(UW_SVS$F2)
hist(UW_SVS$F2)
qqnorm(UW_SVS$F2, pch = 1, frame = FALSE)
qqline(UW_SVS$F2, col = "steelblue", lwd = 2)

UW_SVS_south_ds = subset(UW_SVS, binaryNorthSouth == "generalSouth" & typeTranscription == "deepspeech")
mean_UW_SVS_south_ds = mean(UW_SVS_south_ds$F2)
UW_SVS_south_gt = subset(UW_SVS, binaryNorthSouth == "generalSouth" & typeTranscription == "groundtruth")
mean_UW_SVS_south_gt = mean(UW_SVS_south_gt$F2)
UW_SVS_north_ds = subset(UW_SVS, binaryNorthSouth == "generalNorth" & typeTranscription == "deepspeech")
mean_UW_SVS_north_ds = mean(UW_SVS_north_ds$F2)
UW_SVS_north_gt = subset(UW_SVS, binaryNorthSouth == "generalNorth" & typeTranscription == "groundtruth")
mean_UW_SVS_north_gt = mean(UW_SVS_north_gt$F2)
delta_UW_SVS_ds = abs(mean_UW_SVS_south_ds - mean_UW_SVS_north_ds)
delta_UW_SVS_gt = abs(mean_UW_SVS_south_gt - mean_UW_SVS_north_gt)
delta_UW_SVS_ds
delta_UW_SVS_gt


#====================================================================.
#== Fig10 Overlaps regions                                   =====
#====================================================================.

#== Overlap when y-factor is Labov vowel shift

t = EY_SVS
t$yVal = t$EYshifting
t = subset(t, select = -c(EYshifting) )
overlapRegions = t

t = EH_SVS
t$yVal = t$EHshifting
t = subset(t, select = -c(EHshifting) )
overlapRegions = rbind(overlapRegions, t)

t = AE_SVS
t$yVal = t$AEshifting
t = subset(t, select = -c(AEshifting) )
overlapRegions = rbind(overlapRegions, t)

t = IY_SVS
t$yVal = t$IYshifting
t = subset(t, select = -c(IYshifting) )
overlapRegions = rbind(overlapRegions, t)

t = IH_SVS
t$yVal = t$IHshifting
t = subset(t, select = -c(IHshifting) )
overlapRegions = rbind(overlapRegions, t)

overlapRegions$typeTranscription = as.character(overlapRegions$typeTranscription)
overlapRegions$typeTranscription[overlapRegions$typeTranscription == "deepspeech"] = "DS"
overlapRegions$typeTranscription[overlapRegions$typeTranscription == "groundtruth"] = "GT"
overlapRegions$typeTranscription = as.factor(overlapRegions$typeTranscription)

overlapRegions$binaryNorthSouth = as.character(overlapRegions$binaryNorthSouth)
overlapRegions$binaryNorthSouth[overlapRegions$binaryNorthSouth == "generalNorth"] = "General North"
overlapRegions$binaryNorthSouth[overlapRegions$binaryNorthSouth == "generalSouth"] = "South"
overlapRegions$binaryNorthSouth = as.factor(overlapRegions$binaryNorthSouth)

# Export size: 700x350
ggplot(overlapRegions, aes(x=typeTranscription, y=yVal, fill=binaryNorthSouth)) + 
  facet_grid(cols = vars(vowel))+
  geom_boxplot()+
  labs(title="Region and Transcription Type for Southern Vowel Shift",x="", y = "Vowel Shift", fill="Region")


#== Overlap when y-factor is F2

t = AW_SVS
t$yVal = t$F2
overlapRegions = t

t = OW_SVS
t$yVal = t$F2
overlapRegions = rbind(overlapRegions, t)

t = UW_SVS
t$yVal = t$F2
overlapRegions = rbind(overlapRegions, t)

overlapRegions$typeTranscription = as.character(overlapRegions$typeTranscription)
overlapRegions$typeTranscription[overlapRegions$typeTranscription == "deepspeech"] = "DS"
overlapRegions$typeTranscription[overlapRegions$typeTranscription == "groundtruth"] = "GT"
overlapRegions$typeTranscription = as.factor(overlapRegions$typeTranscription)

overlapRegions$binaryNorthSouth = as.character(overlapRegions$binaryNorthSouth)
overlapRegions$binaryNorthSouth[overlapRegions$binaryNorthSouth == "generalNorth"] = "General North"
overlapRegions$binaryNorthSouth[overlapRegions$binaryNorthSouth == "generalSouth"] = "South"
overlapRegions$binaryNorthSouth = as.factor(overlapRegions$binaryNorthSouth)

# Export size: 700x350
ggplot(overlapRegions, aes(x=typeTranscription, y=yVal, fill=binaryNorthSouth)) + 
  facet_grid(cols = vars(vowel))+
  geom_boxplot()+
  labs(title="Transcriptions and region",x="", y = "F2 (Lobanov-normalized)", fill="Region")


with(overlapRegions, tapply(yVal, list(typeTranscription, vowel, binaryNorthSouth), mean))



#====================================================================.
#== Table2 South LMERs                                          =====
#====================================================================.

#== EY_SVS

t = EY_SVS
t$yVal = t$EYshifting
t$yVal = (t$yVal - min(t$yVal)) / (max(t$yVal) - min(t$yVal)) #norm [0,1]
t$yVal = asin(sqrt(t$yVal))

# takes 45 minutes in a high-performance computer (HPC)
# (any larger models won't run; there isn't enough data)
# estimated age is contributing much less (1x10-9, compared to 1x10-4 for the other factors)
print(Sys.time())
m.EY_SVS.withInteractionsInRandom <- lmer(yVal ~ binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type + (1 + binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type | speaker) + (1 + binaryNorthSouth*typeTranscription + fol_type |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
print(Sys.time())
print(summary(m.EY_SVS.withInteractionsInRandom))

# takes 15 hrs and 10 mins in the HPC
step(m.EY_SVS.withInteractionsInRandom,reduce_fixed=FALSE)

Sys.time()
m.EY_SVS.stepDriven <- lmer(yVal ~ binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type + (1 + binaryNorthSouth + typeTranscription + zEstimatedYOB + fol_type | speaker) + (1 + binaryNorthSouth + typeTranscription + fol_type |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))

Sys.time()
summary(m.EY_SVS.stepDriven)

r.squaredGLMM(m.EY_SVS.stepDriven)
emmeans(m.EY_SVS.stepDriven , specs = pairwise ~ binaryNorthSouth:typeTranscription)

m.EY_SVS.minimal <- lmer(yVal ~ binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type + (1 | speaker) + (1 |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(m.EY_SVS.minimal)

# Model selection
# Data driven model has the lowest AIC. It is significantly different from the minimal model and it
# provides more explanatory power (lower AIC) without being significantly different from them maximal model.

anova(m.EY_SVS.withInteractionsInRandom, m.EY_SVS.stepDriven)
anova(m.EY_SVS.withInteractionsInRandom, m.EY_SVS.minimal)
anova(m.EY_SVS.stepDriven, m.EY_SVS.minimal)

AIC(m.EY_SVS.withInteractionsInRandom)
AIC(m.EY_SVS.stepDriven)
AIC(m.EY_SVS.minimal)

# Assumptions:
# Linearity of residuals: the residuals are distributed around the zero line
# Homoskedascity: the residuals have a blob-like distribution
plot(fitted(m.EY_SVS.stepDriven),resid(m.EY_SVS.stepDriven))
abline(h=0, col="blue")
# Normality of residuals: the distribution is unimodal with similar tails on both sides.
# Moreover, there are 5730 datapoints, so, using the Central Limit Theorem, we can 
# assume that the residuals behave in a manner similar to the Normal distribution.
hist(resid(m.EY_SVS.stepDriven))
qqnorm(resid(m.EY_SVS.stepDriven))
qqline(resid(m.EY_SVS.stepDriven), col = "steelblue", lwd = 2)
length(resid(m.EY_SVS.stepDriven))

plot_model(m.EY_SVS.stepDriven, type='diag') # you can ajust type (see package info: ?plot_model)


#== EH_SVS

t = EH_SVS
t$yVal = t$EHshifting
t$yVal = (t$yVal - min(t$yVal)) / (max(t$yVal) - min(t$yVal)) #norm [0,1]
t$yVal = asin(sqrt(t$yVal))

# takes about 52 minutes in the HPC
print(Sys.time())
m.EH_SVS.withInteractionsInRandom <- lmer(yVal ~ binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type + (1 + binaryNorthSouth*typeTranscription + gender + fol_type | speaker) + (1 + binaryNorthSouth*typeTranscription + fol_type |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
print(Sys.time())
summary(m.EH_SVS.withInteractionsInRandom)

# takes about 10 hours in the HPC
step(m.EH_SVS.withInteractionsInRandom,reduce_fixed=FALSE)

# takes about 14 minutes in the HPC
print(Sys.time())
m.EH_SVS.stepDriven <- lmer(yVal ~ binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type + (1 + binaryNorthSouth + typeTranscription + gender + fol_type | speaker) + (1 + binaryNorthSouth*typeTranscription | word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
print(Sys.time())
summary(m.EH_SVS.stepDriven)


r.squaredGLMM(m.EH_SVS.stepDriven)
emmeans(m.EH_SVS.stepDriven , specs = pairwise ~ binaryNorthSouth:typeTranscription)

m.EH_SVS.minimal <- lmer(yVal ~ binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type + (1 | speaker) + (1 |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
print(summary(m.EH_SVS.minimal))

# Model selection
# Data driven model has the lowest AIC. It is significantly different from the minimal model and it
# provides more explanatory power (lower AIC) without being significantly different from them maximal model.

print(anova(m.EH_SVS.withInteractionsInRandom, m.EH_SVS.stepDriven))
print(anova(m.EH_SVS.withInteractionsInRandom, m.EH_SVS.minimal))
print(anova(m.EH_SVS.stepDriven, m.EH_SVS.minimal))
print(AIC(m.EH_SVS.withInteractionsInRandom))
print(AIC(m.EH_SVS.stepDriven))
print(AIC(m.EH_SVS.minimal))

# Assumptions:
# Linearity of residuals: the residuals are distributed around the zero line
# Homoskedascity: the residuals have a blob-like distribution
plot(fitted(m.EH_SVS.stepDriven),resid(m.EH_SVS.stepDriven))
abline(h=0, col="blue")
# Normality of residuals: the distribution is unimodal with similar tails on both sides.
# Moreover, there are 5730 datapoints, so, using the Central Limit Theorem, we can 
# assume that the residuals behave in a manner similar to the Normal distribution.
hist(resid(m.EH_SVS.stepDriven))
qqnorm(resid(m.EH_SVS.stepDriven))
qqline(resid(m.EH_SVS.stepDriven), col = "steelblue", lwd = 2)
length(resid(m.EH_SVS.stepDriven))

plot_model(m.EH_SVS.stepDriven, type='diag') # you can ajust type (see package info: ?plot_model)


#== IY_SVS

t = IY_SVS
t$yVal = t$IYshifting
t$yVal = (t$yVal - min(t$yVal)) / (max(t$yVal) - min(t$yVal)) #norm [0,1]
t$yVal = asin(sqrt(t$yVal))

# takes about 29 minutes in the HPC
print(Sys.time())
m.IY_SVS.withInteractionsInRandom <- lmer(yVal ~ binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type + (1 + binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type | speaker) + (1 + binaryNorthSouth*typeTranscription + fol_type |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
print(Sys.time())
print(summary(m.IY_SVS.withInteractionsInRandom))

# takes about 8 hrs 15 mins in the HPC
step(m.IY_SVS.withInteractionsInRandom,reduce_fixed=FALSE)

# takes about 14 minutes in the HPC
print(Sys.time())
m.IY_SVS.stepDriven <- lmer(yVal ~ binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type + (1 + typeTranscription + gender + fol_type | speaker) + (1 + binaryNorthSouth * typeTranscription + fol_type | word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
print(Sys.time())
summary(m.IY_SVS.stepDriven)


r.squaredGLMM(m.IY_SVS.stepDriven)
emmeans(m.IY_SVS.stepDriven , specs = pairwise ~ binaryNorthSouth:typeTranscription)

m.IY_SVS.minimal <- lmer(yVal ~ binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type + (1 | speaker) + (1 |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(m.IY_SVS.minimal)

# Model selection
anova(m.IY_SVS.withInteractionsInRandom, m.IY_SVS.stepDriven)
anova(m.IY_SVS.withInteractionsInRandom, m.IY_SVS.minimal)
anova(m.IY_SVS.stepDriven, m.IY_SVS.minimal)
AIC(m.IY_SVS.withInteractionsInRandom)
AIC(m.IY_SVS.stepDriven)
AIC(m.IY_SVS.minimal)

# Assumptions:
# Linearity of residuals: the residuals are distributed around the zero line
# Homoskedascity: the residuals have a blob-like distribution
plot(fitted(m.IY_SVS.stepDriven),resid(m.IY_SVS.stepDriven))
abline(h=0, col="blue")
# Normality of residuals: the distribution is unimodal with similar tails on both sides.
# Moreover, there are 5730 datapoints, so, using the Central Limit Theorem, we can 
# assume that the residuals behave in a manner similar to the Normal distribution.
hist(resid(m.IY_SVS.stepDriven))
qqnorm(resid(m.IY_SVS.stepDriven))
qqline(resid(m.IY_SVS.stepDriven), col = "steelblue", lwd = 2)
length(resid(m.IY_SVS.stepDriven))

plot_model(m.IY_SVS.stepDriven, type='diag') # you can ajust type (see package info: ?plot_model)


#== IH_SVS

t = IH_SVS
t$yVal = t$IHshifting
t$yVal = (t$yVal - min(t$yVal)) / (max(t$yVal) - min(t$yVal)) #norm [0,1]
t$yVal = asin(sqrt(t$yVal))

# takes about 40 minutes in the HPC
Sys.time()
m.IH_SVS.withInteractionsInRandom <- lmer(yVal ~ binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type + (1 + binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type | speaker) + (1 + binaryNorthSouth*typeTranscription + fol_type |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
Sys.time()
print(summary(m.IH_SVS.withInteractionsInRandom))

# takes slightly over 12 hours in the HPC
step(m.IH_SVS.withInteractionsInRandom,reduce_fixed=FALSE)

Sys.time()
m.IH_SVS.stepDriven <- lmer(yVal ~ binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type + (1 + typeTranscription + fol_type | speaker) + (1 + binaryNorthSouth + typeTranscription + fol_type | word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
Sys.time()
summary(m.IH_SVS.stepDriven)


r.squaredGLMM(m.IH_SVS.stepDriven)
emmeans(m.IH_SVS.stepDriven , specs = pairwise ~ binaryNorthSouth:typeTranscription)

m.IH_SVS.minimal <- lmer(yVal ~ binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type + (1 | speaker) + (1 |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(m.IH_SVS.minimal)

# Model selection
# Data driven model has the lowest AIC. It is significantly different from the minimal model and it
# provides more explanatory power (lower AIC) without being significantly different from them maximal model.
anova(m.IH_SVS.withInteractionsInRandom, m.IH_SVS.stepDriven)
anova(m.IH_SVS.withInteractionsInRandom, m.IH_SVS.minimal)
anova(m.IH_SVS.stepDriven, m.IH_SVS.minimal)
AIC(m.IH_SVS.withInteractionsInRandom)
AIC(m.IH_SVS.stepDriven)
AIC(m.IH_SVS.minimal)

# Assumptions:
# Linearity of residuals: the residuals are distributed around the zero line
# Homoskedascity: the residuals have a blob-like distribution
plot(fitted(m.IH_SVS.stepDriven),resid(m.IH_SVS.stepDriven))
abline(h=0, col="blue")
# Normality of residuals: the distribution is unimodal with similar tails on both sides.
# Moreover, there are 5730 datapoints, so, using the Central Limit Theorem, we can 
# assume that the residuals behave in a manner similar to the Normal distribution.
hist(resid(m.IH_SVS.stepDriven))
qqnorm(resid(m.IH_SVS.stepDriven))
qqline(resid(m.IH_SVS.stepDriven), col = "steelblue", lwd = 2)
length(resid(m.IH_SVS.stepDriven))

plot_model(m.IH_SVS.stepDriven, type='diag') # you can ajust type (see package info: ?plot_model)


#== AE_SVS

t = AE_SVS
t$yVal = t$AEshifting
t$yVal = (t$yVal - min(t$yVal)) / (max(t$yVal) - min(t$yVal)) #norm [0,1]
t$yVal = asin(sqrt(t$yVal))

# The m.AE_SVS.withInteractionsInRandom model has too many random effects for the model to run
# The necessarily smaller model takes 48 minutes to run
print(Sys.time())
m.AE_SVS.withInteractionsInRandom <- lmer(yVal ~ binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type + (1 + binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type | speaker) + (1 + binaryNorthSouth + typeTranscription + fol_type |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
print(Sys.time())
print(summary(m.AE_SVS.withInteractionsInRandom))

# Takes about 11 hours to run
print(Sys.time())
step(m.AE_SVS.withInteractionsInRandom,reduce_fixed=FALSE)
print(Sys.time())

# Takes 10 minutes to run
print(Sys.time())
m.AE_SVS.stepDriven <- lmer(yVal ~ binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type + (1 + zEstimatedYOB + fol_type | speaker) + (1 + binaryNorthSouth + typeTranscription + fol_type |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))

print(Sys.time())
summary(m.AE_SVS.stepDriven)

r.squaredGLMM(m.AE_SVS.stepDriven)
emmeans(m.AE_SVS.stepDriven , specs = pairwise ~ binaryNorthSouth:typeTranscription)

m.AE_SVS.minimal <- lmer(yVal ~ binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type + (1 | speaker) + (1 |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
print(summary(m.AE_SVS.minimal))


# Model selection
# Data driven model has the lowest AIC. It is significantly different from the minimal model and it
# provides more explanatory power (lower AIC) without being significantly different from them maximal model.
print(anova(m.AE_SVS.withInteractionsInRandom, m.AE_SVS.stepDriven))
print(anova(m.AE_SVS.withInteractionsInRandom, m.AE_SVS.minimal))
print(anova(m.AE_SVS.stepDriven, m.AE_SVS.minimal))
print(AIC(m.AE_SVS.withInteractionsInRandom))
print(AIC(m.AE_SVS.stepDriven))
print(AIC(m.AE_SVS.minimal))

# Assumptions:
# Linearity of residuals: the residuals are distributed around the zero line
# Homoskedascity: the residuals have a blob-like distribution
plot(fitted(m.AE_SVS.stepDriven),resid(m.AE_SVS.stepDriven))
abline(h=0, col="blue")
# Normality of residuals: the distribution is unimodal with similar tails on both sides.
# Moreover, there are 5730 datapoints, so, using the Central Limit Theorem, we can 
# assume that the residuals behave in a manner similar to the Normal distribution.
hist(resid(m.AE_SVS.stepDriven))
qqnorm(resid(m.AE_SVS.stepDriven))
qqline(resid(m.AE_SVS.stepDriven), col = "steelblue", lwd = 2)
length(resid(m.AE_SVS.stepDriven))

plot_model(m.AE_SVS.stepDriven, type='diag') # you can ajust type (see package info: ?plot_model)


#== AW_SVS

t = AW_SVS
t$yVal = t$F2
t$yVal = (t$yVal - min(t$yVal)) / (max(t$yVal) - min(t$yVal)) #norm [0,1]
t$yVal = asin(sqrt(t$yVal))

# takes about 11 minutes in the HPC
Sys.time()
m.AW_SVS.withInteractionsInRandom <- lmer(yVal ~ binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type + (1 + binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type | speaker) + (1 + binaryNorthSouth*typeTranscription + fol_type |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
Sys.time()
summary(m.AW_SVS.withInteractionsInRandom)

# takes about 2.5 hours in the HPC
step(m.AW_SVS.withInteractionsInRandom,reduce_fixed=FALSE)

# takes about 5 minutes in the HPC
Sys.time()
m.AW_SVS.stepDriven <- lmer(yVal ~ binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type + (1 + binaryNorthSouth + fol_type | speaker) + (1 + binaryNorthSouth + typeTranscription + fol_type |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))

Sys.time()
summary(m.AW_SVS.stepDriven)

r.squaredGLMM(m.AW_SVS.stepDriven)
emmeans(m.AW_SVS.stepDriven , specs = pairwise ~ binaryNorthSouth:typeTranscription)

m.AW_SVS.minimal <- lmer(yVal ~ binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type + (1 | speaker) + (1 |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
summary(m.AW_SVS.minimal)

# Model selection
# Data driven model has the lowest AIC. It is significantly different from the minimal model and it
# provides more explanatory power (lower AIC) without being significantly different from them maximal model.
anova(m.AW_SVS.withInteractionsInRandom, m.AW_SVS.stepDriven)
anova(m.AW_SVS.withInteractionsInRandom, m.AW_SVS.minimal)
anova(m.AW_SVS.stepDriven, m.AW_SVS.minimal)
AIC(m.AW_SVS.withInteractionsInRandom)
AIC(m.AW_SVS.stepDriven)
AIC(m.AW_SVS.minimal)

# Assumptions:
# Linearity of residuals: the residuals are distributed around the zero line
# Homoskedascity: the residuals have a blob-like distribution
plot(fitted(m.AW_SVS.stepDriven),resid(m.AW_SVS.stepDriven))
abline(h=0, col="blue")
# Normality of residuals: the distribution is unimodal with similar tails on both sides.
# Moreover, there are 5730 datapoints, so, using the Central Limit Theorem, we can 
# assume that the residuals behave in a manner similar to the Normal distribution.
hist(resid(m.AW_SVS.stepDriven))
qqnorm(resid(m.AW_SVS.stepDriven))
qqline(resid(m.AW_SVS.stepDriven), col = "steelblue", lwd = 2)
length(resid(m.AW_SVS.stepDriven))

plot_model(m.AW_SVS.stepDriven, type='diag') # you can ajust type (see package info: ?plot_model)

#== OW_SVS

t = OW_SVS
t$yVal = t$F2
t$yVal = (t$yVal - min(t$yVal)) / (max(t$yVal) - min(t$yVal)) #norm [0,1]
t$yVal = asin(sqrt(t$yVal))

# takes about 26 minutes
print(Sys.time())
m.OW_SVS.withInteractionsInRandom <- lmer(yVal ~ binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type + (1 + binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type | speaker) + (1 + binaryNorthSouth*typeTranscription + fol_type |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
print(Sys.time())
print(summary(m.OW_SVS.withInteractionsInRandom))


step(m.OW_SVS.withInteractionsInRandom,reduce_fixed=FALSE)


m.OW_SVS.stepDriven <- lmer(yVal ~ binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type + (1 + typeTranscription + fol_type | speaker) + (1 + binaryNorthSouth + typeTranscription | word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
print(Sys.time())
print(summary(m.OW_SVS.stepDriven))


print(r.squaredGLMM(m.OW_SVS.stepDriven))
emmeans(m.OW_SVS.stepDriven , specs = pairwise ~ binaryNorthSouth:typeTranscription)

m.OW_SVS.minimal <- lmer(yVal ~ binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type + (1 | speaker) + (1 |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
print(summary(m.OW_SVS.minimal))

# Model selection
# Data driven model has the lowest AIC. It is significantly different from the minimal model and it
# provides more explanatory power (lower AIC) without being significantly different from them maximal model.
print(anova(m.OW_SVS.withInteractionsInRandom, m.OW_SVS.stepDriven))
print(anova(m.OW_SVS.withInteractionsInRandom, m.OW_SVS.minimal))
print(anova(m.OW_SVS.stepDriven, m.OW_SVS.minimal))
print(AIC(m.OW_SVS.withInteractionsInRandom))
print(AIC(m.OW_SVS.stepDriven))
print(AIC(m.OW_SVS.minimal))

# Assumptions:
# Linearity of residuals: the residuals are distributed around the zero line
# Homoskedascity: the residuals have a blob-like distribution
plot(fitted(m.OW_SVS.stepDriven),resid(m.OW_SVS.stepDriven))
abline(h=0, col="blue")
# Normality of residuals: the distribution is unimodal with similar tails on both sides.
# Moreover, there are 5730 datapoints, so, using the Central Limit Theorem, we can 
# assume that the residuals behave in a manner similar to the Normal distribution.
hist(resid(m.OW_SVS.stepDriven))
qqnorm(resid(m.OW_SVS.stepDriven))
qqline(resid(m.OW_SVS.stepDriven), col = "steelblue", lwd = 2)
length(resid(m.OW_SVS.stepDriven))

plot_model(m.OW_SVS.stepDriven, type='diag') # you can ajust type (see package info: ?plot_model)

#== UW_SVS

t = UW_SVS
t$yVal = t$F2
t$yVal = (t$yVal - min(t$yVal)) / (max(t$yVal) - min(t$yVal)) #norm [0,1]
t$yVal = asin(sqrt(t$yVal))

# takes about 7 minutes
print(Sys.time())
m.UW_SVS.withInteractionsInRandom <- lmer(yVal ~ binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type + (1 + binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type | speaker) + (1 + binaryNorthSouth*typeTranscription + fol_type |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
print(Sys.time())
print(summary(m.UW_SVS.withInteractionsInRandom))

# takes about 70 minutes
print(Sys.time())
step(m.UW_SVS.withInteractionsInRandom,reduce_fixed=FALSE)
print(Sys.time())

print(Sys.time())
m.UW_SVS.stepDriven <- lmer(yVal ~ binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type + (1 + fol_type | speaker) + (1 + binaryNorthSouth |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
print(Sys.time())
print(summary(m.UW_SVS.stepDriven))


print(r.squaredGLMM(m.UW_SVS.dataDriven))
emmeans(m.UW_SVS.stepDriven , specs = pairwise ~ binaryNorthSouth:typeTranscription)

m.UW_SVS.minimal <- lmer(yVal ~ binaryNorthSouth*typeTranscription + zEstimatedYOB + gender + fol_type + (1 | speaker) + (1 |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
print(summary(m.UW_SVS.minimal))

# Model selection
# Data driven model has the lowest AIC. It is significantly different from the minimal model and it
# provides more explanatory power (lower AIC) without being significantly different from them maximal model.
anova(m.UW_SVS.withInteractionsInRandom, m.UW_SVS.stepDriven)
anova(m.UW_SVS.withInteractionsInRandom, m.UW_SVS.minimal)
anova(m.UW_SVS.stepDriven, m.UW_SVS.minimal)
AIC(m.UW_SVS.withInteractionsInRandom)
AIC(m.UW_SVS.stepDriven)
AIC(m.UW_SVS.minimal)

# Assumptions:
# Linearity of residuals: the residuals are distributed around the zero line
# Homoskedascity: the residuals have a blob-like distribution
plot(fitted(m.UW_SVS.stepDriven),resid(m.UW_SVS.stepDriven))
abline(h=0, col="blue")
# Normality of residuals: the distribution is unimodal with similar tails on both sides.
# Moreover, there are 5730 datapoints, so, using the Central Limit Theorem, we can 
# assume that the residuals behave in a manner similar to the Normal distribution.
hist(resid(m.UW_SVS.stepDriven))
qqnorm(resid(m.UW_SVS.stepDriven))
qqline(resid(m.UW_SVS.stepDriven), col = "steelblue", lwd = 2)
length(resid(m.UW_SVS.stepDriven))

plot_model(m.UW_SVS.stepDriven, type='diag') # you can ajust type (see package info: ?plot_model)



#================================================================================.
#== Fig11 North Vowel Triangle (Ground Truth, GenNorth and InlandNorth)     =====
#====================================================================.

# Inland North versus General North, ground truth

x <- subset(mydata, mydata$compareNorths =="inlandNorth" |
              mydata$compareNorths == "northWithoutInlandNorth")
x <- subset(x, x$typeTranscription == "groundtruth")
nrow(x)

x$compareNorths = as.character(x$compareNorths)
x$compareNorths[x$compareNorths == "inlandNorth"] = "Inland North"
x$compareNorths[x$compareNorths == "northWithoutInlandNorth"] = "General North"
threeColorChoice = c("#0000FF", "#FF0000", "#00FF00", "#FF0000", "#FF3D00", "#FF7A00", "#FF0000", "#FF3D00", "#FF7A00")

to.plot <- c()
to.plot$speaker <- as.character(x$compareNorths)
to.plot$vowel <- as.character(x$vowel)
to.plot$context <- "NA"
to.plot$F1 <- as.numeric(x$F1_LobanovNormed_unscaled)
to.plot$F2 <- as.numeric(x$F2_LobanovNormed_unscaled )
to.plot$F3 <- "NA"
to.plot$gl.F1 <- "NA"
to.plot$gl.F2 <- "NA"
to.plot$gl.F3 <- "NA"
to.plot <- data.frame(to.plot)

# Export size: 600 x 600
vowelplot(compute.means(to.plot, separate=T), color="speakers", color.choice=threeColorChoice, size=0, label="vowels", title="Inland North versus General North\nGround truth only")


#====================================================================.
#== Fig12 North Vowel Triangle (DS, Sphinx)                     =====
#====================================================================.

# Inland North versus General North, Deep Speech

x <- subset(mydata, mydata$compareNorths =="inlandNorth" |
              mydata$compareNorths == "northWithoutInlandNorth")
x <- subset(x, x$typeTranscription == "deepspeech")
nrow(x)

x$compareNorths = as.character(x$compareNorths)
x$compareNorths[x$compareNorths == "inlandNorth"] = "Inland North"
x$compareNorths[x$compareNorths == "northWithoutInlandNorth"] = "General North"
threeColorChoice = c("#0000FF", "#FF0000", "#00FF00", "#FF0000", "#FF3D00", "#FF7A00", "#FF0000", "#FF3D00", "#FF7A00")

to.plot <- c()
to.plot$speaker <- as.character(x$compareNorths)
to.plot$vowel <- as.character(x$vowel)
to.plot$context <- "NA"
to.plot$F1 <- as.numeric(x$F1_LobanovNormed_unscaled)
to.plot$F2 <- as.numeric(x$F2_LobanovNormed_unscaled )
to.plot$F3 <- "NA"
to.plot$gl.F1 <- "NA"
to.plot$gl.F2 <- "NA"
to.plot$gl.F3 <- "NA"
to.plot <- data.frame(to.plot)

# Export size: 600 x 600
vowelplot(compute.means(to.plot, separate=T), color="speakers", color.choice=threeColorChoice, size=0, label="vowels", title="Inland North versus General North\nDeepSpeech only")


# Inland North versus General North, Inhouse

x <- subset(mydata, mydata$compareNorths =="inlandNorth" |
              mydata$compareNorths == "northWithoutInlandNorth")
x <- subset(x, x$typeTranscription == "inhouse")
nrow(x)

x$compareNorths = as.character(x$compareNorths)
x$compareNorths[x$compareNorths == "inlandNorth"] = "Inland North"
x$compareNorths[x$compareNorths == "northWithoutInlandNorth"] = "General North"
threeColorChoice = c("#0000FF", "#FF0000", "#00FF00", "#FF0000", "#FF3D00", "#FF7A00", "#FF0000", "#FF3D00", "#FF7A00")

to.plot <- c()
to.plot$speaker <- as.character(x$compareNorths)
to.plot$vowel <- as.character(x$vowel)
to.plot$context <- "NA"
to.plot$F1 <- as.numeric(x$F1_LobanovNormed_unscaled)
to.plot$F2 <- as.numeric(x$F2_LobanovNormed_unscaled )
to.plot$F3 <- "NA"
to.plot$gl.F1 <- "NA"
to.plot$gl.F2 <- "NA"
to.plot$gl.F3 <- "NA"
to.plot <- data.frame(to.plot)

# Export size: 600 x 600
vowelplot(compute.means(to.plot, separate=T), color="speakers", color.choice=threeColorChoice, size=0, label="vowels", title="Inland North versus General North\nSphinx only")


#====================================================================.
#== Fig13 North Vowel Triangle (All types)                      =====
#====================================================================.

# Choose this for inland!
x <- subset(mydata, mydata$compareNorths == "inlandNorth")
figTitle = "Inland North vowels in\nall three transcription types"
# Choose this for general!
x <- subset(mydata, mydata$compareNorths == "northWithoutInlandNorth")
figTitle = "General North vowels in\nall three transcription types"

nrow(x)

x$typeTranscription = as.character(x$typeTranscription)
x$typeTranscription[x$typeTranscription=="groundtruth"] = "Ground truth"
x$typeTranscription[x$typeTranscription=="inhouse"] = "Sphinx"
x$typeTranscription[x$typeTranscription=="deepspeech"] = "DeepSpeech"
threeColorChoice = c("#0000FF", "#FF0000", "#00FF00", "#FF0000", "#FF3D00", "#FF7A00", "#FF0000", "#FF3D00", "#FF7A00")

to.plot <- c()
to.plot$speaker <- as.character(x$typeTranscription)
to.plot$vowel <- as.character(x$vowel)
to.plot$context <- "NA"
to.plot$F1 <- as.numeric(x$F1_LobanovNormed_unscaled)
to.plot$F2 <- as.numeric(x$F2_LobanovNormed_unscaled )
to.plot$F3 <- "NA"
to.plot$gl.F1 <- "NA"
to.plot$gl.F2 <- "NA"
to.plot$gl.F3 <- "NA"
to.plot <- data.frame(to.plot)

# size is 600x600
vowelplot(compute.means(to.plot, separate=T), color="speakers", color.choice=threeColorChoice, size=0, label="vowels", title=figTitle)


#====================================================================.
#== Fig14 BhatDist North (Bhattacharyya distance)               =====
#====================================================================.

batNorth = subset(batDist, anae=="generalNorth")

tBatNorth = batNorth
tBatNorth$typeTranscription = as.character(tBatNorth$typeTranscription)
tBatNorth$typeTranscription[tBatNorth$typeTranscription=="deepspeech"] = "DeepSpeech"
tBatNorth$typeTranscription[tBatNorth$typeTranscription=="inhouse"] = "Sphinx"
tBatNorth$typeTranscription = as.factor(tBatNorth$typeTranscription)

# size: 1200 x 350
ggplot(tBatNorth, aes(x=vowel, y=bat, fill=typeTranscription)) +
  geom_boxplot()+
  labs(title="Bhattacharyya's Affinity for Inland North vowels by type of transcription",x="Vowel", y = "Bhattacharyya's Affinity", fill="Transcription")+
  ylim(0.6, 1)


#====================================================================.
#== Table3 BhatDist North (Bhattacharyya distance)                =====
#====================================================================.

medBhatDistNorth = with(batNorth, tapply(bat, list(typeTranscription, vowel), median))
medBhatDistNorth = t(medBhatDistNorth)
medBhatDistNorth = data.frame(medBhatDistNorth)
medBhatDistNorth$diff = medBhatDistNorth$deepspeech - medBhatDistNorth$inhouse

medBhatDistNorth <- medBhatDistNorth[order(medBhatDistNorth$diff),]
medBhatDistNorth


#====================================================================.
#== Fig15 Table4 Data for North vowels                                 ====
#== NCS: Northern Cities Shift
#====================================================================.

AE_NCS <- subset(forRbrul, compareNorths=="inlandNorth" | compareNorths=="northWithoutInlandNorth")
AE_NCS <- subset(AE_NCS, vowel=="AE")
AE_NCS <- subset(AE_NCS, typeTranscription != "inhouse")
AE_NCS$typeTranscription <- factor(AE_NCS$typeTranscription)
# calculate AEraising using Labov's formula F2- 2 x F1
AE_NCS$AEraising <- AE_NCS$F2 - 2*AE_NCS$F1
AE_NCS = na.omit(AE_NCS)

AE_NCS_inland_ds = subset(AE_NCS, compareNorths == "inlandNorth" & typeTranscription == "deepspeech")
mean_AE_NCS_inland_ds = mean(AE_NCS_inland_ds$AEraising)
AE_NCS_inland_gt = subset(AE_NCS, compareNorths == "inlandNorth" & typeTranscription == "groundtruth")
mean_AE_NCS_inland_gt = mean(AE_NCS_inland_gt$AEraising)
AE_NCS_north_ds = subset(AE_NCS, compareNorths == "northWithoutInlandNorth" & typeTranscription == "deepspeech")
mean_AE_NCS_north_ds = mean(AE_NCS_north_ds$AEraising)
AE_NCS_north_gt = subset(AE_NCS, compareNorths == "northWithoutInlandNorth" & typeTranscription == "groundtruth")
mean_AE_NCS_north_gt = mean(AE_NCS_north_gt$AEraising)
delta_AE_NCS_ds = abs(mean_AE_NCS_inland_ds - mean_AE_NCS_north_ds)
delta_AE_NCS_gt = abs(mean_AE_NCS_inland_gt - mean_AE_NCS_north_gt)
delta_AE_NCS_ds
delta_AE_NCS_gt




EH_NCS <- subset(forRbrul, compareNorths=="inlandNorth" | compareNorths=="northWithoutInlandNorth")
EH_NCS <- subset(EH_NCS, vowel=="EH")
EH_NCS <- subset(EH_NCS, typeTranscription != "inhouse")
EH_NCS$typeTranscription <- factor(EH_NCS$typeTranscription)
# calculate EHraising using Labov's formula F2- 2 x F1
EH_NCS$EHraising <- EH_NCS$F2 - 2*EH_NCS$F1
EH_NCS = na.omit(EH_NCS)

EH_NCS_inland_ds = subset(EH_NCS, compareNorths == "inlandNorth" & typeTranscription == "deepspeech")
mean_EH_NCS_inland_ds = mean(EH_NCS_inland_ds$EHraising)
EH_NCS_inland_gt = subset(EH_NCS, compareNorths == "inlandNorth" & typeTranscription == "groundtruth")
mean_EH_NCS_inland_gt = mean(EH_NCS_inland_gt$EHraising)
EH_NCS_north_ds = subset(EH_NCS, compareNorths == "northWithoutInlandNorth" & typeTranscription == "deepspeech")
mean_EH_NCS_north_ds = mean(EH_NCS_north_ds$EHraising)
EH_NCS_north_gt = subset(EH_NCS, compareNorths == "northWithoutInlandNorth" & typeTranscription == "groundtruth")
mean_EH_NCS_north_gt = mean(EH_NCS_north_gt$EHraising)
delta_EH_NCS_ds = abs(mean_EH_NCS_inland_ds - mean_EH_NCS_north_ds)
delta_EH_NCS_gt = abs(mean_EH_NCS_inland_gt - mean_EH_NCS_north_gt)
delta_EH_NCS_ds
sd_EH_NCS_general_gt = sd(EH_NCS_north_gt$EHraising)
sd_EH_NCS_general_ds = sd(EH_NCS_north_ds$EHraising)
sd_EH_NCS_general_gt
sd_EH_NCS_general_ds
sd_EH_NCS_inland_gt = sd(EH_NCS_inland_gt$EHraising)
sd_EH_NCS_inland_ds = sd(EH_NCS_inland_ds$EHraising)
sd_EH_NCS_inland_gt
sd_EH_NCS_inland_ds


AO_NCS <- subset(forRbrul, compareNorths=="inlandNorth" | compareNorths=="northWithoutInlandNorth")
AO_NCS <- subset(AO_NCS, vowel=="AO")
AO_NCS <- subset(AO_NCS, typeTranscription != "inhouse")
AO_NCS$typeTranscription <- factor(AO_NCS$typeTranscription)
AO_NCS = na.omit(AO_NCS)

AO_NCS_inland_ds = subset(AO_NCS, compareNorths == "inlandNorth" & typeTranscription == "deepspeech")
mean_AO_NCS_inland_ds = mean(AO_NCS_inland_ds$F2)
AO_NCS_inland_gt = subset(AO_NCS, compareNorths == "inlandNorth" & typeTranscription == "groundtruth")
mean_AO_NCS_inland_gt = mean(AO_NCS_inland_gt$F2)
AO_NCS_north_ds = subset(AO_NCS, compareNorths == "northWithoutInlandNorth" & typeTranscription == "deepspeech")
mean_AO_NCS_north_ds = mean(AO_NCS_north_ds$F2)
AO_NCS_north_gt = subset(AO_NCS, compareNorths == "northWithoutInlandNorth" & typeTranscription == "groundtruth")
mean_AO_NCS_north_gt = mean(AO_NCS_north_gt$F2)
delta_AO_NCS_ds = abs(mean_AO_NCS_inland_ds - mean_AO_NCS_north_ds)
delta_AO_NCS_gt = abs(mean_AO_NCS_inland_gt - mean_AO_NCS_north_gt)
delta_AO_NCS_ds
delta_AO_NCS_gt


AA_NCS <- subset(forRbrul, compareNorths=="inlandNorth" | compareNorths=="northWithoutInlandNorth")
AA_NCS <- subset(AA_NCS, vowel=="AA")
AA_NCS <- subset(AA_NCS, typeTranscription != "inhouse")
AA_NCS$typeTranscription <- factor(AA_NCS$typeTranscription)
AA_NCS = na.omit(AA_NCS)

AA_NCS_inland_ds = subset(AA_NCS, compareNorths == "inlandNorth" & typeTranscription == "deepspeech")
mean_AA_NCS_inland_ds = mean(AA_NCS_inland_ds$F2)
AA_NCS_inland_gt = subset(AA_NCS, compareNorths == "inlandNorth" & typeTranscription == "groundtruth")
mean_AA_NCS_inland_gt = mean(AA_NCS_inland_gt$F2)
AA_NCS_north_ds = subset(AA_NCS, compareNorths == "northWithoutInlandNorth" & typeTranscription == "deepspeech")
mean_AA_NCS_north_ds = mean(AA_NCS_north_ds$F2)
AA_NCS_north_gt = subset(AA_NCS, compareNorths == "northWithoutInlandNorth" & typeTranscription == "groundtruth")
mean_AA_NCS_north_gt = mean(AA_NCS_north_gt$F2)
delta_AA_NCS_ds = abs(mean_AA_NCS_inland_ds - mean_AA_NCS_north_ds)
delta_AA_NCS_gt = abs(mean_AA_NCS_inland_gt - mean_AA_NCS_north_gt)
delta_AA_NCS_ds
delta_AA_NCS_gt


AH_NCS <- subset(forRbrul, compareNorths=="inlandNorth" | compareNorths=="northWithoutInlandNorth")
AH_NCS <- subset(AH_NCS, vowel=="AH")
AH_NCS <- subset(AH_NCS, typeTranscription != "inhouse")
AH_NCS$typeTranscription <- factor(AH_NCS$typeTranscription)
AH_NCS = na.omit(AH_NCS)

AH_NCS_inland_ds = subset(AH_NCS, compareNorths == "inlandNorth" & typeTranscription == "deepspeech")
mean_AH_NCS_inland_ds = mean(AH_NCS_inland_ds$F2)
AH_NCS_inland_gt = subset(AH_NCS, compareNorths == "inlandNorth" & typeTranscription == "groundtruth")
mean_AH_NCS_inland_gt = mean(AH_NCS_inland_gt$F2)
AH_NCS_north_ds = subset(AH_NCS, compareNorths == "northWithoutInlandNorth" & typeTranscription == "deepspeech")
mean_AH_NCS_north_ds = mean(AH_NCS_north_ds$F2)
AH_NCS_north_gt = subset(AH_NCS, compareNorths == "northWithoutInlandNorth" & typeTranscription == "groundtruth")
mean_AH_NCS_north_gt = mean(AH_NCS_north_gt$F2)
delta_AH_NCS_ds = abs(mean_AH_NCS_inland_ds - mean_AH_NCS_north_ds)
delta_AH_NCS_gt = abs(mean_AH_NCS_inland_gt - mean_AH_NCS_north_gt)
delta_AH_NCS_ds
delta_AH_NCS_gt


#====================================================================.
#== Fig15 Overlap regions                                       =====
#====================================================================.

#== Overlap when y-factor is Labov vowel shift

t = AE_NCS
t$yVal = t$AEraising
t = subset(t, select = -c(AEraising) )
overlapRegions = t

t = EH_NCS
t$yVal = t$EHraising
t = subset(t, select = -c(EHraising) )
overlapRegions = rbind(overlapRegions, t)


overlapRegions$typeTranscription = as.character(overlapRegions$typeTranscription)
overlapRegions$typeTranscription[overlapRegions$typeTranscription == "deepspeech"] = "DS"
overlapRegions$typeTranscription[overlapRegions$typeTranscription == "groundtruth"] = "GT"
overlapRegions$typeTranscription = as.factor(overlapRegions$typeTranscription)

overlapRegions$compareNorths = as.character(overlapRegions$compareNorths)
overlapRegions$compareNorths[overlapRegions$compareNorths == "northWithoutInlandNorth"] = "General North"
overlapRegions$compareNorths[overlapRegions$compareNorths == "inlandNorth"] = "Inland North"
overlapRegions$compareNorths = as.factor(overlapRegions$compareNorths)

# Export size: 700x350
ggplot(overlapRegions, aes(x=typeTranscription, y=yVal, fill=compareNorths)) + 
  facet_grid(cols = vars(vowel))+
  geom_boxplot()+
  labs(title="Region and Transcription Type for Northern Cities Shift",x="", y = "Vowel Shift", fill="Region")


#== Overlap when y-factor is F2

t = AA_NCS
t$yVal = t$F2
overlapRegions = t

t = AH_NCS
t$yVal = t$F2
overlapRegions = rbind(overlapRegions, t)

t = AO_NCS
t$yVal = t$F2
overlapRegions = rbind(overlapRegions, t)

overlapRegions$typeTranscription = as.character(overlapRegions$typeTranscription)
overlapRegions$typeTranscription[overlapRegions$typeTranscription == "deepspeech"] = "DS"
overlapRegions$typeTranscription[overlapRegions$typeTranscription == "groundtruth"] = "GT"
overlapRegions$typeTranscription = as.factor(overlapRegions$typeTranscription)

overlapRegions$compareNorths = as.character(overlapRegions$compareNorths)
overlapRegions$compareNorths[overlapRegions$compareNorths == "northWithoutInlandNorth"] = "General North"
overlapRegions$compareNorths[overlapRegions$compareNorths == "inlandNorth"] = "Inland North"
overlapRegions$compareNorths = as.factor(overlapRegions$compareNorths)

# Export size: 700x350
ggplot(overlapRegions, aes(x=typeTranscription, y=yVal, fill=compareNorths)) + 
  facet_grid(cols = vars(vowel))+
  geom_boxplot()+
  labs(title="Region and Transcription Type for Northern Cities Shift",x="", y = "F2 (Lobanov-normalized)", fill="Transcription")




#====================================================================.
#== Table4 North LMERs                                          =====
#====================================================================.

#== AE_NCS

t = AE_NCS
t$yVal = t$AEraising
t$yVal = (t$yVal - min(t$yVal)) / (max(t$yVal) - min(t$yVal)) #norm [0,1]
t$yVal = asin(sqrt(t$yVal))

# takes about 16 minutes
print(Sys.time())
m.AE_NCS.withInteractionsInRandom <- lmer(yVal ~ compareNorths*typeTranscription + zEstimatedYOB + gender + fol_type + (1 + compareNorths*typeTranscription + zEstimatedYOB + gender + fol_type | speaker) + (1 + compareNorths + fol_type |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
print(summary(m.AE_NCS.withInteractionsInRandom))
print(Sys.time())

# takes about 1 hr 55 minutes
step(m.AE_NCS.withInteractionsInRandom,reduce_fixed=FALSE)
print(Sys.time())

print(Sys.time())
m.AE_NCS.stepDriven <- lmer(yVal ~ compareNorths*typeTranscription + zEstimatedYOB + gender + fol_type + (1 + fol_type| speaker) + (1 + fol_type |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
print(Sys.time())
print(summary(m.AE_NCS.stepDriven))

r.squaredGLMM(m.AE_NCS.stepDriven)
emmeans(m.AE_NCS.stepDriven , specs = pairwise ~ compareNorths:typeTranscription)

m.AE_NCS.minimal <- lmer(yVal ~ compareNorths*typeTranscription + zEstimatedYOB + gender + fol_type + (1 | speaker) + (1 |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
print(summary(m.AE_NCS.minimal))

# Model selection
# Data driven model has the lowest AIC. It is significantly different from the minimal model and it
# provides more explanatory power (lower AIC) without being significantly different from them maximal model.
print(anova(m.AE_NCS.withInteractionsInRandom, m.AE_NCS.stepDriven))
print(anova(m.AE_NCS.withInteractionsInRandom, m.AE_NCS.minimal))
print(anova(m.AE_NCS.stepDriven, m.AE_NCS.minimal))

print(AIC(m.AE_NCS.withInteractionsInRandom))
print(AIC(m.AE_NCS.dataDriven))
print(AIC(m.AE_NCS.minimal))

# Assumptions:
# Linearity of residuals: the residuals are distributed around the zero line
# Homoskedascity: the residuals have a blob-like distribution
plot(fitted(m.AE_NCS.dataDriven),resid(m.AE_NCS.dataDriven))
abline(h=0, col="blue")
# Normality of residuals: the distribution is unimodal with similar tails on both sides.
# Moreover, there are 5730 datapoints, so, using the Central Limit Theorem, we can 
# assume that the residuals behave in a manner similar to the Normal distribution.
hist(resid(m.AE_NCS.dataDriven))
qqnorm(resid(m.AE_NCS.dataDriven))
qqline(resid(m.AE_NCS.dataDriven), col = "steelblue", lwd = 2)
length(resid(m.AE_NCS.dataDriven))

plot_model(m.AE_NCS.dataDriven, type='diag') # you can ajust type (see package info: ?plot_model)


#== EH_NCS

t = EH_NCS
t$yVal = t$EHraising
t$yVal = (t$yVal - min(t$yVal)) / (max(t$yVal) - min(t$yVal)) #norm [0,1]
t$yVal = asin(sqrt(t$yVal))

# Takes about 18 minutes
print(Sys.time())
m.EH_NCS.withInteractionsInRandom <- lmer(yVal ~ compareNorths*typeTranscription + zEstimatedYOB + gender + fol_type + (1 + compareNorths*typeTranscription + zEstimatedYOB + gender + fol_type | speaker) + (1 + compareNorths + fol_type |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
print(summary(m.EH_NCS.withInteractionsInRandom))
print(Sys.time())

# takes about 3 hours 10 minutes
print(Sys.time())
step(m.EH_NCS.withInteractionsInRandom,reduce_fixed=FALSE)
print(Sys.time())

print(Sys.time())
m.EH_NCS.stepDriven <- lmer(yVal ~ compareNorths*typeTranscription + zEstimatedYOB + gender + fol_type + (1 + compareNorths + gender + fol_type | speaker) + (1 + fol_type |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
print(Sys.time())
print(summary(m.EH_NCS.stepDriven))

r.squaredGLMM(m.EH_NCS.stepDriven)
emmeans(m.EH_NCS.stepDriven , specs = pairwise ~ compareNorths:typeTranscription)

m.EH_NCS.minimal <- lmer(yVal ~ compareNorths*typeTranscription + zEstimatedYOB + gender + fol_type + (1 | speaker) + (1 |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
print(summary(m.EH_NCS.minimal))

# Model selection
# Data driven model has the lowest AIC. It is significantly different from the minimal model and it
# provides more explanatory power (lower AIC) without being significantly different from them maximal model.
print(anova(m.EH_NCS.withInteractionsInRandom, m.EH_NCS.stepDriven))
print(anova(m.EH_NCS.withInteractionsInRandom, m.EH_NCS.minimal))
print(anova( m.EH_NCS.stepDriven, m.EH_NCS.minimal))

print(AIC(m.EH_NCS.withInteractionsInRandom))
print(AIC(m.EH_NCS.stepDriven))
print(AIC(m.EH_NCS.minimal))

# Assumptions:
# Linearity of residuals: the residuals are distributed around the zero line
# Homoskedascity: the residuals have a blob-like distribution
plot(fitted(m.EH_NCS.stepDriven),resid(m.EH_NCS.stepDriven))
abline(h=0, col="blue")
# Normality of residuals: the distribution is unimodal with similar tails on both sides.
# Moreover, there are 5730 datapoints, so, using the Central Limit Theorem, we can 
# assume that the residuals behave in a manner similar to the Normal distribution.
hist(resid(m.EH_NCS.stepDriven))
qqnorm(resid(m.EH_NCS.stepDriven))
qqline(resid(m.EH_NCS.stepDriven), col = "steelblue", lwd = 2)
length(resid(m.EH_NCS.stepDriven))

plot_model(m.EH_NCS.stepDriven, type='diag') # you can ajust type (see package info: ?plot_model)


#== AO_NCS

t = AO_NCS
t$yVal = t$F2
t$yVal = (t$yVal - min(t$yVal)) / (max(t$yVal) - min(t$yVal)) #norm [0,1]
t$yVal = asin(sqrt(t$yVal))


# The m.AO_NCS.withInteractionsInRandom model has too many random effects for the model to run
# The necessarily smaller model takes about a minute to run
print(Sys.time())
#m.AO_NCS.withInteractionsInRandom <- lmer(yVal ~ compareNorths*typeTranscription + zEstimatedYOB + gender + fol_type + (1 + compareNorths*typeTranscription + zEstimatedYOB + gender + fol_type | speaker) + (1 + compareNorths + fol_type |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
m.AO_NCS.withInteractionsInRandomButSmaller <- lmer(yVal ~ compareNorths*typeTranscription + zEstimatedYOB + gender + fol_type + (1 + compareNorths +typeTranscription + fol_type | speaker) + (1 + compareNorths + fol_type |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
print(summary(m.AO_NCS.withInteractionsInRandomButSmaller))
print(Sys.time())

# Takes about 3 minutes
print(Sys.time())
step(m.AO_NCS.withInteractionsInRandomButSmaller,reduce_fixed=FALSE)
print(Sys.time())

print(Sys.time())
m.AO_NCS.stepDriven <- lmer(yVal ~ compareNorths*typeTranscription + zEstimatedYOB + gender + fol_type + (1 + fol_type | speaker) + (1 + compareNorths + fol_type |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
print(Sys.time())
print(summary(m.AO_NCS.stepDriven))

print(r.squaredGLMM(m.AO_NCS.stepDriven))
emmeans(m.AO_NCS.stepDriven , specs = pairwise ~ compareNorths:typeTranscription)

m.AO_NCS.minimal <- lmer(yVal ~ compareNorths*typeTranscription + zEstimatedYOB + gender + fol_type + (1 | speaker) + (1 |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
print(summary(m.AO_NCS.minimal))

# Model selection
# Data driven model has the lowest AIC. It is significantly different from the minimal model and it
# provides more explanatory power (lower AIC) without being significantly different from them maximal model.
print(anova(m.AO_NCS.withInteractionsInRandomButSmaller, m.AO_NCS.stepDriven))
print(anova(m.AO_NCS.withInteractionsInRandomButSmaller, m.AO_NCS.minimal))
print(anova(m.AO_NCS.stepDriven, m.AO_NCS.minimal))

print(AIC(m.AO_NCS.withInteractionsInRandomButSmaller))
print(AIC(m.AO_NCS.stepDriven))
print(AIC(m.AO_NCS.minimal))

# Assumptions:
# Linearity of residuals: the residuals are distributed around the zero line
# Homoskedascity: the residuals have a blob-like distribution
plot(fitted(m.AO_NCS.stepDriven),resid(m.AO_NCS.stepDriven))
abline(h=0, col="blue")
# Normality of residuals: the distribution is unimodal with similar tails on both sides.
# Moreover, there are 5730 datapoints, so, using the Central Limit Theorem, we can 
# assume that the residuals behave in a manner similar to the Normal distribution.
hist(resid(m.AO_NCS.stepDriven))
qqnorm(resid(m.AO_NCS.stepDriven))
qqline(resid(m.AO_NCS.stepDriven), col = "steelblue", lwd = 2)
length(resid(m.AO_NCS.stepDriven))

plot_model(m.AO_NCS.stepDriven, type='diag') # you can ajust type (see package info: ?plot_model)


#== AA_NCS

t = AA_NCS
t$yVal = t$F2
t$yVal = (t$yVal - min(t$yVal)) / (max(t$yVal) - min(t$yVal)) #norm [0,1]
t$yVal = asin(sqrt(t$yVal))

# takes about 7 minutes
print(Sys.time())
m.AA_NCS.withInteractionsInRandom <- lmer(yVal ~ compareNorths*typeTranscription + zEstimatedYOB + gender + fol_type + (1 + compareNorths*typeTranscription + zEstimatedYOB + gender + fol_type | speaker) + (1 + compareNorths + fol_type |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
print(summary(m.AA_NCS.withInteractionsInRandom))
print(Sys.time())

# takes about 33 minutes
print(Sys.time())
step(m.AA_NCS.withInteractionsInRandom,reduce_fixed=FALSE)
print(Sys.time())

# takes about 1 minute
print(Sys.time())
m.AA_NCS.stepDriven <- lmer(yVal ~ compareNorths*typeTranscription + zEstimatedYOB + gender + fol_type + (1 + typeTranscription + gender + fol_type | speaker) + (1 + fol_type |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
print(Sys.time())
print(summary(m.AA_NCS.stepDriven))

print(r.squaredGLMM(m.AA_NCS.stepDriven))
emmeans(m.AA_NCS.stepDriven , specs = pairwise ~ compareNorths:typeTranscription)

m.AA_NCS.minimal <- lmer(yVal ~ compareNorths*typeTranscription + zEstimatedYOB + gender + fol_type + (1 | speaker) + (1 |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
print(summary(m.AA_NCS.minimal))

emmeans(m.AA_NCS.minimal , specs = pairwise ~ compareNorths:typeTranscription)


# Model selection
# Data driven model has the lowest AIC. It is significantly different from the minimal model and it
# provides more explanatory power (lower AIC) without being significantly different from them maximal model.
print(anova(m.AA_NCS.withInteractionsInRandom, m.AA_NCS.stepDriven))
print(anova(m.AA_NCS.withInteractionsInRandom, m.AA_NCS.minimal))
print(anova(m.AA_NCS.stepDriven, m.AA_NCS.minimal))

print(AIC(m.AA_NCS.withInteractionsInRandom))
print(AIC(m.AA_NCS.stepDriven))
print(AIC(m.AA_NCS.minimal))

# Assumptions:
# Linearity of residuals: the residuals are distributed around the zero line
# Homoskedascity: the residuals have a blob-like distribution
plot(fitted(m.AA_NCS.stepDriven),resid(m.AA_NCS.stepDriven))
abline(h=0, col="blue")
# Normality of residuals: the distribution is unimodal with similar tails on both sides.
# Moreover, there are 5730 datapoints, so, using the Central Limit Theorem, we can 
# assume that the residuals behave in a manner similar to the Normal distribution.
hist(resid(m.AA_NCS.stepDriven))
qqnorm(resid(m.AA_NCS.stepDriven))
qqline(resid(m.AA_NCS.stepDriven), col = "steelblue", lwd = 2)
length(resid(m.AA_NCS.stepDriven))

plot_model(m.AA_NCS.stepDriven, type='diag') # you can ajust type (see package info: ?plot_model)


#== AH_NCS

t = AH_NCS
t$yVal = t$F2
t$yVal = (t$yVal - min(t$yVal)) / (max(t$yVal) - min(t$yVal)) #norm [0,1]
t$yVal = asin(sqrt(t$yVal))

# takes approximately 6 minutes
print(Sys.time())
m.AH_NCS.withInteractionsInRandom <- lmer(yVal ~ compareNorths*typeTranscription + zEstimatedYOB + gender + fol_type + (1 + compareNorths*typeTranscription + zEstimatedYOB + gender + fol_type | speaker) + (1 + compareNorths + fol_type |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
print(summary(m.AH_NCS.withInteractionsInRandom))
print(Sys.time())

# takes approximately 50 minutes
print(Sys.time())
step(m.AH_NCS.withInteractionsInRandom,reduce_fixed=FALSE)
print(Sys.time())

print(Sys.time())
m.AH_NCS.stepDriven <- lmer(yVal ~ compareNorths*typeTranscription + zEstimatedYOB + gender + fol_type + (1+ gender + fol_type | speaker) + (1 + compareNorths +fol_type |word), data=t, REML=F, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
print(Sys.time())
print(summary(m.AH_NCS.stepDriven))

print(r.squaredGLMM(m.AH_NCS.stepDriven))
emmeans(m.AH_NCS.stepDriven , specs = pairwise ~ compareNorths:typeTranscription)

# Model selection
# Data driven model has the lowest AIC. It is significantly different from the minimal model and it
# provides more explanatory power (lower AIC) without being significantly different from them maximal model.

print(anova(m.AH_NCS.withInteractionsInRandom, m.AH_NCS.stepDriven))
print(anova(m.AH_NCS.withInteractionsInRandom, m.AH_NCS.minimal))
print(anova(m.AH_NCS.stepDriven, m.AH_NCS.minimal))

print(AIC(m.AH_NCS.withInteractionsInRandom))
print(AIC(m.AH_NCS.stepDriven))
print(AIC(m.AH_NCS.minimal))

# Assumptions:
# Linearity of residuals: the residuals are distributed around the zero line
# Homoskedascity: the residuals have a blob-like distribution
plot(fitted(m.AH_NCS.stepDriven),resid(m.AH_NCS.stepDriven))
abline(h=0, col="blue")
# Normality of residuals: the distribution is unimodal with similar tails on both sides.
# Moreover, there are 5730 datapoints, so, using the Central Limit Theorem, we can 
# assume that the residuals behave in a manner similar to the Normal distribution.
hist(resid(m.AH_NCS.stepDriven))
qqnorm(resid(m.AH_NCS.stepDriven))
qqline(resid(m.AH_NCS.stepDriven), col = "steelblue", lwd = 2)
length(resid(m.AH_NCS.stepDriven))

plot_model(m.AH_NCS.stepDriven, type='diag') # you can ajust type (see package info: ?plot_model)


#====================================================================.
#== R Citations                                                 =====
#====================================================================.

citation(package = "lme4")
citation(package = "adehabitatHR")  # for function kerneloverlap
citation(package = "stats")
citation(package = "emmeans")
citation()


#====================================================================.
#== Ethnicity and Gender                                        =====
#====================================================================.

t = subset(idea, typeTranscription=="groundtruth")
speakers = unique(t$name)
speakers = as.data.frame(speakers)
speakers$gender = ""

totalSpeakers = nrow(speakers)
for (i in 1:totalSpeakers) {
  temp = subset(idea, typeTranscription == "groundtruth" & name == speakers$speakers[i])
  speakers$gender[i] = as.character(unique(temp$gender))
}

with(speakers, table(gender))


t = subset(idea, typeTranscription=="groundtruth")
speakers = unique(t$name)
speakers = as.data.frame(speakers)
speakers$openEthnicity = ""

totalSpeakers = nrow(speakers)
for (i in 1:totalSpeakers) {
  temp = subset(idea, typeTranscription == "groundtruth" & name == speakers$speakers[i])
  speakers$openEthnicity[i] = as.character(unique(temp$ethnicity))
}

unique(speakers$openEthnicity)

speakers$ethnicity[speakers$openEthnicity == "Caucasian"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "Black"] = "Black"
speakers$ethnicity[speakers$openEthnicity == "African-American"] = "Black"
speakers$ethnicity[speakers$openEthnicity == "Yupâ€™ik"] = "NativeAm"
speakers$ethnicity[speakers$openEthnicity == "Tsimpshean"] = "NativeAm"
speakers$ethnicity[speakers$openEthnicity == "Latina/Chicana"] = "LatAm"
speakers$ethnicity[speakers$openEthnicity == "Hispanic/Chicana"] = "LatAm"
speakers$ethnicity[speakers$openEthnicity == "Navajo"] = "NativeAm"
speakers$ethnicity[speakers$openEthnicity == "white (non-Latina)"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "white"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "white (N. European ancestry)"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "white (half Polish ancestry)"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "white (W. European ancestry)"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "white with Italian & Puerto Rican heritage"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "Hispanic"] = "LatAm"
speakers$ethnicity[speakers$openEthnicity == "Native American"] = "NativeAm"
speakers$ethnicity[speakers$openEthnicity == "Native American (Zuni-Navajo)"] = "NativeAm"
speakers$ethnicity[speakers$openEthnicity == "Caucasian (German ancestry)"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "Puerto Rican/Caucasian"] = "Mixed"
speakers$ethnicity[speakers$openEthnicity == "Native American (Kiowa)"] = "NativeAm"
speakers$ethnicity[speakers$openEthnicity == "Caucasian/Irish-American"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "Native American (Seneca-Iroquois)"] = "NativeAm"
speakers$ethnicity[speakers$openEthnicity == "African American"] = "Black"
speakers$ethnicity[speakers$openEthnicity == "Caucasian (Jewish and Polish roots)"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "Caucasian/European and Cajun"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "Caucasian/European"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "black"] = "Black"
speakers$ethnicity[speakers$openEthnicity == "Caucasian (German"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "Caucasian (Irish/Russian)"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "Caucasian (Slovak ancestry)"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "Asian American"] = "AsianAm"
speakers$ethnicity[speakers$openEthnicity == "Caucasian (Irish and English)"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "Caucasian (Irish descent)"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "Caucasian (Dutch ancestry)"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "Caucasian and Native American"] = "Mixed"
speakers$ethnicity[speakers$openEthnicity == "African American/Jamaican"] = "Black"
speakers$ethnicity[speakers$openEthnicity == "African American/Hispanic"] = "Mixed"
speakers$ethnicity[speakers$openEthnicity == "Caucasian (Italian-American)"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "Jewish"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "Puerto Rican"] = "LatAm"
speakers$ethnicity[speakers$openEthnicity == "Jewish/Caucasian"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "Italian/Indo-Trinidadian"] = "Mixed"
speakers$ethnicity[speakers$openEthnicity == "Puerto Rican and Costa Rican"] = "LatAm"
speakers$ethnicity[speakers$openEthnicity == "Italian"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "Puerto Rican and Ecuadorian"] = "LatAm"
speakers$ethnicity[speakers$openEthnicity == "Dominican"] = "LatAm"
speakers$ethnicity[speakers$openEthnicity == "Haitian"] = "Black"
speakers$ethnicity[speakers$openEthnicity == "Caucasian/Italian-American"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "Caucasian (Greek and Jewish)"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "Polynesian/Caucasian"] = "Mixed"
speakers$ethnicity[speakers$openEthnicity == "Caucasian (Polish"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "Irish-American"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "Latino (Mexican-American)"] = "LatAm"
speakers$ethnicity[speakers$openEthnicity == "Caucasian (European ancestry)"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "Caucasian/Jewish (European ancestry)"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "Hispanic/Mexican-American"] = "LatAm"
speakers$ethnicity[speakers$openEthnicity == "Caucasian/Native American"] = "Mixed"
speakers$ethnicity[speakers$openEthnicity == "Native American (Navajo)"] = "NativeAm"
speakers$ethnicity[speakers$openEthnicity == "Native American (Tohono Oâ€™odham)/black"] = "Mixed"
speakers$ethnicity[speakers$openEthnicity == "Irish/German Caucasian"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "Caucasian/Portuguese"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "Black"] = "Black"
speakers$ethnicity[speakers$openEthnicity == "Cuban American"] = "LatAm"
speakers$ethnicity[speakers$openEthnicity == "Cuban/Caucasian"] = "Mixed"
speakers$ethnicity[speakers$openEthnicity == "Caucasian (with Jewish background)"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "Chicago"] = "NA"
speakers$ethnicity[speakers$openEthnicity == "Armenian"] = "Caucasian"
speakers$ethnicity[speakers$openEthnicity == "Native American/Guatemalan"] = "LatAm"
speakers$ethnicity[speakers$openEthnicity == "black (Creole and Puerto Rican roots)"] = "Black"
speakers$ethnicity[speakers$openEthnicity == ""] = "NA"


# Speakers by ethnic category
with(speakers, table(ethnicity))

# Percentage of speakers by ethnic category
((with(speakers, table(ethnicity)))/totalSpeakers)*100


#====================================================================.
#== Recording durations                                         =====
#====================================================================.

t = subset(idea, typeTranscription=="groundtruth")
speakers = unique(t$name)
speakers = as.data.frame(speakers)
speakers$duration = ""

totalSpeakers = nrow(speakers)
for (i in 1:totalSpeakers) {
  temp = subset(idea, typeTranscription == "groundtruth" & name == speakers$speakers[i])
  speakers$duration[i] = as.character(unique(temp$recordingDuration))
}

speakers$duration = as.double(speakers$duration)

# stats of duration
mean(speakers$duration)
sd(speakers$duration)
min(speakers$duration)
max(speakers$duration)
max(speakers$duration)

totalSeconds = sum(speakers$duration) # total in seconds
durationHours = floor(totalSeconds / 3600)
durationMinutes = floor((totalSeconds - durationHours*60*60)/60)
durationSeconds  = totalSeconds - (durationMinutes*60) - (durationHours * 3600)

totalSeconds
durationHours 
durationMinutes
durationSeconds


hist(speakers$duration, breaks = 50)

nrow(subset(speakers, duration < 60)) / nrow(speakers)  # percentage of speakers with less than a minute
nrow(subset(speakers, duration < 120)) / nrow(speakers)  # percentage of speakers with less than two minutes
