#################################################
### CODE to search and install a specific package
#################################################
# list.of.packages <- c("base","psycho", "tidyverse", "dplyr", "ggplot2", "fs", "lme4", "lmerTest", "sjPlot")
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages)

#######
##Libraries
library(base)
library(psycho)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(fs)
library (lme4)
library(lmerTest)
library(sjPlot)
library(rstatix)
library(Hmisc)
# install.packages("extrafont")
library(extrafont)
library (ggbeeswarm)
library(BayesFactor)
library(bayestestR)
library(glue)
#################################################
### Load Fonts for figures
#################################################
#font_import()
loadfonts(device="win")       #Register fonts for Windows bitmap output
fonts()    
windowsFonts(Times=windowsFont("Times"))

# setwd("C:/Users/cvent/Desktop/RESEARCH/RSA_LondorfLab")

#Select Task Pasive Picture viewing
SCR_RSM<-readRDS("./scripts/SCR_RSM.rds") 
#------------------------
PicEdaAro<-SCR_RSM  %>%
  filter(Task_def==1)
#----
#Matrix should contain VP, Behavioral Variable (e.g., Aro), 
# variable coding the sorting (sortAro), Variables for the vector (variables containing EDA)
PicEdaAro<-PicEdaAro%>%
  select(Vp_ratings,sortAro, contains("EDA"),- contains("EDA_mean")) 

#----------------- DEFINE PARAMETERS FOR RUNING RSMVector for EDA and AROUSAL
#rename to favor the function. 
PicEdaAro<-PicEdaAro%>%
rename(sort= sortAro)

# OriginalData<-PicEdaAro
nVariables<-36 # this variable defines the number of cells e.g., pictures presented
CharContainingNumberFromVector<-4 #Character containing number from variables coding vector.
#e.g., EDA1, EDA2, EDA3... the character is 4: 1, 2, 3...
### Select the variables needed for the loop: VP, Variable sorting, Vector values
VariableType<-"SCR_Vector"
sortType<- "Arousal"


RSMVector(PicEdaAro,nVariables,CharContainingNumberFromVector,VariableType,sortType)


#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#Create Matrix based on single datapoints for EDA SINGLE POINT
PicEdaArosV<-SCR_RSM  %>%
  filter(Task_def==1)%>%
  dplyr::select(Vp_ratings, sortAro,EDA_mean_log)
PicEdaArosV<-PicEdaArosV%>%
  dplyr:: rename(sort= sortAro)
# OriginalData<-PicEdaArosV
MaxValue<- max(PicEdaArosV$EDA_mean_log, na.rm=T)

nVariables<-36 # this variable defines the number of cells e.g., pictures presented
CharContainingNumberFromVector<-4 #Character containing number from variables coding vector.
#e.g., EDA1, EDA2, EDA3... the character is 4: 1, 2, 3...
### Select the variables needed for the loop: VP, Variable sorting, Vector values
VariableType<-"SCR_Single"
sortType<- "Arousal"

RSMSingleValue(PicEdaArosV,nVariables,CharContainingNumberFromVector,VariableType,sortType, MaxValue)



#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#Create Matrix based on single BEHAVIORAL
PicEdaArosV<-SCR_RSM  %>%
  filter(Task_def==1)%>%
  dplyr::select(Vp_ratings, sortAro,Aro)
PicEdaArosV<-PicEdaArosV%>%
  rename(Behavioral= Aro)%>%
  rename(sort= sortAro)
# OriginalData<-PicEdaArosV
MaxValue<- max(PicEdaArosV$Behavioral, na.rm=T)

nVariables<-36 # this variable defines the number of cells e.g., pictures presented
CharContainingNumberFromVector<-3 #Character containing number from variables coding vector.
#e.g., EDA1, EDA2, EDA3... the character is 4: 1, 2, 3...
### Select the variables needed for the loop: VP, Variable sorting, Vector values
VariableType<-"Behav_Single"
sortType<- "Arousal"

RSMSingleValue(PicEdaArosV,nVariables,CharContainingNumberFromVector,VariableType,sortType, MaxValue)

VariableType<-"Behav_Single_AK"
sortType<- "Arousal"

RSMSingleValueAK(PicEdaArosV,nVariables,CharContainingNumberFromVector,VariableType,sortType, MaxValue)

VariableType<-"Behav_Single_invAK"
sortType<- "Arousal"

RSMSingleValueinvAK(PicEdaArosV,nVariables,CharContainingNumberFromVector,VariableType,sortType, MaxValue)

#---------------------------

### RUN PERMUTATION TEST 
# 
#  SCR_Vector_corr_prep_all
# 
#  Behav_Single_corr_prep_all
seedPer<-55
nPermutations<-10000
VariableOriginal<-"SCR_Vector"
VariableModel<- "ArousalRatings"

PermutationTest(SCR_Vector_corr_prep_all,Behav_Single_corr_prep_all,VariableOriginal,VariableModel,nPermutations,seedPer)

#inv AK
seedPer<-55
nPermutations<-10000
VariableOriginal<-"SCR_Vector"
VariableModel<- "ArousalRatings"
PermutationTest(SCR_Vector_corr_prep_all,Behav_Single_invAK_corr_prep_all,VariableOriginal,VariableModel,nPermutations,seedPer)

#---------------------
seedPer<-55
nPermutations<-10000
VariableOriginal<-"SCR_Single"
VariableModel<- "ArousalRatings"
# OriginalData<-SCR_Single_corr_prep_all
# ModelData<-Behav_Single_corr_prep_all
PermutationTest(SCR_Single_corr_prep_all,Behav_Single_corr_prep_all,VariableOriginal,VariableModel,nPermutations,seedPer)

#inv AK
seedPer<-55
nPermutations<-10000
VariableOriginal<-"SCR_Vector"
VariableModel<- "ArousalRatings"
PermutationTest(SCR_Single_corr_prep_all,Behav_Single_invAK_corr_prep_all,VariableOriginal,VariableModel,nPermutations,seedPer)


#---------------------
seedPer<-55
nPermutations<-10000
VariableOriginal<-"SCR_Single"
VariableModel<- "SCR_Vector"
# OriginalData<-SCR_Single_corr_prep_all
# ModelData<-Behav_Single_corr_prep_all
PermutationTest(SCR_Single_corr_prep_all,SCR_Vector_corr_prep_all,VariableOriginal,VariableModel,nPermutations,seedPer)




################ test function RSINdividual 


Model<-"AroRatings"
Data<- "SCRVector"

RSMIndividual(Behav_Single_corr_prep_all,SCR_Vector_corr_prep_all,Model,Data)


Model<-"AroRatings"
Data<- "SCRSingle"
RSMIndividual(Behav_Single_corr_prep_all,SCR_Single_corr_prep_all,Model,Data)



Model<-"AroRatings"
Data<- "SCRSingle"
RSMIndividual(SCR_Single_corr_prep_all,Behav_Single_invAK_corr_prep_all,Model,Data)



Model<-"AroRatings"
Data<- "SCRSingle"
RSMIndividual(SCR_Vector_corr_prep_all,Behav_Single_invAK_corr_prep_all,Model,Data)


Model<-"SCRVector"
Data<- "SCRSingle"
RSMIndividual(SCR_Vector_corr_prep_all,SCR_Single_corr_prep_all,Model,Data)


#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#----------------------- PASSIVE SOUND LISTENING TASK ---------------------------

SoundEdaAro<-SCR_RSM  %>%
  filter(Task_def==2)
#----
#Matrix should contain VP, Behavioral Variable (e.g., Aro), 
# variable coding the sorting (sortAro), Variables for the vector (variables containing EDA)
SoundEdaAro<-SoundEdaAro%>%
  select(Vp_ratings,sortAro, contains("EDA"),- contains("EDA_mean")) 

#----------------- DEFINE PARAMETERS FOR RUNING RSMVector for EDA and AROUSAL
#rename to favor the function. 
SoundEdaAro<-SoundEdaAro%>%
  rename(sort= sortAro)

# OriginalData<-SoundEdaAro
nVariables<-36 # this variable defines the number of cells e.g., Sounds presented
CharContainingNumberFromVector<-4 #Character containing number from variables coding vector.
#e.g., EDA1, EDA2, EDA3... the character is 4: 1, 2, 3...
### Select the variables needed for the loop: VP, Variable sorting, Vector values
VariableType<-"SCR_Vector_PSL"
sortType<- "Arousal"


RSMVector(SoundEdaAro,nVariables,CharContainingNumberFromVector,VariableType,sortType)


#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#Create Matrix based on single datapoints for EDA SINGLE POINT
SoundEdaArosV<-SCR_RSM  %>%
  filter(Task_def==2)%>%
  dplyr::select(Vp_ratings, sortAro,EDA_mean_log)
SoundEdaArosV<-SoundEdaArosV%>%
  dplyr:: rename(sort= sortAro)
# OriginalData<-SoundEdaArosV
MaxValue<- max(SoundEdaArosV$EDA_mean_log, na.rm=T)

nVariables<-36 # this variable defines the number of cells e.g., Sounds presented
CharContainingNumberFromVector<-4 #Character containing number from variables coding vector.
#e.g., EDA1, EDA2, EDA3... the character is 4: 1, 2, 3...
### Select the variables needed for the loop: VP, Variable sorting, Vector values
VariableType<-"SCR_Single_PSL"
sortType<- "Arousal"

RSMSingleValue(SoundEdaArosV,nVariables,CharContainingNumberFromVector,VariableType,sortType, MaxValue)



#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#Create Matrix based on single BEHAVIORAL
SoundEdaArosV<-SCR_RSM  %>%
  filter(Task_def==2)%>%
  dplyr::select(Vp_ratings, sortAro,Aro)
SoundEdaArosV<-SoundEdaArosV%>%
  rename(Behavioral= Aro)%>%
  rename(sort= sortAro)
# OriginalData<-SoundEdaArosV
MaxValue<- max(SoundEdaArosV$Behavioral, na.rm=T)

nVariables<-36 # this variable defines the number of cells e.g., Soundtures presented
CharContainingNumberFromVector<-3 #Character containing number from variables coding vector.
#e.g., EDA1, EDA2, EDA3... the character is 4: 1, 2, 3...
### Select the variables needed for the loop: VP, Variable sorting, Vector values
VariableType<-"Behav_Single_PSL"
sortType<- "Arousal"

RSMSingleValue(SoundEdaArosV,nVariables,CharContainingNumberFromVector,VariableType,sortType, MaxValue)


#AK
VariableType<-"Behav_Single_PSL_AK"
sortType<- "Arousal"

RSMSingleValueAK(SoundEdaArosV,nVariables,CharContainingNumberFromVector,VariableType,sortType, MaxValue)
#invAK
VariableType<-"Behav_Single_PSL_invAK"
sortType<- "Arousal"

RSMSingleValueinvAK(SoundEdaArosV,nVariables,CharContainingNumberFromVector,VariableType,sortType, MaxValue)
#---------------------------

### RUN PERMUTATION TEST 
# 
#  SCR_Vector_corr_prep_all
# 
#  Behav_Single_corr_prep_all
seedPer<-55
nPermutations<-10000
VariableOriginal<-"SCR_Vector"
VariableModel<- "ArousalRatings"

PermutationTest(SCR_Vector_PSL_corr_prep_all,Behav_Single_PSL_corr_prep_all,VariableOriginal,VariableModel,nPermutations,seedPer)

#invAK
seedPer<-55
nPermutations<-10000
VariableOriginal<-"SCR_Vector"
VariableModel<- "ArousalRatings"

PermutationTest(SCR_Vector_PSL_corr_prep_all,Behav_Single_PSL_invAK_corr_prep_all,VariableOriginal,VariableModel,nPermutations,seedPer)


#---------------------
seedPer<-55
nPermutations<-10000
VariableOriginal<-"SCR_Single"
VariableModel<- "ArousalRatings"
# OriginalData<-SCR_Single_corr_prep_all
# ModelData<-Behav_Single_corr_prep_all
PermutationTest(SCR_Single_PSL_corr_prep_all,Behav_Single_PSL_corr_prep_all,VariableOriginal,VariableModel,nPermutations,seedPer)

#invAK
seedPer<-55
nPermutations<-10000
VariableOriginal<-"SCR_Vector"
VariableModel<- "ArousalRatings"

PermutationTest(SCR_Single_PSL_corr_prep_all,Behav_Single_PSL_invAK_corr_prep_all,VariableOriginal,VariableModel,nPermutations,seedPer)



#---------------------
seedPer<-55
nPermutations<-10000
VariableOriginal<-"SCR_Single"
VariableModel<- "SCR_Vector"
# OriginalData<-SCR_Single_corr_prep_all
# ModelData<-Behav_Single_corr_prep_all
PermutationTest(SCR_Single_PSL_corr_prep_all,SCR_Vector_PSL_corr_prep_all,VariableOriginal,VariableModel,nPermutations,seedPer)




################ test function RSINdividual 


Model<-"AroRatings"
Data<- "SCRVector"

RSMIndividual(Behav_Single_PSL_corr_prep_all,SCR_Vector_PSL_corr_prep_all,Model,Data)
#invAK
Model<-"AroRatings"
Data<- "SCRSingle"
RSMIndividual(Behav_Single_PSL_invAK_corr_prep_all,SCR_Vector_PSL_corr_prep_all,Model,Data)


Model<-"AroRatings"
Data<- "SCRSingle"
RSMIndividual(Behav_Single_PSL_corr_prep_all,SCR_Single_PSL_corr_prep_all,Model,Data)

#invAK
Model<-"AroRatings"
Data<- "SCRSingle"
RSMIndividual(Behav_Single_PSL_invAK_corr_prep_all,SCR_Single_PSL_corr_prep_all,Model,Data)


Model<-"SCRVector"
Data<- "SCRSingle"
RSMIndividual(SCR_Vector_PSL_corr_prep_all,SCR_Single_PSL_corr_prep_all,Model,Data)




#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#----------------------- IMAGERY TASK ---------------------------

ImagEdaAro<-SCR_RSM  %>%
  filter(Task_def==3)

ImagEdaAro2<-  ImagEdaAro%>%
  gather(EDA_output, EDA, EDA1:EDA_mean)%>%
  mutate(EDA_output, EDA_num= str_sub(EDA_output, 4))
# ImagEdaAro2$EDA_num<- as.numeric( ImagEdaAro2$EDA_num)

ImagEdaAro2<-ImagEdaAro2%>%
  group_by(Vp_ratings, EDA_num,Code_Stimu)%>%
  summarise(meanAro=(max(sortAro))/2,mean= mean(EDA)  )%>%
  arrange(meanAro)%>%
  arrange(Vp_ratings)%>%
  filter(EDA_num!="_mean")%>%
  select(-Code_Stimu)

ImagEdaAro2$EDA_num<-as.numeric (ImagEdaAro2$EDA_num)
ImagEdaAro2<-ImagEdaAro2%>%
  arrange(EDA_num)%>%
  arrange(meanAro)%>%
arrange(Vp_ratings)


ImagEdaAro2$EDA_num<-paste("EDA",ImagEdaAro2$EDA_num, sep= "")

ImagEdaAro3<- ImagEdaAro2%>%
  spread(EDA_num,mean)%>%
  rename(sortAro= meanAro)
# names(ImagEdaAro3)

#----
#Matrix should contain VP, Behavioral Variable (e.g., Aro), 
# variable coding the sorting (sortAro), Variables for the vector (variables containing EDA)
ImagEdaAro<-ImagEdaAro3%>%
  select(Vp_ratings,sortAro, contains("EDA"),- contains("EDA_mean")) 

#----------------- DEFINE PARAMETERS FOR RUNING RSMVector for EDA and AROUSAL
#rename to favor the function. 
ImagEdaAro<-ImagEdaAro%>%
  rename(sort= sortAro)

# OriginalData<-ImagEdaAro
nVariables<-18 # this variable defines the number of cells e.g., Imagtures presented
CharContainingNumberFromVector<-4 #Character containing number from variables coding vector.
#e.g., EDA1, EDA2, EDA3... the character is 4: 1, 2, 3...
### Select the variables needed for the loop: VP, Variable sorting, Vector values
VariableType<-"SCR_Vector_I"
sortType<- "Arousal"

OriginalData<-ImagEdaAro
RSMVector(ImagEdaAro,nVariables,CharContainingNumberFromVector,VariableType,sortType)


#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#Create Matrix based on single datapoints for EDA SINGLE POINT
ImagEdaArosV<-SCR_RSM  %>%
  filter(Task_def==3)%>%
  dplyr::select(Vp_ratings, Code_Stimu,sortAro,EDA_mean_log)


ImagEdaArosV2<-ImagEdaArosV%>%
  group_by(Vp_ratings,Code_Stimu)%>%
  summarise(meanAro=(max(sortAro))/2, EDA_mean_log= mean(EDA_mean_log)  )%>%
  arrange(meanAro)%>%
  arrange(Vp_ratings)%>%ungroup()%>%
  select(-Code_Stimu)


ImagEdaArosV<-ImagEdaArosV2%>%
  dplyr:: rename(sort= meanAro)
# OriginalData<-ImagEdaArosV
MaxValue<- max(ImagEdaArosV$EDA_mean_log, na.rm=T)

nVariables<-18 # this variable defines the number of cells e.g., Imagtures presented
CharContainingNumberFromVector<-4 #Character containing number from variables coding vector.
#e.g., EDA1, EDA2, EDA3... the character is 4: 1, 2, 3...
### Select the variables needed for the loop: VP, Variable sorting, Vector values
VariableType<-"SCR_Single_I"
sortType<- "Arousal"
OriginalData<-ImagEdaArosV
RSMSingleValue(ImagEdaArosV,nVariables,CharContainingNumberFromVector,VariableType,sortType, MaxValue)



#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#Create Matrix based on single BEHAVIORAL
ImagEdaArosV<-SCR_RSM  %>%
  filter(Task_def==3)%>%
  dplyr::select(Vp_ratings, sortAro,Aro,Code_Stimu)


ImagEdaArosV2<-ImagEdaArosV%>%
  group_by(Vp_ratings,Code_Stimu)%>%
  summarise(meanAro=(max(sortAro))/2, Aro= mean(Aro)  )%>%
  arrange(meanAro)%>%
  arrange(Vp_ratings)%>%ungroup()%>%
  select(-Code_Stimu)


ImagEdaArosV<-ImagEdaArosV2%>%
  rename(Behavioral= Aro)%>%
  rename(sort= meanAro)
# OriginalData<-ImagEdaArosV
MaxValue<- max(ImagEdaArosV$Behavioral, na.rm=T)

nVariables<-18 # this variable defines the number of cells e.g., Imagtures presented
CharContainingNumberFromVector<-3 #Character containing number from variables coding vector.
#e.g., EDA1, EDA2, EDA3... the character is 4: 1, 2, 3...
### Select the variables needed for the loop: VP, Variable sorting, Vector values
VariableType<-"Behav_Single_I"
sortType<- "Arousal"

RSMSingleValue(ImagEdaArosV,nVariables,CharContainingNumberFromVector,VariableType,sortType, MaxValue)
Behav_Single_I_corr_prep_allTest<-Behav_Single_I_corr_prep_all%>%select(var1,var2,mean_all)

#AK
VariableType<-"Behav_Single_I_AK"
sortType<- "Arousal"

RSMSingleValueAK(ImagEdaArosV,nVariables,CharContainingNumberFromVector,VariableType,sortType, MaxValue)
#invAK
VariableType<-"Behav_Single_I_invAK"
sortType<- "Arousal"

RSMSingleValueinvAK(ImagEdaArosV,nVariables,CharContainingNumberFromVector,VariableType,sortType, MaxValue)
#---------------------------

### RUN PERMUTATION TEST 
# 
#  SCR_Vector_corr_prep_all
# 
#  Behav_Single_corr_prep_all
seedPer<-55
nPermutations<-10000
VariableOriginal<-"SCR_Vector"
VariableModel<- "ArousalRatings"

PermutationTest(SCR_Vector_I_corr_prep_all,Behav_Single_I_corr_prep_all,VariableOriginal,VariableModel,nPermutations,seedPer)

#invAK
seedPer<-55
nPermutations<-10000
VariableOriginal<-"SCR_Vector"
VariableModel<- "ArousalRatings"

PermutationTest(SCR_Vector_I_corr_prep_all,Behav_Single_I_invAK_corr_prep_all,VariableOriginal,VariableModel,nPermutations,seedPer)


#---------------------
seedPer<-55
nPermutations<-10000
VariableOriginal<-"SCR_Single"
VariableModel<- "ArousalRatings"
# OriginalData<-SCR_Single_corr_prep_all
# ModelData<-Behav_Single_corr_prep_all
PermutationTest(SCR_Single_I_corr_prep_all,Behav_Single_I_corr_prep_all,VariableOriginal,VariableModel,nPermutations,seedPer)

#invAK
seedPer<-55
nPermutations<-10000
VariableOriginal<-"SCR_Vector"
VariableModel<- "ArousalRatings"

PermutationTest(SCR_Single_I_corr_prep_all,Behav_Single_I_invAK_corr_prep_all,VariableOriginal,VariableModel,nPermutations,seedPer)



#---------------------
seedPer<-55
nPermutations<-10000
VariableOriginal<-"SCR_Single"
VariableModel<- "SCR_Vector"
# OriginalData<-SCR_Single_corr_prep_all
# ModelData<-Behav_Single_corr_prep_all
PermutationTest(SCR_Single_I_corr_prep_all,SCR_Vector_I_corr_prep_all,VariableOriginal,VariableModel,nPermutations,seedPer)




################ test function RSINdividual 


Model<-"AroRatings"
Data<- "SCRVector"

RSMIndividual(Behav_Single_I_corr_prep_all,SCR_Vector_I_corr_prep_all,Model,Data)
#invAK
Model<-"AroRatings"
Data<- "SCRSingle"
RSMIndividual(Behav_Single_I_invAK_corr_prep_all,SCR_Vector_I_corr_prep_all,Model,Data)


Model<-"AroRatings"
Data<- "SCRSingle"
RSMIndividual(Behav_Single_I_corr_prep_all,SCR_Single_I_corr_prep_all,Model,Data)

#invAK
Model<-"AroRatings"
Data<- "SCRSingle"
RSMIndividual(Behav_Single_I_invAK_corr_prep_all,SCR_Single_I_corr_prep_all,Model,Data)


Model<-"SCRVector"
Data<- "SCRSingle"
RSMIndividual(SCR_Vector_I_corr_prep_all,SCR_Single_I_corr_prep_all,Model,Data)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# -------------- ALINAs DATA----------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# RSA based on ASM data
### 1. RSA ASM x SCR
load("./data/asm_scr_proofed.Rdata") 
load("./data/df_rat.Rdata") 
#test with Alina's data set. Please note that not all participants are used for the test.
#Note. how the variable ID is extracted is not applicable to all participants. For this reason, 
#VP201 is filter (N cases for that participant over 1000)

# carlos: 
# scr$ID<- str_sub(scr$Subject,9,-1L)
# scr$picNum<-str_sub(scr$Stimulus,11,-1L)
# df_rat$picNum<-str_sub(df_rat$picture,1,6)

## by alina: 
scr$ID = sub('.*(\\d{3}).*', '\\1', scr$Subject)
scr_sub =scr %>% select(Subject, ID, project, Stimulus,logAmplitudes,NumberOfTimePresented) %>% unique() # to be added  Q_LTE_di, Q_CTQ_di, Q_STAIT_sum, Q_BDI_sum
scr_sub$picNum<-str_sub(scr_sub$Stimulus,11,-1L)
scr_sub$Stimulus = NULL
scr_sub$ID = gsub("[^0-9.-]", "",scr_sub$ID)
scr_sub$ID = gsub("^0", "", scr_sub$ID) 
scr_sub$ID = gsub("^0", "", scr_sub$ID) 

# prep df_rat
df_rat$picNum<-str_sub(df_rat$picture,1,6)
df_rat$picture = NULL
df_rat_scr = merge(scr_sub, df_rat, by = c("ID", "project", "picNum"))
#detach(package:plyr) # if "sanity" is empty
sanity = df_rat_scr %>% group_by(Subject) %>% tally()


### why log amplitudes?
scrAlina<-df_rat_scr%>%select(ID,project, picNum,logAmplitudes,NumberOfTimePresented, arousal_rating)#%>%  
  #filter(ID!="") # not needed anymore

### why is this subject excluded?
#df_ratAlina<-df_rat%>%select(ID,picNum,arousal_rating)%>%
#  filter(ID!=201)

#scr_rat<-inner_join(scrAlina,df_ratAlina, by= c("ID", "picNum"), keep=F)
#N = 252

# get old ID structure back, if needed 
scrAlina$ID <- paste(scrAlina$project ,scrAlina$ID ,  sep = "_")
scrAlina$project = NULL 
scr_rat = scrAlina
# new: n = 494

#correlation between time presented 1 and two
# alina comment: the cor() function works only when the datasets have equal length, which they do, but when removing NA's they become unequal, see code below 
test1 = scr_rat$logAmplitudes[scr_rat$NumberOfTimePresented==1]
test2 = scr_rat$logAmplitudes[scr_rat$NumberOfTimePresented==2]
# table(is.na(test1))
# 
# FALSE  TRUE 
# 17646   138 
# > table(is.na(test2))
# 
# FALSE  TRUE 
# 17570   214 
# --> problem comes from NA differences 
# solution:
test_df = cbind(test1, test2)
test_df = as.data.frame(test_df)
test_df_clean = test_df %>% drop_na()
cor(test_df_clean)

## solve it in dataset scr_rats
# Remove rows with NA values in 'logAmplitudes', considering duplicates in 'ID' and 'picNum'
scr_rat <- scr_rat %>%
  group_by(ID, picNum) %>%
  filter(!any(is.na(logAmplitudes))) %>%
  ungroup()
cor(scr_rat$logAmplitudes[scr_rat$NumberOfTimePresented==1], 
    scr_rat$logAmplitudes[scr_rat$NumberOfTimePresented==2], 
    method= "spearman")

#### end change alina

################################################################
##### new functions added: picture presentation / non-responder 

AverageResponsesacrossruns<- 1# determine whether average across trials if set to zero takes the first run.
if (AverageResponsesacrossruns==1){
  ntrials<-72
  
  scrAlinacounZeros<-scrAlina%>%
    filter(logAmplitudes==0)%>%
    # filter(NumberOfTimePresented==1)%>%
    select(ID,logAmplitudes)
  VPZeros<-as.data.frame(table(scrAlinacounZeros$ID))
  names(VPZeros)<- c("ID", "Freq")
  
}else{
  ntrials<-36
  scrAlinacounZeros<-scrAlina%>%
    filter(logAmplitudes==0)%>%
    filter(NumberOfTimePresented==1)%>%
    select(ID,logAmplitudes)
  VPZeros<-as.data.frame(table(scrAlinacounZeros$ID))
  names(VPZeros)<- c("ID", "Freq")
  
}
removeNonResponders<-0
RespondersRate<- .33
if (removeNonResponders==1){
  VPZerosfilt<-VPZeros%>% filter(Freq> ntrials-ntrials*RespondersRate)
  scrAlina<-anti_join(scrAlina,VPZerosfilt, by = "ID" )
}
# df_ratAlina<-scrAlina%>%select(ID,picNum,arousal_rating) #%>%
# filter(ID!=201)%>%
# filter(ID!= 126 & ID!= 144) #filter nonresponders in any trial
# scr_rat<-inner_join(scrAlina,df_ratAlina, by= c("ID", "picNum"), keep=F)
#View(table(scr_rat$ID)); N = 252 with filtering, N = 166
#correlation between time presented 1 and two 
cor(scr_rat$logAmplitudes[scr_rat$NumberOfTimePresented==1],
    scr_rat$logAmplitudes[scr_rat$NumberOfTimePresented==2], method= "spearman")
#for the test, and given the relatively low correlation take only first presentation
if (AverageResponsesacrossruns==1){
  scr_rat$arousal_rating<-as.numeric(scr_rat$arousal_rating)
  scr_rat_1<-scr_rat%>%
    group_by(ID,picNum)%>%
    dplyr::summarise(logAmplitudes=mean(logAmplitudes, na.rm =T),
                     arousal_rating=mean(arousal_rating,na.rm=T))
  scr_rat_1<-scr_rat_1%>%ungroup()
}else{scr_rat_1<-scr_rat%>%filter(NumberOfTimePresented==1)}



##############################################################

#for the test, and given the relatively low correlation take only first presentation
#scr_rat_1<-scr_rat%>%filter(NumberOfTimePresented==1)

## sort  participant, and arousal 
scr_rat_1<- scr_rat_1%>% 
  arrange_at ("arousal_rating")%>%
  arrange_at("ID")

#create a new variable that "counts the trials based on their arousal rating
j<-1
scr_rat_1$sortAro<-0
for (i in 1:nrow(scr_rat_1)){
  if (i< nrow(scr_rat_1)){
    
    if(scr_rat_1$ID[i]==  scr_rat_1$ID[i+1]){
      scr_rat_1$sortAro[i]<-j
      j <- j+1
    }
    
    else{
      scr_rat_1$sortAro[i]<-j
      j <-1
    }
  }
  else {scr_rat_1$sortAro[i]<-j
  }
}

PicEdaArosVAlina<-scr_rat_1%>%
  dplyr:: rename(sort= sortAro)%>%
  dplyr::select(ID, sort,logAmplitudes)

# OriginalData<-PicEdaArosV
MaxValue<- max(PicEdaArosVAlina$logAmplitudes, na.rm=T)

nVariables<-36 # this variable defines the number of cells e.g., pictures presented
CharContainingNumberFromVector<-4 #Character containing number from variables coding vector.
#e.g., EDA1, EDA2, EDA3... the character is 4: 1, 2, 3...
### Select the variables needed for the loop: VP, Variable sorting, Vector values
VariableType<-"SCR_SingleAlina"
sortType<- "Arousal"
# 1. based in amplitudes
RSMSingleValue(PicEdaArosVAlina,nVariables,CharContainingNumberFromVector,VariableType,sortType, MaxValue)


# RSM for Behavioral DAT
PicBehavArosVAlina<-scr_rat_1%>%
  dplyr:: rename(sort= sortAro)%>%
  dplyr::select(ID, sort,arousal_rating)%>%
  rename(Behavioral= arousal_rating)
PicBehavArosVAlina$Behavioral<-as.numeric(PicBehavArosVAlina$Behavioral)
# OriginalData<-PicEdaArosV
MaxValue<- max(PicBehavArosVAlina$Behavioral, na.rm=T)

nVariables<-36 # this variable defines the number of cells e.g., pictures presented
CharContainingNumberFromVector<-4 #Character containing number from variables coding vector.
#e.g., EDA1, EDA2, EDA3... the character is 4: 1, 2, 3...
### Select the variables needed for the loop: VP, Variable sorting, Vector values
VariableType<-"BehavioralSingleAlina"
sortType<- "Arousal"
OriginalData<-PicBehavArosVAlina

#2. based on arousal ratings 
RSMSingleValue(PicBehavArosVAlina,nVariables,CharContainingNumberFromVector,VariableType,sortType, MaxValue)

VariableType<-"BehavioralSingleAlina_invAK"
sortType<- "Arousal"

# 3.based on arousal invAK 
RSMSingleValueinvAK(PicBehavArosVAlina,nVariables,CharContainingNumberFromVector,VariableType,sortType, MaxValue)


#  Behav_Single_corr_prep_all
seedPer<-55
nPermutations<-10000
VariableOriginal<-"SCR_Vector"
VariableModel<- "ArousalRatings"

# old 
#PermutationTest(SCR_Vector_I_corr_prep_all,Behav_Single_I_corr_prep_all,VariableOriginal,VariableModel,nPermutations,seedPer)
# changed in our zoom call 
PermutationTest(SCR_SingleAlina_corr_prep_all,BehavioralSingleAlina_corr_prep_all,VariableOriginal,VariableModel,nPermutations,seedPer)


#invAK
seedPer<-55
nPermutations<-10000
VariableOriginal<-"SCR_Vector"
VariableModel<- "ArousalRatings"

#the output matrices for the Permutation test are created in the "RSM" functions

PermutationTest(SCR_SingleAlina_corr_prep_all,BehavioralSingleAlina_corr_prep_all,VariableOriginal,VariableModel,nPermutations,seedPer)

PermutationTest(SCR_SingleAlina_corr_prep_all,BehavioralSingleAlina_invAK_corr_prep_all,VariableOriginal,VariableModel,nPermutations,seedPer)


Model<-"SCRSingle"
Data<- "SCRArousal"
RSMIndividual(SCR_SingleAlina_corr_prep_all,BehavioralSingleAlina_corr_prep_all,Model,Data)


Model<-"SCRSingle"
Data<- "SCRArousal"
RSMIndividual(SCR_SingleAlina_corr_prep_all,BehavioralSingleAlina_invAK_corr_prep_all,Model,Data)

# *****
# added by Alina
### 2. RSA ASM x EMG
load("./data/allData_proofed.Rdata") 
load("./data/df_rat.Rdata") 

# CAUTION variable here got an "s" more: NumberOfTimesPresented not NumberOfTimePresented

emg$ID = sub('.*(\\d{3}).*', '\\1', emg$Subject)
emg_sub =emg %>% select(Subject, ID, project, Stimulus,TAmplitudes,NumberOfTimesPresented) %>% unique() # to be added  Q_LTE_di, Q_CTQ_di, Q_STAIT_sum, Q_BDI_sum
emg_sub$picNum<-str_sub(emg_sub$Stimulus,11,-1L)
emg_sub$Stimulus = NULL
emg_sub$ID = gsub("[^0-9.-]", "",emg_sub$ID)
emg_sub$ID = gsub("^0", "", emg_sub$ID) 
emg_sub$ID = gsub("^0", "", emg_sub$ID) 

# prep df_rat
df_rat$picNum<-str_sub(df_rat$picture,1,6)
df_rat$picture = NULL
df_rat_emg = merge(emg_sub, df_rat, by = c("ID", "project", "picNum"))
#detach(package:plyr) # if "sanity" is empty
sanity = df_rat_emg %>% group_by(Subject) %>% tally()


### why log amplitudes?
emgAlina<-df_rat_emg%>%select(ID,project, picNum,TAmplitudes,NumberOfTimesPresented, valence_rating)

# get old ID structure back, if needed 
emgAlina$ID <- paste(emgAlina$project ,emgAlina$ID ,  sep = "_")
emgAlina$project = NULL 
emg_rat = emgAlina
# new: n = 401

#### end change alina

## sort  participant, and valence 
emg_rat_1<- emg_rat%>% 
  arrange_at ("valence_rating")%>%
  arrange_at("ID")

#create a new variable that "counts the trials based on their arousal rating
j<-1
emg_rat_1$sortAro<-0
for (i in 1:nrow(emg_rat_1)){
  if (i< nrow(emg_rat_1)){
    
    if(emg_rat_1$ID[i]==  emg_rat_1$ID[i+1]){
      emg_rat_1$sortAro[i]<-j
      j <- j+1
    }
    
    else{
      emg_rat_1$sortAro[i]<-j
      j <-1
    }
  }
  else {emg_rat_1$sortAro[i]<-j
  }
}

PicEdaArosVAlina<-emg_rat_1%>%
  dplyr:: rename(sort= sortAro)%>%
  dplyr::select(ID, sort,TAmplitudes)

# OriginalData<-PicEdaArosV
MaxValue<- max(PicEdaArosVAlina$TAmplitudes, na.rm=T)

nVariables<-36 # this variable defines the number of cells e.g., pictures presented
CharContainingNumberFromVector<-4 #Character containing number from variables coding vector.
#e.g., EDA1, EDA2, EDA3... the character is 4: 1, 2, 3...
### Select the variables needed for the loop: VP, Variable sorting, Vector values

# 1.
VariableType<-"emg_SingleAlina"
sortType<- "Valence"
RSMSingleValue(PicEdaArosVAlina,nVariables,CharContainingNumberFromVector,VariableType,sortType, MaxValue)

#############################
# RSM for Behavioral DAT
PicBehavArosVAlina<-emg_rat_1%>%
  dplyr:: rename(sort= sortAro)%>%
  dplyr::select(ID, sort,valence_rating)%>%
  rename(Behavioral= valence_rating)
PicBehavArosVAlina$Behavioral<-as.numeric(PicBehavArosVAlina$Behavioral)
# OriginalData<-PicEdaArosV
MaxValue<- max(PicBehavArosVAlina$Behavioral, na.rm=T)

nVariables<-36 # this variable defines the number of cells e.g., pictures presented
CharContainingNumberFromVector<-4 #Character containing number from variables coding vector.
#e.g., EDA1, EDA2, EDA3... the character is 4: 1, 2, 3...
### Select the variables needed for the loop: VP, Variable sorting, Vector values

#2.
VariableType<-"BehavioralSingleAlina"
sortType<- "Valence"
OriginalData<-PicBehavArosVAlina
RSMSingleValue(PicBehavArosVAlina,nVariables,CharContainingNumberFromVector,VariableType,sortType, MaxValue)

# 3. 
VariableType<-"BehavioralSingleAlina_invAK"
sortType<- "Valence"
RSMSingleValueAK(PicBehavArosVAlina,nVariables,CharContainingNumberFromVector,VariableType,sortType, MaxValue)


#  Behav_Single_corr_prep_all
seedPer<-55
nPermutations<-10000
VariableOriginal<-"SCR_Vector" # also EMG?
VariableModel<- "ValenceRatings"

PermutationTest(emg_SingleAlina_corr_prep_all,BehavioralSingleAlina_corr_prep_all,VariableOriginal,VariableModel,nPermutations,seedPer)

#AK
seedPer<-55
nPermutations<-10000
VariableOriginal<-"SCR_Vector"# also EMG?
VariableModel<- "ValenceRatings"

#the output matrices for the Permutation test are created in the "RSM" functions

#PermutationTest(emg_SingleAlina_corr_prep_all,BehavioralSingleAlina_corr_prep_all,VariableOriginal,VariableModel,nPermutations,seedPer)

PermutationTest(emg_SingleAlina_corr_prep_all,BehavioralSingleAlina_invAK_corr_prep_all,VariableOriginal,VariableModel,nPermutations,seedPer)


Model<-"SCRSingle"
Data<- "SCRValence"
RSMIndividual(emg_SingleAlina_corr_prep_all,BehavioralSingleAlina_corr_prep_all,Model,Data)


Model<-"SCRSingle"
Data<- "SCRValence"
RSMIndividual(emg_SingleAlina_corr_prep_all,BehavioralSingleAlina_invAK_corr_prep_all,Model,Data)
