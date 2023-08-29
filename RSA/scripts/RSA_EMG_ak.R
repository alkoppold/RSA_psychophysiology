### RSA ASM x EMG
# *****
# AK 23.08.2023

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


# first run
# RSAFunctions_only_withRegression.R file

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




###########################################################
###########################################################
# Regression analyses, 29.08.2023
###########################################################
###########################################################
# line 874
OriginalData<- emg_SingleAlina_corr_prep_all
ModelData1<- BehavioralSingleAlina_invAK_corr_prep_all #BehavioralSingleAlina_SV_invAK_corr_prep_all
ModelData2<- BehavioralSingleAlina_corr_prep_all #BehavioralSingleAlina_SV_NN_corr_prep_all
whichModel<- c("ModelData1")#,"ModelData2")
VariableOriginal<- "SCRinvAKAroSingleValue" #"SCRinvAKAroSingleValue"
VariableModel2<- "NNVal"
VariableModel1<- "invAKVal"
seedPer<-49
nPermutations<-10000
PlotRSMReg<-1
Category = "VP"
PermutationTestRegModels(OriginalData,ModelData1,ModelData2,VariableOriginal,VariableModel1,VariableModel2, Category,nPermutations,seedPer,PlotRSMReg)
