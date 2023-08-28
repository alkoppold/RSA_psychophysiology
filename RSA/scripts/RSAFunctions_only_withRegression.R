
#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
# Functions to round the floor and ceiling function to a decimal.
floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)
ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)


#Function to extract the last n characters: 
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#FUNCTION CREATE TO COMPARE AVERAGED RSMs based on a vector
#-----------testing data
# nVariables<-nrow(SCR_SingleAlina_VP_sortbyCTQSumm_4ISC) # this variable defines the number of cells e.g., pictures presented
# CharContainingNumberFromVector<-4 #Character containing number from variables coding vector.
# #e.g., EDA1, EDA2, EDA3... the character is 4: 1, 2, 3...
# ### Select the variables needed for the loop: VP, Variable sorting, Vector values
# VariableType<-"ISC_SCR_RSMSingleValue"
# sortType<- "CTQScores"
# OriginalData<-SCRRaw_SingleAlina_VP_sortbyCTQSumm_4ISC
# 
#-----------testing data end

RSMVector<-function(OriginalData,nVariables,CharContainingNumberFromVector,VariableType,sortType){
  
  
  sp<-OriginalData%>%
    dplyr::select(1,2, 3:ncol(OriginalData))
  
  ##### Matrices necessary.
  Vp_count<- table(OriginalData[1])
  Vp_count<-as.data.frame(Vp_count)
  names(Vp_count)[1]<-"VP"
  Vp_count<-as_tibble(Vp_count)
  
  ##Create a 3d matrix where the correlations will be saved. 
  corrMatrixVp<- array(0, dim = c(nVariables,nVariables,nrow(Vp_count)))
  
  
  for(i in 1: length(Vp_count$VP)) {
    #################################################
    #Create Matrix for Physiological variable for each participant
    #################################################
    # Single<-sp%>%
    #   filter(.[[1]] == Vp_count$VP[i])
    Single<-sp%>%
      filter_at(c(1) , all_vars(.==Vp_count$VP[i]))
    
    gaSingle<- Single%>%
      gather(output, value, 3:ncol(Single))%>%
      mutate(output, num= str_sub(output, CharContainingNumberFromVector))%>% #change the name to numeric to play around with the vector 
      dplyr::select(-output)%>%
      dplyr::arrange (.[[1]],"num")
    gaSingle$num<- as.numeric( gaSingle$num)
    spSingle<- gaSingle%>%
      spread(sort,value)
    spSingle<-spSingle%>%
      arrange_at ("num")%>%ungroup()%>%
      dplyr::select(3:ncol(spSingle))
    corrspSingle<- rcorr(as.matrix(spSingle), type="pearson")
    #put all correlations in a matrix
    corrMatrixVp[,,i]<-corrspSingle$r
    corr_prep<- corrspSingle$r%>% #corr_matrix_pic_EDA_aro_min$r %>% #corr_matrix_pic_EDA_aro
      as_tibble(rownames = "var1") %>%
      gather(var2, value, -var1) 
    corr_prep$var1<-as.numeric(corr_prep$var1)
    corr_prep$var2<-as.numeric(corr_prep$var2)
    VP_name<- Vp_count$VP[i]
    Vp_name<-as.character(Vp_count$VP[i])
    new_name<- paste("VP",Vp_name, sep= "")
    corr_prep<-corr_prep%>%
      arrange_at("var1")%>%
      arrange_at("var2")
    names(corr_prep)[3]<- new_name
    if(i==1){
      corr_prep_all<- corr_prep
    }else{corr_prep_all<- inner_join(corr_prep_all, corr_prep, by= c("var1"= "var1", "var2"= "var2"))}
  }
  
  corr_prep_all$mean_all<- rowMeans(subset(corr_prep_all,select= c(3:ncol(corr_prep_all)), na.rm = FALSE))
  
  corrMatrixVpName<-paste(VariableType,"corrMatrixVp",sep= "_")
  corr_prep_allName<-paste(VariableType,"corr_prep_all",sep= "_")
  
  assign(corr_prep_allName, corr_prep_all,envir = globalenv())
  assign(corrMatrixVpName, corrMatrixVp,envir = globalenv())
  ################################################################
  #-------------- Plot RSM Matrix
  mypalette<- colors()[c(73, 30, 28, 124, 635, 86,384,144, 148)]
  
  corr_prepPlot<- corr_prep_all%>%
    dplyr::select(var1, var2, mean_all)
  min<- floor_dec(min(corr_prepPlot[,3]),5)
  max<-ceiling_dec( max(corr_prepPlot[,3][corr_prepPlot[,3] != max(corr_prepPlot[,3])]),5) 
  cols <- rev(rainbow(20)[-20])
  
  plotRSM<-ggplot(data=corr_prepPlot, aes(x = var1, y = var2, fill = mean_all)) + #value2
    geom_tile() +
    labs(x = paste("Trials sorted by", sortType), y = paste("Trials sorted by", sortType), fill = "", title = "Similarity Matrix") +
    coord_fixed() +    
    theme_minimal()+
    scale_y_continuous(trans = "reverse")+
    scale_fill_gradientn(colours=rev(mypalette),na.value = colors() [c(73)], #heat.colors(12)
                         breaks=c(min, max),labels=c(min, max),
                         limits=c(min,max))+ theme(legend.position="bottom")+ 
    theme(text=element_text(family="Times", size=12))
  
  return(plotRSM)
  
  
  
}

#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#FUNCTION CREATE TO COMPARE AVERAGED RSMs based on single values.
# OriginalData<-BehavCTQAlina

RSMSingleValue<-function(OriginalData,nVariables,CharContainingNumberFromVector,VariableType,sortType,MaxValue){
  
  
  sp<-OriginalData%>%
    dplyr::select(1,2, 3:ncol(OriginalData))
  # OriginalData[1]<-as.character(OriginalData[1])
  ##### Matrices necessary.
  Vp_count<- table(OriginalData[1])
  Vp_count<-as.data.frame(Vp_count)
  names(Vp_count)[1]<-"VP"
  Vp_count<-as_tibble(Vp_count)
  
  ##Create a 3d matrix where the correlations will be saved. 
  corrMatrixVp<- array(0, dim = c(nVariables,nVariables,nrow(Vp_count)))
  corrMatrixVpSingle<- as_tibble(data.frame(matrix(data= 0, nrow=nVariables, ncol=nVariables))) 
  # class(corrMatrixVp)
  for(i in 1: length(Vp_count$VP)) {
    #################################################
    #Create Matrix for Physiological variable for each participant
    #################################################
    Single<-sp%>%
      dplyr::filter(.[[1]] == Vp_count$VP[i])
    gaSingle<- Single%>%
      spread(sort,ncol(Single))
    for(ii in 1:(length(gaSingle)-1)) {
      for (iii in 1:(length(gaSingle)-1)) {
        corrMatrixVpSingle[ii,iii]= MaxValue-(abs(gaSingle[1,ii+1]-gaSingle[1,iii+1]))
      }
    }
    corrMatrixVpSingle2<-as.matrix(corrMatrixVpSingle)
    corrMatrixVp[,,i]<-corrMatrixVpSingle2
    corr_prep<- corrMatrixVpSingle%>% 
      as_tibble(rownames = "var1") %>%
      gather(var2, value, -var1) 
    corr_prep$var2<-str_sub(corr_prep$var2,2,-1L)
    corr_prep$var1<-as.numeric(corr_prep$var1)
    corr_prep$var2<-as.numeric(corr_prep$var2)
    VP_name<- Vp_count$VP[i]
    Vp_name<-as.character(Vp_count$VP[i])
    new_name<- paste("VP",Vp_name, sep= "")
    corr_prep<-corr_prep%>%
      arrange_at("var1")%>%
      arrange_at("var2")
    names(corr_prep)[3]<- new_name
    if(i==1){
      corr_prep_all<- corr_prep
    }else{corr_prep_all<- inner_join(corr_prep_all, corr_prep, by= c("var1"= "var1", "var2"= "var2"))}
  }
  corr_prep_all<-corr_prep_all[ , colSums(is.na(corr_prep_all)) == 0] #Remove columns with NA N = 37
  corr_prep_all$mean_all<- rowMeans(subset(corr_prep_all,select= c(3:ncol(corr_prep_all)), na.rm = T))
  
  corrMatrixVpName<-paste(VariableType,"corrMatrixVp",sep= "_")
  corr_prep_allName<-paste(VariableType,"corr_prep_all",sep= "_")
  
  assign(corr_prep_allName, corr_prep_all,envir = globalenv())
  assign(corrMatrixVpName, corrMatrixVp,envir = globalenv())
  ################################################################
  #-------------- Plot RSM Matrix
  mypalette<- colors()[c(73, 30, 28, 124, 635, 86,384,144, 148)]
  # n_colors<-15
  # mypalette <- viridis_pal(option = "D")(n_colors)               # Apply viridis_pal function
  corr_prepPlot<- corr_prep_all%>%
    dplyr::select(var1, var2, mean_all)
  min<- floor_dec(min(corr_prepPlot[,3]),5)
  max<-ceiling_dec( max(corr_prepPlot[,3][corr_prepPlot[,3] != max(corr_prepPlot[,3])]),5) 
  cols <- rev(rainbow(20)[-20])
  
  plotRSM<-ggplot(data=corr_prepPlot, aes(x = var1, y = var2, fill = mean_all)) + #value2
    geom_tile() +
    labs(x = paste("Trials sorted by", sortType), y = paste("Trials sorted by", sortType), fill = "", title = "Similarity Matrix") +
    coord_fixed() +    
    theme_minimal()+
    scale_y_continuous(trans = "reverse")+
    scale_fill_gradientn(colours=rev(mypalette),na.value = colors() [c(73)], #heat.colors(12)
                         breaks=c(min, max),labels=c(min, max),
                         limits=c(min,max))+ theme(legend.position="bottom")+ 
    theme(text=element_text(family="Times", size=12))
  
  return(plotRSM)
  
  
  
}

#----------------------------------------------------------------------------

RSMSingleValueAK<-function(OriginalData,nVariables,CharContainingNumberFromVector,VariableType,sortType,MaxValue){
  
  
  sp<-OriginalData%>%
    dplyr::select(1,2, 3:ncol(OriginalData))
  
  ##### Matrices necessary.
  Vp_count<- table(OriginalData[1])
  Vp_count<-as.data.frame(Vp_count)
  names(Vp_count)[1]<-"VP"
  Vp_count<-as_tibble(Vp_count)
  
  ##Create a 3d matrix where the correlations will be saved. 
  corrMatrixVp<- array(0, dim = c(nVariables,nVariables,nrow(Vp_count)))
  corrMatrixVpSingle<- as_tibble(data.frame(matrix(data= 0, nrow=nVariables, ncol=nVariables))) 
  # class(corrMatrixVp)
  for(i in 1: length(Vp_count$VP)) {
    #################################################
    #Create Matrix for Physiological variable for each participant
    #################################################
    Single<-sp%>%
      filter(.[[1]] == Vp_count$VP[i])
    gaSingle<- Single%>%
      spread(sort,ncol(Single))
    for(ii in 1:(length(gaSingle)-1)) {
      for (iii in 1:(length(gaSingle)-1)) {
        corrMatrixVpSingle[ii,iii]= ((gaSingle[1,ii+1]+gaSingle[1,iii+1])/2)/nVariables
      }
    }
    corrMatrixVpSingle2<-as.matrix(corrMatrixVpSingle)
    corrMatrixVp[,,i]<-corrMatrixVpSingle2
    corr_prep<- corrMatrixVpSingle%>% 
      as_tibble(rownames = "var1") %>%
      gather(var2, value, -var1) 
    corr_prep$var2<-str_sub(corr_prep$var2,2,-1L)
    corr_prep$var1<-as.numeric(corr_prep$var1)
    corr_prep$var2<-as.numeric(corr_prep$var2)
    VP_name<- Vp_count$VP[i]
    Vp_name<-as.character(Vp_count$VP[i])
    new_name<- paste("VP",Vp_name, sep= "")
    corr_prep<-corr_prep%>%
      arrange_at("var1")%>%
      arrange_at("var2")
    names(corr_prep)[3]<- new_name
    if(i==1){
      corr_prep_all<- corr_prep
    }else{corr_prep_all<- inner_join(corr_prep_all, corr_prep, by= c("var1"= "var1", "var2"= "var2"))}
  }
  corr_prep_all<-corr_prep_all[ , colSums(is.na(corr_prep_all)) == 0] #Remove columns with NA N = 37
  corr_prep_all$mean_all<- rowMeans(subset(corr_prep_all,select= c(3:ncol(corr_prep_all)), na.rm = T))
  
  corrMatrixVpName<-paste(VariableType,"corrMatrixVp",sep= "_")
  corr_prep_allName<-paste(VariableType,"corr_prep_all",sep= "_")
  
  assign(corr_prep_allName, corr_prep_all,envir = globalenv())
  assign(corrMatrixVpName, corrMatrixVp,envir = globalenv())
  ################################################################
  #-------------- Plot RSM Matrix
  mypalette<- colors()[c(73, 30, 28, 124, 635, 86,384,144, 148)]
  
  corr_prepPlot<- corr_prep_all%>%
    dplyr::select(var1, var2, mean_all)
  min<- floor_dec(min(corr_prepPlot[,3]),5)
  max<-ceiling_dec( max(corr_prepPlot[,3][corr_prepPlot[,3] != max(corr_prepPlot[,3])]),4) 
  cols <- rev(rainbow(20)[-20])
  
  plotRSM<-ggplot(data=corr_prepPlot, aes(x = var1, y = var2, fill = mean_all)) + #value2
    geom_tile() +
    labs(x = paste("Trials sorted by", sortType), y = paste("Trials sorted by", sortType), fill = "", title = "Similarity Matrix") +
    coord_fixed() +    
    theme_minimal()+
    scale_y_continuous(trans = "reverse")+
    scale_fill_gradientn(colours=rev(mypalette),na.value = colors() [c(73)], #heat.colors(12)
                         breaks=c(min, max),labels=c(min, max),
                         limits=c(min,max))+ theme(legend.position="bottom")+ 
    theme(text=element_text(family="Times", size=12))
  
  return(plotRSM)
  
  
  
}

#----------------------------------------------------------------------------
# OriginalData<-BehavCTQAlina
RSMSingleValueinvAK<-function(OriginalData,nVariables,CharContainingNumberFromVector,VariableType,sortType,MaxValue){
  
  
  sp<-OriginalData%>%
    dplyr::select(1,2, 3:ncol(OriginalData))
  
  ##### Matrices necessary.
  Vp_count<- table(OriginalData[1])
  Vp_count<-as.data.frame(Vp_count)
  names(Vp_count)[1]<-"VP"
  Vp_count<-as_tibble(Vp_count)
  
  ##Create a 3d matrix where the correlations will be saved. 
  corrMatrixVp<- array(0, dim = c(nVariables,nVariables,nrow(Vp_count)))
  corrMatrixVpSingle<- as_tibble(data.frame(matrix(data= 0, nrow=nVariables, ncol=nVariables))) 
  # class(corrMatrixVp)
  for(i in 1: length(Vp_count$VP)) {
    #################################################
    #Create Matrix for Physiological variable for each participant
    #################################################
    Single<-sp%>%
      filter(.[[1]] == Vp_count$VP[i])
    gaSingle<- Single%>%
      spread(sort,ncol(Single))
    for(ii in 1:(length(gaSingle)-1)) {
      for (iii in 1:(length(gaSingle)-1)) {
        corrMatrixVpSingle[ii,iii]= (1/((gaSingle[1,ii+1]+gaSingle[1,iii+1])/2))/nVariables
      }
    }
    corrMatrixVpSingle2<-as.matrix(corrMatrixVpSingle)
    corrMatrixVp[,,i]<-corrMatrixVpSingle2
    corr_prep<- corrMatrixVpSingle%>% 
      as_tibble(rownames = "var1") %>%
      gather(var2, value, -var1) 
    corr_prep$var2<-str_sub(corr_prep$var2,2,-1L)
    corr_prep$var1<-as.numeric(corr_prep$var1)
    corr_prep$var2<-as.numeric(corr_prep$var2)
    VP_name<- Vp_count$VP[i]
    Vp_name<-as.character(Vp_count$VP[i])
    new_name<- paste("VP",Vp_name, sep= "")
    corr_prep<-corr_prep%>%
      arrange_at("var1")%>%
      arrange_at("var2")
    names(corr_prep)[3]<- new_name
    if(i==1){
      corr_prep_all<- corr_prep
    }else{corr_prep_all<- inner_join(corr_prep_all, corr_prep, by= c("var1"= "var1", "var2"= "var2"))}
  }
  corr_prep_all<-corr_prep_all[ , colSums(is.na(corr_prep_all)) == 0] #Remove columns with NA N = 37
  # corr_prep_alltest<-corr_prep_all%>%
  #   filter(is.na(VP1)) #column3
  corr_prep_all$mean_all<- rowMeans(subset(corr_prep_all,select= c(3:ncol(corr_prep_all)), na.rm = T))
  
  corrMatrixVpName<-paste(VariableType,"corrMatrixVp",sep= "_")
  corr_prep_allName<-paste(VariableType,"corr_prep_all",sep= "_")
  
  assign(corr_prep_allName, corr_prep_all,envir = globalenv())
  assign(corrMatrixVpName, corrMatrixVp,envir = globalenv())
  ################################################################
  #-------------- Plot RSM Matrix
  mypalette<- colors()[c(73, 30, 28, 124, 635, 86,384,144, 148)]
  
  corr_prepPlot<- corr_prep_all%>%
    dplyr::select(var1, var2, mean_all)
  min<- floor_dec(min(corr_prepPlot[,3]),6)
  max<-ceiling_dec( max(corr_prepPlot[,3][corr_prepPlot[,3] != max(corr_prepPlot[,3])]),6) 
  cols <- rev(rainbow(20)[-20])
  
  plotRSM<-ggplot(data=corr_prepPlot, aes(x = var1, y = var2, fill = mean_all)) + #value2
    geom_tile() +
    labs(x = paste("Trials sorted by", sortType), y = paste("Trials sorted by", sortType), fill = "", title = "Similarity Matrix") +
    coord_fixed() +    
    theme_minimal()+
    scale_y_continuous(trans = "reverse")+
    scale_fill_gradientn(colours=rev(mypalette),na.value = colors() [c(73)], #heat.colors(12)
                         breaks=c(min, max),labels=c(min, max),
                         limits=c(min,max))+ theme(legend.position="bottom")+ 
    theme(text=element_text(family="Times", size=12))
  
  return(plotRSM)
  
  
  
}
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
#-------- FUNCTION CREATED TO CALCULTE RSM VALUES BASES ON FACTOR OR CHARACTER VARIABLE
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
# OriginalData<-ImagBehavArosCat

RSMSingleValueFactor<-function(OriginalData,nVariables,CharContainingNumberFromVector,VariableType,sortType,MaxValue){
  
  
  sp<-OriginalData%>%
    dplyr::select(1,2, 3:ncol(OriginalData))
  # OriginalData[1]<-as.character(OriginalData[1])
  ##### Matrices necessary.
  Vp_count<- table(OriginalData[1])
  Vp_count<-as.data.frame(Vp_count)
  names(Vp_count)[1]<-"VP"
  Vp_count<-as_tibble(Vp_count)
  
  ##Create a 3d matrix where the correlations will be saved. 
  corrMatrixVp<- array(0, dim = c(nVariables,nVariables,nrow(Vp_count)))
  corrMatrixVpSingle<- as_tibble(data.frame(matrix(data= 0, nrow=nVariables, ncol=nVariables))) 
  # class(corrMatrixVp)
  for(i in 1: length(Vp_count$VP)) {
    #################################################
    #Create Matrix for Physiological variable for each participant
    #################################################
    Single<-sp%>%
      dplyr::filter(.[[1]] == Vp_count$VP[i])
    gaSingle<- Single%>%
      spread(sort,ncol(Single))
    for(ii in 1:(length(gaSingle)-1)) {
      for (iii in 1:(length(gaSingle)-1)) {
        corrMatrixVpSingle[ii,iii]= ifelse(gaSingle[1,ii+1]==gaSingle[1,iii+1],1,0)
      }
    }
    corrMatrixVpSingle2<-as.matrix(corrMatrixVpSingle)
    corrMatrixVp[,,i]<-corrMatrixVpSingle2
    corr_prep<- corrMatrixVpSingle%>% 
      as_tibble(rownames = "var1") %>%
      gather(var2, value, -var1) 
    corr_prep$var2<-str_sub(corr_prep$var2,2,-1L)
    corr_prep$var1<-as.numeric(corr_prep$var1)
    corr_prep$var2<-as.numeric(corr_prep$var2)
    VP_name<- Vp_count$VP[i]
    Vp_name<-as.character(Vp_count$VP[i])
    new_name<- paste("VP",Vp_name, sep= "")
    corr_prep<-corr_prep%>%
      arrange_at("var1")%>%
      arrange_at("var2")
    names(corr_prep)[3]<- new_name
    if(i==1){
      corr_prep_all<- corr_prep
    }else{corr_prep_all<- inner_join(corr_prep_all, corr_prep, by= c("var1"= "var1", "var2"= "var2"))}
  }
  corr_prep_all<-corr_prep_all[ , colSums(is.na(corr_prep_all)) == 0] #Remove columns with NA N = 37
  corr_prep_all$mean_all<- rowMeans(subset(corr_prep_all,select= c(3:ncol(corr_prep_all)), na.rm = T))
  
  corrMatrixVpName<-paste(VariableType,"corrMatrixVp",sep= "_")
  corr_prep_allName<-paste(VariableType,"corr_prep_all",sep= "_")
  
  assign(corr_prep_allName, corr_prep_all,envir = globalenv())
  assign(corrMatrixVpName, corrMatrixVp,envir = globalenv())
  ################################################################
  #-------------- Plot RSM Matrix
  mypalette<- colors()[c(73, 30, 28, 124, 635, 86,384,144, 148)]
  # n_colors<-15
  # mypalette <- viridis_pal(option = "D")(n_colors)               # Apply viridis_pal function
  corr_prepPlot<- corr_prep_all%>%
    dplyr::select(var1, var2, mean_all)
  min<- floor_dec(min(corr_prepPlot[,3]),5)
  max<-ceiling_dec( max(corr_prepPlot[,3][corr_prepPlot[,3] != max(corr_prepPlot[,3])]),5) 
  cols <- rev(rainbow(20)[-20])
  
  plotRSM<-ggplot(data=corr_prepPlot, aes(x = var1, y = var2, fill = mean_all)) + #value2
    geom_tile() +
    labs(x = paste("Trials sorted by", sortType), y = paste("Trials sorted by", sortType), fill = "", title = "Similarity Matrix") +
    coord_fixed() +    
    theme_minimal()+
    scale_y_continuous(trans = "reverse")+
    scale_fill_gradientn(colours=rev(mypalette),na.value = colors() [c(73)], #heat.colors(12)
                         breaks=c(min, max),labels=c(min, max),
                         limits=c(min,max))+ theme(legend.position="bottom")+ 
    theme(text=element_text(family="Times", size=12))
  
  return(plotRSM)
  
  
  
}

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
################################################################
#-------------- Plot RSM Matrix
# dataOriginal#data coming up from RSM functions it should contain the variable mean_all
# mypalette#colors
# minimum#scale
# maximum# scale
# Category# VP or trials
# sortType #variable to sort on

PlotMatrix<-function(dataOriginal, minimum,maximum,Category,sortType, mypalette){
 if (length(mypalette)==0){
   
mypalette<- colors()[c(73, 30, 28, 124, 635, 86,384,144, 148)]
 } 
corr_prepPlot<- dataOriginal%>%
  dplyr::select(var1, var2, mean_all)
if (length(minimum)==0){
minimum<- floor_dec(min(corr_prepPlot[,3]),5)
}
if (length(maximum)==0){
maximum<-ceiling_dec( max(corr_prepPlot[,3][corr_prepPlot[,3] != max(corr_prepPlot[,3])]),5) 
}
cols <- rev(rainbow(20)[-20])

plotRSM<-ggplot(data=corr_prepPlot, aes(x = var1, y = var2, fill = mean_all)) + #value2
  geom_tile() +
  labs(x = paste( Category, "sorted by", sortType), y = paste(Category, "sorted by", sortType), fill = "", title = "Similarity Matrix") +
  coord_fixed() +    
  theme_minimal()+
  scale_y_continuous(trans = "reverse")+
  scale_fill_gradientn(colours=rev(mypalette),na.value = colors() [c(73)], #heat.colors(12)
                       breaks=c(minimum, maximum),labels=c(minimum, maximum),
                       limits=c(minimum,maximum))+ theme(legend.position="bottom")+ 
  theme(text=element_text(family="Times", size=12))

return(plotRSM)
}

#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
PermutationTest<-function(OriginalData,ModelData,VariableOriginal,VariableModel,nPermutations,seedPer){
  OriginalData<-OriginalData%>%dplyr::select(var1,var2,mean_all)%>%
    filter(OriginalData$var1>OriginalData$var2)%>%
    filter(var1 > 0)
  OriginalData$RSA<-VariableOriginal
  ModelData<-ModelData%>%dplyr::select(var1,var2,mean_all)%>%
    filter(ModelData$var1>ModelData$var2)%>%
    filter(var1 > 0)
  ModelData$RSA<- VariableModel
  
  
  CombinedOriginalModel<- rbind(OriginalData, ModelData)
  
  
  
  test.cor<- cor(CombinedOriginalModel$mean_all[CombinedOriginalModel$RSA==VariableOriginal],
                 CombinedOriginalModel$mean_all[CombinedOriginalModel$RSA==VariableModel], method= "spearman")
  # test.cor.bf = correlationBF(CombinedOriginalModel$mean_all[CombinedOriginalModel$RSA==VariableOriginal],
  #                    CombinedOriginalModel$mean_all[CombinedOriginalModel$RSA==VariableModel"])
  # 
  # samples = correlationBF(y = CombinedOriginalModel$mean_all[CombinedOriginalModel$RSA==VariableOriginal"],
  #                         x =  CombinedOriginalModel$mean_all[CombinedOriginalModel$RSA==VariableModel],
  #                         posterior = TRUE, iterations = 10000)
  # plot(samples[,"rho"])
  
  set.seed(seedPer)  
  
  Perm.test.corr <- rep(0, nPermutations) #matrix were to save the permutations
  n <- length(CombinedOriginalModel$RSA)  
  # the variable we will resample from 
  #     (note, could use the labels(RSA) too, and "shuffle this")
  variable <- CombinedOriginalModel$RSA
  # initialize a matrix to store the permutation of the name
  PermSamplesOther <- matrix(0, nrow=n, ncol=nPermutations)
  # now, get those permutation samples from the labels, using a loop
  for(i in 1:nPermutations){
    PermSamplesOther[,i] <- sample(variable, size= n, replace=FALSE)
  }
  
  #Loop to perform the nPermutations permutations
  for (i in 1:nPermutations){
    
    # calculate the perm-test-stat1 and save it
    Perm_test_corr_sort<- CombinedOriginalModel%>%
      select ( - RSA)
    Perm<- PermSamplesOther[,i]
    Perm_test_corr_sort <- data.frame(Perm_test_corr_sort, Perm)
    Perm_test_corr_sort<-as_tibble(Perm_test_corr_sort)%>%
      arrange(var1)%>%
      arrange(var2)%>%
      arrange(Perm)
    
    Perm.test.corr[i] <- cor(Perm_test_corr_sort$mean_all[Perm_test_corr_sort$Perm==VariableOriginal],
                             Perm_test_corr_sort$mean_all[Perm_test_corr_sort$Perm==VariableModel],method = c("spearman"))
  }
  
  
  
  
  Perm.test.corr<- as_tibble(Perm.test.corr)
  Perm.test.corr<- Perm.test.corr%>%
    arrange_at(1)
  

  lim<-Perm.test.corr[nPermutations/100*95,1]
  lim<-as.numeric(lim)
  plot.perm<-
    ggplot(Perm.test.corr, aes(x=value)) +
    geom_histogram(
        bins= 20,
       fill=I("blue"), 
        col=I("black"), 
        alpha=I(.2))+
    geom_vline(xintercept=lim, linetype="dashed", 
               color = "black", size=2 )+
    geom_vline(xintercept=test.cor,
               color = "red", size=2 )+ theme(text=element_text(family="Times", size=12),
                                              axis.text=element_text(size=25),axis.line.y= element_line(),
                                              axis.line.x= element_line())+
    
    theme(legend.position="bottom")+
    theme(#panel.grid.minor = element_blank(),
      panel.background = element_blank())
  significanceTest<- which(Perm.test.corr>test.cor)
  
  permPValue<-(nPermutations-significanceTest[1])/nPermutations
  
  
  bf =correlationBF(CombinedOriginalModel$mean_all[CombinedOriginalModel$RSA==VariableOriginal],
                    CombinedOriginalModel$mean_all[CombinedOriginalModel$RSA==VariableModel],nullInterval = c(Perm.test.corr[1,], Perm.test.corr[nPermutations/100*95,]))
  
  # bf
  bfevidence<- bf[2] / bf[1]
  # View(bfevidence)
  
  template <- "A permutation Test was was conducted to compare {VariableOriginal} and  {VariableModel}.
The relation between both RMS was: {test.cor}. The permutation test with percentil 95 was {lim}, the relationship had a p-value of {permPValue}.
the BF10 index was: {bfevidence} "
  
  
  template2<-glue::glue(
    template,
    
    VariableOriginal   = VariableOriginal,
    VariableModel     = VariableModel,
    test.cor   =round( test.cor,2),
    lim     = round(lim,2),
    permPValue     = permPValue,
    bfevidence= round(bfevidence@bayesFactor[["bf"]],2)
  )
  
  mylist<-list (plot.perm, template2)
    # return(mylist)
  return(mylist)
  
}


# 
# OriginalData<- Corr_SCRinvAKAlinaSingleValue_corr_prep_all
# ModelData1<-CTQ_SumScoresinvAK_corr_prep_all
# ModelData2<- CTQ_SumScores_corr_prep_all
# whichModel<- c("ModelData1","ModelData2")
# VariableOriginal<-"SCRinvAKAroSingleValueCorrelationinterindividualCTQ"
# VariableModel2<- "CTQ"
# VariableModel1<- "invAKCTQ"
# Category = "VP"

#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#----------------- FUNCTION CREATED TO run compare RSMs with each other for each individual separately. 

# tableData<-SCR_SingleAlina_corr_prep_all
# tableModel<-BehavioralSingleAlina_invAK_corr_prep_all
RSMIndividual<-function(tableModel,tableData,Model,Data){
  
  #Select the biggest table
  if (ncol(tableModel)>ncol(tableData)){
    tableMother<-tableData
    tableDaugther<-tableModel
  }else {tableMother<-tableModel
  tableDaugther<-tableData}
  # tableCorrelation<-as_tibble(vector(mode='numeric',length=ncol(tableMother)))
  tableCorrelation<-as_tibble(matrix(ncol=2, nrow=ncol(tableMother)-3))
  names(tableCorrelation)<-c("VP","Correlation")
  tableDaugther2<-tableDaugther%>%
    filter(tableDaugther$var1>tableDaugther$var2)
  tableMother2<-tableMother%>%
    filter(tableMother$var1>tableMother$var2)
  for(i in 3:(ncol(tableMother)-1)) { #starts in column three
    
    single<-tableMother2[,i]
    single[,2]<-tableDaugther2[,which( colnames(tableDaugther2)==names(tableMother2[i]))]
    corrspSingle<- rcorr(as.matrix(single), type="spearman")
    tableCorrelation[i-2,2]<-corrspSingle$r[2]
    tableCorrelation[i-2,1]<-names(tableMother2[i])
  }
  
  tableCorrelation<-tableCorrelation%>%drop_na()
  corrMatrixVpName<-paste(Model,Data,"tableCorrelation",sep= "_")
  assign(corrMatrixVpName, tableCorrelation,envir = globalenv())
  
  tableCorrelationPlot<-tableCorrelation$Correlation
  tableCorrelationPlot<-as_tibble(tableCorrelationPlot)
  tableCorrelationPlot$x<-1
  onewayTTest<- t.test(tableCorrelation$Correlation, mu = 0, alternative = "two.sided")
  
  plotDistribution <- ggplot(tableCorrelationPlot, aes(y = value,x = x)) +
    geom_violin(position = position_dodge(width = 0.9)) +
    geom_quasirandom(dodge.width = 0.9, varwidth = TRUE)+
    stat_summary(fun="mean",color="red", geom="point", size=5)+
    # annotate("text", x=.75, y=max(tableCorrelation$value)+mean(tableCorrelation$value), label= paste("p-value:",onewayTTest[["p.value"]], sep = " "),
    #          col="blue", size=5, parse=TRUE)+
    
    annotate("text", x=.75, y=max(tableCorrelationPlot$value), label= paste("p-value:",onewayTTest[["p.value"]], sep = " "),
             col="blue", size=5, parse=TRUE)+
    
    annotate("text", x=.75, y=mean(tableCorrelationPlot$value),#*3,
             label= paste("t-value:",onewayTTest[["statistic"]], sep = " "),
             col="blue", size=5, parse=TRUE)+
    annotate("text", x=.75, y=min(tableCorrelationPlot$value),#*4, 
             label= paste("Mean:",mean(tableCorrelationPlot$value), sep = " "),
             col="blue", size=5, parse=TRUE)
  plotDistribution
}

#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#----------------- FUNCTION CREATED TO run compare RSMs with each other for each individual separately, coloring dots based on a third variable
# # 
# tableData<-SCR_SingleAlina_corr_prep_all
# tableModel<-BehavioralSingleAlina_SV_invAK_corr_prep_all
# tableColorData<-SCR_SingleAlina_corr_prep_allfilt_sp
# VariableColor<- "Q_CTQ_sum"
RSMIndividualCol<-function(tableModel,tableData,Model,Data,tableColorData,VariableColor){
  
  #Select the biggest table
  if (ncol(tableModel)>ncol(tableData)){
    tableMother<-tableData
    tableDaugther<-tableModel
  }else {tableMother<-tableModel
  tableDaugther<-tableData}
  # tableCorrelation<-as_tibble(vector(mode='numeric',length=ncol(tableMother)))
  tableCorrelation<-as_tibble(matrix(ncol=2, nrow=ncol(tableMother)-3))
  names(tableCorrelation)<-c("VP","Correlation")
  tableDaugther2<-tableDaugther%>%
    filter(tableDaugther$var1>tableDaugther$var2)
  tableMother2<-tableMother%>%
    filter(tableMother$var1>tableMother$var2)
  for(i in 3:(ncol(tableMother)-1)) { #starts in column three
    
    single<-tableMother2[,i]
    single[,2]<-tableDaugther2[,which( colnames(tableDaugther2)==names(tableMother2[i]))]
    corrspSingle<- rcorr(as.matrix(single), type="spearman")
    tableCorrelation[i-2,2]<-corrspSingle$r[2]
    tableCorrelation[i-2,1]<-names(tableMother2[i])
  }
  
  tableCorrelation<-tableCorrelation%>%drop_na()
  tableColorData2<-tableColorData%>%select(names(tableCorrelation[1]),VariableColor)
  names(tableColorData2)= c("VP","Q")
  tableCorrelation<-inner_join(tableCorrelation,tableColorData2, by = names(tableCorrelation[1]))
  corrMatrixVpName<-paste(Model,Data,"coloredby",VariableColor,"tableCorrelation",sep= "_")
  assign(corrMatrixVpName, tableCorrelation,envir = globalenv())
  
  tableCorrelationPlot<-tableCorrelation%>% select(-VP)
  tableCorrelationPlot$x<-1
  
  onewayTTest<- t.test(tableCorrelation$Correlation, mu = 0, alternative = "two.sided")
  
  plotDistribution <- ggplot(tableCorrelationPlot, aes(y = Correlation,x = x)) +
    geom_violin(position = position_dodge(width = 0.9)) +
    geom_quasirandom(
      aes(colour = Q,
          size = 3 ,
          # dodge.width = 0.9, varwidth = TRUE
      ),
    )+
    stat_summary(fun="mean",color="red", geom="point", size=5)+
    # geom_point(aes(colour = cut(Q_CTQ_sum, c(-Inf, 30, 45, Inf))),
    #            size = 3) +
    # geom_tile(aes(fill = Q))+
    scale_colour_gradientn(colours = terrain.colors(10))+
    # annotate("text", x=.75, y=max(tableCorrelation$value)+mean(tableCorrelation$value), label= paste("p-value:",onewayTTest[["p.value"]], sep = " "),
    #          col="blue", size=5, parse=TRUE)+
    
    annotate("text", x=.75, y=max(tableCorrelationPlot$Correlation), label= paste("p-value:",onewayTTest[["p.value"]], sep = " "),
             col="blue", size=5, parse=TRUE)+
    
    annotate("text", x=.75, y=mean(tableCorrelationPlot$Correlation),#*3,
             label= paste("t-value:",onewayTTest[["statistic"]], sep = " "),
             col="blue", size=5, parse=TRUE)+
    annotate("text", x=.75, y=min(tableCorrelationPlot$Correlation),#*4, 
             label= paste("Mean:",mean(tableCorrelationPlot$Correlation), sep = " "),
             col="blue", size=5, parse=TRUE)
  plotDistribution
}

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#---------------- RSM Individual regression
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# tableData<-SCR_Single_corr_prep_all
#   tableModel1<-Behav_Single_invAK_corr_prep_all
#   tableModel2<-Behav_Single_corr_prep_all
#   variableOriginal<-"SCRSingle"
#   VariableModel1<-"AroRatingsinvAK"
#   VariableModel2<-"AroRatingsNN"
RSMIndividualRegression<-function(OriginalData,tableModel1,tableModel2,VariableOriginal,VariableModel1,VariableModel2){ 
  #Select the biggest table
  if (ncol(tableModel1)>ncol(tableData)){
    tableMother<-tableData
    tableDaugther<-tableModel1
  }else {tableMother<-tableModel1
  tableDaugther<-tableData}
  tableCousin<-tableModel2
  # tableCorrelation<-as_tibble(vector(mode='numeric',length=ncol(tableMother)))
  tableCorrelation<-as_tibble(matrix(ncol=3, nrow=ncol(tableMother)-3))
  names(tableCorrelation)<-c("VP",paste("Similarity",VariableOriginal,VariableModel1,sep="_"),
                             paste("Similarity",VariableOriginal,VariableModel,sep="_"))
  tableDaugther2<-tableDaugther%>%
    filter(tableDaugther$var1>tableDaugther$var2)
  tableMother2<-tableMother%>%
    filter(tableMother$var1>tableMother$var2)
  tableCousin2<-tableCousin%>%
    filter(tableCousin$var1>tableCousin$var2)
  for(i in 3:(ncol(tableMother)-1)) { #starts in column three
    single<-tableMother2[,i]
    single[,2]<-tableDaugther2[,which( colnames(tableDaugther2)==names(tableMother2[i]))]
    single[,3]<-tableCousin2[,which( colnames(tableCousin2)==names(tableMother2[i]))]
    names(single)<-c("Original","M1","M2")
    modelM1<- lm(M1~M2, data=single )
    modelM2<- lm(M2~M1, data=single )
    corrspSingleM1<-as_tibble(matrix(ncol=3, nrow=nrow(single)))
    corrspSingleM1[,1]<-single$Original
    corrspSingleM1[,2]<-modelM1[["residuals"]]
    corrspSingleM1[,3]<-modelM2[["residuals"]]
    corrOriModelsRes<-rcorr(as.matrix(corrspSingleM1), type="spearman")
    
    tableCorrelation[i-2,2]<-corrOriModelsRes$r[2]
    tableCorrelation[i-2,3]<-corrOriModelsRes$r[3]
    tableCorrelation[i-2,1]<-names(tableMother2[i])
  }
  
  tableCorrelation<-tableCorrelation%>%drop_na()
  corrMatrixVpName<-paste(variableOriginal,VariableModel1,VariableModel2,"tableRegression",sep= "_")
  assign(corrMatrixVpName, tableCorrelation,envir = globalenv())
  
  tableCorrelationPlot<-tableCorrelation%>% select(-VP)
  tableCorrelationPlot$x<-1
  tableCorrelationPlotM<-colMeans(tableCorrelationPlot)
  #Plot regression1
  onewayTTestM1<- t.test(tableCorrelation[,2], mu = 0, alternative = "two.sided")
  
  plotDistributionM1 <- ggplot(tableCorrelationPlot, aes_string(y = names(tableCorrelationPlot[1]),x = "x")) +
    geom_violin(position = position_dodge(width = 0.9)) +
    geom_quasirandom(dodge.width = 0.9, varwidth = TRUE)+
    stat_summary(fun="mean",color="red", geom="point", size=5)+
  
    scale_colour_gradientn(colours = terrain.colors(10))+
 
    annotate("text", x=.75, y=max(tableCorrelationPlot[,1]), label= paste("p-value:",onewayTTestM1[["p.value"]], sep = " "),
             col="blue", size=5, parse=TRUE)+
    
    annotate("text", x=.75, y=tableCorrelationPlotM[1],#*3,
             label= paste("t-value:",onewayTTestM1[["statistic"]], sep = " "),
             col="blue", size=5, parse=TRUE)+
    annotate("text", x=.75, y=min(tableCorrelationPlot[,1]),#*4, 
             label= paste("Mean:",tableCorrelationPlotM[1], sep = " "),
             col="blue", size=5, parse=TRUE)+
    ggtitle(names(tableCorrelationPlot[1]))
  
  
  #Plot regression1
  onewayTTestM2<- t.test(tableCorrelation[,3], mu = 0, alternative = "two.sided")
  
  plotDistributionM2 <- ggplot(tableCorrelationPlot, aes_string(y = names(tableCorrelationPlot[2]),x = "x")) +
    geom_violin(position = position_dodge(width = 0.9)) +
    geom_quasirandom(dodge.width = 0.9, varwidth = TRUE)+
    stat_summary(fun="mean",color="red", geom="point", size=5)+
    
    scale_colour_gradientn(colours = terrain.colors(10))+
    
    annotate("text", x=.75, y=max(tableCorrelationPlot[,2]), label= paste("p-value:",onewayTTestM2[["p.value"]], sep = " "),
             col="blue", size=5, parse=TRUE)+
    
    annotate("text", x=.75, y=tableCorrelationPlotM[2],#*3,
             label= paste("t-value:",onewayTTestM2[["statistic"]], sep = " "),
             col="blue", size=5, parse=TRUE)+
    annotate("text", x=.75, y=min(tableCorrelationPlot[,2]),#*4, 
             label= paste("Mean:",tableCorrelationPlotM[2], sep = " "),
             col="blue", size=5, parse=TRUE)+
    ggtitle(names(tableCorrelationPlot[2]))
  
  ggarrange(plotDistributionM1,plotDistributionM2 , 
            labels = c("A", "B"),
            ncol = 2, nrow = 1)
  
  
}
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#-----------------------------Permutation test with regression-----------------------------------------------
#----------------------------------------------------------------------------
# # #PermutationwithRegression
# OriginalData<- SCR_SingleAlina_corr_prep_all
# ModelData1<-BehavioralSingleAlina_SV_invAK_corr_prep_all
# ModelData2<- BehavioralSingleAlina_SV_NN_corr_prep_all
# whichModel<- c("ModelData1")#,"ModelData2")
# VariableOriginal<-"SCRinvAKAroSingleValue"
# VariableModel2<- "NNAro"
# VariableModel1<- "invAKAro"
# seedPer<-49
# nPermutations<-100
# PlotRSMReg<-1
# Category = "VP"

PermutationTestRegModels<-function(OriginalData,ModelData1,ModelData2,VariableOriginal,VariableModel1,VariableModel2,
                                   Category,nPermutations,seedPer,PlotRSMReg){ 
  #variable whichModel would define which model to operate. because that if does not work, both models are operated
  OriginalData<-OriginalData%>%dplyr::select(var1,var2,mean_all)%>%
    filter(OriginalData$var1>OriginalData$var2)%>%
    filter(var1 > 0)
  OriginalData$RSA<-VariableOriginal
  
  Models<- inner_join(ModelData1,ModelData2, by=c('var1', 'var2'),suffix= c("M1","M2"))
  modelM1<- lm(mean_allM1~mean_allM2, data=Models )
  modelM2<- lm(mean_allM2~mean_allM1, data=Models )
  for (ii in 1:length(whichModel)) {
    # if (whichModel[ii]== deparse(substitute(ModelData1))){ # does not work
    #Determine which permutation to run, depending on whichModel variable 
    #EachModelSeparately
    ModelData1_2Plot<-ModelData1
    ModelData1_2Plot$mean_all<-modelM1[["residuals"]]
    
    ModelData1$mean_all<-modelM1[["residuals"]]
    
    ModelData1s<-ModelData1%>%dplyr::select(var1,var2,mean_all)
    ModelData1s<-ModelData1s%>%
      filter(ModelData1s$var1>ModelData1s$var2)%>%
      filter(var1 > 0)
    ModelData1s$RSA<-VariableModel1
    
    
    CombinedOriginalModel1<- rbind(OriginalData, ModelData1s)
    test.cor<- cor(CombinedOriginalModel1$mean_all[CombinedOriginalModel1$RSA==VariableOriginal],
                   CombinedOriginalModel1$mean_all[CombinedOriginalModel1$RSA==VariableModel1], method= "spearman")
    
    
    
    set.seed(seedPer)  
    
    Perm.test.corr <- rep(0, nPermutations) #matrix were to save the permutations
    n <- length(CombinedOriginalModel1$RSA)  
    # the variable we will resample from 
    #     (note, could use the labels(RSA) too, and "shuffle this")
    variable <- CombinedOriginalModel1$RSA
    # initialize a matrix to store the permutation of the name
    PermSamplesOther <- matrix(0, nrow=n, ncol=nPermutations)
    # now, get those permutation samples from the labels, using a loop
    for(i in 1:nPermutations){
      PermSamplesOther[,i] <- sample(variable, size= n, replace=FALSE)
    }
    
    #Loop to perform the nPermutations permutations
    for (i in 1:nPermutations){
      # calculate the perm-test-stat1 and save it
      Perm_test_corr_sort<- CombinedOriginalModel1%>%
        select ( - RSA)
      
      Perm<- PermSamplesOther[,i]
      Perm_test_corr_sort <- data.frame(Perm_test_corr_sort, Perm)
      Perm_test_corr_sort<-as_tibble(Perm_test_corr_sort)%>%
        arrange(var1)%>%
        arrange(var2)%>%
        arrange(Perm)
      
      Perm.test.corr[i] <- cor(Perm_test_corr_sort$mean_all[Perm_test_corr_sort$Perm==VariableOriginal],
                               Perm_test_corr_sort$mean_all[Perm_test_corr_sort$Perm==VariableModel1],method = c("spearman"))
    }  
    
    
    
    
    Perm.test.corr<- as_tibble(Perm.test.corr)
    Perm.test.corr<- Perm.test.corr%>%
      arrange_at(1)
    
    
    lim<-Perm.test.corr[nPermutations/100*95,1]
    lim<-as.numeric(lim)
    plot.perm<-
      ggplot(Perm.test.corr, aes(x=value)) +
      geom_histogram(
        bins= 20,
        fill=I("blue"), 
        col=I("black"), 
        alpha=I(.2))+
      geom_vline(xintercept=lim, linetype="dashed", 
                 color = "black", size=2 )+
      geom_vline(xintercept=test.cor,
                 color = "red", size=2 )+ theme(text=element_text(family="Times", size=12),
                                                axis.text=element_text(size=25),axis.line.y= element_line(),
                                                axis.line.x= element_line())+
      
      theme(legend.position="bottom")+
      theme(#panel.grid.minor = element_blank(),
        panel.background = element_blank())
    significanceTest<- which(Perm.test.corr>test.cor)
    
    permPValue<-(nPermutations-significanceTest[1])/nPermutations
    
    
    bf =correlationBF(CombinedOriginalModel1$mean_all[CombinedOriginalModel1$RSA==VariableOriginal],
                      CombinedOriginalModel1$mean_all[CombinedOriginalModel1$RSA==VariableModel1],nullInterval = c(Perm.test.corr[1,], Perm.test.corr[nPermutations/100*95,]))
    
    # bf
    bfevidence<- bf[2] / bf[1]
    # View(bfevidence)
    
    template <- "A permutation Test was was conducted to compare {VariableOriginal} and  {VariableModel1}.
The relation between both RMS was: {test.cor}. The permutation test with percentil 95 was {lim}, the relationship had a p-value of {permPValue}.
the BF10 index was: {bfevidence} "
    
    
    template2<-glue::glue(
      template,
      
      VariableOriginal   = VariableOriginal,
      VariableModel     = VariableModel1,
      test.cor   =round( test.cor,2),
      lim     = round(lim,2),
      permPValue     = permPValue,
      bfevidence= round(bfevidence@bayesFactor[["bf"]],2)
    )
    
    
    if (PlotRSMReg==1){
      #-------------- Plot RSM Matrix
      mypalette<- colors()[c(73, 30, 28, 124, 635, 86,384,144, 148)]
      
      corr_prepPlot<- ModelData1_2Plot%>%
        dplyr::select(var1, var2, mean_all)
      min<- floor_dec(min(corr_prepPlot[,3]),5)
      max<-ceiling_dec( max(corr_prepPlot[,3][corr_prepPlot[,3] != max(corr_prepPlot[,3])]),5)
      cols <- rev(rainbow(20)[-20])
      
      plotRSM1<-ggplot(data=corr_prepPlot, aes(x = var1, y = var2, fill = mean_all)) + #value2
        geom_tile() +
        labs(x = paste(Category, "sorted by", VariableModel1), y = paste(Category, "sorted by", VariableModel1), fill = "", title = "Similarity Matrix") +
        coord_fixed() +
        theme_minimal()+
        scale_y_continuous(trans = "reverse")+
        scale_fill_gradientn(colours=rev(mypalette),na.value = colors() [c(73)], #heat.colors(12)
                             breaks=c(min, max),labels=c(min, max),
                             limits=c(min,max))+ theme(legend.position="bottom")+
        theme(text=element_text(family="Times", size=12))
    }
    
    # }  else if (whichModel[ii]== deparse(substitute(ModelData2))){ #does not work....
    
    #EachModelSeparately
    ModelData2_2Plot<-ModelData2
    ModelData2_2Plot$mean_all<-modelM2[["residuals"]]
    
    ModelData2$mean_all<-modelM2[["residuals"]]
    ModelData2<-ModelData2%>%dplyr::select(var1,var2,mean_all)%>%
      filter(ModelData2$var1>ModelData2$var2)%>%
      filter(var1 > 0)
    ModelData2$RSA<-VariableModel2
    
    
    CombinedOriginalModel2<- rbind(OriginalData, ModelData2)
    test.cor<- cor(CombinedOriginalModel2$mean_all[CombinedOriginalModel2$RSA==VariableOriginal],
                   CombinedOriginalModel2$mean_all[CombinedOriginalModel2$RSA==VariableModel2], method= "spearman")
    
    
    
    set.seed(seedPer)  
    
    Perm.test.corr <- rep(0, nPermutations) #matrix were to save the permutations
    n <- length(CombinedOriginalModel2$RSA)  
    # the variable we will resample from 
    #     (note, could use the labels(RSA) too, and "shuffle this")
    variable <- CombinedOriginalModel2$RSA
    # initialize a matrix to store the permutation of the name
    PermSamplesOther <- matrix(0, nrow=n, ncol=nPermutations)
    # now, get those permutation samples from the labels, using a loop
    for(i in 1:nPermutations){
      PermSamplesOther[,i] <- sample(variable, size= n, replace=FALSE)
    }
    
    #Loop to perform the nPermutations permutations
    for (i in 1:nPermutations){
      
      # calculate the perm-test-stat1 and save it
      Perm_test_corr_sort<- CombinedOriginalModel2%>%
        select ( - RSA)
      Perm<- PermSamplesOther[,i]
      Perm_test_corr_sort <- data.frame(Perm_test_corr_sort, Perm)
      Perm_test_corr_sort<-as_tibble(Perm_test_corr_sort)%>%
        arrange(var1)%>%
        arrange(var2)%>%
        arrange(Perm)
      
      Perm.test.corr[i] <- cor(Perm_test_corr_sort$mean_all[Perm_test_corr_sort$Perm==VariableOriginal],
                               Perm_test_corr_sort$mean_all[Perm_test_corr_sort$Perm==VariableModel2],method = c("spearman"))
    }
    
    
    
    
    Perm.test.corr<- as_tibble(Perm.test.corr)
    Perm.test.corr<- Perm.test.corr%>%
      arrange_at(1)
    
    
    lim<-Perm.test.corr[nPermutations/100*95,1]
    lim<-as.numeric(lim)
    plot.perm2<-
      ggplot(Perm.test.corr, aes(x=value)) +
      geom_histogram(
        bins= 20,
        fill=I("blue"), 
        col=I("black"), 
        alpha=I(.2))+
      geom_vline(xintercept=lim, linetype="dashed", 
                 color = "black", size=2 )+
      geom_vline(xintercept=test.cor,
                 color = "red", size=2 )+ theme(text=element_text(family="Times", size=12),
                                                axis.text=element_text(size=25),axis.line.y= element_line(),
                                                axis.line.x= element_line())+
      
      theme(legend.position="bottom")+
      theme(#panel.grid.minor = element_blank(),
        panel.background = element_blank())
    significanceTest<- which(Perm.test.corr>test.cor)
    
    permPValue<-(nPermutations-significanceTest[1])/nPermutations
    
    
    bf =correlationBF(CombinedOriginalModel2$mean_all[CombinedOriginalModel2$RSA==VariableOriginal],
                      CombinedOriginalModel2$mean_all[CombinedOriginalModel2$RSA==VariableModel2],nullInterval = c(Perm.test.corr[1,], Perm.test.corr[nPermutations/100*95,]))
    
    # bf
    bfevidence<- bf[2] / bf[1]
    # View(bfevidence)
    
    template <- "A permutation Test was was conducted to compare {VariableOriginal} and  {VariableModel2}.
The relation between both RMS was: {test.cor}. The permutation test with percentil 95 was {lim}, the relationship had a p-value of {permPValue}.
the BF10 index was: {bfevidence} "
    
    
    template3<-glue::glue(
      template,
      
      VariableOriginal   = VariableOriginal,
      VariableModel     = VariableModel2,
      test.cor   =round( test.cor,2),
      lim     = round(lim,2),
      permPValue     = permPValue,
      bfevidence= round(bfevidence@bayesFactor[["bf"]],2)
    )
    
    if (PlotRSMReg==1){
      #-------------- Plot RSM Matrix
      mypalette<- colors()[c(73, 30, 28, 124, 635, 86,384,144, 148)]
      
      corr_prepPlot<- ModelData2_2Plot%>%
        dplyr::select(var1, var2, mean_all)
      min<- floor_dec(min(corr_prepPlot[,3]),5)
      max<-ceiling_dec( max(corr_prepPlot[,3][corr_prepPlot[,3] != max(corr_prepPlot[,3])]),5) 
      cols <- rev(rainbow(20)[-20])
      
      plotRSM2<-ggplot(data=corr_prepPlot, aes(x = var1, y = var2, fill = mean_all)) + #value2
        geom_tile() +
        labs(x = paste(Category, "sorted by", VariableModel2), y = paste(Category, "sorted by", VariableModel2), fill = "", title = "Similarity Matrix") +
        coord_fixed() +    
        theme_minimal()+
        scale_y_continuous(trans = "reverse")+
        scale_fill_gradientn(colours=rev(mypalette),na.value = colors() [c(73)], #heat.colors(12)
                             breaks=c(min, max),labels=c(min, max),
                             limits=c(min,max))+ theme(legend.position="bottom")+ 
        theme(text=element_text(family="Times", size=12))
      
    }
    #}  if that does not work whichModel
  }
  if (PlotRSMReg==1){
    mylist<-list (plot.perm,plotRSM1, template2, plot.perm2, plotRSM2, template3)
  }  else{  mylist<-list (plot.perm, template2, plot.perm2, template3)
  }
  
  return(mylist)
}

