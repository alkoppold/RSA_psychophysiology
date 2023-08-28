
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
  min<- floor_dec(min(corr_prepPlot[,3]),2)
  max<-ceiling_dec( max(corr_prepPlot[,3][corr_prepPlot[,3] != max(corr_prepPlot[,3])]),2) 
  cols <- rev(rainbow(20)[-20])
  
  plotRSM<-ggplot(data=corr_prepPlot, aes(x = var1, y = var2, fill = mean_all)) + #value2
    geom_tile() +
    labs(x = paste("Trials sorted by", sortType), y = paste("Trials sorted by", sortType), fill = "", title = "Correlation Matrix") +
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


RSMSingleValue<-function(OriginalData,nVariables,CharContainingNumberFromVector,VariableType,sortType,MaxValue){
  
  
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
  
  corr_prepPlot<- corr_prep_all%>%
    dplyr::select(var1, var2, mean_all)
  min<- floor_dec(min(corr_prepPlot[,3]),2)
  max<-ceiling_dec( max(corr_prepPlot[,3][corr_prepPlot[,3] != max(corr_prepPlot[,3])]),2) 
  cols <- rev(rainbow(20)[-20])
  
  plotRSM<-ggplot(data=corr_prepPlot, aes(x = var1, y = var2, fill = mean_all)) + #value2
    geom_tile() +
    labs(x = paste("Trials sorted by", sortType), y = paste("Trials sorted by", sortType), fill = "", title = "Correlation Matrix") +
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
  min<- floor_dec(min(corr_prepPlot[,3]),2)
  max<-ceiling_dec( max(corr_prepPlot[,3][corr_prepPlot[,3] != max(corr_prepPlot[,3])]),2) 
  cols <- rev(rainbow(20)[-20])
  
  plotRSM<-ggplot(data=corr_prepPlot, aes(x = var1, y = var2, fill = mean_all)) + #value2
    geom_tile() +
    labs(x = paste("Trials sorted by", sortType), y = paste("Trials sorted by", sortType), fill = "", title = "Correlation Matrix") +
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
  min<- floor_dec(min(corr_prepPlot[,3]),2)
  max<-ceiling_dec( max(corr_prepPlot[,3][corr_prepPlot[,3] != max(corr_prepPlot[,3])]),2) 
  cols <- rev(rainbow(20)[-20])
  
  plotRSM<-ggplot(data=corr_prepPlot, aes(x = var1, y = var2, fill = mean_all)) + #value2
    geom_tile() +
    labs(x = paste("Trials sorted by", sortType), y = paste("Trials sorted by", sortType), fill = "", title = "Correlation Matrix") +
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






#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#----------------- FUNCTION CREATED TO run compare RSMs with each other for each individual separately. 

tableData<-SCR_SingleAlina_corr_prep_all
tableModel<-BehavioralSingleAlina_invAK_corr_prep_all
RSMIndividual<-function(tableModel,tableData,Model,Data){
  
  #Select the biggest table
  if (ncol(tableModel)>ncol(tableData)){
    tableMother<-tableData
    tableDaugther<-tableModel
  }else {tableMother<-tableModel
  tableDaugther<-tableData}
  tableCorrelation<-as_tibble(vector(mode='numeric',length=ncol(tableMother)))
  tableDaugther2<-tableDaugther%>%
    filter(tableDaugther$var1>tableDaugther$var2)
  tableMother2<-tableMother%>%
    filter(tableMother$var1>tableMother$var2)
  for(i in 3:(ncol(tableMother)-1)) { #starts in column three
    
    single<-tableMother2[,i]
    single[,2]<-tableDaugther2[,which( colnames(tableDaugther2)==names(tableMother2[i]))]
    corrspSingle<- rcorr(as.matrix(single), type="spearman")
    tableCorrelation[i-2,]<-corrspSingle$r[2]
  }
  tableCorrelation<-tableCorrelation%>%drop_na()
  corrMatrixVpName<-paste(Model,Data,"tableCorrelation",sep= "_")
  assign(corrMatrixVpName, tableCorrelation,envir = globalenv())
  
  tableCorrelationPlot<-tableCorrelation
  tableCorrelationPlot$x<-1
  
  onewayTTest<- t.test(tableCorrelation, mu = 0, alternative = "two.sided")
  
  plotDistribution <- ggplot(tableCorrelationPlot, aes(y = value,x = x)) +
    geom_violin(position = position_dodge(width = 0.9)) +
    geom_quasirandom(dodge.width = 0.9, varwidth = TRUE)+
    stat_summary(fun="mean",color="red", geom="point", size=5)+
    # annotate("text", x=.75, y=max(tableCorrelation$value)+mean(tableCorrelation$value), label= paste("p-value:",onewayTTest[["p.value"]], sep = " "),
    #          col="blue", size=5, parse=TRUE)+
    
    annotate("text", x=.75, y=max(tableCorrelation$value), label= paste("p-value:",onewayTTest[["p.value"]], sep = " "),
             col="blue", size=5, parse=TRUE)+
    
    annotate("text", x=.75, y=mean(tableCorrelation$value),#*3,
             label= paste("t-value:",onewayTTest[["statistic"]], sep = " "),
             col="blue", size=5, parse=TRUE)+
    annotate("text", x=.75, y=min(tableCorrelation$value),#*4, 
             label= paste("Mean:",mean(tableCorrelationPlot$value), sep = " "),
             col="blue", size=5, parse=TRUE)
  plotDistribution
}
