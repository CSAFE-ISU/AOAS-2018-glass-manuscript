#functions

#take the pairwise differences between panes
Data_take_diff<-function(log_glass_data, k){
  Data<-NULL
  D<-NULL
  
  for (s in 1:(nrow(log_glass_data)-1)){
    for (d in (s+1):nrow(log_glass_data)){
      Info<-log_glass_data[s,-c(1,2)]-log_glass_data[d,-c(1,2)]
      D<-data.frame(s,d,Info)
      Data<-rbind(Data,D)
    }
  }
  Data$type<-k
  return(Data)
}


# test ROC by the model on given test data
test_roc <- function(model, data) {
  
  pROC::roc(data$type,
            predict(model, data, type = "prob")[, "X0"])
  
}

# estimate TPR at desired FPR value (p) on given test data
get_TPR_fixed_FPR<-function(data, p){
  finaldata<-NULL
  
  for (i in 1:length(levels(data$Method))){
    
    D1<-subset(data, (p-0.005) > one_minus_Spec & Method==levels(data$Method)[i])
    D2<-subset(data, one_minus_Spec > (p+0.005) & Method==levels(data$Method)[i])
    
    x1<-D1[which.max(D1$one_minus_Spec),"one_minus_Spec"]
    x2<-D2[which.min(D2$one_minus_Spec),"one_minus_Spec"]
    
    y1<-D1[which.max(D1$one_minus_Spec),"Sen"]
    y2<-D2[which.min(D2$one_minus_Spec),"Sen"]
    
    lx<-x2-x1
    lx1<-p-x1
    lx2<-x2-p
    
    fixed_TPR<-y1*(lx1/lx)+y2*(lx2/lx)
    Method<-paste0(levels(data$Method)[i], ', p=', p) 
    dd<-data.frame(Method, fixed_TPR)
    finaldata<-rbind(finaldata,dd)
    finaldata<-finaldata[order(-finaldata$fixed_TPR),]
  }
  return(finaldata)
}


#prediction error table by random forest (FPR, FNR)
err.table<-function(table){
  A<-table[1,1]
  B<-table[1,2]
  C<-table[2,1]
  D<-table[2,2]
  
  one_minus_sen<-1-A/(A+C)
  one_minus_spec<-1-D/(B+D)
  overall.err<-round(1-sum(diag(table))/sum(table), 4)
  
  Rate<-c(one_minus_sen,one_minus_spec,overall.err)
  return(data.frame(Rate))
}

#prediction error table by BART (FPR, FNR)
err.table.bart<-function(table){
  A<-table[1,1]
  B<-table[2,1]
  C<-table[1,2]
  D<-table[2,2]
  
  one_minus_sen<-1-A/(A+C)
  one_minus_spec<-1-D/(B+D)
  overall.err<-c(table[3,3])
  
  Rate<-c(one_minus_sen,one_minus_spec,overall.err)
  return(data.frame(Rate))
}




# produce FPR, FNR, TPR, TNR by standard or mnsc methods
stad_err<-function(same, diff){
  
  a<-sum(same[,3]>4) # pred diff when same is true - type 2 error
  b<-sum(diff[,3]<4) # pred same when diff is true - type 1 error
  na<-nrow(same)
  nb<-nrow(diff)
  Sen<-(1-a/na)
  one_minus_spec<-b/nb
  one_minus_sen<-a/na
  Spec<-1-b/nb
  overall_err<-(a+b)/(na+nb)
  
  return(data.frame(Sen, one_minus_spec, Spec, one_minus_sen, overall_err))
  
}


# From the given pane, store the information of desired number of control fragments as list
# The information of control fragments is in S.
control_samples<-function(data, S, i, k){
  
  X_c<-list()  
  c<-NULL
  c<-S[k,2:length(S)]
  
  for ( j in 1:length(c)){
    X_c[[j]]<-c(data[data$Piece == levels(data$Piece)[c[,j]],i])[1:5]}
  
  return(X_c)
}



# store the differnece between Q (one recovered) and three K's (three controls) as data frame
# The difference between Q and three K's, when they are from the same pane
# Along the S, the number of controls can be changed 

rf_test_same_ncontrol<-function(originaldata,S){
  D<-NULL
  Data<-NULL
  diff<-NULL
  diff<-c(0,0)
  originaldata$Piece<-as.factor(originaldata$Piece)
  
  for (k in 1:nrow(S)){
    r<-S[k,1]
    c1<-S[k,2]
    c2<-S[k,3]
    c3<-S[k,4]
    
    for (i in 3:20){
      X_r1<-c(originaldata[originaldata$Piece==levels(originaldata$Piece)[r],i])[1:5]
      X_c<-control_samples(originaldata, S, i, k)
      
      log_X_r1<-log(X_r1)
      log_X_c<-lapply(X_c,log)
      
      mean_log_X_r<-mean(log_X_r1, na.rm=T)
      mean_log_X_c<-mean(unlist(log_X_c), na.rm=T)
      
      diff[i]<-mean_log_X_r-mean_log_X_c
      
    }
    
    Info<-c(diff[3:20])
    Data<-rbind(Data,diff[3:20])
    Data<-data.frame(Data)
    Data$type<-"X0"
    names(Data)<-c("Li7", "Na23","Mg25","Al27","K39","Ca42","Ti49","Mn55","Fe57","Rb85",
                   "Sr88", "Zr90","Ba137","La139","Ce140","Nd146","Hf178", "Pb208","type")
  }
  return(Data) 
}

# store the differnece between Q (one recovered) and three K's (three controls) as data frame
# The difference between Q and three K's, when they are from two different panes
# Along the S, the number of controls can be changed 


rf_test_diff_ncontrol<-function(data1, data2, nfrag1, S2){
  D<-NULL
  Data<-NULL
  diff<-c(0,0)
  
  data1$Piece<-as.factor(data1$Piece)
  data2$Piece<-as.factor(data2$Piece)
  
  for (r in 1:nfrag1){
    
    for (k in 1:30){
      c1<-S2[k,2]
      c2<-S2[k,3]
      c3<-S2[k,4]
      
      for (i in 3:20){
        
        X_r1<-c(data1[data1$Piece==levels(data1$Piece)[r],i])[1:5]
        X_c<-control_samples(data2, S2, i, k)
        
        log_X_r1<-log(X_r1)
        log_X_c<-lapply(X_c,log)
        
        mean_log_X_r<-mean(log_X_r1, na.rm=T)
        mean_log_X_c<-mean(unlist(log_X_c), na.rm=T)
        
        diff[i]<-mean_log_X_r-mean_log_X_c
        
      }
      
      Info<-c(diff[3:20])
      Data<-rbind(Data,diff[3:20])
    }
  }
  
  Data<-data.frame(Data)
  Data$type<-"X1"
  names(Data)<-c("Li7", "Na23","Mg25","Al27","K39","Ca42","Ti49","Mn55","Fe57","Rb85",
                 "Sr88", "Zr90","Ba137","La139","Ce140","Nd146","Hf178", "Pb208","type")
  return(Data)
}

# version 2 - 5 times repeating, instead of 30 times

rf_test_diff_ncontrol_ver2<-function(data1, data2, nfrag1, S2){
  D<-NULL
  Data<-NULL
  diff<-c(0,0)
  
  data1$Piece<-as.factor(data1$Piece)
  data2$Piece<-as.factor(data2$Piece)
  
  for (r in 1:nfrag1){
    
    for (k in 1:5){
      c1<-S2[k,2]
      c2<-S2[k,3]
      c3<-S2[k,4]
      
      for (i in 3:20){
        
        X_r1<-c(data1[data1$Piece==levels(data1$Piece)[r],i])[1:5]
        X_c<-control_samples(data2, S2, i, k)
        
        log_X_r1<-log(X_r1)
        log_X_c<-lapply(X_c,log)
        
        mean_log_X_r<-mean(log_X_r1, na.rm=T)
        mean_log_X_c<-mean(unlist(log_X_c), na.rm=T)
        
        diff[i]<-mean_log_X_r-mean_log_X_c
        
      }
      
      Info<-c(diff[3:20])
      Data<-rbind(Data,diff[3:20])
    }
  }
  
  Data<-data.frame(Data)
  Data$type<-"X1"
  names(Data)<-c("Li7", "Na23","Mg25","Al27","K39","Ca42","Ti49","Mn55","Fe57","Rb85",
                 "Sr88", "Zr90","Ba137","La139","Ce140","Nd146","Hf178", "Pb208","type")
  return(Data)
}

# From the classfication table by random forest, it produces a table of PPV and NPV

rf_ppv_npv<-function(table){
  
  ppv<-table[1,1]/(table[1,1]+table[2,1])
  npv<-table[2,2]/(table[1,2]+table[2,2])
  
  return(data.frame(ppv, npv))
}

# From the classfication table by BART, it produces a table of PPV and NPV
bart_ppv_npv<-function(table){
  
  ppv<-table[1,1]/(table[1,1]+table[1,2])
  npv<-table[2,2]/(table[2,1]+table[2,2])
  
  return(data.frame(ppv, npv))
}


# Generate the differnece of Q-K1, Q-K2, Q-K3 separately and then get the scores of three comparisons.
# Among three scores, get average, maximum and minimum of them
# Q and three K's are from the same pane

rf_test_same_3indep<-function(originaldata,S){
  D<-NULL
  Data<-NULL
  Data2<-NULL
  Data3<-NULL
  
  
  diff1<-NULL
  diff1<-c(0,0)
  diff2<-NULL
  diff2<-c(0,0)
  diff3<-NULL
  diff3<-c(0,0)
  
  originaldata$Piece<-as.factor(originaldata$Piece)
  
  for (k in 1:nrow(S)){
    r<-S[k,1]
    c1<-S[k,2]
    c2<-S[k,3]
    c3<-S[k,4]
    
    for (i in 3:20){
      X_r1<-c(originaldata[originaldata$Piece==levels(originaldata$Piece)[r],i])[1:5]
      
      X_c1<-c(originaldata[originaldata$Piece == levels(originaldata$Piece)[c1],i])[1:5]
      X_c2<-c(originaldata[originaldata$Piece == levels(originaldata$Piece)[c2],i])[1:5]
      X_c3<-c(originaldata[originaldata$Piece == levels(originaldata$Piece)[c3],i])[1:5]
      
      log_X_r1<-log(X_r1)
      log_X_c1<-log(X_c1)
      log_X_c2<-log(X_c2)
      log_X_c3<-log(X_c3)
      
      
      mean_log_X_r<-mean(log_X_r1, na.rm=T)
      mean_log_X_c1<-mean(log_X_c1, na.rm=T)
      mean_log_X_c2<-mean(log_X_c2, na.rm=T)
      mean_log_X_c3<-mean(log_X_c3, na.rm=T)
      
      diff1[i]<-mean_log_X_r-mean_log_X_c1
      diff2[i]<-mean_log_X_r-mean_log_X_c2
      diff3[i]<-mean_log_X_r-mean_log_X_c3
    }
    
    Data<-rbind(diff1[3:20],diff2[3:20],diff3[3:20])
    Data<-data.frame(Data)
    Data$type<-"X0"
    names(Data)<-c("Li7", "Na23","Mg25","Al27","K39","Ca42","Ti49","Mn55","Fe57","Rb85",
                   "Sr88", "Zr90","Ba137","La139","Ce140","Nd146","Hf178", "Pb208","type")
    
    RF_score<-predict(plan14_down_fit,  Data,type="prob")[,"X0"]
    mean_score<-mean(RF_score)
    max_score<-max(RF_score)
    min_score<-min(RF_score)
    
    Data2<-data.frame(k,mean_score,max_score,min_score)
    Data3<-rbind(Data3, Data2)
  }
  Data3$type<-"X0"
  return(Data3) 
}


# Generate the differnece of Q-K1, Q-K2, Q-K3 separately and then get the scores by RF of three comparisons.
# Among three scores, get average, maximum and minimum of them
# Q and three K's are from two different panes

rf_test_diff_3indep<-function(data1, data2, nfrag1, S2){
  D<-NULL
  Data<-NULL
  Data2<-NULL
  Data3<-NULL
  
  
  diff1<-NULL
  diff1<-c(0,0)
  diff2<-NULL
  diff2<-c(0,0)
  diff3<-NULL
  diff3<-c(0,0)
  
  data1$Piece<-as.factor(data1$Piece)
  data2$Piece<-as.factor(data2$Piece)
  
  for (r in 1:nfrag1){
    
    for (k in 1:30){
      c1<-S2[k,2]
      c2<-S2[k,3]
      c3<-S2[k,4]
      
      for (i in 3:20){
        X_r1<-c(data1[data1$Piece==levels(data1$Piece)[r],i])[1:5]
        
        X_c1<-c(data2[data2$Piece == levels(data2$Piece)[c1],i])[1:5]
        X_c2<-c(data2[data2$Piece == levels(data2$Piece)[c2],i])[1:5]
        X_c3<-c(data2[data2$Piece == levels(data2$Piece)[c3],i])[1:5]
        
        
        log_X_r1<-log(X_r1)
        log_X_c1<-log(X_c1)
        log_X_c2<-log(X_c2)
        log_X_c3<-log(X_c3)
        
        
        mean_log_X_r<-mean(log_X_r1, na.rm=T)
        mean_log_X_c1<-mean(log_X_c1, na.rm=T)
        mean_log_X_c2<-mean(log_X_c2, na.rm=T)
        mean_log_X_c3<-mean(log_X_c3, na.rm=T)
        
        diff1[i]<-mean_log_X_r-mean_log_X_c1
        diff2[i]<-mean_log_X_r-mean_log_X_c2
        diff3[i]<-mean_log_X_r-mean_log_X_c3
      }
      
      Data<-rbind(diff1[3:20],diff2[3:20],diff3[3:20])
      Data<-data.frame(Data)
      Data$type<-"X1"
      names(Data)<-c("Li7", "Na23","Mg25","Al27","K39","Ca42","Ti49","Mn55","Fe57","Rb85",
                     "Sr88", "Zr90","Ba137","La139","Ce140","Nd146","Hf178", "Pb208","type")
      
      RF_score<-predict(plan14_down_fit,  Data, type="prob")[,"X0"]
      mean_score<-mean(RF_score)
      max_score<-max(RF_score)
      min_score<-min(RF_score)
      
      Data2<-data.frame(k,mean_score,max_score,min_score)
      Data3<-rbind(Data3, Data2)
    }
  }
  Data3$type<-"X1"
  return(Data3) 
}


# Using threshold 0.5, get FPR, TPR, FNR, FPR of the prediction results by random forest
RF_err_col<-function(same, diff){
  
  a<-sum(same<0.5) # pred diff when same is true - type 2 error
  b<-sum(diff>0.5) # pred same when diff is true - type 1 error
  na<-length(same)
  nb<-length(diff)
  Sen<-(1-a/na)
  one_minus_spec<-b/nb
  one_minus_sen<-a/na
  Spec<-1-b/nb
  overall_err<-(a+b)/(na+nb)
  
  return(data.frame(Sen, one_minus_spec, Spec, one_minus_sen, overall_err))
  
}


# Generate the differnece of Q-K1, Q-K2, Q-K3 separately and then get the scores by BART of three comparisons.
# Among three scores, get average, maximum and minimum of them
# Q and three K's are from the same pane

bart_test_same_3indep<-function(originaldata,S){
  D<-NULL
  Data<-NULL
  Data2<-NULL
  Data3<-NULL
  Data_X<-NULL
  
  diff1<-NULL
  diff1<-c(0,0)
  diff2<-NULL
  diff2<-c(0,0)
  diff3<-NULL
  diff3<-c(0,0)
  
  originaldata$Piece<-as.factor(originaldata$Piece)
  
  for (k in 1:nrow(S)){
    r<-S[k,1]
    c1<-S[k,2]
    c2<-S[k,3]
    c3<-S[k,4]
    
    for (i in 3:20){
      X_r1<-c(originaldata[originaldata$Piece==levels(originaldata$Piece)[r],i])[1:5]
      
      X_c1<-c(originaldata[originaldata$Piece == levels(originaldata$Piece)[c1],i])[1:5]
      X_c2<-c(originaldata[originaldata$Piece == levels(originaldata$Piece)[c2],i])[1:5]
      X_c3<-c(originaldata[originaldata$Piece == levels(originaldata$Piece)[c3],i])[1:5]
      
      log_X_r1<-log(X_r1)
      log_X_c1<-log(X_c1)
      log_X_c2<-log(X_c2)
      log_X_c3<-log(X_c3)
      
      
      mean_log_X_r<-mean(log_X_r1, na.rm=T)
      mean_log_X_c1<-mean(log_X_c1, na.rm=T)
      mean_log_X_c2<-mean(log_X_c2, na.rm=T)
      mean_log_X_c3<-mean(log_X_c3, na.rm=T)
      
      diff1[i]<-mean_log_X_r-mean_log_X_c1
      diff2[i]<-mean_log_X_r-mean_log_X_c2
      diff3[i]<-mean_log_X_r-mean_log_X_c3
    }
    
    Data<-rbind(diff1[3:20],diff2[3:20],diff3[3:20])
    Data<-data.frame(Data)
    Data$type<-"X0"
    names(Data)<-c("Li7", "Na23","Mg25","Al27","K39","Ca42","Ti49","Mn55","Fe57","Rb85",
                   "Sr88", "Zr90","Ba137","La139","Ce140","Nd146","Hf178", "Pb208","type")
    
    Data_X<-Data[,1:18]
    
    bart_score<-predict(fit_plan14_down_bart, Data_X, type="prob")  
    mean_score<-mean(bart_score)
    max_score<-max(bart_score)
    min_score<-min(bart_score)
    
    Data2<-data.frame(k,mean_score,max_score,min_score)
    Data3<-rbind(Data3, Data2)
  }
  Data3$type<-"X0"
  return(Data3) 
}


# Generate the differnece of Q-K1, Q-K2, Q-K3 separately and then get the scores by BART of three comparisons.
# Among three scores, get average, maximum and minimum of them
# Q and three K's are from two different panes

bart_test_diff_3indep<-function(data1, data2, nfrag1, S2){
  D<-NULL
  Data<-NULL
  Data2<-NULL
  Data3<-NULL
  
  
  diff1<-NULL
  diff1<-c(0,0)
  diff2<-NULL
  diff2<-c(0,0)
  diff3<-NULL
  diff3<-c(0,0)
  
  data1$Piece<-as.factor(data1$Piece)
  data2$Piece<-as.factor(data2$Piece)
  
  for (r in 1:nfrag1){
    
    for (k in 1:30){
      c1<-S2[k,2]
      c2<-S2[k,3]
      c3<-S2[k,4]
      
      for (i in 3:20){
        X_r1<-c(data1[data1$Piece==levels(data1$Piece)[r],i])[1:5]
        
        X_c1<-c(data2[data2$Piece == levels(data2$Piece)[c1],i])[1:5]
        X_c2<-c(data2[data2$Piece == levels(data2$Piece)[c2],i])[1:5]
        X_c3<-c(data2[data2$Piece == levels(data2$Piece)[c3],i])[1:5]
        
        
        log_X_r1<-log(X_r1)
        log_X_c1<-log(X_c1)
        log_X_c2<-log(X_c2)
        log_X_c3<-log(X_c3)
        
        
        mean_log_X_r<-mean(log_X_r1, na.rm=T)
        mean_log_X_c1<-mean(log_X_c1, na.rm=T)
        mean_log_X_c2<-mean(log_X_c2, na.rm=T)
        mean_log_X_c3<-mean(log_X_c3, na.rm=T)
        
        diff1[i]<-mean_log_X_r-mean_log_X_c1
        diff2[i]<-mean_log_X_r-mean_log_X_c2
        diff3[i]<-mean_log_X_r-mean_log_X_c3
      }
      
      Data<-rbind(diff1[3:20],diff2[3:20],diff3[3:20])
      Data<-data.frame(Data)
      Data$type<-"X1"
      names(Data)<-c("Li7", "Na23","Mg25","Al27","K39","Ca42","Ti49","Mn55","Fe57","Rb85",
                     "Sr88", "Zr90","Ba137","La139","Ce140","Nd146","Hf178", "Pb208","type")
      
      Data_X<-Data[,1:18]
      
      bart_score<-predict(fit_plan14_down_bart, Data_X, type="prob")  
      mean_score<-mean(bart_score)
      max_score<-max(bart_score)
      min_score<-min(bart_score)
      
      Data2<-data.frame(k,mean_score,max_score,min_score)
      Data3<-rbind(Data3, Data2)
    }
  }
  Data3$type<-"X1"
  return(Data3) 
}


# Generate score from standard E2330 from the mean of Q and mean fo three K's, when they are from the same pane

stad_test_same_ncontrol<-function(originaldata,S){
  D<-NULL
  Data<-NULL
  crit<-NULL
  crit<-c(0,0)
  originaldata$Piece<-as.factor(originaldata$Piece)
  
  for (k in 1:nrow(S)){
    r<-S[k,1]
    c1<-S[k,2]
    c2<-S[k,3]
    c3<-S[k,4]
    
    for (i in 3:20){
      X_r1<-c(originaldata[originaldata$Piece==levels(originaldata$Piece)[r],i])[1:5]
      X_c<-control_samples(originaldata, S, i, k)
      
      mean_X_r<-mean(X_r1,na.rm=T)
      mean_X_c<-mean(unlist(X_c), na.rm=T)
      
      sc1<-sqrt(var(unlist(X_c), na.rm=T))
      #sc1<-sqrt(var(originaldata[originaldata$Piece==levels(originaldata$Piece)[d],i]))
      sc2<-0.03*mean_X_c
      sc<-max(sc1,sc2)
      crit[i]<-abs((mean_X_r-mean_X_c)/sc)
    }
    score_stad<-max(crit)
    D<-data.frame(k,k,score_stad)
    Data<-rbind(Data,D)
  }
  return(Data) 
}

# Generate score from standard E2330 from the mean of Q and mean fo three K's, when they are from two different panes

stad_test_diff_ncontrol<-function(data1, data2, nfrag1, S2){
  D<-NULL
  Data<-NULL
  crit<-NULL
  crit<-c(0,0)
  data1$Piece<-as.factor(data1$Piece)
  data2$Piece<-as.factor(data2$Piece)
  
  for (r in 1:nfrag1){
    
    for (k in 1:30){
      c1<-S2[k,2]
      c2<-S2[k,3]
      c3<-S2[k,4]
      
      for (i in 3:20){
        
        X_r1<-c(data1[data1$Piece==levels(data1$Piece)[r],i])[1:5]
        X_c<-control_samples(data2, S2, i, k)
        
        mean_X_r<-mean(X_r1,na.rm=T)
        mean_X_c<-mean(unlist(X_c), na.rm=T)
        
        sc1<-sqrt(var(unlist(X_c), na.rm=T))
        #sc1<-sqrt(var(originaldata[originaldata$Piece==levels(originaldata$Piece)[d],i]))
        sc2<-0.03*mean_X_c
        sc<-max(sc1,sc2)
        crit[i]<-abs((mean_X_r-mean_X_c)/sc)
      }
      score_stad<-max(crit)
      D<-data.frame(r,k,score_stad)
      Data<-rbind(Data,D)
    }}
  return(Data)
}



stad_test_diff_ncontrol_ver2<-function(data1, data2, nfrag1, S2){
  D<-NULL
  Data<-NULL
  crit<-NULL
  crit<-c(0,0)
  data1$Piece<-as.factor(data1$Piece)
  data2$Piece<-as.factor(data2$Piece)
  
  
  for (r in 1:nfrag1){
    
    for (k in 1:5){
      c1<-S2[k,2]
      c2<-S2[k,3]
      c3<-S2[k,4]
      
      for (i in 3:20){
        
        X_r1<-c(data1[data1$Piece==levels(data1$Piece)[r],i])[1:5]
        X_c<-control_samples(data2, S2, i, k)
        
        mean_X_r<-mean(X_r1,na.rm=T)
        mean_X_c<-mean(unlist(X_c), na.rm=T)
        
        sc1<-sqrt(var(unlist(X_c), na.rm=T))
        #sc1<-sqrt(var(originaldata[originaldata$Piece==levels(originaldata$Piece)[d],i]))
        sc2<-0.03*mean_X_c
        sc<-max(sc1,sc2)
        crit[i]<-abs((mean_X_r-mean_X_c)/sc)
      }
      score_stad<-max(crit)
      D<-data.frame(r,k,score_stad)
      Data<-rbind(Data,D)
    }
    
  }
  return(Data)
}


# FPR, FNR, TPR, TNR by the score of standard - E2330
stad_err<-function(same, diff){
  
  a<-sum(same[,3]>4) # pred diff when same is true - type 2 error
  b<-sum(diff[,3]<4) # pred same when diff is true - type 1 error
  na<-nrow(same)
  nb<-nrow(diff)
  Sen<-(1-a/na)
  one_minus_spec<-b/nb
  one_minus_sen<-a/na
  Spec<-1-b/nb
  overall_err<-(a+b)/(na+nb)
  
  return(data.frame(Sen, one_minus_spec, Spec, one_minus_sen, overall_err))
  
}

# PPV and NPV by the score of standard - E2330
stad_ppv_npv<-function(same, diff){
  
  a<-sum(same[,3]>4) # pred diff when same is true - type 2 error
  b<-sum(diff[,3]<4) # pred same when diff is true - type 1 error
  na<-nrow(same)
  nb<-nrow(diff)
  
  ppv<-(na-a)/na
  npv<-(nb-b)/nb
  
  return(data.frame(ppv, npv))
  
}

# Generate score from modified n sigma criterion (Weis, 2011) from the mean of Q and mean fo three K's, when they are from the same pane

mnsc_test_same_ncontrol<-function(originaldata,S){
  D<-NULL
  Data<-NULL
  crit<-NULL
  crit<-c(0,0)
  originaldata$Piece<-as.factor(originaldata$Piece)
  
  for (k in 1:nrow(S)){
    r<-S[k,1]
    c1<-S[k,2]
    c2<-S[k,3]
    c3<-S[k,4]
    
    for (i in 3:20){
      X_r1<-c(originaldata[originaldata$Piece==levels(originaldata$Piece)[r],i])[1:5]
      X_c<-control_samples(originaldata,S, i, k)
      
      mean_X_r<-mean(X_r1,na.rm=T)
      mean_X_c<-mean(unlist(X_c), na.rm=T)
      
      num<-exp(abs(log(mean_X_r)-log(mean_X_c)))-1
      deno<-FRSD[i]
      
      crit[i]<-num/deno
    }
    score_mnsc<-max(crit)
    D<-data.frame(k,k,score_mnsc)
    Data<-rbind(Data,D)
  }
  return(Data) 
}


# Generate score from modified n sigma criterion (Weis, 2011) from the mean of Q and mean fo three K's, when they are from two different panes

mnsc_test_diff_ncontrol<-function(data1, data2, nfrag1, S2){
  D<-NULL
  Data<-NULL
  crit<-NULL
  crit<-c(0,0)
  data1$Piece<-as.factor(data1$Piece)
  data2$Piece<-as.factor(data2$Piece)
  
  
  for (r in 1:nfrag1){
    
    for (k in 1:30){
      c1<-S2[k,2]
      c2<-S2[k,3]
      c3<-S2[k,4]
      
      for (i in 3:20){
        
        X_r1<-c(data1[data1$Piece==levels(data1$Piece)[r],i])[1:5]
        X_c<-control_samples(data2,S2, i, k)
        
        mean_X_r<-mean(X_r1,na.rm=T)
        mean_X_c<-mean(unlist(X_c), na.rm=T)
        
        num<-exp(abs(log(mean_X_r)-log(mean_X_c)))-1
        deno<-FRSD[i]
        
        crit[i]<-num/deno
      }
      score_mnsc<-max(crit)
      D<-data.frame(k,k,score_mnsc)
      Data<-rbind(Data,D)
    }
  }
  return(Data)
}


# Generate score as hotelling T^2 statistic using shrinkage covariance estimation (Campbell and Curran, 2009)
# from the mean of Q and mean fo three K's, when they are from the same pane


hot_test_same_ncontrol<-function(originaldata,S){
  D<-NULL
  Data<-NULL
  originaldata$Piece<-as.factor(originaldata$Piece)
  
  for (k in 1:nrow(S)){
    r<-S[k,1]
    
    X_r1<-as.matrix(originaldata[originaldata$Piece==levels(originaldata$Piece)[r],c(3:20)])
    X_c<-control_matrix(originaldata,S, k)
    
    if(nrow(X_r1)>5) X_r1<-X_r1[1:5,]
    
    for ( h in 1:length(X_c)){
      if(nrow(X_c[[h]])>5) X_c[[h]]<-X_c[[h]][1:5,]
    }
    
    combine_X_c<-do.call(rbind, X_c)
    
    hot<-hotelling.test(X_r1, combine_X_c, shrinkage = TRUE, perm = TRUE, B=100)
    hot_stat<-hot$stats$statistic
    hot_pval<-hot$pval
    hot_pval_res<-ifelse(hot_pval<0.05, "D", "S")
    D<-data.frame(k,k,hot_stat,hot_pval_res)
    Data<-rbind(Data,D)
  }
  return(Data) 
}


# Generate score as hotelling T^2 statistic using shrinkage covariance estimation (Campbell and Curran, 2009)
# from the mean of Q and mean fo three K's, when they are from two different panes


hot_test_diff_ncontrol<-function(data1, data2, nfrag1, S2){
  D<-NULL
  Data<-NULL
  crit<-NULL
  crit<-c(0,0)
  data1$Piece<-as.factor(data1$Piece)
  data2$Piece<-as.factor(data2$Piece)
  
  for (r in 1:nfrag1){
    
    for (k in 1:30){
      X_r1<-as.matrix(data1[data1$Piece==levels(data1$Piece)[r],c(3:20)])
      X_c<-control_matrix(data2,S2, k)
      
      if(nrow(X_r1)>5) X_r1<-X_r1[1:5,]
      
      for ( h in 1:length(X_c)){
        if(nrow(X_c[[h]])>5) X_c[[h]]<-X_c[[h]][1:5,]
      }
      
      combine_X_c<-do.call(rbind, X_c)
      
      hot<-hotelling.test(X_r1, combine_X_c, shrinkage = TRUE, perm = TRUE, B=100)
      hot_stat<-hot$stats$statistic
      hot_pval<-hot$pval
      hot_pval_res<-ifelse(hot_pval<0.05, "D", "S")
      D<-data.frame(k,k,hot_stat,hot_pval_res)
      Data<-rbind(Data,D)
    }
  }
  return(Data)
}


# FPR, FNR, TPR, TNR by the prediction results from hotelling T^2 statistic
hot_err<-function(same, diff){
  
  a<-sum(same[,4]=="D") # pred diff when same is true - type 2 error
  b<-sum(diff[,4]=="S") # pred same when diff is true - type 1 error
  na<-nrow(same)
  nb<-nrow(diff)
  Sen<-(1-(a/na))
  one_minus_spec<-b/nb
  one_minus_sen<-a/na
  Spec<-(1-(b/nb))
  overall_err<-(a+b)/(na+nb)
  
  return(data.frame(Sen, one_minus_spec, Spec, one_minus_sen, overall_err))
  
}

# PPV and NPV by the prediction results from hotelling T^2 statistic
hot_ppv_npv<-function(same, diff){
  
  a<-sum(same[,4]=="D") # pred diff when same is true - type 2 error
  b<-sum(diff[,4]=="S") # pred same when diff is true - type 1 error
  na<-nrow(same)
  nb<-nrow(diff)
  
  ppv<-(na-a)/na
  npv<-(nb-b)/nb
  
  return(data.frame(ppv, npv))
  
}




# Generate score as Optimum test statistic (Parker and Holford, 1968)
# from the mean of Q and mean fo three K's, when they are from the same pane

optim_test_same_ncontrol<-function(originaldata,S){
  D<-NULL
  Data<-NULL
  originaldata$Piece<-as.factor(originaldata$Piece)
  d2i<-c(0,0)
  
  for (k in 1:nrow(S)){
    r<-S[k,1]
    c1<-S[k,2]
    c2<-S[k,3]
    c3<-S[k,4]
    
    for (i in 3:20){
      X_r1<-c(originaldata[originaldata$Piece==levels(originaldata$Piece)[r],i])
      X_c<-control_samples(originaldata,S, i, k)
      
      if(length(X_r1)>5) X_r1<-X_r1[1:5]
      for ( h in 1:length(X_c)){
        if(length(X_c[[h]])>5) X_c[[h]]<-X_c[[h]][1:5]
      }
      
      log_X_c<-lapply(X_c,log)
      
      log_X_r1<-log(X_r1)
      
      ui<-var(log_X_r1,na.rm=T)
      vi<-var(c(unlist(log_X_c)),na.rm=T)
      
      nr<-length(log_X_r1)
      nc<-length(c(unlist(log_X_c)))
      
      mean_log_X_r<-mean(log_X_r1, na.rm=T)
      mean_log_X_c<-mean(unlist(log_X_c), na.rm=T)
      
      li<-sqrt((ui*(nr-1)+vi*(nc-1))/(nr+nc-2)*(1/nr+1/nc))
      
      d2i[i]<-((mean_log_X_r-mean_log_X_c)/li)^2
    }
    
    nu<-(nr+nc-2)
    H=sum(nu*log(1+d2i/nu))
    Final<-ifelse(H>36.81, "D","S")
    D<-data.frame(k,k, H, Final)
    Data<-rbind(Data,D)
  }
  return(Data) 
}


# Generate score as Optimum test statistic (Parker and Holford, 1968)
# from the mean of Q and mean fo three K's, when they are from two different panes

optim_test_diff_ncontrol<-function(data1, data2, nfrag1, S2){
  D<-NULL
  Data<-NULL
  d2i<-c(0,0)
  
  data1$Piece<-as.factor(data1$Piece)
  data2$Piece<-as.factor(data2$Piece)
  
  for (r in 1:nfrag1){
    for (k in 1:30){
      c1<-S2[k,2]
      c2<-S2[k,3]
      c3<-S2[k,4]
      
      for (i in 3:20){
        X_r1<-c(data1[data1$Piece==levels(data1$Piece)[r],i])
        X_c<-control_samples(data2,S2, i, k)
        
        if(length(X_r1)>5) X_r1<-X_r1[1:5]
        for ( h in 1:length(X_c)){
          if(length(X_c[[h]])>5) X_c[[h]]<-X_c[[h]][1:5]
        }
        
        log_X_c<-lapply(X_c,log)
        
        log_X_r1<-log(X_r1)
        
        
        ui<-var(log_X_r1,na.rm=T)
        vi<-var(c(unlist(log_X_c)),na.rm=T)
        
        
        nr<-length(log_X_r1)
        nc<-length(c(unlist(log_X_c)))
        
        mean_log_X_r<-mean(log_X_r1, na.rm=T)
        mean_log_X_c<-mean(unlist(log_X_c), na.rm=T)
        
        li<-sqrt((ui*(nr-1)+vi*(nc-1))/(nr+nc-2)*(1/nr+1/nc))
        
        d2i[i]<-((mean_log_X_r-mean_log_X_c)/li)^2
      }
      
      nu<-(nr+nc-2)
      H=sum(nu*log(1+d2i/nu))
      Final<-ifelse(H>36.81, "D","S")
      D<-data.frame(k,k, H, Final)
      Data<-rbind(Data,D)
    }}
  return(Data) 
}

# FPR, FNR, TPR, TNR by the prediction results from Optimum test statistic
optim_err<-function(same, diff){
  
  a<-sum(same[,3]>30.03) # pred diff when same is true - type 2 error
  b<-sum(diff[,3]<30.03) # pred same when diff is true - type 1 error
  na<-nrow(same)
  nb<-nrow(diff)
  Sen<-(1-a/na)
  one_minus_spec<-b/nb
  one_minus_sen<-a/na
  Spec<-1-b/nb
  overall_err<-(a+b)/(na+nb)
  
  return(data.frame(Sen, one_minus_spec, Spec, one_minus_sen, overall_err))
  
}

# PPV and NPV by the prediction results from Optimum test statistic
optim_ppv_npv<-function(same, diff){
  
  a<-sum(same[,3]>30.03) # pred diff when same is true - type 2 error
  b<-sum(diff[,3]<30.03) # pred same when diff is true - type 1 error
  na<-nrow(same)
  nb<-nrow(diff)
  ppv<-(na-a)/na
  npv<-(nb-b)/nb
  
  return(data.frame(ppv, npv))
  
}

# Produce TPR, TNR, FPR, FNR with desired threshold
ROC.data<-function(Data, Test.set, thres, Method){
  
  nsame<-sum(Test.set$type=="X0") 
  ntotal<-nrow(Test.set)
  
  Sen<-NULL
  Spec<-NULL
  one_minus_Spec<-NULL
  one_minus_Sen<-NULL
  
  
  for (i in 1:length(thres)){
    A<-sum(Data$X0[1:nsame]>=thres[i])   # true same - pred same : A
    C<-sum(Data$X0[1:nsame]<thres[i])   # true same - pred diff : C 
    B<-sum(Data$X0[(nsame+1):ntotal]>=thres[i]) # true diff - pred same : B
    D<-sum(Data$X0[(nsame+1):ntotal]<thres[i]) # true diff - pred diff : D
    
    Sen[i]<-A/(A+C)
    Spec[i]<-D/(B+D)
    one_minus_Spec[i]<-(1-Spec[i])
    one_minus_Sen[i]<-(1-Sen[i])
    
  }
  
  D=data.frame(Sen=Sen, one_minus_Spec=one_minus_Spec, Spec=Spec, one_minus_Sen, Method=Method)
  return(D)
}

# Generate TPR, TNR, FPR, FNR by a set of threshold when prediction methods are RF, BART
ROC.data_col<-function(pred_prob, Test.set, thres, Method){
  
  nsame<-15300
  ndiff<-150060
  ntotal<-165360
  
  Sen<-NULL
  Spec<-NULL
  one_minus_Spec<-NULL
  one_minus_Sen<-NULL
  
  
  for (i in 1:length(thres)){
    A<-sum(pred_prob[1:nsame]>=thres[i])   # true same - pred same : A
    C<-sum(pred_prob[1:nsame]<thres[i])   # true same - pred diff : C 
    B<-sum(pred_prob[(nsame+1):ntotal]>=thres[i]) # true diff - pred same : B
    D<-sum(pred_prob[(nsame+1):ntotal]<thres[i]) # true diff - pred diff : D
    
    Sen[i]<-A/(A+C)
    Spec[i]<-D/(B+D)
    one_minus_Spec[i]<-(1-Spec[i])
    one_minus_Sen[i]<-(1-Sen[i])
    
  }
  
  D=data.frame(thres, Sen=Sen, one_minus_Spec=one_minus_Spec, Spec=Spec, one_minus_Sen, Method=Method)
  return(D)
}


# Generate TPR, TNR, FPR, FNR by a set of threshold when prediction methods are standard, modifired s sigma criteria, 
# Optimum test statistic and Hotelling T^2 statistic
ROC.stad<-function(Data, nsame, ntotal, thres, Method){
  Sen<-NULL
  Spec<-NULL
  one_minus_Spec<-NULL
  one_minus_Sen<-NULL
  
  for (i in 1:length(thres)){
    A<-sum(Data[1:nsame]<=thres[i])  # true same - pred same : A
    C<-sum(Data[1:nsame]>thres[i])   # true same - pred diff : C
    B<-sum(Data[(nsame+1):ntotal]<thres[i]) # true diff - pred same : B
    D<-sum(Data[(nsame+1):ntotal]>=thres[i]) # true diff - pred diff : D
    
    Sen[i]<-A/(A+C)
    Spec[i]<-D/(B+D)
    one_minus_Spec[i]<-(1-Spec[i])
    one_minus_Sen[i]<-(1-Sen[i])
  }
  
  D=data.frame(thres, Sen=Sen, one_minus_Spec=one_minus_Spec, Spec=Spec, one_minus_Sen, Method=Method)
  return(D)
}



# Using the information of TPR, TNR, FPR, FNR by several methods (in the list), 
# 10000 thresholds for BART and RF are generated from 0 to 1
# 10000 thresholds for Standard, Modified s sigma criteria, Hotelling T^2 statistic, optimum test statistic are generated
# frm minimum and maximum of predicted score by each method

ROC_data_list2<-function(ROC_list){
  
  thres<-seq(1,0, length.out=10000)
  thres_stad<-seq(min(ROC_list[[5]]), max(ROC_list[[5]]), length.out=10000)
  thres_new_mnsc<-seq(min(ROC_list[[6]]), max(ROC_list[[6]]), length.out = 10000)
  thres_hot<-seq(min(ROC_list[[7]]), max(ROC_list[[7]]), length.out = 10000)
  thres_optim<-seq(min(ROC_list[[8]]), max(ROC_list[[8]]), length.out = 10000)
  
  ROC_data<-list()
  
  ROC_data[[1]]<-ROC.data_col(ROC_list[[1]],plan16_test.set,thres, "Down-RF-Mean")
  ROC_data[[2]]<-ROC.data_col(ROC_list[[2]],plan16_test.set,thres, "Down-BART-Mean")
  ROC_data[[3]]<-ROC.data_col(ROC_list[[3]],plan16_test.set,thres, "Down-RF-Min")
  ROC_data[[4]]<-ROC.data_col(ROC_list[[4]],plan16_test.set,thres, "Down-BART-Min")
  
  nsame<-15300
  ndiff<-150060
  ntotal<-165360
  
  ROC_data[[5]]<-ROC.stad(ROC_list[[5]], nsame, ntotal, thres_stad, "Standard-ASTM")
  ROC_data[[6]]<-ROC.stad(ROC_list[[6]], nsame, ntotal, thres_new_mnsc, "Modified-s-sigma")
  ROC_data[[7]]<-ROC.stad(ROC_list[[7]], nsame, ntotal, thres_hot, "Hot.shirinkage")
  ROC_data[[8]]<-ROC.stad(ROC_list[[8]], nsame, ntotal, thres_optim, "Opt.test")
  
  
  return(ROC_data)
  
} 

# Calculate AUC with TPR and FPR values
simple_auc<-function(TPR, FPR){
  # inputs already sorted, best scores first 
  dFPR <- c(diff(FPR), 0)
  dTPR <- c(diff(TPR), 0)
  sum(TPR * dFPR) + sum(dTPR * dFPR)/2
}


# Calculate EER of the desired ROC curve 
eer.ROC<-function(ROCdata){
  
  ROC<-ROCdata
  missrates <-ROC$one_minus_Sen 
  farates <- ROC$one_minus_Spec 
  
  # Find the point on the ROC with miss slightly >= fa, and the point
  # next to it with miss slightly < fa.
  
  dists <- missrates - farates;
  idx1 <- which( dists == min( dists[ dists >= 0 ] ) );
  idx2 <- which( dists == max( dists[ dists < 0 ] ) );
  stopifnot( length( idx1 ) == 1 );
  stopifnot( length( idx2 ) == 1 );
  stopifnot( abs( idx1 - idx2 ) == 1 );
  
  
  
  # Extract the two points as (x) and (y), and find the point on the
  # line between x and y where the first and second elements of the
  # vector are equal.  Specifically, the line through x and y is:
  #   x + a*(y-x) for all a, and we want a such that
  #   x[1] + a*(y[1]-x[1]) = x[2] + a*(y[2]-x[2]) so
  #   a = (x[1] - x[2]) / (y[2]-x[2]-y[1]+x[1])
  
  x <- c( missrates[idx1], farates[idx1] );
  y <- c( missrates[idx2], farates[idx2] );
  a <- ( x[1] - x[2] ) / ( y[2] - x[2] - y[1] + x[1] );
  eer <- x[1] + a * ( y[1] - x[1] );
  
  return(eer)}


# AUC, EER are caculated and the optimum thresihold are found in ROC curves in ROC_list.
auc_eer_opt.t_generation<-function(ROC_list){
  
  auc<-list()
  eer<-list()
  opt.t<-list()
  one_minus_Spec<-list()
  one_minus_Sen<-list()
  
  
  for (i in 1:length(ROC_list)){
    auc[[i]]<-with(ROC_list[[i]], simple_auc(Sen, one_minus_Spec))
    eer[[i]]<-c(eer.ROC(unique(ROC_list[[i]][, -1])))
    opt.t[[i]]<-ROC_list[[i]]$thres[which.max(ROC_list[[i]]$Sen + ROC_list[[i]]$Spec)]
    
    one_minus_Spec[[i]]<-ROC_list[[i]][which(ROC_list[[i]]$thres == opt.t[[i]]),]$one_minus_Spec
    one_minus_Sen[[i]]<-ROC_list[[i]][which(ROC_list[[i]]$thres == opt.t[[i]]),]$one_minus_Sen
  }
  
  auc_combined<-c(do.call(rbind, auc))
  eer_combined<-c(do.call(rbind, eer))
  opt.t_combined<-c(do.call(rbind, opt.t))
  one_minus_Spec_c<-c(do.call(rbind, one_minus_Spec))
  one_minus_Sen_c<-c(do.call(rbind, one_minus_Sen))
  
  name<-c('Down-RF','Down-BART','Down-RF-Min','Down-BART-Min',"Standard-ASTM", "Modified-n-sigma", "Hot.shirinkage", "Opt.test")
  auc_eer_data_summary<-data.frame(name, auc_combined,eer_combined,opt.t_combined,one_minus_Spec_c,one_minus_Sen_c)
  auc_eer_data_summary2<-auc_eer_data_summary[order(auc_eer_data_summary$auc_combined, decreasing = TRUE), ]
  
  return(auc_eer_data_summary2)
  
}


# ROC_list is row-bined and then calculate TPR at fixed FPR (P)

get_TPR_fixed_FPR_combined<-function(data, P){
  finaldata<-NULL
  fixed_TPR<-NULL
  
  for (i in 1:length(levels(data$Method))){
    for ( j in 1:length(P)){
      p<-P[j]
      D1<-subset(data, (p-0.01) > one_minus_Spec & Method==levels(data$Method)[i])
      D2<-subset(data, one_minus_Spec > (p+0.01) & Method==levels(data$Method)[i])
      
      x1<-D1[which.max(D1$one_minus_Spec),"one_minus_Spec"]
      x2<-D2[which.min(D2$one_minus_Spec),"one_minus_Spec"]
      
      y1<-D1[which.max(D1$one_minus_Spec),"Sen"]
      y2<-D2[which.min(D2$one_minus_Spec),"Sen"]
      
      lx<-x2-x1
      lx1<-p-x1
      lx2<-x2-p
      ly<-y2-y1
      
      fixed_TPR[j]<-y1+ly*(lx1/lx)
      
    }
    
    Method<-paste0(levels(data$Method)[i]) 
    dd<-data.frame(Method, fixed_TPR[1], fixed_TPR[2], fixed_TPR[3])
    names(dd)[c(2,3,4)]<-c(paste("fixed TPR, p=",P[1]),paste("fixed TPR, p=",P[2]),paste("fixed TPR, p=",P[3]) )
    finaldata<-rbind(finaldata,dd)
  }
  
  finaldata<-finaldata[order(-finaldata[,2]),]
  return(finaldata)
}


# get the score of Standard - E2330 when one Q and one K are compared when Q and K are from the same pane

stad_test_same<-function(originaldata, nfrag){
  D<-NULL
  Data<-NULL
  crit<-NULL
  crit<-c(0,0)
  originaldata$Piece<-as.factor(originaldata$Piece)
  for (s in 1:(nfrag-1)){
    for (d in (s+1):nfrag){
      for (i in 3:20){
        X_r<-mean(originaldata[originaldata$Piece==levels(originaldata$Piece)[s],i])
        X_c<-mean(originaldata[originaldata$Piece==levels(originaldata$Piece)[d],i])
        sc1<-sqrt(var(originaldata[originaldata$Piece==levels(originaldata$Piece)[d],i]))
        sc2<-0.03*X_c
        sc<-max(sc1,sc2)
        crit[i]<-abs((X_r-X_c)/sc)
      }
      score_stad<-max(crit)
      D<-data.frame(s,d,score_stad)
      Data<-rbind(Data,D)
    }}
  return(Data) 
}

# get the score of Standard - E2330 when one Q and one K are compared when Q and K are from two different panes

stad_test_diff<-function(data1, data2, frag1, frag2){
  D<-NULL
  Data<-NULL
  crit<-NULL
  crit<-c(0,0)
  data1$Piece<-as.factor(data1$Piece)
  data2$Piece<-as.factor(data2$Piece)
  for (s in 1:frag1){
    for (d in 1:frag2){
      for (i in 3:20){
        X_r<-mean(data1[data1$Piece==levels(data1$Piece)[s],i])
        X_c<-mean(data2[data2$Piece==levels(data2$Piece)[d],i])
        sc1<-sqrt(var(data2[data2$Piece==levels(data2$Piece)[d],i]))
        sc2<-0.03*X_c
        sc<-max(sc1,sc2)
        crit[i]<-abs((X_r-X_c)/sc)
      }
      score_stad<-max(crit)
      D<-data.frame(s,d,score_stad)
      Data<-rbind(Data,D)
    }}
  return(Data)
}


# Get density estimation of mated and non-mated scores using KDE with optimal bandwidth.
# For the test scores, calculate the ratio of densities at desired test score
# For NA values with no density, 0.000000000001 is assgiend

SLR_by_density<-function(mated_score, nonmated_score, test_score){
  
  f1<-approxfun(density(mated_score))
  f2<-approxfun(density( nonmated_score))
  d1<-f1(test_score)
  d2<-f2(test_score)
  
  d1[which(is.na(d1))]<-0.000000000001
  d2[which(is.na(d2))]<-0.000000000001
  
  SLR<-d1/d2
  
  return(SLR)
}

# simulation study of random forest and modified n sigma criteria
# Two fragments from the data are selected and simulated using multivariate log-normal distribution.
# The prediction scores of two set of simulated fragments are obtained by RF and MSSC.
# Repeat entire simulation and comparison 1000 times and produce error rate by two methods.

simulation_mvn_mnsc_rf<-function(frag1, frag2, trueclass){
  
  wrong_mnsc<-NULL  
  Final_mnsc<-NULL
  D_mnsc<-NULL
  Data_mnsc<-NULL
  wrong_rf<-NULL
  Final_rf<-NULL
  D_rf<-NULL
  Data_rf<-NULL
  R<-NULL
  data<-NULL
  
  log_frag1<-log(frag1[,-c(1,2)])
  mean_log_frag1<-colMeans(log_frag1)
  cov_log_frag1<-cov(log_frag1)
  
  log_frag2<-log(frag2[,-c(1,2)])
  mean_log_frag2<-colMeans(log_frag2)
  cov_log_frag2<-cov(log_frag2)
  
  for ( j in 1:1000){
    
    X_recovered<-mvrnorm(n = 5, mu=mean_log_frag1, Sigma=cov_log_frag1, tol = 1e-6)
    X_control<-mvrnorm(n = 5, mu=mean_log_frag2, Sigma=cov_log_frag2,  tol = 1e-6)
    X_E_recovered<-exp(X_recovered)
    X_E_control<-exp(X_control)
    
    for (i in 1:18){
      X_r<-mean(X_E_recovered[,i])
      X_c<-mean(X_E_control[,i])
      Lower<-X_r/(1+4*FRSD2[i])
      Upper<-X_r*(1+4*FRSD2[i])
      R[i]<-ifelse(Lower<= X_c & X_c<= Upper, "S", "D")
    }
    
    Final_mnsc<-ifelse(sum(R=="D")>0, "D","S")
    D_mnsc<-data.frame(Final_mnsc)
    Data_mnsc<-rbind(Data_mnsc,D_mnsc)
    
    mean_log_X_r<-colMeans(X_recovered)
    mean_log_X_c<-colMeans(X_control)
    log_diff<-mean_log_X_r-mean_log_X_c
    log_diff<-matrix(log_diff,nrow = 1, byrow =TRUE)
    log_diff<-data.frame(log_diff)
    names(log_diff)<-c("Li7","Na23","Mg25","Al27","K39","Ca42","Ti49","Mn55","Fe57","Rb85","Sr88","Zr90",
                       "Ba137","La139","Ce140","Nd146","Hf178","Pb208")
    
    Final_rf<-predict(plan14_down_fit, log_diff)  
    D_rf<-data.frame(Final_rf)
    Data_rf<-rbind(Data_rf,D_rf)
    
  }
  
  Result_data<-cbind(Data_mnsc, Data_rf)
  
  wrong_mnsc<-ifelse(trueclass=="SP","D", "S")
  wrong_rf<-ifelse(trueclass=="SP", "X1","X0")
  
  error_mnsc<-nrow(subset(Result_data,Final_mnsc==wrong_mnsc))/1000 #
  error_rf<-nrow(subset(Result_data,Final_rf==wrong_rf))/1000 #0
  
  data<-data.frame(error_mnsc, error_rf)
  return(data)
  
}

# Var of AUC by Delong's method 

delong_var_AUC_RF_smote<-function(test.set){
  
  data<-plan14_test.set
  data_pred<-predict(plan14_smote_fit, data, type = "prob")[, "X0"]
  
  data_combine<- data.frame(data$type, data_pred)
  datanewY<-data_combine$data_pred[1:5990]
  datanewX<-data_combine$data_pred[5991:129795]
  rateX<-NULL
  Vx<-NULL
  var_AUC<-NULL
  for ( j in 1:length(datanewX)){
    for ( i in 1:length(datanewY)){
      if (datanewX[j]<datanewY[i]) {
        rateX[i]<-1
      } else if (datanewX[j]==datanewY[i]) {
        rateX[i]<-0.5
      } else rateX[i]<-0
    }
    Vx[j]<-mean(rateX)
  }
  
  rateY<-NULL
  Vy<-NULL
  
  for ( j in 1:length(datanewY)){
    for ( i in 1:length(datanewX)){
      if (datanewX[i]<datanewY[j]) {
        rateY[i]<-1
      } else if (datanewX[i]==datanewY[j]) {
        rateY[i]<-0.5
      } else rateY[i]<-0
    }
    Vy[j]<-mean(rateY)
    
  }
  
  var_AUC<-var(Vx)/length(Vx)+var(Vy)/length(Vy) 
  return(var_AUC)
  
}