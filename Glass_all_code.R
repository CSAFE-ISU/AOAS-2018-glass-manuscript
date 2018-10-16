########################################### Load the packages ###############################################
library(plyr)
library(dplyr)
library(ggplot2)
library(caret)
library(DMwR) 
library(ROSE)
library(randomForest)
library(reshape2)
library(corrplot)
library(e1071)
library(readr)
library(MASS)
library(ranger)
library(tidyr)
library(ROCR)
library(Boruta)
library(GGally)
library(RColorBrewer)
options(java.parameters = "-Xmx20g")
library(bartMachine)
library(ggpubr)
library(pROC)
library(matrixcalc)
library(Hotelling)
library(purrr)

######################################## Load glass data ########################################################
#1-31 panes are from CompanyA, 
# Get the files names
setwd("C:/Users/drsoy/Dropbox/Forensic Science/Glass-AOAS1211-materials/glassdata/CompanyA-glass") 
files = list.files(pattern="*.csv")
length(files) #31 panes

glass_compA<-NULL
glass_compA <- lapply(files,function(x) read.csv(x, stringsAsFactors = FALSE))

head(glass_compA[[1]])
#17 panes are from CompanyB
setwd("C:/Users/drsoy/Dropbox/Forensic Science/Glass-AOAS1211-materials/glassdata/CompanyB-glass") 
files = list.files(pattern="*.csv")
length(files) #17 panes

glass_compB<-NULL
glass_compB <- lapply(files,function(x) read.csv(x, stringsAsFactors = FALSE))

head(glass_compB[[1]])

#entire panes
all_glass<-append(glass_compA, glass_compB)
length(all_glass) #48
## Dataset 2 : Load Pane X with 34 fragments from Weis et al (2011) paper 
## Dataset 1 and 2 by Weis et al (2011) are not allowed to be in the public
#getwd()
#setwd("/home/sypark/glass")
#all_glass[[49]]<-read.csv("glass_same_NoSi_NoSn.csv")



#Log transformation first and then take the mean of replicates at each fragment
log_all_glass<-list()
mean_log_glass<-list()
for ( i in 1:length(all_glass)){
  log_all_glass[[i]]<-cbind(all_glass[[i]][,c(1:4)],log(all_glass[[i]][,-c(1:4)]))
  mean_log_glass[[i]]<-ddply(log_all_glass[[i]][,-c(1:2)], "fragment", colwise(mean))
}


## Dataset 1 : Load data with 62 fragments from all different pane from Weis(2011) paper 
#weis_diff<-read.csv('glass_diff_NoSi_NoSn.csv')
#log_weisdiff<-cbind(weis_diff[,c(1,2)],log(weis_diff[,-c(1,2)]))
#mean_log_weisdiff<-ddply(log_weisdiff, "fragment", colwise(mean))

################################################# train and test panes ####################################################
setwd("/home/sypark/glass")
#set 1 - train panes : AA-AAQ, BA-BP : 28 panes + Data1 (62 panes/fragments) from (Weis et al, 2011)
#set 2 - test pane AB-AAR, BB-BR, X : 21 panes  
plan14_set2_test<-list()
plan14_set1<-mean_log_glass[-c(2,4,6,8,10,12,17,20,22,26,28,31,33,35,37,39,41,43,46,48,49)]  #28 panes
plan14_set2_test<-mean_log_glass[c(2,4,6,8,10,12,17,20,22,26,28,31,33,35,37,39,41,43,46,48,49)] #21 panes


# Load Data1 as list by pane
#list_weis_diff<-list()
#for ( i in 1:62){
#  list_weis_diff[[i]]<-subset(mean_log_weisdiff, fragment==levels(weis_diff$fragment)[i])
#}

#add Data1 into train set
#plan14_set1_diff<-plan14_set1
#k<-length(plan14_set1_diff) #28
#plan14_set1_diff[k:(k+62-1)]<-list_weis_diff
#length(plan14_set1_diff) #90=28+62

############################################ train set - pairs from the same pane ###############################################

# Take pairwise difference between fragments within the same pane (using only train panes)
# Assign "X0" to variable type, meaning same pane 
plan14_same_type<-list()
for (i in 1: (length(plan14_set1))){
  R<- Data_take_diff(plan14_set1[[i]],"X0")
  pane1<-paste0("same",i)
  frag1<-R[,1]
  pane2<-paste0("same",i)
  frag2<-R[,2]
  plan14_same_type[[i]]<-data.frame(pane1, frag1, pane2, frag2, R[,-c(1:2)])
}

#Combine into one data frame
plan14_same_btw_iowa<-do.call(rbind,plan14_same_type)
plan14_same_btw_iowa$pane1<-as.factor(plan14_same_btw_iowa$pane1)
plan14_same_btw_iowa$pane2<-as.factor(plan14_same_btw_iowa$pane2)
plan14_same_btw_iowa$frag1<-as.factor(plan14_same_btw_iowa$frag1)
plan14_same_btw_iowa$frag2<-as.factor(plan14_same_btw_iowa$frag2)
plan14_same_btw_iowa$type<-as.factor(plan14_same_btw_iowa$type)

############################################# test set - pairs from the same pane ###############################################

# Take pairwise difference between fragments within the same pane (using only test panes)
# Assign "X0" to variable type, meaning same pane 
plan14_same_type_test<-list()
for (i in 1: (length(plan14_set2_test))){
  R<- Data_take_diff(plan14_set2_test[[i]],"X0")
  pane1<-paste0("same",i)
  frag1<-R[,1]
  pane2<-paste0("same",i)
  frag2<-R[,2]
  plan14_same_type_test[[i]]<-data.frame(pane1, frag1, pane2, frag2, R[,-c(1:2)])
}

plan14_same_btw_iowa_test<-do.call(rbind,plan14_same_type_test)
plan14_same_btw_iowa_test$pane1<-as.factor(plan14_same_btw_iowa_test$pane1)
plan14_same_btw_iowa_test$pane2<-as.factor(plan14_same_btw_iowa_test$pane2)
plan14_same_btw_iowa_test$frag1<-as.factor(plan14_same_btw_iowa_test$frag1)
plan14_same_btw_iowa_test$frag2<-as.factor(plan14_same_btw_iowa_test$frag2)
plan14_same_btw_iowa_test$type<-as.factor(plan14_same_btw_iowa_test$type)

########################################### train set - pairs from two different panes #############################################
# Take pairwise difference between fragments between training panes
# Assign "X1" to variable type, meaning different pane 

Data<-NULL
for (i in 1:(length(plan14_set1_diff)-1)){
  M1<-plan14_set1_diff[[i]][,3:20]
  
  for (j in (i+1):length(plan14_set1_diff)){
    M2<-plan14_set1_diff[[j]][,3:20]
    
    for (s in 1: nrow(M1)){
      for (d in 1: nrow(M2)){
        Info<-M1[s,]-M2[d,]
        pane1<-i
        frag1<-s
        pane2<-j
        frag2<-d
        Result<-data.frame(pane1, frag1, pane2, frag2,Info)
        Data<-rbind(Data,Result)
      }}
  }
}

Data$type<-"X1"
plan14_diff_btw_iowa<-Data
plan14_diff_btw_iowa$pane1<-as.factor(plan14_diff_btw_iowa$pane1)
plan14_diff_btw_iowa$pane2<-as.factor(plan14_diff_btw_iowa$pane2)
plan14_diff_btw_iowa$frag1<-as.factor(plan14_diff_btw_iowa$frag1)
plan14_diff_btw_iowa$frag2<-as.factor(plan14_diff_btw_iowa$frag2)
plan14_diff_btw_iowa$type<-as.factor(plan14_diff_btw_iowa$type)

############################################## test set - pairs from two different panes ###################################################
# Take pairwise difference between fragments between test panes
# Assign "X1" to variable type, meaning different pane 

Data<-NULL
for (i in 1:(length(plan14_set2_test)-1)){
  M1<-plan14_set2_test[[i]][,3:20]
  
  for (j in (i+1):length(plan14_set2_test)){
    M2<-plan14_set2_test[[j]][,3:20]
    
    for (s in 1: nrow(M1)){
      for (d in 1: nrow(M2)){
        Info<-M1[s,]-M2[d,]
        pane1<-i
        frag1<-s
        pane2<-j
        frag2<-d
        Result<-data.frame(pane1, frag1, pane2, frag2,Info)
        Data<-rbind(Data,Result)
      }
    }
  }
}  

Data$type<-"X1"
plan14_diff_btw_iowa_test<-Data
plan14_diff_btw_iowa_test$pane1<-as.factor(plan14_diff_btw_iowa_test$pane1)
plan14_diff_btw_iowa_test$pane2<-as.factor(plan14_diff_btw_iowa_test$pane2)
plan14_diff_btw_iowa_test$frag1<-as.factor(plan14_diff_btw_iowa_test$frag1)
plan14_diff_btw_iowa_test$frag2<-as.factor(plan14_diff_btw_iowa_test$frag2)
plan14_diff_btw_iowa_test$type<-as.factor(plan14_diff_btw_iowa_test$type)


#combine mated and non-mated pairs in train panes
plan14_total_iowa_data<-NULL
plan14_total_iowa_data<-rbind(plan14_same_btw_iowa,plan14_diff_btw_iowa)


#combine mated and non-mated pairs in test panes
plan14_total_iowa_data_test<-NULL
plan14_total_iowa_data_test<-rbind(plan14_same_btw_iowa_test,plan14_diff_btw_iowa_test)

plan14_total_glass<-plan14_total_iowa_data[,-c(1:4)]
plan14_total_glass_test<-plan14_total_iowa_data_test[,-c(1:4)]

plan14_train.set<-plan14_total_glass
plan14_test.set<-plan14_total_glass_test

#setwd("/home/sypark/glass")
#write.csv(plan14_train.set, "plan14_train.set.csv")
#write.csv(plan14_test.set, "plan14_test.set.csv")

#nrow(plan14_train.set) #268278
#nrow(plan14_test.set) #129795

################################### Exploratory analysis with graphs ##################################################
#FIGURE 3 : box plots
rbind_mean_log_glass<-do.call(rbind, mean_log_glass)
Manufacturer<-c(rep("CompanyA", 742), rep("CompanyB", 405), rep("X", 34))
Pane<-c(rep(1:19,each = 24), rep(20,22), rep(21:31,each = 24), rep(1,24), rep(2,23), rep(3:8,24), rep(9:10, 23), 
        rep(11:17,24), rep(1,34)) # without Data2, the last rep(1,34) should be removed.

rbind_mean_log_glass$Pane<-Pane
rbind_mean_log_glass$Manufacturer<-Manufacturer

data<-rbind_mean_log_glass
data$Pane<-as.factor(data$Pane)
data$Manufacturer<-as.factor(data$Manufacturer)

graphlist<-list()

graphlist[[1]]<-qplot(Pane, Na23, data = data, geom = "boxplot", outlier.size = 0.5)+
  facet_grid(.~Manufacturer, scales="free",space = "free")+
  labs(x="Pane", y="Na")+scale_fill_grey(start = 0, end = 1)+
  theme(legend.position="none",axis.text.x=element_blank(),text = element_text(size=18), axis.title.x=element_blank())+
  scale_y_continuous(breaks=c(11.5,11.6), limits = c(min(data$Na23),max(data$Na23)))


graphlist[[2]]<-qplot(Pane, Ti49, data = data, geom = "boxplot",outlier.size = 0.5)+
  facet_grid(.~Manufacturer, scales="free",space = "free")+
  labs(x="Pane", y="Ti")+scale_fill_grey(start = 0, end = 1)+
  theme(legend.position="none",axis.text.x=element_blank(),text = element_text(size=18), axis.title.x=element_blank())+
  scale_y_continuous(breaks=c(4.8, 5.1, 5.4), limits = c(min(data$Ti49),max(data$Ti49)))

graphlist[[3]]<-qplot(Pane, Zr90, data = data, geom = "boxplot",outlier.size = 0.5)+
  facet_grid(.~Manufacturer, scales="free",space = "free")+
  labs(x="Pane", y="Zr")+scale_fill_grey(start = 0, end = 1)+
  theme(legend.position="none",axis.text.x=element_blank(),text = element_text(size=18), axis.title.x=element_blank())+
  scale_y_continuous(breaks=c(3.3, 3.7, 4.1), limits = c(min(data$Zr90),max(data$Zr90)))

graphlist[[4]]<-qplot(Pane, Hf178, data = data, geom = "boxplot",outlier.size = 0.5)+
  facet_grid(.~Manufacturer, scales="free",space = "free")+
  labs(x="Pane", y="Hf")+scale_fill_grey(start = 0, end = 1)+
  theme(legend.position="none",axis.text.x=element_blank(),text = element_text(size=18), axis.title.x=element_blank())+
  scale_y_continuous(breaks=c(-0.5, 0, 0.5), limits = c(min(data$Hf178),max(data$Hf178)))

ggarrange(graphlist[[1]],graphlist[[2]],graphlist[[3]],graphlist[[4]], ncol=1, nrow=4)


# newer version
graphlist<-list()

graphlist[[1]]<-qplot(Pane, Na23, data = data, geom = "boxplot", outlier.size = 0.5)+
  facet_grid(.~Manufacturer, scales="free",space = "free")+
  labs(x="Pane", y="Na")+scale_fill_grey(start = 0, end = 1)+
  theme(legend.position="none",axis.text.x=element_blank(),text = element_text(size=20))+scale_y_continuous(breaks=c(11.4,11.5,11.6), limits = c(min(data$Na23),max(data$Na23)))

graphlist[[2]]<-qplot(Pane, Ti49, data = data, geom = "boxplot",outlier.size = 0.5)+
  facet_grid(.~Manufacturer, scales="free",space = "free")+
  labs(x="Pane", y="Ti")+scale_fill_grey(start = 0, end = 1)+
  theme(legend.position="none",axis.text.x=element_blank(),text = element_text(size=20))+scale_y_continuous(limits = c(min(data$Ti49),max(data$Ti49)))

graphlist[[3]]<-qplot(Pane, Zr90, data = data, geom = "boxplot",outlier.size = 0.5)+
  facet_grid(.~Manufacturer, scales="free",space = "free")+
  labs(x="Pane", y="Zr")+scale_fill_grey(start = 0, end = 1)+
  theme(legend.position="none",axis.text.x=element_blank(),text = element_text(size=20))+scale_y_continuous(limits = c(min(data$Zr90),max(data$Zr90)))

graphlist[[4]]<-qplot(Pane, Hf178, data = data, geom = "boxplot",outlier.size = 0.5)+
  facet_grid(.~Manufacturer, scales="free",space = "free")+
  labs(x="Pane", y="Hf")+scale_fill_grey(start = 0, end = 1)+
  theme(legend.position="none",axis.text.x=element_blank(),text = element_text(size=20))+scale_y_continuous(limits = c(min(data$Hf178),max(data$Hf178)))

ggarrange(graphlist[[1]],graphlist[[2]],graphlist[[3]],graphlist[[4]], ncol=1, nrow=4)


#FIGURE 4 : correlation plot

# pane CAR
ggcorr(mean_log_glass[[31]][,3:20], geom = "blank", label = TRUE, hjust = 0.75) +
  geom_point(size = 10, aes(color = coefficient > 0, alpha = abs(coefficient) > 0.5)) +
  scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) +
  guides(color = FALSE, alpha = FALSE)+
  scale_color_manual(values=c("black", "dark gray"))

# pane X (Data 2)
ggcorr(mean_log_glass[[49]][,3:20], geom = "blank", label = TRUE, hjust = 0.75) +
  geom_point(size = 10, aes(color = coefficient > 0, alpha = abs(coefficient) > 0.5)) +
  scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) +
  guides(color = FALSE, alpha = FALSE)+
  scale_color_manual(values=c("black", "dark gray"))



#FIGURE 5 : histogram 
# class "X0" means same pane (Known Match), "X1" means different panes (Known Non-Match). 
glass_data<-plan14_total_glass
glass_data$class<-glass_data$type
glass_data$class <- revalue(glass_data$class, c("X0"="KM", "X1"="KNM"))

graphlist<-list()
graphlist[[1]]<-ggplot(glass_data,aes(x=Zr90,group=class,fill=class))+
  geom_histogram(position="dodge")+theme_bw()+  labs(x = "Zr")+ geom_vline(xintercept = 0)+
  scale_fill_manual(values=c("black", "dark gray"))+theme_bw(base_size = 15)
graphlist[[2]]<-ggplot(glass_data,aes(x=Ti49,group=class,fill=class))+
  geom_histogram(position="dodge")+theme_bw()+  labs(x = "Ti")+ geom_vline(xintercept = 0)+
  scale_fill_manual(values=c("black", "dark gray"))+theme_bw(base_size = 15)
graphlist[[3]]<-ggplot(glass_data,aes(x=Hf178,group=class,fill=class))+
  geom_histogram(position="dodge")+theme_bw()+  labs(x = "Hf")+ geom_vline(xintercept = 0)+
  scale_fill_manual(values=c("black", "dark gray"))+theme_bw(base_size = 15)
graphlist[[4]]<-ggplot(glass_data,aes(x=Nd146,group=class,fill=class))+
  geom_histogram(position="dodge")+theme_bw()+  labs(x = "Nd")+ geom_vline(xintercept = 0)+
  scale_fill_manual(values=c("black", "dark gray"))+theme_bw(base_size = 15)


#figure 5
ggarrange(graphlist[[1]],graphlist[[2]],graphlist[[3]],graphlist[[4]],
          ncol=2, nrow=2, common.legend = TRUE, legend="bottom")



# Figure 13 within variance and between variance of 18 elements of 31 panes from CompanyA 
M1<-list()
for ( i in 1:31){
  M1[[i]]<-ddply(all_glass[[i]], "fragment", colwise(mean))
}

diag(cov(M1[,3:20])) #within pane cov of mean of frags

#between pane - cov of means of panes
mean_31_panes<-lapply(M1, colMeans)
rbind_mean_31_panes<-do.call(rbind,mean_31_panes)


C1<-list()
R<-NULL
D<-NULL
for ( i in 1:31){
  C1[[i]]<-ddply(all_glass[[i]], "fragment", colwise(mean))
  R<-diag(cov(C1[[i]][,3:20])) #within pane cov
  D<-rbind(D,R)}

within_pane_cov<-D
between_pane_cov<-diag(cov(rbind_mean_31_panes))[-c(1,2)] #between pane cov
within_pane_cov<-data.frame(as.matrix(within_pane_cov))


pane<-c(1:31)
within_pane_cov$pane<-pane
p1<-qplot(pane, Li7, data=within_pane_cov)+  
  geom_hline(aes(yintercept=between_pane_cov[1]), color="dark gray")+
  ylab("Li")+xlab("Pane")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 20),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
p2<-qplot(pane, Na23, data=within_pane_cov)+ 
  geom_hline(aes(yintercept=between_pane_cov[2]), color="dark gray")+
  ylab("Na")+xlab("Pane")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 20),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
p3<-qplot(pane, Mg25, data=within_pane_cov)+ 
  geom_hline(aes(yintercept=between_pane_cov[3]), color="dark gray")+
  ylab("Mg")+xlab("Pane")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 20),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
p4<-qplot(pane, Al27, data=within_pane_cov)+
  geom_hline(aes(yintercept=between_pane_cov[4]), color="dark gray")+
  ylab("Al")+xlab("Pane")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 20),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
p5<-qplot(pane, K39, data=within_pane_cov)+ 
  geom_hline(aes(yintercept=between_pane_cov[5]), color="dark gray")+
  ylab("K")+xlab("Pane")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 20),
                                panel.background = element_blank(), axis.line = element_line(colour = "black"))
p6<-qplot(pane, Ca42, data=within_pane_cov)+ 
  geom_hline(aes(yintercept=between_pane_cov[6]), color="dark gray")+
  ylab("Ca")+xlab("Pane")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 20),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
p7<-qplot(pane, Ti49, data=within_pane_cov)+ 
  geom_hline(aes(yintercept=between_pane_cov[7]), color="dark gray")+
  ylab("Ti")+xlab("Pane")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 20),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
p8<-qplot(pane, Mn55, data=within_pane_cov)+ 
  geom_hline(aes(yintercept=between_pane_cov[8]), color="dark gray")+
  ylab("Mn")+xlab("Pane")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 20),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
p9<-qplot(pane, Fe57, data=within_pane_cov)+ 
  geom_hline(aes(yintercept=between_pane_cov[9]), color="dark gray")+
  ylab("Fe")+xlab("Pane")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 20),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
p10<-qplot(pane, Rb85, data=within_pane_cov)+ 
  geom_hline(aes(yintercept=between_pane_cov[10]), color="dark gray")+
  ylab("Rb")+xlab("Pane")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 20),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
p11<-qplot(pane, Sr88, data=within_pane_cov)+ 
  geom_hline(aes(yintercept=between_pane_cov[11]), color="dark gray")+
  ylab("Sr")+xlab("Pane")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 20),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
p12<-qplot(pane, Zr90, data=within_pane_cov)+ 
  geom_hline(aes(yintercept=between_pane_cov[12]), color="dark gray")+
  ylab("Zr")+xlab("Pane")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 20),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
p13<-qplot(pane, Ba137, data=within_pane_cov)+ 
  geom_hline(aes(yintercept=between_pane_cov[13]), color="dark gray")+
  ylab("Ba")+xlab("Pane")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 20),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
p14<-qplot(pane, La139, data=within_pane_cov)+ 
  geom_hline(aes(yintercept=between_pane_cov[14]), color="dark gray")+
  ylab("La")+xlab("Pane")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 20),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
p15<-qplot(pane, Ce140, data=within_pane_cov)+ 
  geom_hline(aes(yintercept=between_pane_cov[15]), color="dark gray")+
  ylab("Ce")+xlab("Pane")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 20),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
p16<-qplot(pane, Nd146, data=within_pane_cov)+ 
  geom_hline(aes(yintercept=between_pane_cov[16]), color="dark gray")+
  ylab("Nd")+xlab("Pane")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 20),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
p17<-qplot(pane, Hf178, data=within_pane_cov)+ 
  geom_hline(aes(yintercept=between_pane_cov[17]), color="dark gray")+
  ylab("Hf")+xlab("Pane")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 20),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
p18<-qplot(pane, Pb208, data=within_pane_cov)+ 
  geom_hline(aes(yintercept=between_pane_cov[18]), color="dark gray")+
  ylab("Pb")+xlab("Pane")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 20),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))

#ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,
#          p10,p11,p12,p13,p14,p15,p16,p17,p18, ncol=3, nrow=6)


ggarrange(p13, p5, p7, p14, p3, p6, ncol=3, nrow=2)



##################################### Random forest with 5 sampling tehnique ################################################
ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 3,
                     savePredictions="final",
                     summaryFunction = twoClassSummary,
                     classProbs = TRUE
)


#random forest without any sampling method (original version)
# mtry is recommended to use sqrt(# of variables)=sqrt(18) so 3 and 4 are tried to find the optimal one
fit_plan14_rf <- train (type~.,
                        plan14_train.set,
                        method = "rf",
                        tuneGrid = expand.grid(
                          .mtry = c(3,4)),
                        metric="ROC",
                        preProc=c("center", "scale"),
                        trControl = ctrl)

pred_plan14_rf <- predict(fit_plan14_rf, plan14_test.set)                           
table_plan14.rf<-confusionMatrix(pred_plan14_rf, plan14_test.set$type)$table
rf.caret.er_plan14<- round(1-sum(diag(table_plan14.rf))/sum(table_plan14.rf), 4)

#class probability for the same pane is used as similarity score for RF
pred_plan14_prob_rf <- predict(fit_plan14_rf, plan14_test.set,type="prob")  


###### Build down-sampled model

# Use the same seed to ensure same cross-validation splits
ctrl$seeds <- fit_plan14_rf$control$seeds
ctrl$sampling <- "down"

plan14_down_fit <- train(type ~ .,
                         data = plan14_train.set,
                         method = "rf",
                         tuneGrid = expand.grid(
                           .mtry = c(3,4)),
                         verbose = FALSE,
                         #weights = model_weights,
                         preProc=c("center", "scale"),
                         metric = "ROC",
                         trControl = ctrl)

#class probability for the same pane is used as similarity score for RF - down sampling
pred_plan14_prob_down_rf <- predict(plan14_down_fit, plan14_test.set,type="prob")
pred_plan14_down_rf <- predict(plan14_down_fit, plan14_test.set)                           
table_plan14_down.rf<-confusionMatrix(pred_plan14_down_rf, plan14_test.set$type)$table
rf.caret.er_down_plan14<- round(1-sum(diag(table_plan14_down.rf))/sum(table_plan14_down.rf), 4)

# Build smote model

ctrl$sampling <- "smote"
plan14_smote_fit <- train(type ~ .,
                          data = plan14_train.set,
                          method = "rf",
                          tuneGrid = expand.grid(
                            .mtry = c(3,4)),
                          verbose = FALSE,
                          #weights = model_weights,
                          preProc=c("center", "scale"),
                          metric = "ROC",
                          trControl = ctrl)

#class probability for the same pane is used as similarity score for RF-SMOTE sampling
pred_plan14_prob_smote_rf <- predict(plan14_smote_fit, plan14_test.set,type="prob")
pred_plan14_smote_rf <- predict(plan14_smote_fit, plan14_test.set)                           
table_plan14_smote.rf<-confusionMatrix(pred_plan14_smote_rf, plan14_test.set$type)$table
rf.caret.er_smote_plan14<- round(1-sum(diag(table_plan14_smote.rf))/sum(table_plan14_smote.rf), 4)

# Build weighted model

# Create model weights (they sum to one)
model_weights <- ifelse(plan14_train.set$type == "X0",
                        (1/table(plan14_train.set$type)[1])*0.5,
                        (1/table(plan14_train.set$type)[2])*0.5)

sum(model_weights) #1
plan14_weighted_fit <- train(type ~ .,
                             data = plan14_train.set,
                             method = "rf",
                             tuneGrid = expand.grid(
                               .mtry = c(3,4)),
                             verbose = FALSE,
                             weights = model_weights,
                             preProc=c("center", "scale"),
                             metric = "ROC",
                             trControl = ctrl)

#class probability for the same pane is used as similarity score for RF - weighted 
pred_plan14_prob_weighted_rf <- predict(plan14_weighted_fit, plan14_test.set,type="prob")

# Build up-sampled model

ctrl$sampling <- "up"
plan14_up_fit <- train(type ~ .,
                       data = plan14_train.set,
                       method = "rf",
                       tuneGrid = expand.grid(
                         .mtry = c(3,4)),
                       verbose = FALSE,
                       #weights = model_weights,
                       preProc=c("center", "scale"),
                       metric = "ROC",
                       trControl = ctrl)

#class probability for the same pane is used as similarity score for RF- up sampling
pred_plan14_prob_up_rf <- predict(plan14_up_fit, plan14_test.set,type="prob")

# Build rose model

ctrl$sampling <- "rose"
plan14_rose_fit <- train(type ~ .,
                         data = plan14_train.set,
                         method = "rf",
                         tuneGrid = expand.grid(
                           .mtry = c(3,4)),
                         verbose = FALSE,
                         #weights = model_weights,
                         preProc=c("center", "scale"),
                         metric = "ROC",
                         trControl = ctrl)

#class probability for the same pane is used as similarity score for RF - ROSE sampling
pred_plan14_prob_rose_rf <- predict(plan14_rose_fit, plan14_test.set,type="prob")

# Create ROC dataset of 5 random forest with different sampling techniques
plan14_rf_list <- list(original = fit_plan14_rf,
                       weighted = plan14_weighted_fit,
                       down = plan14_down_fit,
                       up = plan14_up_fit,
                       SMOTE = plan14_smote_fit,
                       ROSE= plan14_rose_fit)

plan14_rf_list_roc <- plan14_rf_list %>%
  map(test_roc, data = plan14_test.set)

results_rf_list_roc_14 <- list(NA)
num_mod <- 1

for(the_roc in plan14_rf_list_roc){
  
  results_rf_list_roc_14[[num_mod]] <- 
    data_frame(tpr = the_roc$sensitivities,
               fpr = 1 - the_roc$specificities,
               model = names(plan14_rf_list)[num_mod])
  
  num_mod <- num_mod + 1
  
}
plan14_results_df_roc <- do.call(rbind,results_rf_list_roc_14)


#### FIGURE 6 : AUC, Sensitivity and specificity of RF with 5 sampling technique 
resampling<-resamples(plan14_rf_list)
bwtheme  <- canonical.theme(color = FALSE)
bwplot(resampling, par.settings = bwtheme)
bwplot(resampling)

### AUC of the ROC by DeLong's method
data<-plan14_test.set
roc_rf_down<-pROC::roc(data$type, predict(plan14_down_fit, data, type = "prob")[, "X0"], levels=c("X1","X0"),ci=TRUE)
# 95% CI: 0.976-0.9781 (DeLong)
roc_rf_SMOTE<-pROC::roc(data$type,predict( plan14_smote_fit, data, type = "prob")[, "X0"], levels=c("X1","X0"),ci=TRUE)
# 95% CI: 0.9773-0.9793 (DeLong)

### var and cov of ROC 
pROC::cov(roc_rf_down, roc_rf_SMOTE, method="delong") #2.654748e-07
pROC::var(roc_rf_down, method="delong") #2.853215e-07
pROC::var(roc_rf_SMOTE, method="delong") #2.630052e-07
delong_var_AUC_RF_smote(plan14_test.set) #2.630052e-07

# For getting a 95 % conservative (wider) CI of AUC using number of independent panes, the function delong_var_AUC_RF_smote is revised. 
# Instead of using number of pairwise comparisons of mated and non-mated pairs, the denominator is changed into 21 as number of test panes.


##### FIGURE 7 : ROC curves for RF by sampling techniques
plan14_results_df_roc <- do.call(rbind,results_rf_list_roc_14 )
custom_col <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")

graphlist<-list()

graphlist[[1]]<-ggplot(aes(x = fpr,  y = tpr, group = model), data =plan14_results_df_roc) +
  geom_line(aes(color = model, linetype=model), size = 1) +
  scale_color_manual(values = custom_col) +
  geom_abline(intercept = 0, slope = 1, color = "gray", size = 1) +
  theme_bw(base_size = 12)+xlab("1-Specificity (FPR)") +ylab("Sensitivity (TPR)")

graphlist[[2]]<-ggplot(aes(x = fpr,  y = tpr, group = model), data = plan14_results_df_roc) +
  geom_line(aes(color = model, linetype=model), size = 1) +
  scale_color_manual(values = custom_col) +
  geom_abline(intercept = 0, slope = 1, color = "gray", size = 1) +
  theme_bw(base_size = 12)+xlab("1-Specificity (FPR) < 0.3") +ylab("Sensitivity (TPR) > 0.9")+
  xlim(0.05, 0.3)+ylim(0.9, 1)

ggarrange(graphlist[[1]],graphlist[[2]], ncol=1, nrow=2, common.legend = TRUE, legend="bottom")


ROCdata_plan14_weighted_rf<-ROC.data(pred_plan14_prob_weighted_rf,plan14_test.set,thres, "weighted_RF")
ROCdata_plan14_down_rf<-ROC.data(pred_plan14_prob_down_rf,plan14_test.set,thres, "down_RF")
ROCdata_plan14_up_rf<-ROC.data(pred_plan14_prob_up_rf,plan14_test.set,thres, "up_RF")
ROCdata_plan14_smote_rf<-ROC.data(pred_plan14_prob_smote_rf,plan14_test.set,thres, "smote_RF")
ROCdata_plan14_rose_rf<-ROC.data(pred_plan14_prob_rose_rf,plan14_test.set,thres, "rose_RF")


# TABLE 2 : TPR at fixed 5%, 10%, 15% FPR
ROCdata_RF<-rbind(ROCdata_plan14_rf,ROCdata_plan14_weighted_rf,#ROCdata_plan14_rose_rf,
                  ROCdata_plan14_down_rf,ROCdata_plan14_up_rf,ROCdata_plan14_smote_rf)
P<-c(0.05, 0.1, 0.15)
TPRs_RF<-get_TPR_fixed_FPR_combined(ROCdata_RF,P)
TPRs_RF

##################################### variable importance ######################################
rf_down_imp<-varImp(plan14_down_fit, scale = TRUE)
plot(rf_down_imp, top=18, col="black")

rf_smote_imp<-varImp(plan14_smote_fit, scale = TRUE)
plot(rf_smote_imp, top=18, col="black")
################################### BART #################################################

#Fit BART with down sampled data

plan14_train_X<-plan14_train.set[,1:18]
plan14_train_Y<-plan14_train.set[,19]

plan14_test_X<-plan14_test.set[,1:18]
plan14_test_Y<-plan14_test.set[,19]

### down sample majority class into the same number of minority class
down_plan14_train.set<-downSample(plan14_train_X, plan14_train_Y)

names(down_plan14_train.set)[19]<-"type"
table(down_plan14_train.set$type)
# X0   X1 
#7705 7705 

plan14_train_down_X<-down_plan14_train.set[,1:18]
plan14_train_down_Y<-down_plan14_train.set[,19]

#10 fold cross-validation BART
fit_plan14_down_bart <- bartMachineCV(plan14_train_down_X, plan14_train_down_Y, k_folds = 10, num_tree_cvs = 100, serialize = TRUE)

#variable importance 
vs_bart<-var_selection_by_permute(fit_plan14_down_bart,num_permute_samples = 10, bottom_margin = 10)
vs_bart$important_vars_local_names
vs_bart$important_vars_global_max_names


#explore performance on test data
plan14_down_oos_perf = bart_predict_for_test_data(fit_plan14_down_bart, plan14_test_X, plan14_test_Y)
table_plan14_down.bart<-plan14_down_oos_perf$confusion_matrix
bart.caret.er_down_plan14<- round(table_plan14_down.bart[3,3], 4)

#class probability for the same pane is used as similarity score for BART
pred_plan14_prob_down_bart <- predict(fit_plan14_down_bart, plan14_test_X, type="prob")  

#### FIGURE 8 : class probability by BART of 100 random KM and KNM 
table(plan14_test_Y) # In test set, #X0 (mated pairs) :5990 # X1 (non-mated pairs) :123805 

# Randomly select 100 from 1 to 5990, also select 100 from 5991 to 129795
s<-sample(1:5990, 100)
d<-sample(5991:129795, 100)

inter<-calc_credible_intervals(fit_plan14_down_bart, plan14_test_X[s, ])
ind<-1:100
pred<-predict(fit_plan14_down_bart, plan14_test_X[s, ], type="prob")
credible_KM<-data.frame(ind,pred,inter )

# KM 
data<-credible_KM
newdata <- data[order(-pred),] 
newdata$ind  <- factor(newdata$ind, levels = newdata$ind[order(-newdata$pred)])
p11<-ggplot(newdata, aes(x=ind, y=pred)) +
  geom_point(aes(x=ind, y=pred)) +
  geom_errorbar(aes(ymin=ci_lower_bd, ymax=ci_upper_bd), colour="gray60", width=.5, position=pd) +
  ylab("Class Probability")+xlab("KM")+theme_bw(base_size = 15)+
  theme(axis.ticks.x=element_blank(),axis.text.x=element_blank())


# KNM
inter2<-calc_credible_intervals(fit_plan14_down_bart, plan14_test_X[d, ])
pred2<-predict(fit_plan14_down_bart, plan14_test_X[d, ], type="prob")
ind<-1:100
credible_KNM<-data.frame(ind,pred2,inter2)
head(credible_KNM)

data2<-credible_KNM
newdata2 <- data2[order(pred),] 
newdata2$ind  <- factor(newdata2$ind, levels = newdata2$ind[order(-newdata2$pred)])

p22<-ggplot(newdata2, aes(x=ind, y=pred2)) +
  geom_point(aes(x=ind, y=pred2)) +
  geom_errorbar(aes(ymin=ci_lower_bd, ymax=ci_upper_bd), colour="gray60", width=.5, position=pd) +
  ylab("Class Probability")+xlab("KNM")+theme_bw(base_size = 15)+
  theme(axis.ticks.x=element_blank(),axis.text.x=element_blank())

#FIGURE 8
ggarrange(p11, p22)

######################################### testing with 3 controls (K) and 1 recovered (Q) ###############################
# For each pane, pick Q first and three random fragments (K1,K2,K3) are selected in the same pane, which is not the same
# with Q. Generate the list containing the information which number of fragments are selected in each pane. 
# In each list, first column is fragment Q and rest of three fragments are K1, K2, K3. 

test_panes<-all_glass[c(2,4,6,8,10,12,17,20,22,26,28,31,33,35,37,39,41,43,46,48,49)] 
list_sample_ind_3control<-list()

for (i in 1:(length(test_panes))){
  originaldata<-test_panes[[i]]
  originaldata$fragment<-as.factor(originaldata$fragment)
  nfrag<-length(levels(originaldata$fragment))
  #level.frag<-as.numeric(levels(originaldata$fragment))
  sample_ind<-NULL
  
  for ( j in 1:nfrag){
    row<-1:nfrag  
    row1<-row[-j]
    Info<-NULL
    for (h in 1:30){
      control<-sample(row1, 3)
      Info<-rbind(Info, c(j, control))}
    sample_ind<-rbind(sample_ind, Info )
  }
  
  list_sample_ind_3control[[i]]<-data.frame(sample_ind)
}


## Generate differences between Q and three Ks as data frame for getting prediction by RF
## Q and three K's are from the same pane
rf_test_data_same_3control<-list()

for (i in 1:(length(test_panes))){
  originaldata<-test_panes[[i]]
  originaldata$fragment<-as.factor(originaldata$fragment)
  S<-list_sample_ind_3control[[i]]
  rf_test_data_same_3control[[i]]<-rf_test_same_ncontrol(originaldata,S)
}
plan16_rf_test_data_same_3control<-do.call(rbind,rf_test_data_same_3control)


# Generate differences between Q and three Ks as data frame for getting prediction by RF 
## Q and three K's are from two different panes
D<-NULL
Data<-NULL
#plan16_rf_test_data_diff_3control<-NULL
for (j in 1:(length(test_panes)-1)){
  
  Data1<-test_panes[[j]]
  Data1$fragment<-as.factor(Data1$fragment)
  nfrag1<-length(levels(Data1$fragment))
  S1<-list_sample_ind_3control[[j]]
  
  for ( k in (j+1): length(test_panes)){
    
    Data2<-test_panes[[k]]
    Data2$fragment<-as.factor(Data2$fragment)
    nfrag2<-length(levels(Data2$fragment))
    S2<-list_sample_ind_3control[[k]]
    D<-rf_test_diff_ncontrol(Data1, Data2, nfrag1, S2)
    Data<-rbind(Data,D)}
  
}

#plan16_rf_test_data_diff_3control<-NULL
plan16_rf_test_data_diff_3control<-Data

#Create difference of log transformed values bewteen Q and mean of K. 
plan16_rf_test_data_combined_3control<-rbind(plan16_rf_test_data_same_3control,plan16_rf_test_data_diff_3control)
plan16_test.set<-plan16_rf_test_data_combined_3control

### Get class probability from random forest with down sampling 
pred_plan16_prob_rf_3control<- predict(plan14_down_fit, plan16_rf_test_data_combined_3control, type="prob")

### Get binary prediction of class from random forest with down sampling
pred_plan16_rf_3control <- predict(plan14_down_fit, plan16_rf_test_data_combined_3control)                           
table_plan16_rf_3control<-confusionMatrix(pred_plan16_rf_3control, plan16_rf_test_data_combined_3control$type)$table

table<-table_plan16_rf_3control
error_FPR_3control_rf<-table[1,2]/(table[1,2]+table[2,2]) #FPR 0.09636812
error_FNR_3control_rf<-table[2,1]/(table[1,1]+table[2,1]) #FNR 0.02352941


###### Predict the class between one Q and three K's by BART
plan16_3control_test_X<-plan16_rf_test_data_combined_3control[,1:18]
plan16_3control_test_Y<-as.factor(plan16_rf_test_data_combined_3control[,19])

#explore performance on test data
plan16_3control_oos_perf = bart_predict_for_test_data(fit_plan14_down_bart, plan16_3control_test_X, plan16_3control_test_Y)
table_plan16_3control.bart<-plan16_3control_oos_perf$confusion_matrix

### Get class probability from BART
pred_plan16_prob_bart_3control<- predict(fit_plan14_down_bart, plan16_3control_test_X, type="prob")  

table<-table_plan16_3control.bart
error_FPR_3control_bart<-table[2,1]/(table[2,2]+table[2,1]) #FPR 0.09538851
error_FNR_3control_bart<-table[1,2]/(table[1,1]+table[1,2]) #FNR 0.02196078


# Get random forest score from the three different comparison between Q-K1, Q-K2, Q-K3.
# Then get average of three scores, minimum and maximum of scores. 
# Q-K1, Q-K2, Q-K3 comparisons and average or minimum 
# Q and three K's are from the same pane
rf_test_score_same_3indep<-list()

for (i in 1:(length(test_panes))){
  originaldata<-test_panes[[i]]
  originaldata$fragment<-as.factor(originaldata$fragment)
  #nfrag<-length(levels(originaldata$fragment))
  S<-list_sample_ind[[i]]
  rf_test_score_same_3indep[[i]]<-rf_test_same_3indep(originaldata,S)
}

#plan16_rf_test_score_same_3indep<-NULL
plan16_rf_test_score_same_3indep<-do.call(rbind,rf_test_score_same_3indep)



# Get random forest score from the three different comparison between Q-K1, Q-K2, Q-K3.
# Then get average of three scores, minimum and maximum of scores. 
# Q-K1, Q-K2, Q-K3 comparisons and average or minimum 
# Q and three K's are from two different panes

D<-NULL
Data<-NULL
#plan16_rf_test_score_diff_3indep<-NULL
for (j in 1:(length(test_panes)-1)){
  
  Data1<-test_panes[[j]]
  Data1$fragment<-as.factor(Data1$fragment)
  nfrag1<-length(levels(Data1$fragment))
  S1<-list_sample_ind[[j]]
  
  for ( k in (j+1): length(test_panes)){
    
    Data2<-test_panes[[k]]
    Data2$fragment<-as.factor(Data2$fragment)
    nfrag2<-length(levels(Data2$fragment))
    S2<-list_sample_ind[[k]]
    
    D<-rf_test_diff_3indep(Data1, Data2, nfrag1, S2)
    Data<-rbind(Data,D)}
  
}
plan16_rf_test_score_diff_3indep<-Data
plan16_rf_test_score_3indep_combined<-rbind(plan16_rf_test_score_same_3indep,plan16_rf_test_score_diff_3indep)

### TPR, FPR, TNR, FNR by mean scores of three comparisons (Q-K1, Q-K2, Q-K3)
same<-plan16_rf_test_score_same_3indep$mean_score
diff<-plan16_rf_test_score_diff_3indep$mean_score

oob_err_RF_mean<-RF_err_col(same, diff)
#        Sen one_minus_spec     Spec one_minus_sen overall_err
# 0.9752941     0.09087698 0.909123    0.02470588  0.08475448

### TPR, FPR, TNR, FNR by max score of three comparisons (Q-K1, Q-K2, Q-K3)
same<-plan16_rf_test_score_same_3indep$max_score
diff<-plan16_rf_test_score_diff_3indep$max_score

oob_err_RF_max<-RF_err_col(same, diff)
#        Sen one_minus_spec      Spec one_minus_sen overall_err
# 0.9922876      0.1272291 0.8727709   0.007712418   0.1161708

### TPR, FPR, TNR, FNR by min score of three comparisons (Q-K1, Q-K2, Q-K3)
same<-plan16_rf_test_score_same_3indep$min_score
diff<-plan16_rf_test_score_diff_3indep$min_score

oob_err_RF_min<-RF_err_col(same, diff)
#Sen one_minus_spec      Spec one_minus_sen overall_err
# 0.8710458     0.05639744 0.9436026     0.1289542  0.06311079






# Get BART score from the three different comparison between Q-K1, Q-K2, Q-K3.
# Then get average of three scores, minimum and maximum of scores. 
# Q and three K's are from the same pane

bart_test_score_same_3indep<-list()

for (i in 1:(length(test_panes))){
  originaldata<-test_panes[[i]]
  originaldata$fragment<-as.factor(originaldata$fragment)
  #nfrag<-length(levels(originaldata$fragment))
  S<-list_sample_ind[[i]]
  bart_test_score_same_3indep[[i]]<-bart_test_same_3indep(originaldata,S)
}

plan16_bart_test_score_same_3indep<-NULL
plan16_bart_test_score_same_3indep<-do.call(rbind,bart_test_score_same_3indep)


# Get BART score from the three different comparison between Q-K1, Q-K2, Q-K3.
# Then get average of three scores, minimum and maximum of scores. 
# Q and three K's are from two different panes
D<-NULL
Data<-NULL
#plan16_bart_test_score_diff_3indep<-NULL
for (j in 1:(length(test_panes)-1)){
  
  Data1<-test_panes[[j]]
  Data1$fragment<-as.factor(Data1$fragment)
  nfrag1<-length(levels(Data1$fragment))
  S1<-list_sample_ind[[j]]
  
  for ( k in (j+1): length(test_panes)){
    
    Data2<-test_panes[[k]]
    Data2$fragment<-as.factor(Data2$fragment)
    nfrag2<-length(levels(Data2$fragment))
    S2<-list_sample_ind[[k]]
    
    D<-bart_test_diff_3indep(Data1, Data2, nfrag1, S2)
    Data<-rbind(Data,D)}
  
}
plan16_bart_test_score_diff_3indep<-Data

#plan16_bart_test_score_3indep_combined<-NULL
plan16_bart_test_score_3indep_combined<-rbind(plan16_bart_test_score_same_3indep,plan16_bart_test_score_diff_3indep)



### TPR, FPR, TNR, FNR by mean score of three comparisons (Q-K1, Q-K2, Q-K3)
same<-plan16_bart_test_score_same_3indep$mean_score
diff<-plan16_bart_test_score_diff_3indep$mean_score
oob_err_bart_mean<-RF_err_col(same, diff)
#Sen         one_minus_spec      Spec one_minus_sen overall_err
# 0.9736601     0.08645209 0.9135479    0.02633987  0.08089018


### TPR, FPR, TNR, FNR by max score of three comparisons (Q-K1, Q-K2, Q-K3)
same<-plan16_bart_test_score_same_3indep$max_score
diff<-plan16_bart_test_score_diff_3indep$max_score
oob_err_bart_max<-RF_err_col(same, diff)
# Sen one_minus_spec      Spec one_minus_sen overall_err
# 0.9949673      0.1324204 0.8675796    0.00503268   0.1206338

### TPR, FPR, TNR, FNR by min score of three comparisons (Q-K1, Q-K2, Q-K3)
same<-plan16_bart_test_score_same_3indep$min_score
diff<-plan16_bart_test_score_diff_3indep$min_score
oob_err_bart_min<-RF_err_col(same, diff)
#       Sen    one_minus_spec      Spec one_minus_sen overall_err
# 0.8658824     0.04904038 0.9509596     0.1341176  0.05691219





# Get Standard - E2330 score from the three different comparison between Q-K1, Q-K2, Q-K3.
# Then get average of three scores, minimum and maximum of scores. 
# Q and three K's are from the same pane

stad_same_score_3control<-list()

for (i in 1:(length(test_panes))){
  originaldata<-test_panes[[i]]
  originaldata$fragment<-as.factor(originaldata$fragment)
  #nfrag<-length(levels(originaldata$fragment))
  S<-list_sample_ind_3control[[i]]
  stad_same_score_3control[[i]]<-stad_test_same_ncontrol(originaldata,S)
}
plan16_stad_same_test_score_3control<-do.call(rbind,stad_same_score_3control)



# Get Standard - E2330 score from the three different comparison between Q-K1, Q-K2, Q-K3.
# Then get average of three scores, minimum and maximum of scores. 
# Q and three K's are from two different panes

D<-NULL
Data<-NULL

for (j in 1:(length(test_panes)-1)){
  
  Data1<-test_panes[[j]]
  Data1$fragment<-as.factor(Data1$fragment)
  nfrag1<-length(levels(Data1$fragment))
  S1<-list_sample_ind_3control[[j]]
  
  for ( k in (j+1): length(test_panes)){
    
    Data2<-test_panes[[k]]
    Data2$fragment<-as.factor(Data2$fragment)
    nfrag2<-length(levels(Data2$fragment))
    S2<-list_sample_ind_3control[[k]]
    
    D<-stad_test_diff_ncontrol(Data1, Data2, nfrag1, S2)
    Data<-rbind(Data,D)}
  
}

#plan16_stad_diff_test_score_3control<-NULL
plan16_stad_diff_test_score_3control<-Data

plan16_stad_test_score_3control<-c(plan16_stad_same_test_score_3control[,3],plan16_stad_diff_test_score_3control[,3])

# TPR, FPR, TNR, FNR by standard method
oob_err_stad_3control<-stad_err(plan16_stad_same_test_score_3control,plan16_stad_diff_test_score_3control)
#       Sen one_minus_spec      Spec one_minus_sen overall_err
# 0.9441176      0.1866453 0.8133547    0.05588235   0.1745464
error_FPR_3control_stad<-oob_err_stad_3control[,2]
error_FNR_3control_stad<-oob_err_stad_3control[,4]


# Fixed relative standard deviation for 18 elements (FRSD - Weis et al, 2011)
# first two are just zeros (meaning nothing)
FRSD<-c(0.000, 0.000, 0.030, 0.030, 0.030, 0.055, 0.030, 0.037, 0.036, 0.030, 0.030,
        0.030, 0.039, 0.081, 0.033, 0.060, 0.030, 0.057, 0.089, 0.038)

# Get Modified n sigma criteria score (Weis et al, 2011) from the three different comparison between Q-K1, Q-K2, Q-K3.
# Then get average of three scores, minimum and maximum of scores. 
# Q and three K's are from the same pane

mnsc_same_score_3control<-list()

for (i in 1:(length(test_panes))){
  originaldata<-test_panes[[i]]
  originaldata$fragment<-as.factor(originaldata$fragment)
  #nfrag<-length(levels(originaldata$fragment))
  S<-list_sample_ind_3control[[i]]
  mnsc_same_score_3control[[i]]<-mnsc_test_same_ncontrol(originaldata,S)
}


plan16_mnsc_same_test_score_3control<-do.call(rbind,mnsc_same_score_3control)
nrow(plan16_mnsc_same_test_score_3control) #15300
head(plan16_mnsc_same_test_score_3control)

# Get Modified n sigma criteria score (Weis et al, 2011) from the three different comparison between Q-K1, Q-K2, Q-K3.
# Then get average of three scores, minimum and maximum of scores. 
# Q and three K's are from two different panes

plan16_mnsc_diff_test_score<-NULL
D<-NULL
Data<-NULL

for (j in 1:(length(test_panes)-1)){
  
  Data1<-test_panes[[j]]
  Data1$fragment<-as.factor(Data1$fragment)
  nfrag1<-length(levels(Data1$fragment))
  S1<-list_sample_ind_3control[[j]]
  
  for ( k in (j+1): length(test_panes)){
    
    Data2<-test_panes[[k]]
    Data2$fragment<-as.factor(Data2$fragment)
    nfrag2<-length(levels(Data2$fragment))
    S2<-list_sample_ind_3control[[k]]
    
    D<-mnsc_test_diff_ncontrol(Data1, Data2, nfrag1, S2)
    Data<-rbind(Data,D)}
  
}

#plan16_mnsc_diff_test_score_3control<-NULL
plan16_mnsc_diff_test_score_3control<-Data
plan16_mnsc_test_score_3control<-c(plan16_mnsc_same_test_score_3control[,3],plan16_mnsc_diff_test_score_3control[,3])

# TPR, FPR, TNR, FNR by MNSC method
oob_err_mnsc_3control<-stad_err(plan16_mnsc_same_test_score_3control,plan16_mnsc_diff_test_score_3control)
#        Sen one_minus_spec      Spec one_minus_sen overall_err
# 0.5517647      0.0627549 0.9372451     0.4482353  0.09842163
error_FPR_3control_mnsc<-oob_err_mnsc_3control[,2]
error_FNR_3control_mnsc<-oob_err_mnsc_3control[,4]



# Get Hotelling T^2 statistic with shrinkage covariance estimator (Campbell and Curran, 2009) as score 
# from the three different comparison between Q-K1, Q-K2, Q-K3.
# Then get average of three scores, minimum and maximum of scores. 
# Q and three K's are from the same pane
hot_same_score_3control<-list()

for (i in 1:(length(test_panes))){
  originaldata<-test_panes[[i]]
  originaldata$fragment<-as.factor(originaldata$fragment)
  #nfrag<-length(levels(originaldata$fragment))
  S<-list_sample_ind_3control[[i]]
  hot_same_score_3control[[i]]<-hot_test_same_ncontrol(originaldata,S)
}

plan16_hot_same_test_score_3control<-do.call(rbind,hot_same_score_3control)

# Get Hotelling T^2 statistic with shrinkage covariance estimator (Campbell and Curran, 2009) as score 
# from the three different comparison between Q-K1, Q-K2, Q-K3.
# Then get average of three scores, minimum and maximum of scores. 
# Q and three K's are from two different panes

plan16_hot_diff_test_score<-NULL
D<-NULL
Data<-NULL

for (j in 1:(length(test_panes)-1)){
  
  Data1<-test_panes[[j]]
  Data1$fragment<-as.factor(Data1$fragment)
  nfrag1<-length(levels(Data1$fragment))
  S1<-list_sample_ind_3control[[j]]
  
  for ( k in (j+1): length(test_panes)){
    
    Data2<-test_panes[[k]]
    Data2$fragment<-as.factor(Data2$fragment)
    nfrag2<-length(levels(Data2$fragment))
    S2<-list_sample_ind_3control[[k]]
    
    D<-hot_test_diff_ncontrol(Data1, Data2, nfrag1, S2)
    Data<-rbind(Data,D)}
  
}

#plan16_hot_diff_test_score_3control<-NULL
plan16_hot_diff_test_score_3control<-Data

plan16_hot_test_score_3control<-c(plan16_hot_same_test_score_3control[,3],plan16_hot_diff_test_score_3control[,3])

# TPR, FPR, TNR, FNR by Hotelling T^2 statistic
oob_err_hot_3control<-hot_err(plan16_hot_same_test_score_3control,plan16_hot_diff_test_score_3control)
#         Sen one_minus_spec      Spec one_minus_sen overall_err
# 0.05411765    0.000193256 0.9998067     0.9458824  0.08769352

error_FPR_3control_hot<-oob_err_hot_3control[,2]
error_FNR_3control_hot<-oob_err_hot_3control[,4]



# Get Optimum test statistic (Parker and Holford, 1989) as score from the three different comparison between Q-K1, Q-K2, Q-K3.
# Then get average of three scores, minimum and maximum of scores. 
# Q and three K's are from the same pane

optim_same_score_3control<-list()

for (i in 1:(length(test_panes))){
  originaldata<-test_panes[[i]]
  originaldata$fragment<-as.factor(originaldata$fragment)
  #nfrag<-length(levels(originaldata$fragment))
  S<-list_sample_ind_3control[[i]]
  optim_same_score_3control[[i]]<-optim_test_same_ncontrol(originaldata,S)
}

plan16_optim_same_test_score_3control<-do.call(rbind,optim_same_score_3control)

# Get Optimum test statistic (Parker and Holford, 1989) as score from the three different comparison between Q-K1, Q-K2, Q-K3.
# Then get average of three scores, minimum and maximum of scores. 
# Q and three K's are from the same pane

plan16_optim_diff_test_score_3control<-NULL
D<-NULL
Data<-NULL

for (j in 1:(length(test_panes)-1)){
  
  Data1<-test_panes[[j]]
  Data1$fragment<-as.factor(Data1$fragment)
  nfrag1<-length(levels(Data1$fragment))
  S1<-list_sample_ind_3control[[j]]
  
  for ( k in (j+1): length(test_panes)){
    
    Data2<-test_panes[[k]]
    Data2$fragment<-as.factor(Data2$fragment)
    nfrag2<-length(levels(Data2$fragment))
    S2<-list_sample_ind_3control[[k]]
    
    D<-optim_test_diff_ncontrol(Data1, Data2, nfrag1, S2)
    Data<-rbind(Data,D)}
  
}

#plan16_optim_diff_test_score_3control<-NULL
plan16_optim_diff_test_score_3control<-Data

# TPR, FPR, TNR, FNR by Optimum test statistic
optim_err(plan16_optim_same_test_score_3control,plan16_optim_diff_test_score_3control)
oob_err_optim_3control<-optim_err(plan16_optim_same_test_score_3control,plan16_optim_diff_test_score_3control)
#        Sen one_minus_spec     Spec one_minus_sen overall_err
# 0.2343137    0.004165001    0.995835     0.7656863  0.07462506

plan16_optim_test_score_3control<-c(plan16_optim_same_test_score_3control[,3],plan16_optim_diff_test_score_3control[,3])
optim_err(plan16_optim_same_test_score_3control,plan16_optim_diff_test_score_3control)

plan16_optim_test_score_3control<-c(plan16_optim_same_test_score_3control[,3],plan16_optim_diff_test_score_3control[,3])
error_FPR_3control_optim<-oob_err_optim_3control[,2]
error_FNR_3control_optim<-oob_err_optim_3control[,4]

################################################ ROC curve #############################################################
ROC_list_3control_ver2<-list()

ROC_list_3control_ver2[[1]]<-pred_plan16_prob_rf_3control[,"X0"]
ROC_list_3control_ver2[[2]]<-pred_plan16_prob_bart_3control
ROC_list_3control_ver2[[3]]<-plan16_rf_test_score_3indep_combined[,4]
ROC_list_3control_ver2[[4]]<-plan16_bart_test_score_3indep_combined[,4]
ROC_list_3control_ver2[[5]]<-plan16_stad_test_score_3control
ROC_list_3control_ver2[[6]]<-plan16_mnsc_test_score_3control
ROC_list_3control_ver2[[7]]<-plan16_hot_test_score_3control
ROC_list_3control_ver2[[8]]<-plan16_optim_test_score_3control

ROC_3control_ver2<-ROC_data_list2(ROC_list_3control_ver2)
ROC_combined_3control_ver2<-do.call(rbind, ROC_3control_ver2)

ROC_list<-ROC_3control_ver2
auc_eer_3control_ver2<-auc_eer_opt.t_generation(ROC_list)


# TPR at fixed FPR at 5% ,10%, 15% in ROC curves
data<-ROC_combined_3control_ver2
P<-c(0.05, 0.1, 0.15)
TPR_3control<-get_TPR_fixed_FPR_combined(data,P)
TPR_3control

# ROC curves with 10,000 thresholds 
ROC<-ROC_combined_3control_ver2
custom_col<-brewer.pal(8, "Dark2")

graphlist[[1]]<-ggplot(aes(x = one_minus_Spec,  y = Sen, group = Method), data = ROC) +
  geom_line(aes(color = Method, linetype=Method), size = 1) +
  scale_color_manual(values = custom_col) +
  geom_abline(intercept = 0, slope = 1, color = "gray", size = 1) +
  theme_bw(base_size = 15)+xlab("1-Specificity (FPR)") +ylab("Sensitivity (TPR)")


graphlist[[2]]<-ggplot(ROC, aes(x = one_minus_Spec, y = Sen, shape = Method)) +
  geom_line(aes(color = Method, linetype=Method), size = 1) +
  scale_color_manual(values = custom_col) +
  geom_abline(intercept = 0, slope = 1, color = "gray", size = 1) +
  theme_bw(base_size = 15)+xlab("1-Specificity (FPR) < 0.3") +ylab("Sensitivity (TPR) > 0.9")+
  xlim(0.05, 0.3)+ylim(0.9,1)

#figure 8
ggarrange(graphlist[[1]],graphlist[[2]], ncol=1, nrow=2, common.legend = TRUE, legend="bottom")



error_FPR_3control<-c(error_FPR_3control_rf, error_FPR_3control_bart, error_FPR_3control_rf8, error_FPR_3control_bart8,
                      error_FPR_3control_stad, error_FPR_3control_mnsc, error_FPR_3control_hot, error_FPR_3control_optim)

error_FNR_3control<-c(error_FNR_3control_rf, error_FNR_3control_bart, error_FNR_3control_rf8, error_FNR_3control_bart8,
                      error_FNR_3control_stad, error_FNR_3control_mnsc, error_FNR_3control_hot, error_FNR_3control_optim)

############################################# SLR ############################################################
# Using training panes (train set), construct the distribution of KM scores and KNM scores by random forest
# and Standard - E2330 

# To create database using training panes, we also picked random three K's and one Q to generate mated and non-mated scores.

# For known match pairs, generate pairwise differences within each pane in training panes.
# For known non-match pairs, generate pairwise differences between panes of CompanyA and CompanyB,
# also between fragments, one from CompanyA or CompanyB and the other from 62 fragments from all different panes,
# and between 62 fragments from all different panes. (from Weis 2011 paper)

#Scores by random forest with down sampling on train set
# training panes
train_panes<-all_glass[-c(2,4,6,8,10,12,17,20,22,26,28,31,33,35,37,39,41,43,46,48,49)] 
length(train_panes) #28

# create index information of three random K's in each training pane
#list_train_ind_3control<-list()

for (i in 1:(length(train_panes))){
  originaldata<-train_panes[[i]]
  originaldata$fragment<-as.factor(originaldata$fragment)
  nfrag<-length(levels(originaldata$fragment))
  #level.frag<-as.numeric(levels(originaldata$fragment))
  sample_ind<-NULL
  
  for ( j in 1:nfrag){
    row<-1:nfrag  
    row1<-row[-j]
    Info<-NULL
    for (h in 1:10){
      control<-sample(row1, 3)
      Info<-rbind(Info, c(j, control))}
    sample_ind<-rbind(sample_ind, Info )
  }
  
  list_train_ind_3control[[i]]<-data.frame(sample_ind)
}


# Create difference of three Ks and Q as data frame for random forest prediction

rf_train_data_same_3control<-list()

for (i in 1:(length(train_panes))){
  originaldata<-train_panes[[i]]
  originaldata$fragment<-as.factor(originaldata$fragment)
  S<-list_train_ind_3control[[i]]
  rf_train_data_same_3control[[i]]<-rf_test_same_ncontrol(originaldata,S)
}
plan16_rf_train_data_same_3control<-do.call(rbind,rf_train_data_same_3control)

# Generate differences between Q and three Ks as data frame for getting prediction by RF 
## Q and three K's are from two different panes
D<-NULL
Data<-NULL
#plan16_rf_test_data_diff_3control<-NULL

for (j in 1:(length(train_panes)-1)){
  
  Data1<-train_panes[[j]]
  Data1$fragment<-as.factor(Data1$fragment)
  nfrag1<-length(levels(Data1$fragment))
  S1<-list_train_ind_3control[[j]]
  
  for ( k in (j+1): length(train_panes)){
    
    Data2<-train_panes[[k]]
    Data2$fragment<-as.factor(Data2$fragment)
    nfrag2<-length(levels(Data2$fragment))
    S2<-list_train_ind_3control[[k]]
    D<-rf_test_diff_ncontrol_ver2(Data1, Data2, nfrag1, S2)
    Data<-rbind(Data,D)}
  
}
#plan16_rf_train_data_diff_3control<-NULL
#plan16_rf_train_data_diff_3control<-Data


# difference with 62 fragments from Weis et al, 2011 (1-Q) and other panes (3-Ks)
D<-NULL
Data<-NULL
#plan16_rf_test_data_diff_3control<-NULL

Data1<-weis_diff
Data1$fragment<-as.factor(Data1$fragment)
nfrag1<-length(levels(Data1$fragment))

for ( k in 1: length(train_panes)){
  
  Data2<-train_panes[[k]]
  Data2$fragment<-as.factor(Data2$fragment)
  nfrag2<-length(levels(Data2$fragment))
  S2<-list_train_ind_3control[[k]]
  D<-rf_test_diff_ncontrol_ver2(Data1, Data2, nfrag1, S2)
  Data<-rbind(Data,D)
}

#plan16_rf_train_data_diff_3control_with_weis<-NULL
plan16_rf_train_data_diff_3control_with_weis<-Data
nrow(plan16_rf_train_data_diff_3control_with_weis) #1736

#Create difference of log transformed values bewteen Q and mean of K. 
plan16_rf_train_same_combined_3control<-rbind(plan16_rf_train_data_same_3control)
plan16_rf_train_diff_combined_3control<-rbind(plan16_rf_train_data_diff_3control, plan16_rf_train_data_diff_3control_with_weis)

### Get class probability from random forest with down sampling 
pred_plan16_train_rf_3control_same<- predict(plan14_down_fit, plan16_rf_train_same_combined_3control, type="prob")[,"X0"]
pred_plan16_train_rf_3control_diff<- predict(plan14_down_fit, plan16_rf_train_diff_combined_3control, type="prob")[,"X0"]

length(pred_plan16_train_rf_3control_same) #6710
length(pred_plan16_train_rf_3control_diff) #54020

# combine the predicted score
RF_score_3control<-c(pred_plan16_train_rf_3control_same, pred_plan16_train_rf_3control_diff)
class<-c(rep("X0", 6710), rep("X1", 54020))
train_score_rf_3control<-data.frame(RF_score_3control, class)

ggplot(train_score_rf_3control, aes(x=RF_score_3control, color=class)) +
  geom_histogram(aes(y = ..density.., color = class, fill = class),  alpha = 0.5, position = "identity") + 
  geom_density(aes(color = class), size =1)+
  geom_vline(xintercept=0.5) +
  theme(legend.position = "bottom")+xlab("Score by RF")+
  scale_color_manual(values=c("black", "dark gray"))+
  scale_fill_manual(values=c("black", "dark gray"))

#Test Standard method using training panes

#Mated pairs within train panes
#stad_same_train_3control<-list()
for (i in 1:(length(train_panes))){
  originaldata<-train_panes[[i]]
  originaldata$fragment<-as.factor(originaldata$fragment)
  #nfrag<-length(levels(originaldata$fragment))
  S<-list_train_ind_3control[[i]]
  stad_same_train_3control[[i]]<-stad_test_same_ncontrol(originaldata,S)
}
plan16_stad_same_train_score_3control<-do.call(rbind,stad_same_train_3control)

# Non-mated pairs between train panes
D<-NULL
Data<-NULL

for (j in 1:(length(train_panes)-1)){
  
  Data1<-train_panes[[j]]
  Data1$fragment<-as.factor(Data1$fragment)
  nfrag1<-length(levels(Data1$fragment))
  S1<-list_train_ind_3control[[j]]
  
  for ( k in (j+1): length(train_panes)){
    
    Data2<-train_panes[[k]]
    Data2$fragment<-as.factor(Data2$fragment)
    nfrag2<-length(levels(Data2$fragment))
    S2<-list_train_ind_3control[[k]]
    
    D<-stad_test_diff_ncontrol_ver2(Data1, Data2, nfrag1, S2)
    Data<-rbind(Data,D)}
}

#plan16_stad_diff_train_score_3control<-NULL
plan16_stad_diff_train_score_3control<-Data
nrow(plan16_stad_diff_train_score_3control) # 45340



## difference with 62 fragments from Weis et al, 2011 (1-Q) and other panes (3-Ks)
D<-NULL
Data<-NULL
#plan16_rf_test_data_diff_3control<-NULL

Data1<-weis_diff
Data1$fragment<-as.factor(Data1$fragment)
nfrag1<-length(levels(Data1$fragment))
#S1<-list_train_ind_3control[[j]]

for ( k in 1: length(train_panes)){
  
  Data2<-train_panes[[k]]
  Data2$fragment<-as.factor(Data2$fragment)
  nfrag2<-length(levels(Data2$fragment))
  S2<-list_train_ind_3control[[k]]
  D<-stad_test_diff_ncontrol_ver2(Data1, Data2, nfrag1, S2)
  Data<-rbind(Data,D)}

#plan16_stad_train_data_diff_3control_with_weis<-NULL
plan16_stad_train_data_diff_3control_with_weis<-Data
nrow(plan16_stad_train_data_diff_3control_with_weis) #1736

plan16_stad_train_score_3control<-c(plan16_stad_same_train_score_3control[,3],
                                    plan16_stad_diff_train_score_3control[,3],
                                    plan16_stad_train_data_diff_3control_with_weis[,3])

stad_score_3control<-c(plan16_stad_same_train_score_3control[,3],
                       plan16_stad_diff_train_score_3control[,3],
                       plan16_stad_train_data_diff_3control_with_weis[,3])

train_score_stad_3control<-data.frame(stad_score_3control, class)

ggplot(train_score_stad_3control, aes(x=stad_score_3control, color=class)) +
  geom_histogram(aes(y = ..density.., color = class, fill = class),  alpha = 0.5, position = "identity") + 
  geom_density(aes(color = class), size =1)+
  geom_vline(xintercept=4)+scale_x_continuous() +
  theme(legend.position = "bottom")+xlab("Score by Standard")+
  scale_color_manual(values=c("black", "dark gray"))+
  scale_fill_manual(values=c("black", "dark gray"))+
  scale_x_continuous(limits = c(0, 100),breaks = c(0, 4, 25, 50, 75, 100)) 

####### Score data frame
Score_rf_stad2<-data.frame(RF_score_3control, stad_score_3control, class )
names(Score_rf_stad)<-c("Score_RF", "Score_stad", "class")
Score_data<-Score_rf_stad2
Score_data$class <- revalue(Score_data$class, c("X0"="KM", "X1"="KNM"))

head(Score_data)

graphlist<-list()
#down sampling
graphlist[[1]]<-ggplot(Score_data, aes(x=RF_score_3control, color=class)) +
  geom_histogram(aes(y = ..density.., color = class, fill = class),  alpha = 0.5, position = "identity") + 
  geom_density(aes(color = class), size =1)+
  geom_vline(xintercept=0.5) +
  theme(legend.position = "bottom")+xlab("Score by RF")+
  scale_color_manual(values=c("black", "dark gray"))+
  scale_fill_manual(values=c("black", "dark gray"))

graphlist[[2]]<-ggplot(Score_data, aes(x=stad_score_3control, color=class)) +
  geom_histogram(aes(y = ..density.., color = class, fill = class),  alpha = 0.5, position = "identity") + 
  geom_density(aes(color = class), size =1)+
  geom_vline(xintercept=4)+scale_x_continuous() +
  theme(legend.position = "bottom")+xlab("Score by Standard")+
  scale_color_manual(values=c("black", "dark gray"))+
  scale_fill_manual(values=c("black", "dark gray"))+
  scale_x_continuous(limits = c(0, 100),breaks = c(0, 4, 25, 50, 75, 100)) 

ggarrange(graphlist[[1]],graphlist[[2]], ncol=1, nrow=2, common.legend = TRUE, legend="bottom")

# Summaries 
head(subset(down_train_score_rf_3control, class=="X1"))
sum(subset(down_train_score_rf_3control, class=="X1")$RF_score_3control>0.5)/nrow(subset(down_train_score_rf_3control, class="X1")) # 4.17%
sum(subset(down_train_score_rf_3control, class=="X1")$RF_score_3control>0.2)/nrow(subset(down_train_score_rf_3control, class="X1")) # 8.06%

head(subset(train_score_stad_3control, class=="X1"))
summary(subset(train_score_stad_3control, class=="X1")$stad_score_3control)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.494    4.522   30.530   61.982   39.490 4878.593 
sum(subset(train_score_stad_3control, class=="X1")$stad_score_3control>300)/nrow(subset(train_score_stad_3control, class=="X1")) # 4%
sum(subset(train_score_stad_3control, class=="X1")$stad_score_3control>500)/nrow(subset(train_score_stad_3control, class=="X1")) # 4%
sum(subset(train_score_stad_3control, class=="X1")$stad_score_3control<25)/nrow(subset(train_score_stad_3control, class=="X1")) # 45.77%



nrow(mean_log_glass[[33]])

# test score - TABLE 6 
frag_CB2_CB24<-mean_log_glass[[2]][2, -c(1,2)]-mean_log_glass[[2]][24, -c(1,2)] # AB-2 - AB-24    : SS 
frag_GB1_GB5<-mean_log_glass[[33]][1, -c(1,2)]-mean_log_glass[[33]][5, -c(1,2)] # BB-1 - BB-5  : SS 
frag_CB2_GB2<-mean_log_glass[[2]][2, -c(1,2)]-mean_log_glass[[33]][2, -c(1,2)]  # AB-2 - BB-2    : DS 
frag_GB2_CF2<-mean_log_glass[[2]][2, -c(1,2)]-mean_log_glass[[6]][2, -c(1,2)]  # AB-2 - AF-2     : DS 
frag_GB14_GD2<-mean_log_glass[[33]][14, -c(1,2)]-mean_log_glass[[35]][2, -c(1,2)] # BB-14 - BD-2    : DS 
score_test<-rbind(frag_CB2_CB24,frag_GB1_GB5,frag_CB2_GB2,frag_GB2_CF2,frag_GB14_GD2)


#score of test set by rf
score_test_by_rf <- predict(plan14_down_fit, score_test ,type="prob")[, "X0"]
score_test_stad<-all_glass[-c(2,4,6,8,10,12,17,20,22,26,28,31,33,35,37,39,41,43,46,48,49)] 

Data<-NULL
D<-NULL
Pane1<-c(2,33,2,2,33)
frag1<-c(2,1,2,2,14)
Pane2<-c(2,33,33,6,35)
frag2<-c(24,5,2,2,2)

for (j in 1:5){
  crit<-c(0,0)
  for (i in 3:20){
    X_r<-mean(all_glass[[Pane1[j]]][all_glass[[Pane1[j]]]$fragment==frag1[j],i])
    X_c<-mean(all_glass[[Pane2[j]]][all_glass[[Pane2[j]]]$fragment==frag2[j],i])
    sc1<-sqrt(var(all_glass[[Pane1[j]]][all_glass[[Pane1[j]]]$fragment==frag1[j],i]))
    sc2<-0.03*X_c
    sc<-max(sc1,sc2)
    crit[i]<-abs((X_r-X_c)/sc)
  }
  #score_stad<-max(crit)
  score_stad<-max(crit)
  D<-score_stad
  Data<-c(Data,D)
}
Data
score_test_by_stad<-Data

#score stad
Match<-c("SS", "SS", "DS", "DS", "DS")
data.frame(score_test_by_rf, score_test_by_stad, Match)


#density plot for rf scores
s_rf_KM2<-filter(Score_data, class=="KM")$RF_score_3control
s_rf_KNM2<-filter(Score_data, class=="KNM")$RF_score_3control
SLR_rf2<-SLR_by_density(s_rf_KM2, s_rf_KNM2, score_test_by_rf)

#stad
s_stad_KM2<-filter(Score_data, class=="KM")$stad_score_3control
s_stad_KNM2<-filter(Score_data, class=="KNM")$stad_score_3control

SLR_stad2<-SLR_by_density(s_stad_KM2, s_stad_KNM2, score_test_by_stad)
SLR_test_result2<-data.frame(score_test_by_rf, SLR_rf2, score_test_by_stad, SLR_stad2, Match)
SLR_test_result2

#score_test_by_rf      SLR_rf2 score_test_by_stad    SLR_stad2 Match
#            0.986 1.758091e+02           1.595799 1.439588e+01    SS
#            0.946 3.233160e+01           3.201886 2.358419e+00    SS
#            0.000 4.296193e-14          35.934699 3.588261e-11    DS
#            0.160 9.138466e-03           5.229163 4.714212e-01    DS
#            0.556 4.178134e-01           3.784561 1.424079e+00    DS



# simulation 

frag_CB_24<-filter(all_glass[[2]], fragment==24) #AB
frag_CAR_24<-filter(all_glass[[31]], fragment==24) #AAR
frag_CL_14<-filter(all_glass[[12]], fragment==14) #AB
frag_CL_24<-filter(all_glass[[12]], fragment==24) #AL
frag_GJ_14<-filter(all_glass[[41]], fragment==14) #AB
frag_GJ_24<-filter(all_glass[[41]], fragment==24) #BB
frag_CB_14<-filter(all_glass[[2]], fragment==14) #AB
frag_CAR_14<-filter(all_glass[[31]], fragment==14) #AAR
frag_CAR_24<-filter(all_glass[[31]], fragment==24) #AAR
frag_weis_104G<-filter(all_glass[[49]], fragment=="104G")


# simulation by multivariate log normal distribution
# performance test of RF and modified n sigma criterion 

FRSD2<-FRSD[-c(1,2)]

simulation_mvn_mnsc_rf(frag_CB_14, frag_CB_14, "SP")

simulation_mvn_mnsc_rf(frag_CB_24, frag_CAR_14, "DP")

simulation_mvn_mnsc_rf(frag_CB_14, frag_CB_24, "SP")





################## Out of bag test of Standard and Modified n sigma criteria using 6, 9, 12 controls ##########################

################## Standard E2330 

#  6 controls (K1, K2, K3, ..., K6) are randomly selected
#list_sample_ind_6control<-list()
for (i in 1:(length(test_panes))){
  originaldata<-test_panes[[i]]
  originaldata$fragment<-as.factor(originaldata$fragment)
  nfrag<-length(levels(originaldata$fragment))
  #level.frag<-as.numeric(levels(originaldata$fragment))
  sample_ind<-NULL
  
  for ( j in 1:nfrag){
    row<-1:nfrag  
    row1<-row[-j]
    Info<-NULL
    for (h in 1:30){
      control<-sample(row1, 6)
      Info<-rbind(Info, c(j, control))}
    
    sample_ind<-rbind(sample_ind, Info )
  }
  
  list_sample_ind_6control[[i]]<-data.frame(sample_ind)
  
}

error_FPR_6control<-NULL
error_FNR_6control<-NULL

### Mated pairs 
#stad_same_score_6control<-list()
for (i in 1:(length(test_panes))){
  originaldata<-test_panes[[i]]
  originaldata$fragment<-as.factor(originaldata$fragment)
  #nfrag<-length(levels(originaldata$fragment))
  S<-list_sample_ind_6control[[i]]
  stad_same_score_6control[[i]]<-stad_test_same_ncontrol(originaldata,S)
}

plan16_stad_same_test_score_6control<-do.call(rbind,stad_same_score_6control)

### Non-Mated pairs
D<-NULL
Data<-NULL
for (j in 1:(length(test_panes)-1)){
  
  Data1<-test_panes[[j]]
  Data1$fragment<-as.factor(Data1$fragment)
  nfrag1<-length(levels(Data1$fragment))
  S1<-list_sample_ind_6control[[j]]
  
  for ( k in (j+1): length(test_panes)){
    
    Data2<-test_panes[[k]]
    Data2$fragment<-as.factor(Data2$fragment)
    nfrag2<-length(levels(Data2$fragment))
    S2<-list_sample_ind_6control[[k]]
    D<-stad_test_diff_ncontrol(Data1, Data2, nfrag1, S2)
    Data<-rbind(Data,D)}
}
#plan16_stad_diff_test_score_6control<-NULL
#plan16_stad_diff_test_score_6control<-Data
plan16_stad_test_score_6control<-c(plan16_stad_same_test_score_6control[,3],plan16_stad_diff_test_score_6control[,3])

oob_err_stad_6control<-stad_err(plan16_stad_same_test_score_6control,plan16_stad_diff_test_score_6control)
#Sen one_minus_spec      Spec one_minus_sen overall_err
#0.9824183      0.1948421 0.8051579     0.0175817    0.178441

error_FPR_6control_stad<-oob_err_stad_6control[,2]
error_FNR_6control_stad<-oob_err_stad_6control[,4]


#  9 controls (K1, K2, K3, ..., K9) are randomly selected 
#list_sample_ind_9control<-list()
for (i in 1:(length(test_panes))){
  originaldata<-test_panes[[i]]
  originaldata$fragment<-as.factor(originaldata$fragment)
  nfrag<-length(levels(originaldata$fragment))
  #level.frag<-as.numeric(levels(originaldata$fragment))
  sample_ind<-NULL
  for ( j in 1:nfrag){
    row<-1:nfrag  
    row1<-row[-j]
    Info<-NULL
    for (h in 1:30){
      control<-sample(row1, 9)
      Info<-rbind(Info, c(j, control))}
    
    sample_ind<-rbind(sample_ind, Info )
  }
  
  list_sample_ind_9control[[i]]<-data.frame(sample_ind)
}

### Mated pairs 
#stad_same_score_9control<-list()
for (i in 1:(length(test_panes))){
  originaldata<-test_panes[[i]]
  originaldata$fragment<-as.factor(originaldata$fragment)
  #nfrag<-length(levels(originaldata$fragment))
  S<-list_sample_ind_9control[[i]]
  stad_same_score_9control[[i]]<-stad_test_same_ncontrol(originaldata,S)
}

plan16_stad_same_test_score_9control<-do.call(rbind,stad_same_score_9control)

### Non-Mated pairs 
D<-NULL
Data<-NULL
for (j in 1:(length(test_panes)-1)){
  
  Data1<-test_panes[[j]]
  Data1$fragment<-as.factor(Data1$fragment)
  nfrag1<-length(levels(Data1$fragment))
  S1<-list_sample_ind_9control[[j]]
  
  for ( k in (j+1): length(test_panes)){
    
    Data2<-test_panes[[k]]
    Data2$fragment<-as.factor(Data2$fragment)
    nfrag2<-length(levels(Data2$fragment))
    S2<-list_sample_ind_9control[[k]]
    
    D<-stad_test_diff_ncontrol(Data1, Data2, nfrag1, S2)
    Data<-rbind(Data,D)}
  
}

#plan16_stad_diff_test_score_9control<-NULL
plan16_stad_diff_test_score_9control<-Data

plan16_stad_test_score_9control<-c(plan16_stad_same_test_score_9control[,3],plan16_stad_diff_test_score_9control[,3])
oob_err_stad_9control<-stad_err(plan16_stad_same_test_score_9control,plan16_stad_diff_test_score_9control)
#       Sen one_minus_spec      Spec one_minus_sen overall_err
# 0.993268      0.2016927 0.7983073   0.006732026   0.1836538

error_FPR_9control_stad<-oob_err_stad_9control[,2]
error_FNR_9control_stad<-oob_err_stad_9control[,4]


#  12 controls (K1, K2, K3, ..., K12) are randomly selected 
#list_sample_ind_12control<-list()
for (i in 1:(length(test_panes))){
  originaldata<-test_panes[[i]]
  originaldata$fragment<-as.factor(originaldata$fragment)
  nfrag<-length(levels(originaldata$fragment))
  #level.frag<-as.numeric(levels(originaldata$fragment))
  
  sample_ind<-NULL
  
  for ( j in 1:nfrag){
    row<-1:nfrag  
    row1<-row[-j]
    Info<-NULL
    for (h in 1:30){
      control<-sample(row1, 12)
      Info<-rbind(Info, c(j, control))}
    
    sample_ind<-rbind(sample_ind, Info )
  }
  list_sample_ind_12control[[i]]<-data.frame(sample_ind)
}

# Mated pairs

#stad_same_score_12control<-list()
for (i in 1:(length(test_panes))){
  originaldata<-test_panes[[i]]
  originaldata$fragment<-as.factor(originaldata$fragment)
  #nfrag<-length(levels(originaldata$fragment))
  S<-list_sample_ind_12control[[i]]
  stad_same_score_12control[[i]]<-stad_test_same_ncontrol(originaldata,S)
}

plan16_stad_same_test_score_12control<-do.call(rbind,stad_same_score_12control)

# Non-Mated pairs
D<-NULL
Data<-NULL

for (j in 1:(length(test_panes)-1)){
  
  Data1<-test_panes[[j]]
  Data1$fragment<-as.factor(Data1$fragment)
  nfrag1<-length(levels(Data1$fragment))
  S1<-list_sample_ind_12control[[j]]
  
  for ( k in (j+1): length(test_panes)){
    
    Data2<-test_panes[[k]]
    Data2$fragment<-as.factor(Data2$fragment)
    nfrag2<-length(levels(Data2$fragment))
    S2<-list_sample_ind_12control[[k]]
    
    D<-stad_test_diff_ncontrol(Data1, Data2, nfrag1, S2)
    Data<-rbind(Data,D)}
}

#plan16_stad_diff_test_score_12control<-NULL
plan16_stad_diff_test_score_12control<-Data

plan16_stad_test_score_12control<-c(plan16_stad_same_test_score_12control[,3],plan16_stad_diff_test_score_12control[,3])


oob_err_stad_12control<-stad_err(plan16_stad_same_test_score_12control,plan16_stad_diff_test_score_12control)
#       Sen one_minus_spec      Spec one_minus_sen overall_err
# 0.993268      0.2016927 0.7983073   0.006732026   0.1836538

error_FPR_12control_stad<-oob_err_stad_12control[,2]
error_FNR_12control_stad<-oob_err_stad_12control[,4]


# Modified 4 sigma criteria

#  6 controls (K1, K2, K3, ..., K6) are randomly selected 

# Mated pairs
#mnsc_same_score_6control<-list()
for (i in 1:(length(test_panes))){
  originaldata<-test_panes[[i]]
  originaldata$fragment<-as.factor(originaldata$fragment)
  #nfrag<-length(levels(originaldata$fragment))
  S<-list_sample_ind_6control[[i]]
  mnsc_same_score_6control[[i]]<-mnsc_test_same_ncontrol(originaldata,S)
}

plan16_mnsc_same_test_score_6control<-do.call(rbind,mnsc_same_score_6control)

# Non-Mated pairs
#plan16_mnsc_diff_test_score<-NULL
D<-NULL
Data<-NULL

for (j in 1:(length(test_panes)-1)){
  
  Data1<-test_panes[[j]]
  Data1$fragment<-as.factor(Data1$fragment)
  nfrag1<-length(levels(Data1$fragment))
  S1<-list_sample_ind_6control[[j]]
  
  for ( k in (j+1): length(test_panes)){
    
    Data2<-test_panes[[k]]
    Data2$fragment<-as.factor(Data2$fragment)
    nfrag2<-length(levels(Data2$fragment))
    S2<-list_sample_ind_6control[[k]]
    
    D<-mnsc_test_diff_ncontrol(Data1, Data2, nfrag1, S2)
    Data<-rbind(Data,D)}
  
}

#plan16_mnsc_diff_test_score_6control<-NULL
plan16_mnsc_diff_test_score_6control<-Data
plan16_mnsc_test_score_6control<-c(plan16_mnsc_same_test_score_6control[,3],plan16_mnsc_diff_test_score_6control[,3])

oob_err_mnsc_6control<-stad_err(plan16_mnsc_same_test_score_6control,plan16_mnsc_diff_test_score_6control)
#Sen one_minus_spec      Spec one_minus_sen overall_err
# 0.5697386     0.06457417 0.9354258     0.4302614  0.09840953

error_FPR_6control_mnsc<-oob_err_mnsc_6control[,2]
error_FNR_6control_mnsc<-oob_err_mnsc_6control[,4]

#  9 controls (K1, K2, K3, ..., K9) are randomly selected 
# Mated pairs
# mnsc_same_score_9control<-list()
for (i in 1:(length(test_panes))){
  originaldata<-test_panes[[i]]
  originaldata$fragment<-as.factor(originaldata$fragment)
  #nfrag<-length(levels(originaldata$fragment))
  S<-list_sample_ind_9control[[i]]
  mnsc_same_score_9control[[i]]<-mnsc_test_same_ncontrol(originaldata,S)
}
plan16_mnsc_same_test_score_9control<-do.call(rbind,mnsc_same_score_9control)


# Non-Mated pairs
plan16_mnsc_diff_test_score<-NULL
D<-NULL
Data<-NULL

for (j in 1:(length(test_panes)-1)){
  
  Data1<-test_panes[[j]]
  Data1$fragment<-as.factor(Data1$fragment)
  nfrag1<-length(levels(Data1$fragment))
  S1<-list_sample_ind_9control[[j]]
  
  for ( k in (j+1): length(test_panes)){
    
    Data2<-test_panes[[k]]
    Data2$fragment<-as.factor(Data2$fragment)
    nfrag2<-length(levels(Data2$fragment))
    S2<-list_sample_ind_9control[[k]]
    
    D<-mnsc_test_diff_ncontrol(Data1, Data2, nfrag1, S2)
    Data<-rbind(Data,D)}
  
}

#plan16_mnsc_diff_test_score_9control<-NULL
plan16_mnsc_diff_test_score_9control<-Data
plan16_mnsc_test_score_9control<-c(plan16_mnsc_same_test_score_9control[,3],plan16_mnsc_diff_test_score_9control[,3])

oob_err_mnsc_9control<-stad_err(plan16_mnsc_same_test_score_9control,plan16_mnsc_diff_test_score_9control)
#Sen      one_minus_spec      Spec    one_minus_sen overall_err
# 0.5815686     0.06616687 0.9338331     0.4184314  0.09876028

error_FPR_9control_mnsc<-oob_err_mnsc_9control[,2]
error_FNR_9control_mnsc<-oob_err_mnsc_9control[,4]

##  12 controls (K1, K2, K3, ..., K12) are randomly selected 

# Mated pairs
#mnsc_same_score_12control<-list()
for (i in 1:(length(test_panes))){
  originaldata<-test_panes[[i]]
  originaldata$fragment<-as.factor(originaldata$fragment)
  #nfrag<-length(levels(originaldata$fragment))
  S<-list_sample_ind_12control[[i]]
  mnsc_same_score_12control[[i]]<-mnsc_test_same_ncontrol(originaldata,S)
}

plan16_mnsc_same_test_score_12control<-do.call(rbind,mnsc_same_score_12control)

# Non-Mated pairs
#plan16_mnsc_diff_test_score<-NULL
D<-NULL
Data<-NULL

for (j in 1:(length(test_panes)-1)){
  
  Data1<-test_panes[[j]]
  Data1$fragment<-as.factor(Data1$fragment)
  nfrag1<-length(levels(Data1$fragment))
  S1<-list_sample_ind_12control[[j]]
  
  for ( k in (j+1): length(test_panes)){
    
    Data2<-test_panes[[k]]
    Data2$fragment<-as.factor(Data2$fragment)
    nfrag2<-length(levels(Data2$fragment))
    S2<-list_sample_ind_12control[[k]]
    
    D<-mnsc_test_diff_ncontrol(Data1, Data2, nfrag1, S2)
    Data<-rbind(Data,D)}
}

#plan16_mnsc_diff_test_score_12control<-NULL
plan16_mnsc_diff_test_score_12control<-Data

plan16_mnsc_test_score_12control<-c(plan16_mnsc_same_test_score_12control[,3],plan16_mnsc_diff_test_score_12control[,3])

oob_err_mnsc_12control<-stad_err(plan16_mnsc_same_test_score_12control,plan16_mnsc_diff_test_score_12control)
#Sen      one_minus_spec      Spec    one_minus_sen overall_err
# 0.5815686     0.06616687 0.9338331     0.4184314  0.09876028

error_FPR_12control_mnsc<-oob_err_mnsc_12control[,2]
error_FNR_12control_mnsc<-oob_err_mnsc_12control[,4]


#Classification result using 3, 6, 9, 12 controls

stad_FPR<-c(error_FPR_3control_stad, error_FPR_6control_stad,error_FPR_9control_stad,error_FPR_12control_stad)
stad_FNR<-c(error_FNR_3control_stad, error_FNR_6control_stad, error_FNR_9control_stad, error_FNR_12control_stad) 

mnsc_FPR<-c(error_FPR_3control_mnsc, error_FPR_6control_mnsc,error_FPR_9control_mnsc,error_FPR_12control_mnsc)
mnsc_FNR<-c(error_FNR_3control_mnsc, error_FNR_6control_mnsc, error_FNR_9control_mnsc, error_FNR_12control_mnsc) 


class_result<-data.frame(rbind(stad_FNR, stad_FPR, mnsc_FNR, mnsc_FPR))
names(class_result)<-c("3controls","6controls","9controls","12controls")
class_result

#         3controls  6controls   9controls  12controls
#stad_FNR 0.05588235 0.01758170 0.006732026 0.004248366
#stad_FPR 0.18664534 0.19484206 0.201692656 0.204338265
#mnsc_FNR 0.44823529 0.43026144 0.418431373 0.420261438
#mnsc_FPR 0.06275490 0.06457417 0.066166867 0.067413035





