#libraryd packages####
library(ggplot2)
library(RColorBrewer)
library(plotrix)
library(plyr)
library(dplyr)
library(vegan)
library(ggExtra)
library(gridExtra)
library(cowplot)
library(tidyr)

#function for calculating abundance of species in the model#####
SIH<-function(species=9,dispersal=0.01, patches=30){
  N<- matrix(10,ncol=species,nrow=patches) # Community x Species abundance matrix
  R<-rep(10*(species/10),patches) #Initial resources in each community
  
  rInput<-150 #resource input
  rLoss<-10 #resource loss 
  eff<-0.2 #conversion efficiency
  mort<-0.2 #mortality
  Ext<- 0.1 #extinction Threshold
  
  ePeriod<-40000 #period of env sinusoidal fluctuations
  eAMP<-1 #amplitude of envrionment sinusoidal fluctuations
  
  Tmax<-140000 #number of time steps in Sim
  DT<- 0.08 # % size of discrete "time steps"

  eOptimum<-1-seq(0,eAMP, by=eAMP/(species-1)) #species environmental optima
  
  #dispersal conditions#
  dispersal_matrix<-matrix(1/(patches-1),patches,patches)
  diag(dispersal_matrix)<-0
  
  calc.immigration <- function(N,a,dispersal_matrix) dispersal_matrix%*%N*rep(a,each=patches)
  
  disp<-rep(dispersal,species)
  dispM<-matrix(rep(disp,each=patches),patches,species)
  
  Prod<-matrix(NA,species*patches,40000)
  Abund<-Prod
  
  for(TS in 1:Tmax){
    Immigrants<-calc.immigration(N,disp,dispersal_matrix) # calculates inmmigration
    envt.v<-0.5*eAMP*(sin((2*pi/ePeriod)*TS+1+(1:patches)*2*pi/patches)+1) #calculates current environment
    consume <- 0.1*(1.5-abs(sapply(eOptimum,'-',envt.v))) #calculates current comsumption
    Nt <- N*(1+DT*(eff*R*consume - dispM - mort)) + DT*Immigrants #abundance step
    Rt <- DT*rInput+R*(1-DT*(rLoss + rowSums(consume*N))) #resource step    
    N <- Nt * (Nt>Ext) # set to 0 if below extinction threshold
    R <- Rt
    
    if(TS>=100000){ #samples data after time step 100 000
      Prod[,(TS-100000)] <- c(t(eff*consume*R*N)) #Productivity
      Abund[,(TS-100000)] <- c(t(N)) #Abundance
    }
  } 
  
  Prod<-array(t(Prod),dim=c(40000,species,patches))
  Prod<-Prod[seq(1,40000,100),,]
  Abund<-array(t(Abund),dim=c(40000,species,patches))
  Abund<-Abund[seq(1,40000,100),,]
  return(list(Prod=Prod, Abund=Abund))
}

#function for generating functional traits####
trait_distribute<-function(species,functions, Prop.Cont){
  repeat{
    MTraits<-t(matrix(runif(species*functions),species,functions))
    MTraits[rbinom(species*functions, 1, Prop.Cont)==0]<-0 
    MTraits<-decostand(MTraits,"total",2)
    if(sum(colSums(MTraits))==species){break}
  }
  return(MTraits)
}

#functions for calculating multifunctionality at various thresholds####
FuncMaxfun<-function(adf,thresh=0.7,Local.max=NA,Reg.max=NA){
  l.max<-Local.max*thresh
  r.max<-Reg.max*thresh
  LfuncMaxed<-mean(rowMeans(apply(adf>=l.max,3,rowSums)))
  RfuncMaxed<-mean(rowSums(apply(adf,2,rowSums)>=r.max))
  ret <- data.frame(LfuncMaxed,RfuncMaxed)
  names(ret) <- c("LfuncMaxed","RfuncMaxed")
  ret$nFunc <- dim(adf)[2]
  ret
}

FuncMaxfuns<-function(adf,Local.max=NA,Reg.max=NA, threshmin=0.1,threshmax=2,threshstep=0.1){
  ret <- ddply(data.frame(thresholds = seq(threshmin, threshmax, threshstep)), .variables = .(thresholds), function(x) {
    FuncMaxfun(adf, thresh = x[1, 1],Local.max=Local.max,Reg.max=Reg.max)
  })
  ret
}

FuncMaxSDfun<-function(adf,thresh=0.7,Local.max=NA,Reg.max=NA){
  l.max<-Local.max*thresh
  r.max<-Reg.max*thresh
  LfuncMaxedCV<-mean(apply(apply(adf>=l.max,3,rowSums),2,sd))#/mean(rowMeans(apply(adf>=l.max,3,rowSums)))
  RfuncMaxedCV<-sd(rowSums(apply(adf,2,rowSums)>=r.max))#/mean(rowSums(apply(adf,2,rowSums)>=r.max))
  ret <- data.frame(LfuncMaxedCV,RfuncMaxedCV)
  names(ret) <- c("LfuncMaxedCV","RfuncMaxedCV")
  ret$nFunc <- dim(adf)[2]
  ret
}

FuncMaxSDfuns<-function(adf,Local.max=NA,Reg.max=NA, threshmin=0.1,threshmax=2,threshstep=0.1){
  ret <- ddply(data.frame(thresholds = seq(threshmin, threshmax, threshstep)), .variables = .(thresholds), function(x) {
    FuncMaxSDfun(adf, thresh = x[1, 1],Local.max=Local.max,Reg.max=Reg.max)
  })
  ret
}

####

#Metacommunity multifunctionality simulation####
runs<-50 #number of replicates
species<-9 #number of species
patches<-30 #number of patches
functions<-7 #number of functions
Prop.Cont<-c(0.1,0.25,0.5,0.75,1) #proportion of species that contribute to each species on average
DispV<-c(0.0001,0.0005,0.001,0.0015,0.005,0.01,0.05,0.1,0.5,1) #dispersal rates

#calculate abundances
SIH_data<-sapply(DispV,SIH,species=species,patches=patches)

#calculate maximum level of function that can be simultaneously produced for all species
Equal_Traits<-matrix(1/functions,functions,species)
Equal.max<-matrix(NA,length(DispV),functions)
Reg.max<-NA
for(i in 1:length(DispV)){
  LFunc_rate1<-apply(SIH_data[["Abund",i]],1,function(M){Equal_Traits%*%M})
  LFunc_rate<-array(t(LFunc_rate1),dim=c(nrow(SIH_data[["Abund",1]]),functions,patches))
  Equal.max[i,]<-apply(LFunc_rate,2,mean)
  Reg.max[i]<-mean(apply(LFunc_rate,2,rowSums))
}
Local.max<-max(Equal.max)
Reg.max<-max(Reg.max)

#generate thresholds for production
threshmin<-0.25
threshmax<-1.25
threshstep<-0.25
threshseq<-seq(threshmin,threshmax,by=threshstep)


#data storage
FuncMaxed.DF<-data.frame(run=rep(1:runs,each=(length(Prop.Cont)*length(DispV)*length(threshseq))),prop.cont=rep(Prop.Cont,each=(length(DispV)*length(threshseq))),dispersal=rep(DispV,each=length(threshseq)),thresholds=threshseq,LfuncMaxed=NA,RfuncMaxed=NA,nFunc=NA)
FuncMaxedSD.DF<-FuncMaxed.DF
Prop.patch.df<-data.frame(Proportion_patches=NA,Threshold=threshseq,Functions=rep(1:functions,each=length(threshseq)),Rep=rep(1:runs,each=functions*length(threshseq)),prop.cont=rep(Prop.Cont,each=runs*functions*length(threshseq)),Dispersal=rep(DispV,each=runs*length(Prop.Cont)*functions*length(threshseq)))
CV1<-data.frame(run=rep(1:runs,each=(length(DispV)*length(Prop.Cont))),prop.cont=rep(Prop.Cont, each=length(DispV)),dispersal=DispV,Local=NA,Regional=NA)

#simulate multifunctionality for number of reps
pb <- txtProgressBar(min = 0, max = runs, style = 3) #progress bar
for(r in 1:runs){
  for(p in 1:length(Prop.Cont)){
    MTraits<-trait_distribute(species=species,functions=functions,Prop.Cont=Prop.Cont[p]) #generate functional traits
    rownames(MTraits)<-1:functions 
    colnames(MTraits)<-1:species
    for(i in 1:length(DispV)){
      LFunc_rate1<-apply(SIH_data[["Abund",i]],1,function(M){MTraits%*%M}) #calculate function production based on abundances
      LFunc_rate<-array(t(LFunc_rate1),dim=c(nrow(SIH_data[["Abund",1]]),functions,30)) #reformat function into an array
      
      Prop.patch.df$Proportion_patches[Prop.patch.df$Rep == r & Prop.patch.df$prop.cont==Prop.Cont[p] & Prop.patch.df$Dispersal==DispV[i]]<-sapply(1:functions,function(x){sapply(threshseq,function(y){mean(colMeans(apply(LFunc_rate>Local.max*y,3,rowSums)>=x))})}) #calculate the proportion of patches that produce each level of multifunction
      #calculate the temporal variability of each individual function
      CV1$Local[CV1$run==r &  CV1$dispersal == DispV[i] & CV1$prop.cont == Prop.Cont[p]]<-mean(apply(LFunc_rate,3,function(x){apply(x,2,sd,na.rm=T)/apply(x,2,mean,na.rm=T)}),na.rm=T) 
      CV1$Regional[CV1$run==r &  CV1$dispersal == DispV[i] & CV1$prop.cont == Prop.Cont[p]]<-mean(apply(apply(LFunc_rate,2,rowSums),2,function(x){sd(x,na.rm=T)/mean(x,na.rm=T)}),na.rm=T)     
      
      #calculate multifunctionality
      FuncMaxed.DF[FuncMaxed.DF$run==r & FuncMaxed.DF$prop.cont==Prop.Cont[p]& FuncMaxed.DF$dispersal == DispV[i],-(1:3)]<-FuncMaxfuns(LFunc_rate,Local.max=Local.max,Reg.max=Reg.max,threshmin = threshmin,threshmax = threshmax,threshstep = threshstep)
    }
  }
  Sys.sleep(0.1)
  setTxtProgressBar(pb, r) #progress bar
}

options(scipen=999)

#Data processing#####
#calculate means and SD of multifunction
Reg_Multifunc_means<-aggregate(FuncMaxed.DF$RfuncMaxed,by=list(Prop.Cont=FuncMaxed.DF$prop.cont,Threshold=FuncMaxed.DF$threshold,Dispersal=FuncMaxed.DF$dispersal),mean)
Reg_Multifunc_means$sd<-aggregate(FuncMaxed.DF$RfuncMaxed,by=list(Prop.Cont=FuncMaxed.DF$prop.cont,Threshold=FuncMaxed.DF$threshold,Dispersal=FuncMaxed.DF$dispersal),sd)$x
Local_Multifunc_means<-aggregate(FuncMaxed.DF$LfuncMaxed,by=list(Prop.Cont=FuncMaxed.DF$prop.cont,Threshold=FuncMaxed.DF$threshold,Dispersal=FuncMaxed.DF$dispersal),mean)
Local_Multifunc_means$sd<-aggregate(FuncMaxed.DF$LfuncMaxed,by=list(Prop.Cont=FuncMaxed.DF$prop.cont,Threshold=FuncMaxed.DF$threshold,Dispersal=FuncMaxed.DF$dispersal),sd)$x

Thresh.df_All<-rbind(Reg_Multifunc_means,Local_Multifunc_means)
Thresh.df_All$Scale<-c(rep("Regional",nrow(Reg_Multifunc_means)),rep("Local",nrow(Reg_Multifunc_means)))
Thresh.df_All<-Thresh.df_All[order(Thresh.df_All$Prop.Cont,decreasing=T),]
Thresh.df_All$Threshold<-as.factor(Thresh.df_All$Threshold)

Thresh.df_sub<-Thresh.df_All[Thresh.df_All$Prop.Cont==0.1|Thresh.df_All$Prop.Cont==0.5|Thresh.df_All$Prop.Cont==1,]

#calculate means and SD of temporal variability of individual functions
CV1.df<-aggregate(CV1$Local,by = list(Dispersal=CV1$dispersal,Prop.Cont=CV1$prop.cont),mean, na.rm=T)
names(CV1.df)<-c("Dispersal","Prop.Cont","local")
CV1.df$local.sd<-aggregate(CV1$Local,by = list(Dispersal=CV1$dispersal,Prop.Cont=CV1$prop.cont),sd, na.rm=T)$x
CV1.df$regional<-aggregate(CV1$Regional,by = list(Dispersal=CV1$dispersal,Prop.Cont=CV1$prop.cont),mean, na.rm=T)$x
CV1.df$regional.sd<-aggregate(CV1$Regional,by = list(Dispersal=CV1$dispersal,Prop.Cont=CV1$prop.cont),sd, na.rm=T)$x

CV.df<-gather(CV1.df,key = Scale,value=CV,local,regional)[,-c(3:4)]
CV.df$sd<-gather(CV1.df,key = Scale,value=SD,local.sd,regional.sd)[,6]

#calculate species diversity and biomass
for(i in 1:length(DispV)){
  hold.df<-data.frame(Diversity=c(mean(exp(apply(SIH_data[["Abund",i]],3,diversity))),mean(exp(diversity(apply(SIH_data[["Abund",i]],1,rowSums),MARGIN = 2)))),
             Function1=mean(apply(SIH_data[["Abund",i]],1,colSums))/functions,
             Scale=c("Local","Regional"),Dispersal=DispV[i])
if(i==1){
  Diversity.df<-hold.df
} else{Diversity.df<-rbind(Diversity.df,hold.df)}
  }

#Figures####
#Figure 1####
p1<-ggplot(Diversity.df,aes(x=Dispersal,y=Diversity,group=Scale, linetype=Scale))+
  geom_line(size=1)+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"))+
  theme_bw(base_size = 15)+
  theme(legend.justification=c(0,0), legend.position=c(0,0.75))+
  removeGrid()+
  ylab("Effective species diversity")

p2<-ggplot(Diversity.df,aes(x=Dispersal,y=Function1,group=Scale, linetype=Scale))+
  geom_line(size=1)+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"))+
  theme_bw(base_size = 15)+
  theme(legend.position="none")+
  removeGrid()+
  ylab("Local community biomass")

p3<-ggplot(CV.df,aes(x=Dispersal,y=CV,group=interaction(Scale,Prop.Cont), linetype=Scale, color=factor(Prop.Cont,ordered = T)))+
  geom_line(size=1)+
  geom_errorbar(aes(ymax=CV+sd, ymin=CV-sd))+
  scale_color_manual(values = grey(seq(0.8,0,length=5)),name="Functional\noverlap")+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"))+
  theme_bw(base_size = 15)+
  guides(linetype=F)+
  removeGrid()+
  theme(legend.justification=c(0,0), legend.position=c(0,0.5))+
  ylab("Individual function CV")

pdf(file = "Figure 1.pdf",height = 4.5, width=15.5)
plot_grid(p1,p2,p3, labels=c("(a)", "(b)","(c)"), ncol = 3, nrow = 1)
dev.off()



#Figure 2####
ggplot(Diversity.df,aes(x=Diversity,y=Function1,group=Scale, linetype=Scale))+
  geom_path(aes(linetype=Scale))+
  scale_linetype_manual(values=c(2,1))+
  geom_point(aes(shape=Scale,color=as.factor(Dispersal)))+
  scale_color_manual(values = rev(brewer.pal(10,"RdYlBu")),name="Dispersal")+
  scale_shape_manual(values=c(17,19))
  theme_bw(base_size = 15)+
  removeGrid()+
  ylab("Effective species diversity")

par(mfrow=c(1,2),las=1, pty='s',mar=c(5,5,2,2))
plot((RFun1*functions)~R.Div, type='b', pch="", xlab="Effective species diversity", ylab="Local community biomass", lwd=2,cex.lab=cex.labV)
points(RFun1*functions~R.Div,col=rev(brewer.pal(length(DispV),name = "RdYlBu")),pch=19)
points(RFun1*functions~L.Div, pch="", col=1, lwd=2, lty=2, type='b')
points(RFun1*functions~L.Div,col=rev(brewer.pal(length(DispV),name = "RdYlBu")),pch=17)
mtext("(a)",adj=-0.22, cex=1.2)

plot(Reg_Multifunc_means$x[Reg_Multifunc_means$Prop.Cont==0.5 & Reg_Multifunc_means$Threshold==0.5]~R.Div, type='b', pch="", xlab="Effective species diversity", ylab="Multifunctionality", lwd=2,cex.lab=cex.labV, ylim=c(2.5,6.1))
points(Reg_Multifunc_means$x[Reg_Multifunc_means$Prop.Cont==0.5 & Reg_Multifunc_means$Threshold==0.5]~R.Div,col=rev(brewer.pal(length(DispV),name = "RdYlBu")),pch=19)
points(Local_Multifunc_means$x[Local_Multifunc_means$Prop.Cont==0.5 & Local_Multifunc_means$Threshold==0.5]~L.Div, pch="", col=1, lwd=2, lty=2, type='b')
points(Local_Multifunc_means$x[Local_Multifunc_means$Prop.Cont==0.5 & Local_Multifunc_means$Threshold==0.5]~L.Div,col=rev(brewer.pal(length(DispV),name = "RdYlBu")),pch=17)
legend("topleft",legend=c("Local", "Regional"), lwd=2, pch=c(17,19), lty=c(2,1), bty='n', title="Scale",cex=0.8)
legend("bottomright", legend=DispV,col=rev(brewer.pal(length(DispV),name = "RdYlBu")), pch=19, bty='n', title= "Dispersal", cex=0.8)
mtext("(b)",adj=-0.22, cex=1.2)

#Figure 3####
ggplot(summarize(group_by(Prop.patch.df,prop.cont,Dispersal,Functions,Threshold),Mean=mean(Proportion_patches),SD=sd(Proportion_patches)),aes(x=Functions,y=Mean,group=as.factor(Dispersal),color=as.factor(Dispersal)))+
  geom_errorbar(aes(ymin=Mean-SD,ymax=Mean+SD),width=0)+
  scale_color_manual(values = rev(brewer.pal(10,"RdYlBu")),name="Dispersal")+
  geom_point(size=2)+
  geom_line(size=1.5)+
  ylab("Proportion of patches")+
  facet_grid(Threshold~prop.cont)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(breaks=seq(0,1,by=0.25))

#Figure 4####
ggplot(data=Thresh.df_sub,aes(y=x,x = Dispersal,group=Threshold, color=Threshold))+
  geom_errorbar(aes(ymax=x+sd,ymin=x-sd, size=1.5),width=0.1)+
  geom_point()+
  geom_line(size=1.5)+
  scale_size(range=c(0.2, 2), guide=FALSE)+
  facet_grid(Prop.Cont~Scale)+
  scale_color_manual(values = brewer.pal(n = length(threshseq)+1,name = "Blues")[-1])+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"))+
  #scale_y_continuous(limits=c(0,7))+
  ylab("Multifunctionality") +
  xlab("Dispersal")+
  theme_bw(base_size = 16)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(vjust=-0.5),axis.title.y=element_text(vjust=1.5))

#Supplementary figures
Thresh.df_Both<-Reg_Multifunc_means
Thresh.df_Both$Local<-Local_Multifunc_means$x
Thresh.df_Both$Local_SD<-Local_Multifunc_means$SD
Thresh.df_Both$Prop.Cont<-as.factor(Thresh.df_Both$Prop.Cont)
Thresh.df_Both$Dispersal<-as.factor(Thresh.df_Both$Dispersal)
Thresh.df_Both<-Thresh.df_Both[Thresh.df_Both$Dispersal!=0.0015 & Thresh.df_Both$Dispersal!=1,]

Thresh.df_Both<-Thresh.df_Both[order(Thresh.df_Both$Prop.Cont,decreasing=T),]
Thresh.df_Both<-Thresh.df_Both[Thresh.df_Both$Dispersal!=0,]

ColV<-rev(brewer.pal(8,name = "RdYlBu"))

ggplot(data=Thresh.df_Both,aes(y=Local,x = x, fill=Dispersal))+
  geom_point(aes(size=Prop.Cont),pch=21)+
  facet_wrap(~Threshold,nrow = 2)+
  scale_fill_manual(values = ColV[])+
  scale_size_discrete(range = c(1.8,6.5),name="Functional\noverlap")+
  xlab("Regional multifunction") +
  ylab("Local multifunction")+
  theme_bw(base_size = 16)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title.x=element_text(vjust=-0.5),axis.title.y=element_text(vjust=1.8),legend.title.align=0.5, strip.text.x=element_text(face="bold"))+
  #theme(panel.background = element_rect(fill = 'grey60'))+
  geom_abline(intercept=0,slope=1,linetype=2,color="black")+
  guides(fill=guide_legend(override.aes = list(size=5)))+
  theme(legend.key = element_blank())+
  ylim(0,7)+
  xlim(0,7)+
  coord_fixed(ratio=1)

Thresh.df_All$Div<-rep(L.Div,each=length(threshseq))
Thresh.df_All$Div[Thresh.df_All$Scale=="Regional"]<-rep(rep(R.Div,each=length(threshseq)),length(Prop.Cont))

ggplot(data=Thresh.df_All,aes(x=Div,y=x, group=Scale))+
  coord_fixed()+
  geom_path(aes(linetype=Scale))+
  scale_linetype_manual(values=c(2,1))+
  geom_point(aes(shape=Scale,color=as.factor(Dispersal)))+
  scale_color_manual(values = rev(brewer.pal(10,"RdYlBu")),name="Dispersal")+
  scale_shape_manual(values=c(17,19))+
  facet_grid(Threshold~Prop.Cont)+
  #scale_fill_manual(values = c("#377eb8","#e41a1c"))+
  xlab("Effective species diversity") +
  ylab("Multifunctionality")+
  theme_bw(base_size = 16)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title.x=element_text(vjust=-0.5),axis.title.y=element_text(vjust=1.8),legend.title.align=0.5, strip.text.x=element_text(face="bold"))+
  #theme(panel.background = element_rect(fill = 'grey60'))+
  theme(legend.key = element_blank())+
  scale_x_continuous(breaks=c(0,3,6,9))+
  coord_cartesian()+
  theme(aspect.ratio=1)