library(haven)
library(ggplot2)

df=read_sas("kidpan_data.sas7bdat")

#Restrict to Kidney Organ
dat=df[df$WL_ORG=="KI",]

#Sort By CTR ID
dat$CTR_CD=factor(dat$LISTING_CTR_CODE,labels=1:317)
dat=dat[order(dat$CTR_CD),]
dat=dat[!is.na(dat$CTR_CD),]

#Restrict to Time Window
dat$END_DT=as.Date(dat$END_DATE,format="%m/%d/%y")
dat$LISTING_DT=as.Date(dat$INIT_DATE,format="%m/%d/%y")
dat$TX_DT=as.Date(dat$TX_DATE,format="%m/%d/%y")
dat=dat[dat$LISTING_DT <= as.Date("2016-08-31") & dat$LISTING_DT > as.Date("2011-08-31"),]

#Create Outcome Variable
dat$TIME=as.numeric(dat$END_DT-dat$LISTING_DT)
dat=dat[dat$TIME > 0 & !is.na(dat$TIME),]
dat$STATUS=as.numeric(dat$REM_CD %in% c(2,4,14,15,18,19) & dat$TIME < 5*365.25)

#Create Predictor Variables
dat$Blood_Type="A"
dat$Blood_Type[dat$ABO=="B"]="B"
dat$Blood_Type[dat$ABO %in% c("A1B","A2B","AB")]="AB"
dat$Blood_Type[dat$ABO=="O"]="O"
dat$Blood_Type[is.na(dat$ABO)]=NA

dat$BMI="Underweight"
dat$BMI[dat$BMI_TCR >= 18.5 & dat$BMI_TCR < 24.9]="Normal"
dat$BMI[dat$BMI_TCR >= 25.0 & dat$BMI_TCR < 29.9]="Overweight"
dat$BMI[dat$BMI_TCR >= 30]="Obese"
dat$BMI[is.na(dat$BMI_TCR)]=NA

dat$Diabetes=1
dat$Diabetes[dat$DIAB==1]=0
dat$Diabetes[dat$DIAB==998 | dat$DIAB==5]=NA
dat$Diabetes[is.na(dat$DIAB)]=NA

dat$CAN_PREV_TX=as.numeric(dat$PREV_KI_TX=="Y")

#Restrict to White and Black
dat=dat[dat$ETHCAT %in% c(1,2),]
f1=function(x){return(sum(x==1))}
f2=function(x){return(sum(x==2))}
size=sapply(split(dat$CTR_CD,dat$CTR_CD),length)
size1=sapply(split(dat$ETHCAT,dat$CTR_CD),f1)
size2=sapply(split(dat$ETHCAT,dat$CTR_CD),f2)
dat$size=rep(size,size)
dat$size1=rep(size1,size)
dat$size2=rep(size2,size)
dat=dat[dat$size1 > 25 & dat$size2 > 25,]

m=length(unique(dat$CTR_CD))
dat$CTR_CD=factor(dat$CTR_CD,labels=1:m)
for(i in 1:(m-1)) {
  dat[,paste("CTR",i,sep="")]=ifelse(dat$CTR_CD==i,1,ifelse(dat$CTR_CD==m,-1,0))
  dat[,paste("WHITE_","CTR",i,sep="")]=dat[,paste("CTR",i,sep="")]*(dat$ETHCAT==1)
  dat[,paste("BLACK_","CTR",i,sep="")]=dat[,paste("CTR",i,sep="")]*(dat$ETHCAT==2)
}

dat$RACE=as.numeric(dat$ETHCAT==1)

#Fit Model
vars=c("TIME","STATUS","CTR_CD","RACE",
       paste0("WHITE_","CTR",1:(m-1)),
       paste0("BLACK_","CTR",1:(m-1)),
       "INIT_AGE",
       "Blood_Type",
       "BMI",
       "Diabetes",
       "CAN_PREV_TX")
dat=dat[,vars]
dat=dat[complete.cases(dat),]
write.csv(dat,"//Users//nicholashartman//Downloads//OPTN//dat_full.csv")

library(lme4)
g=glmer(STATUS~(RACE|CTR_CD)+RACE+INIT_AGE+factor(Blood_Type)+factor(BMI)+Diabetes+CAN_PREV_TX,data=dat,family="binomial")
summary(g)

dat=read.csv("//Users//nicholashartman//Downloads//OPTN//dat.csv")
m=length(unique(dat$CTR_CD))
g=glm(as.formula(paste("STATUS~RACE+",
                       paste0("WHITE_","CTR",1:(m-1),collapse="+"),
                       "+",paste0("BLACK_","CTR",1:(m-1),collapse="+"),
                       "+","INIT_AGE",
                       "+","factor(Blood_Type)",
                       "+","factor(BMI)",
                       "+","Diabetes",
                       "+","CAN_PREV_TX")),data=dat,family="binomial")

#Calculate Measures
num_equality=g$coefficients[substr(names(g$coefficients),1,5)=="WHITE"]-g$coefficients[substr(names(g$coefficients),1,5)=="BLACK"]
num_equity=g$coefficients[2]+g$coefficients[substr(names(g$coefficients),1,5)=="WHITE"]-g$coefficients[substr(names(g$coefficients),1,5)=="BLACK"]
ind_white=which(substr(names(g$coefficients),1,5)=="WHITE")
ind_black=which(substr(names(g$coefficients),1,5)=="BLACK")
var_num_equality=diag(vcov(g))[ind_white]+diag(vcov(g))[ind_black]-2*diag(vcov(g)[ind_white,ind_black])
var_num_equity=diag(vcov(g))[2]+diag(vcov(g))[ind_white]+diag(vcov(g))[ind_black]+2*vcov(g)[2,ind_white]-2*vcov(g)[2,ind_black]-2*diag(vcov(g)[ind_white,ind_black])
Equality=num_equality/sqrt(var_num_equality)
Equity=num_equity/sqrt(var_num_equity)

setwd("//Users//nicholashartman//Downloads//OPTN")
set.seed(3)
library(ggplot2)
x=rnorm(50,1)
y=x+rnorm(50,1,0.25)
dat=data.frame(x,y)
dat$shape="1"
dat$shape[dat$x > 1 & dat$y > 1]="2"
dat$shape[dat$x < 1 & dat$y > 1]="3"
g=ggplot(dat,aes(x,y))+geom_point(aes(shape=shape,colour=shape),size=2)
g=g+scale_shape_manual(values=c("1"=15,"2"=16,"3"=17))
g=g+scale_colour_manual(values=c("1"="gray","2"="black","3"="gray40"))
g=g+theme_classic()
g=g+xlab("Standardized Equality Statistic")+ylab("Standardized Equity Statistic")
g=g+theme(text=element_text(size=18))
g=g+geom_vline(aes(xintercept=1),colour="red",linetype="dashed",size=1)+geom_hline(aes(yintercept=1),colour="red",linetype="dashed",size=1)
g=g+theme(legend.position="none")
g=g+annotate("text",1.85,0.25,label="Control Limits",colour="black",size=6,hjust=0)
g=g+annotate("text",2.25,1.75,label="Significantly\nDiscriminatory and Inequitable",colour="black",size=6,hjust=0.5)
g=g+annotate("text",-0.7,0.25,label="Not Significantly\nDiscriminatory or Inequitable",colour="black",size=6,hjust=0.5)
g=g+annotate("text",-0.7,1.75,label="Significantly\nInequitable Only",colour="black",size=6,hjust=0.5)
g=g+geom_segment(aes(x = 1.8, y = 0.25, xend = 1, yend = 0.25),
                 arrow = arrow(length = unit(0.5, "cm")))
g=g+geom_segment(aes(x = 2.25, y = 0.35, xend = 2.25, yend = 1),
                 arrow = arrow(length = unit(0.5, "cm")))
g=g+theme(axis.ticks=element_blank(),axis.text=element_blank())
g
ggsave(plot=g,filename="Disparity_Assessment_Chart.eps")

library(ggplot2)
dat=data.frame(Equality,Equity)
dat$shape="1"
dat$shape[dat$Equality > 1.96 & dat$Equity > 1.96]="2"
dat$shape[dat$Equality < 1.96 & dat$Equity > 1.96]="3"
g=ggplot(dat,aes(Equality,Equity))+geom_point(aes(shape=shape,colour=shape),size=2)
g=g+scale_shape_manual(values=c("1"=15,"2"=16,"3"=17))
g=g+scale_colour_manual(values=c("1"="gray","2"="black","3"="gray40"))
g=g+theme_classic()
g=g+xlab("Standardized Equality Statistic")+ylab("Standardized Equity Statistic")
g=g+theme(text=element_text(size=18))
g=g+geom_vline(aes(xintercept=1.96),colour="red",linetype="dashed",size=1)+geom_hline(aes(yintercept=1.96),colour="red",linetype="dashed",size=1)
g=g+theme(legend.position="none")
g
ggsave(plot=g,filename="Disparity_Assessment_Chart_Access.eps")

mean(Equality > 1.96 & Equity > 1.96)
mean(Equality < 1.96 & Equity > 1.96)
1-mean(Equality > 1.96 & Equity > 1.96)-mean(Equality < 1.96 & Equity > 1.96)