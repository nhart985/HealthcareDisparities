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

m=length(unique(dat$CTR_CD))
g=glm(as.formula(paste("STATUS~RACE+",
                       paste0("WHITE_","CTR",1:(m-1),collapse="+"),
                       "+",paste0("BLACK_","CTR",1:(m-1),collapse="+"),
                       "+","INIT_AGE",
                       "+","factor(Blood_Type)",
                       "+","factor(BMI)",
                       "+","Diabetes",
                       "+","CAN_PREV_TX")),data=dat,family="binomial")
summary(g)

#Calculate Measures
beta_equality=g$coefficients[c(1,2,357:length(g$coefficients))]
X_equality=model.matrix(g)[,c(1,2,357:length(g$coefficients))]
beta_equity=g$coefficients[c(357:length(g$coefficients))]
X_equity=model.matrix(g)[,c(357:length(g$coefficients))]
dat$Xbeta_equality=as.numeric(X_equality%*%beta_equality)
dat$Xbeta_equity=as.numeric(X_equity%*%beta_equity)
beta=g$coefficients
dat$Xbeta=as.numeric(model.matrix(g)%*%beta)

O_w=sapply(split(dat$STATUS[dat$RACE==1],
                 dat$CTR_CD[dat$RACE==1]),sum)
O_b=sapply(split(dat$STATUS[dat$RACE==0],
                 dat$CTR_CD[dat$RACE==0]),sum)

dat$exp_equality=plogis(dat$Xbeta_equality)
dat$exp_equity=plogis(dat$Xbeta_equity)
dat$n_temp_equality=plogis(dat$Xbeta_equality)*(1-plogis(dat$Xbeta_equality))
dat$n_temp_equity=plogis(dat$Xbeta_equity)*(1-plogis(dat$Xbeta_equity))
dat$o_var=plogis(dat$Xbeta)*(1-plogis(dat$Xbeta))

E_equality_w=sapply(split(dat$exp_equality[dat$RACE==1],
                          dat$CTR_CD[dat$RACE==1]),sum)
E_equity_w=sapply(split(dat$exp_equity[dat$RACE==1],
                        dat$CTR_CD[dat$RACE==1]),sum)
E_equality_b=sapply(split(dat$exp_equality[dat$RACE==0],
                          dat$CTR_CD[dat$RACE==0]),sum)
E_equity_b=sapply(split(dat$exp_equity[dat$RACE==0],
                        dat$CTR_CD[dat$RACE==0]),sum)
n_equality_w=sapply(split(dat$n_temp_equality[dat$RACE==1],
                          dat$CTR_CD[dat$RACE==1]),sum)
n_equity_w=sapply(split(dat$n_temp_equity[dat$RACE==1],
                        dat$CTR_CD[dat$RACE==1]),sum)
n_equality_b=sapply(split(dat$n_temp_equality[dat$RACE==0],
                          dat$CTR_CD[dat$RACE==0]),sum)
n_equity_b=sapply(split(dat$n_temp_equity[dat$RACE==0],
                        dat$CTR_CD[dat$RACE==0]),sum)
o_var_w=sapply(split(dat$o_var[dat$RACE==1],
                     dat$CTR_CD[dat$RACE==1]),sum)
o_var_b=sapply(split(dat$o_var[dat$RACE==0],
                     dat$CTR_CD[dat$RACE==0]),sum)

T_equality_w=(O_w-E_equality_w)/n_equality_w
T_equity_w=(O_w-E_equity_w)/n_equity_w
T_equality_b=(O_b-E_equality_b)/n_equality_b
T_equity_b=(O_b-E_equity_b)/n_equity_b
Den_Equality=sqrt(o_var_w/(n_equality_w^2)+o_var_b/(n_equality_b^2))
Den_Equity=sqrt(o_var_w/(n_equity_w^2)+o_var_b/(n_equity_b^2))
Equality=(T_equality_w-T_equality_b)/Den_Equality
Equity=(T_equity_w-T_equity_b)/Den_Equity

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
g=g+xlab("Standardized Equality Measure")+ylab("Standardized Equity Measure")
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
g=g+xlab("Standardized Equality Measure")+ylab("Standardized Equity Measure")
g=g+theme(text=element_text(size=18))
g=g+geom_vline(aes(xintercept=1.96),colour="red",linetype="dashed",size=1)+geom_hline(aes(yintercept=1.96),colour="red",linetype="dashed",size=1)
g=g+theme(legend.position="none")
g
ggsave(plot=g,filename="Disparity_Assessment_Chart_Access.eps")

mean(Equality > 1.96 & Equity > 1.96)
mean(Equality < 1.96 & Equity > 1.96)











