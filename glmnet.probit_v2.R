####Data#######################################################
library(maxLik)
library(ROCR)
library(glmnet)
library(mice)
library("pROC")
###Initialization###
setwd("C:/Users/xliu203/Desktop/new Vector")
fileName <- dir()
datalist <- vector("list",94)
j=1;
for(i in 1:94){
  if(!is.element(i, c(80:81,88,90:94))){
    datalist[[j]]=read.table(fileName[i],header=T)
    j=j+1
  }
}

patient=read.csv("C:/Users/xliu203/Desktop/Combine data.csv")
Acute=patient$Acute.GI.Toxicity
y=R=Acute
x1=patient$age.at.IMRT
x2=patient$Hormones
x3=patient$T.stage.greater.than.2a...1.
x4=patient$gleason.score
x5=patient$psa.prior.to.IMRT
x6=patient$Diabetes
x7=patient$Prostate.vol.cc.
x7=(x7-mean(x7))/sd(x7)
x8=patient$Statins.2
x9=patient$Neoadjuvant
x10=patient$Concurrent
x11=patient$Adjuvant
X=cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11)
colnames(X)=c("age","Hormones","T.stage","gleason.score",
              "psa","Diabetes","Prostate","Statins", 
              "Neoadjuvant","Concurrent","Adjuvant")
tempData=mice(X, m=5, maxit=50, meth='pmm', seed=500)
X=as.matrix(complete(tempData, 5))
X=X[-c(80:81,88,90:94),];
y=y[-c(80:81,88,90:94)]
################################################################################################

####Algorithm####
eps=0.001
lambda.max=100
lambda.min=eps*lambda.max
lambda=seq(from=log(lambda.max), to=log(lambda.min), length.out=100)
lambda=exp(lambda)
tau=0.192;
likelihood=function(X, y, beta){
  
  w1=tau/mean(y)
  w0=(1-tau)/(1-mean(y))
  w=rep(0,length(y))
  for(i in 1:length(y)){
    w[i]=w1*y[i]+w0*(1-y[i])   ###weight###
  }
  X1=X
  wllh=rep(0, nrow(X))
  for(i in 1:nrow(X)){
    if(pnorm(X1[i,]%*%beta)>1-1e-7){wllh[i]=w[i]*(1-y[i])*(-16.12)}
    else if(pnorm(X1[i,]%*%beta)<1e-7){wllh[i]=w[i]*y[i]*(-16.12)}
    else{wllh[i]=(y[i]*log(pnorm(X1[i,]%*%beta)/(1-pnorm(X1[i,]%*%beta)))+log(1-pnorm(X1[i,]%*%beta)))*w[i]}
  }
  
  return(sum(wllh))
}

glmnet_probit=function(X, y, penalty.factor,lambda){
  X=cbind(1,X)
  beta.initial=rep(0, ncol(X))
  beta0=beta.initial
  eps=0.001
  penalty.factor=c(0,penalty.factor)
  beta=matrix(0,ncol=length(lambda), nrow=ncol(X))
  deviance=rep(0,length(lambda))
  BIC=rep(0,length(lambda))
  AIC=rep(0,length(lambda))
  AICC=rep(0,length(lambda))
  w1=tau/mean(y)
  w0=(1-tau)/(1-mean(y))
  omega=rep(0,length(y))
  for(i in 1:length(y)){
    omega[i]=w1*y[i]+w0*(1-y[i])   ###weight###
  }
  for(i in 1:length(lambda)){
    s=lambda[i]
    beta.hat=beta.initial
    beta.final=c(1, rep(0.5, ncol(X)-1))
    while(sum(abs(beta.hat-beta.final))>eps){
      beta.final=beta.hat
      z=X%*%beta.hat+(y-pnorm(X%*%beta.hat))/dnorm(X%*%beta.hat)
      w=(dnorm(X%*%beta.hat)^2/pnorm(X%*%beta.hat)/(1-pnorm(X%*%beta.hat)))*omega
      for(j in 1:nrow(X)){
        if(dnorm(X[j,]%*%beta.hat)<=1e-8){
          z[j]=X[j,]%*%beta.hat+1e-5
          w[j]=1e-5*omega[j]
        }
      }
      
      beta.temp=rep(1,length(beta.final))
      while(sum(abs(beta.hat-beta.temp))>eps){
        beta.temp=beta.hat
        for(index in 1:ncol(X)){
          beta.lasso=sum(w*X[,index]*(z-X[,-index]%*%beta.hat[-index]))
          if(penalty.factor[index]==1){
            if(beta.lasso>0&&s<abs(beta.lasso)) beta.hat[index]=(beta.lasso-s)/sum(w*X[,index]^2)
            if(beta.lasso<0&&s<abs(beta.lasso)) beta.hat[index]=(beta.lasso+s)/sum(w*X[,index]^2)
            if(s>=abs(beta.lasso)) beta.hat[index]=0
            
          }
          if(penalty.factor[index]==0){
            beta.hat[index]=beta.lasso/sum(w*X[,index]^2)
          }
        }##end for
        
      }##end while
    }            
    beta.initial=beta.final
    beta[,i]=beta.final
    deviance[i]=-2*(likelihood(X,y,beta0)-likelihood(X,y, beta.final))
    AIC=-2*likelihood(X,y,beta.final)+sum(beta.final!=0)*2
    BIC[i]=-2*likelihood(X,y,beta.final)+sum(beta.final!=0)*log(nrow(X))
    AICC[i]=-2*likelihood(X,y,beta.final)+sum(beta.final!=0)*2+2*sum(beta.final!=0)*(sum(beta.final!=0)+1)/(nrow(X)-1-sum(beta.final!=0))
  }
  result=list(lambda=lambda, beta=beta, deviance=deviance, BIC=BIC, AICC=AICC)
  return(result)
}
object=glmnet_probit(X_sample,y_sample,c(0,rep(1,14)),lambda)

summary(glm(y~as.matrix(X), family=binomial(link="probit")))

cv_glmnet_probit=function(X,y,penalty.factor,lambda){
  folds=sample(1:5,size=nrow(X),replace=TRUE)
  cv=matrix(0,ncol=length(lambda),nrow=5)
  for(i in 1:5){
    Xc=X[folds!=i,];yc=y[folds!=i]
    beta=glmnet_probit(Xc,yc,penalty.factor,lambda)[[2]]
    for(j in 1:length(lambda)){
      cv[i,j]=-2*(likelihood(cbind(1,X[folds==i,]),y[folds==i], beta[,j]))
      
    }
  }
  cv.mean=apply(cv,2,median)
  result=list(cv=cv.mean, lambda=lambda, beta=beta)
  return(result)
}
object1=(cv_glmnet_probit(X,y,c(0,rep(1,11)),lambda))
plot(log(object1[[2]]), object1[[1]])

predict_glmnet_probit=function(X, y,newx, penalty.factor, lambda){
  
  object1=glmnet_probit(X,y,penalty.factor,lambda)
  object2=cv_glmnet_probit(X,y,penalty.factor,lambda)
  lambda=object2[[2]]
  cv=object2[[1]]
  k=which(cv==min(cv))[1]
  s=lambda[k]
  beta=object1[[2]][,k]
  predict=pnorm(cbind(1,newx)%*%beta)
  return(list(predict,cv[k]))
}

predict_glmnet_probit(X,y,newx,rep(1,11),lambda)

###################################################
################modified LKB#######################
nc=seq(0, 1, length.out=40)
geud=rep(0,86)
BIC=rep(0,length(nc))
AIC=rep(0,length(nc))
AICC=rep(0,length(nc))
Beta_AICC=matrix(0,ncol=50, nrow=13)
Beta_BIC=matrix(0,ncol=50, nrow=13)
Beta_AIC=matrix(0,ncol=50, nrow=13)

for(k in 1:length(nc)){
  
  n=nc[k]
  for(i in 1:86){
    d=datalist[[i]][[1]]
    geud[i]=(sum(d^(1/n))/length(d))^(n)  ##calculate gEUD
  }
  
  Xa=cbind(X,geud);
  pre=glmnet_probit(Xa, y,c(1,1,1,1,1,1,1,1,1,1,1,0), lambda)
  j1=which.min(pre$AICC)[1]
  AICC[k]=pre$AICC[j1]
  Beta_AICC[,k]=pre$beta[,j1]
  j2=which.min(pre$AIC)[1]
  AIC[k]=pre$AICC[j2]
  Beta_AIC[,k]=pre$beta[,j2]
  j3=which.min(pre$BIC)[1]
  BIC[k]=pre$BIC[j1]
  Beta_BIC[,k]=pre$beta[,j3]
}


k=which.min(AIC)
k=which.min(AICC)
k=which.min(BIC)
Beta_AIC[,k]
Beta_AICC[,k]
Beta_BIC[,k]

n=nc[k]
for(i in 1:86){
  d=datalist[[i]][[1]]
  geud[i]=(sum(d^(1/n))/length(d))^(n)  ##calculate gEUD
}
Xa=cbind(X,geud);
X1=Xa[,c(2,3,5,6,7,8,12)]
Data=data.frame(X1,y)
w1=tau/mean(y)
w0=(1-tau)/(1-mean(y))

predict=rep(0,86)
for(j in 1:86){
  Data1=Data[-j,]
  prob1=glm(y~., family=binomial(link="probit"), data=Data1)
  beta=coef(prob1)
  X2=as.matrix(cbind(1,Data1[,1:ncol(X1)]))
  eta=X2%*%beta
  w=dnorm(eta)^2/pnorm(eta)/(1-pnorm(eta))*(w0*pnorm(eta)+w1*(1-pnorm(eta)))
  Hes=matrix(0,nrow=ncol(X2), ncol=ncol(X2))
  for(i in 1:85){
    Hes=Hes+w[i]*X2[i,]%*%t(X2[i,])
  }
  Q=(X2)%*%solve(Hes)%*%t(X2)
  ksi=eta/2*diag(Q)
  b=solve(Hes)%*%t(X2)%*%diag(as.vector(w))%*%ksi ####bias
  beta.hat=beta-b
  predict[j]=pnorm(as.matrix(cbind(1,Data[j,1:ncol(X1)]))%*%beta.hat)
}


pred=prediction(predict, y)
perf1=performance(pred,"tpr","fpr")
auc.tmp=performance(pred,"auc")
auc1=as.numeric(auc.tmp@y.values)
auc1
ci.auc(y,predict)

prob=glm(y~., family=binomial(link="probit"), data=Data)
beta=coef(prob)
X2=as.matrix(cbind(1,Data[,1:4]))
eta=X2%*%beta
w=dnorm(eta)^2/pnorm(eta)/(1-pnorm(eta))*(w0*pnorm(eta)+w1*(1-pnorm(eta)))
Hes=matrix(0,nrow=5, ncol=5)
for(i in 1:85){
  Hes=Hes+w[i]*X2[i,]%*%t(X2[i,])
}
Q=(X2)%*%solve(Hes)%*%t(X2)
ksi=eta/2*diag(Q)
b=solve(Hes)%*%t(X2)%*%diag(as.vector(w))%*%ksi ####bias
beta.hat=beta-b

######################################################
##########Calculate the bias##########################
n=nc[7]
geud=rep(0,86)
for(i in 1:86){
  d=datalist[[i]][[1]]
  geud[i]=(sum(d^(1/n))/length(d))^(n)  ##calculate gEUD
}
Xa=cbind(X,geud);
object=cv_glmnet_probit(Xa, y ,c(1,1,1,1,1,1,1,1,1,1,1,0), lambda)
beta=object$beta[,which.min(object$cv)]
Data=data.frame(Xa[,c(1,4,5,7,8,12)],y)
prob1=glm(y~., family=binomial(link="probit"), data=Data)
beta=coef(prob1)
X2=as.matrix(cbind(1,Data[,1:6]))
eta=X2%*%beta
w=dnorm(eta)^2/pnorm(eta)/(1-pnorm(eta))*(w0*pnorm(eta)+w1*(1-pnorm(eta)))
Hes=matrix(0,nrow=6, ncol=6)
for(i in 1:79){
  Hes=Hes+w[i]*X1[i,]%*%t(X1[i,])
}
Q=(X1)%*%solve(Hes)%*%t(X1)
ksi=eta/2*diag(Q)
solve(Hes)%*%t(X1)%*%diag(as.vector(w))%*%ksi
##############################################################
#################LKB model####################################
pred_LKB=rep(0,86)
for(j in 1:86){
  dntcpl <- function(param){
    
    TD_50=param[1];
    m0=param[2];
    n=param[3];
    lntcp=0;
    if(TD_50>50&&TD_50<100&&m0>0&&m0<1&&n>0&&n<1){
      
      for(i in 1:86){
        if(i!=j){
          d=datalist[[i]][[1]]
          eud=(sum(d^(1/n))/length(d))^(n)
          t=(eud-TD_50)/(m0*TD_50)
          lntcp=lntcp+R[i]*log(pnorm(t))+(1-R[i])*log(1-pnorm(t))
        }
        
      }
    }else {lntcp=-100}
    return(lntcp)
  }
  model=maxLik(dntcpl,start=c(65,0.2,0.3),method="nm")
  summary(model)
  coef=coef(model)
  TD_50=coef[1]
  m0=coef[2]
  n=coef[3]
  d.pre=datalist[[j]][[1]]
  eud=(sum(d.pre^(1/n))/length(d.pre))^(n)
  t=(eud-TD_50)/(m0*TD_50)
  pred_LKB[j]=pnorm(t)
}
pred=prediction(pred_LKB, y)
perf2=performance(pred,"tpr","fpr")
auc.tmp=performance(pred,"auc")
auc2=as.numeric(auc.tmp@y.values)
auc2
ci.auc(y,pred_LKB)
###############################################################################
###############################################################################

est=function(nc, X, y){
  cvlist=list(100); temp=NULL; Xlist=list()
  for(l in 1:length(nc)){
    n=nc[l]
    geud=rep(0,94)
    for(i in 1:94){
      d=datalist[[i]][[1]]
      geud[i]=(sum(d^(1/n))/length(d))^(n)  ##calculate gEUD
    }
    Xa=cbind(X,geud);Xlist[[l]]=Xa
    cv.lasso=cv.glmnet(Xa,y,alpha=1,penalty.factor=c(1,1,1,1,1,1,1,1,1,1,1,0)
                       ,family=binomial(link="probit"), nfolds=10)
    if(coef(cv.lasso,s="lambda.min")[13]>=0){
      lambda=cv.lasso$lambda.min
      dev=cv.lasso$glmnet.fit$dev[which(cv.lasso$lambda==lambda)]
      temp=c(temp, dev)
      cvlist[[l]]=cv.lasso
    }else{
      cv.lasso=cv.glmnet(Xa,y,alpha=1,penalty.factor=c(1,1,1,1,1,1,1,1,1,1,1,0),
                         exclude=12,family=binomial, nfolds=10)
      lambda=cv.lasso$lambda.min
      dev=cv.lasso$glmnet.fit$dev[which(cv.lasso$lambda==lambda)]
      temp=c(temp, dev)
      cvlist[[l]]=cv.lasso
      
    }
    
  }
  k=which(temp==min(temp))
  Q=rep(0,79)
  for(i in 1:94){
    cv.lasso=cv.glmnet(Xlist[[k]][-i,],y[-i],alpha=1,penalty.factor=c(1,1,1,1,1,1,1,1,1,1,1,0)
                       ,family="binomial", nfolds=10)
    Q[i]=predict(cv.lasso, Xlist[[k]], s="lambda.min")[i]
  }
  pred=prediction(Q, y)
  perf5=performance(pred,"tpr","fpr")
  auc.tmp=performance(pred,"auc")
  auc=as.numeric(auc.tmp@y.values)
  return(c(nc[k], auc))
}
MLKB=NULL
for(i in 1:50){
  MLKB=rbind(MLKB,est(nc, X, R))
}
colnames(MLKB)=c("n","auc")
#############################################################################################
#############################statistical model###############################################
z0=matrix(0,86,17)
for(i in 1:86){
  for(j in 1:17){
    d=datalist[[i]][[1]]
    z0[i,j]=quantile(datalist[[i]][[1]],0.05*j+0.05)
  }
}
colnames(z0)=c("D90","D85","D80","D75","D70","D65","D60","D55","D50",
               "D45","D40","D35","D30","D25","D20","D15","D10")

Xc=cbind(X, z0)
pred_stat=rep(0,86)
for(i in 1:86){
  lasso1=glmnet(Xc[-i,], y[-i], family="binomial")
  AICC=(1-lasso1$dev.ratio)*lasso1$nulldev+2*colSums(lasso1$beta!= 0)+
    2*colSums(lasso1$beta!= 0)*(colSums(lasso1$beta!= 0)+1)/(nrow(X)-1-colSums(lasso1$beta!= 0))
  j=which.min(AICC)
  pred_stat[i]=predict(lasso1, Xc, s=lasso1$lambda[j], type="response")[i]
}
pred3=prediction(pred_stat, y)
perf3=performance(pred,"tpr","fpr")
auc.tmp3=performance(pred3,"auc")
auc3=as.numeric(auc.tmp3@y.values)
ci.auc(y,pred_stat)

lasso0=glmnet(Xc, y, family="binomial")
preee=predict(lasso0, Xc, s=lasso1$lambda[j], type="response")
AICC=(1-lasso0$dev.ratio)*lasso0$nulldev+2*colSums(lasso0$beta!= 0)+
  2*colSums(lasso0$beta!= 0)*(colSums(lasso0$beta!= 0)+1)/(nrow(X)-1-colSums(lasso0$beta!= 0))
j=which.min(AICC)
lasso0$beta[,j]

#######################################333
##########################################
roc.test(y, data.frame(predict, pred_stat),method="specificity", specificity=0.8)
#################################################################################################
###############################ROC curve#########################################################
plot(perf2,col="1")
par(new=TRUE)
plot(perf3,col="2")
par(new=TRUE)
plot(perf1,col="3")
legend(0.73,0.7,c("LKB","statistical model","mLKB"),
       lty=c(1,1,1),lwd=c(1.5,1.5,1.5), col=c(1,2,3))

##################################################################################################
######################Find the possible variables in statistical and modified model###############
z0=matrix(0,94,18)
for(i in 1:94){
  for(j in 1:18){
    d=datalist[[i]][[1]]
    z0[i,j]=quantile(datalist[[i]][[1]],0.05*j)
  }
}
colnames(z0)=c("D90","D85","D80","D75","D70","D65","D60","D55","D50",
               "D45","D40","D35","D30","D25","D20","D15","D10", "D5")

standard=apply(z0[1:79,1:18],2,mean)
error=rep(0,15)
for(i in 1:15){
  error[i]=sum(z0[79+i,1:18]-standard)^2/sum(standard^2)
}
error
which(error<=0.7)
