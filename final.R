
library(glmnet)
library(MASS)

#data preperation

###orthogonal data
set.seed(100)
n=500

p=15
X=scale(matrix(rnorm(n*p),nrow = n))
b=c(1,1.5,2,-1,2.5,rep(0,p-5))
y=X %*% b+rnorm(n,0,9) 
sdata=data.frame(y,X)


##for lasso and adaptive lasso
ols_coef=coef(lm(sdata$y~.,data = sdata))   #ols estimate
gamma=2
weight=1/(abs(ols_coef[-1]))^gamma

#residualbootstrap
#for residual
cvfit=cv.glmnet(X,y,alpha=1,penalty.factor=weight)
adlassomodel=glmnet(X,y,alpha=1,family="gaussian",penalty.factor=weight)
adlassoest=coef(adlassomodel,s=cvfit$lambda.min) # estimates
res=y-predict(cvfit,s=cvfit$lambda.1se,newx=X) #residual
res=res-mean(res)


#for bootstrap
indicator=rep(0,p)
indicator[b!=0]=1
B=500
beta_B=matrix(0,nrow=p+1,ncol = B)
count=0
for(i in 1:B){
  s=sample(1:nrow(X),nrow(X),replace=TRUE)
  bootres=res[s]
  y_new=X%*%adlassoest[-1]+bootres
  data_B=data.frame(y_new,X)
  ols_coef=coef(lm(data_B$y_new~.,data = data_B))   #ols estimate
  weight=1/(abs(ols_coef[-1]))^gamma
  cvfit=cv.glmnet(X,y_new,alpha=1,penalty.factor=weight)
  lambda.min=cvfit$lambda.min
  beta_B[,i]=as.numeric(coef(cvfit,s="lambda.min"))
  indicator_B=rep(0,p)
  indicator_B[beta_B[-1,i]!=0]=1
  if(setequal(indicator_B,indicator)==T){
    count=count+1
  }
}





se=apply(beta_B,1, sd)

conf=apply(beta_B,1, quantile,probs=c(0.025,0.975))

lasso_equal0=numeric(p+1)
for(i in 1:p+1){
  lasso_equal0[i]=sum(beta_B[i,]==0)/B
}

adalasso_equal0=numeric(p+1)
for(i in 1:p+1){
  adalasso_equal0[i]=sum(beta_B[i,]==0)/B
}


#p_value=numeric(p+1)
#for(i in 1:p+1){
#  p_value[i]=min(sum(beta_B[i,]<lassoest[i]),sum(beta_B[i,]>lassoest[i]))/B
#}
boxplot(t(beta_B[2:16,]),col="cyan",xlab="estimated coefficients",main="lasso")
abline(h=0,col="black")
y0s=c(b[1:15])
x0s=1:15-0.4
x1s=1:15+0.4
segments(x0=x0s,x1=x1s,y0=y0s,col="red",lwd=5)
segments(x0=x0s+0.2,x1=x1s-0.2,y0=conf[1,2:16],col="green",lwd=3,lty=3)
segments(x0=x0s+0.2,x1=x1s-0.2,y0=conf[2,2:16],col="green",lwd=3,lty=3)

#legend("topright",legend = c("line passing through 0","actual parameters given\n in simulation","95% confidence interval"),
#       col=c("black","red","green"),lty = c(1,1,3),lwd=c(1,5,3))

cbind(1:15,b,as.numeric(adlassoest)[-1],apply(beta_B, 1, mean)[-1],se[-1],t(conf[,-1]))




##modified bootstrap
set.seed(100)
n=500

p=15
X=scale(matrix(rnorm(n*p),nrow = n))
b=c(1,1.5,2,-1,2.5,rep(0,p-5))
y=X %*% b+rnorm(n,0,9) 
sdata=data.frame(y,X)
#lasso fit
cvfit=cv.glmnet(X,y,alpha=1)
lassomodel=glmnet(X,y,alpha=1,family="gaussian")
pred=predict(lassomodel,s=cvfit$lambda.min,newx=X)
est=coef(cvfit,s=cvfit$lambda.min)[,1]
est=est[-1]
#10 fold cross validation selection of delta
library(caret)
fold=createFolds(y,k=10)
sdel=function(delta){
  rss=numeric(10)
  for(i in 1:10){
    y_test=y[fold[[i]]]
    X_test=X[fold[[i]],]
    y_train=y[-fold[[i]]]
    X_train=X[-fold[[i]],]
    model=cv.glmnet(X_train,y_train,alpha=1)
    lmodel=glmnet(X_train,y_train,alpha=1,family="gaussian")
    lest=coef(model,s=model$lambda.min)[,1]
    lest=lest[-1]
    lseq=1/(n^delta)
    mest=numeric()
    for(j in 1:p){
      if(lest[j]<lseq && lest[j]>-lseq){
        mest[j]=0
      }else{
        mest[j]=lest[j]
      }
    }
    lpred=X_test%*%mest
    rss[i]=mean((y_test-lpred)^2)
  }
  return(mean(rss))
}
del=seq(0.01,0.49,by=0.01)
delta_mse=numeric(length(del))
for(i in 1:length(del)){
  delta_mse[i]=sdel(del[i])
}
par(mfrow=c(1,1))
plot(del,delta_mse,type="l",xlab = "delta",ylab = "CV-mse average")
for(i in 1:length(del)){
  if(delta_mse[i]==min(delta_mse)){
    del_min=del[i]
    break
  }
}
del_min
abline(v=del_min)
cbind(del,delta_mse)



#Estimate Modification
seq=1/(n^del_min)
est1=numeric()
for(i in 1:p){
  if(est[i]<seq && est[i]>-seq){
    est1[i]=0
  }else{
    est1[i]=est[i]
  }
}

#modified residual bootstrap
indicator=rep(0,p)
indicator[b!=0]=1
count=0
B=500
beta_B=matrix(0,nrow=p+1,ncol = B)
modres=y-X%*%est1
modres=modres-mean(modres)
for(i in 1:B){
  s=sample(1:n,n,replace = T)
  bootres=modres[s]
  y_B=X[s,] %*% est1+bootres
  X_B=X[s,]
  cvfit=cv.glmnet(X_B,y_B,alpha=1)
  lambda.min=cvfit$lambda.min
  
  beta_B[,i]=as.numeric(coef(cvfit,s="lambda.min"))
  indicator_B=rep(0,p)
  indicator_B[beta_B[-1,i]!=0]=1
  if(setequal(indicator_B,indicator)==T){
    count=count+1
  }
}
prob_B=count/B

prob_B




se=apply(beta_B,1, sd)
conf=apply(beta_B,1, quantile,probs=c(0.025,0.975))
bootest=apply(beta_B,1,mean)

boxplot(t(beta_B[2:16,]),col="cyan",xlab="estimated coefficients",main="modified bootstrap")
abline(h=0,col="black")
y0s=c(b[1:15])
x0s=1:15-0.4
x1s=1:15+0.4
segments(x0=x0s,x1=x1s,y0=y0s,col="red",lwd=5)
segments(x0=x0s+0.2,x1=x1s-0.2,y0=conf[1,2:16],col="green",lwd=3,lty=3)
segments(x0=x0s+0.2,x1=x1s-0.2,y0=conf[2,2:16],col="green",lwd=3,lty=3)
abline(h=0)
#legend("topright",legend = c("line passing through 0","actual parameters given\n in simulation","95% confidence interval"),
#       col=c("black","red","green"),lty = c(1,1,3),lty = c(1,1,3),lwd=c(1,5,3))



mod_equal0=numeric(p+1)
for(i in 1:p+1){
  mod_equal0[i]=sum(beta_B[i,]==0)/B
}



barplot(lasso_equal0[2:21],ylim = c(0,1)
        ,main="lasso",names=1:20,col = "cyan")
abline(h=0.5)
barplot(adalasso_equal0[2:21],ylim = c(0,1)
        ,main="adaptive lasso (gamma=2)",names=1:20,col = "cyan")
abline(h=0.5)
barplot(mod_equal0[2:21],ylim = c(0,1)
        ,main="modified bootstrap",names=1:20,col = "cyan")
abline(h=0.5)


#__________________________________________________________________
#setup 2
set.seed(100)
n=300
p=400
X=matrix(rnorm(n*p),nrow=n)
b=c(1,0.5,0.64,-0.4,0.6,0.5,0.7,1.2,0.61,0.8,rep(0,p-10))
y=X %*% b+rnorm(n,0,1) 
sdata=scale(data.frame(y,X))



cvfit_W=cv.glmnet(X,y,alpha=0)
ridge_W=glmnet(X,y,alpha=0,family="gaussian")
coef_w=coef(ridge_W,s=cvfit_W$lambda.min) # estimates
gamma=0
weight=1/(abs(coef_w[-1]))^gamma

#residualbootstrap
#for residual
cvfit=cv.glmnet(X,y,alpha=1,penalty.factor=weight)
adlassomodel=glmnet(X,y,alpha=1,family="gaussian",penalty.factor=weight)
adlassoest=coef(adlassomodel,s=cvfit$lambda.min) # estimates
res=y-predict(cvfit,s=cvfit$lambda.1se,newx=X) #residual
res=res-mean(res)

#for bootstrap
indicator=rep(0,p)
indicator[b!=0]=1
B=500
beta_B=matrix(0,nrow=p+1,ncol = B)
count=0
for(i in 1:B){
  s=sample(1:nrow(X),nrow(X),replace=TRUE)
  bootres=res[s]
  y_new=X%*%adlassoest[-1]+bootres
  
  cvfit_B=cv.glmnet(X,y_new,alpha=0)
  ridge_B=glmnet(X,y_new,alpha=0,family="gaussian")
  coef_B=coef(ridge_B,s=cvfit_B$lambda.min) # estimates
  weight=1/(abs(coef_B[-1]))^gamma
  
  cvfit=cv.glmnet(X,y_new,alpha=1,penalty.factor=weight)
  lambda.min=cvfit$lambda.min
  beta_B[,i]=as.numeric(coef(cvfit,s="lambda.min"))
  indicator_B=rep(0,p)
  indicator_B[beta_B[-1,i]!=0]=1
  if(setequal(indicator_B,indicator)==T){
    count=count+1
  }
}



se=apply(beta_B,1, sd)

conf=apply(beta_B,1, quantile,probs=c(0.025,0.975))

lasso_equal0=numeric(p+1)
for(i in 1:p+1){
  lasso_equal0[i]=sum(beta_B[i,]==0)/B
}

adalasso_equal0=numeric(p+1)
for(i in 1:p+1){
  adalasso_equal0[i]=sum(beta_B[i,]==0)/B
}


#p_value=numeric(p+1)
#for(i in 1:p+1){
#  p_value[i]=min(sum(beta_B[i,]<lassoest[i]),sum(beta_B[i,]>lassoest[i]))/B
#}
boxplot(t(beta_B[2:16,]),col="cyan",xlab="estimated coefficients",main="lasso")
abline(h=0,col="black")
y0s=c(b[1:15])
x0s=1:15-0.4
x1s=1:15+0.4
segments(x0=x0s,x1=x1s,y0=y0s,col="red",lwd=5)
segments(x0=x0s+0.2,x1=x1s-0.2,y0=conf[1,2:16],col="green",lwd=3,lty=3)
segments(x0=x0s+0.2,x1=x1s-0.2,y0=conf[2,2:16],col="green",lwd=3,lty=3)

#legend("topright",legend = c("line passing through 0","actual parameters given\n in simulation","95% confidence interval"),
#       col=c("black","red","green"),lty = c(1,1,3),lwd=c(1,5,3))

cbind(1:15,b,as.numeric(adlassoest)[-1],apply(beta_B, 1, mean)[-1],se[-1],t(conf[,-1]))



