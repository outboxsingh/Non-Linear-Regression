#Part-1
data=read.table(file.choose())
library(polynom)
data
colnames(data)=c("t","y")
plot(data$y,xlab="t",ylab="y",main = "Scatterplot of y against t",col="black",pch=20)


#Part-2

n=25
k=2
y=data$y
t=data$t
z=matrix(0,ncol=(k+1),nrow=(n-k))
for (i in 1:(n-k))
{
  
  for (j in 1:(k+1))
  {
    z[i,j]=y[i+j-1]
  }
}
v=t(z)%*%(z)
e=eigen(v)
e$vectors
e$values
p=as.polynomial(c(e$vectors[1,3],e$vectors[2,3],e$vectors[3,3]))
pz=solve(p)
pz
bhat=log(pz)
x=matrix(0,nrow=25,ncol=2)
for (i in 1:25)
{
  for (j in 1:2)
  {
    x[i,j]=exp(bhat[j]*i)
  }
}
ahat=solve(t(x)%*%x)%*%t(x)%*%y
yhat=array(0)
for (t in 1:25)
{
yhat[t]=(ahat[1]*exp(t*bhat[1])+ahat[2]*exp(t*bhat[2]))
}
plot(yhat,type="l",ylab="Predicted values",xlab="t",main ="Predicted values from Prony's Method")
plot(y-yhat,ylab="Error",xlab="t",main ="Error from Prony's Method")


#Part-3

library(ggplot2)
library(MASS)
a=seq(-3,0.1,0.1)
b=seq(-3,0.1,0.1)
y1=matrix(0,nrow=length(a),ncol=length(b))
expand=expand.grid(a,b)
t=1:25
for(i in 1:nrow(expand))
{
  x=matrix(c(expand[i,1]*t,expand[i,2]*t),ncol=2)
  alp=ginv(t(x)%*%x)%*%t(x)%*%y
  yhat=alp[1]*exp(t*expand[i,1])+alp[2]*exp(t*expand[i,2])
  expand[i,3]=sum((y-yhat)^2)
}
colnames(expand)=c("x","y","z")
plot1 <- ggplot(expand, aes(x = x, y = y, z = z)) + stat_contour()+labs(title = "Contour plot",x="beta1",y="beta2")
plot1

#part-4

a1a=ahat[1]
a2a=ahat[2]
b1a=bhat[1]
b2a=bhat[2]
theta_app=c(a1a,b1a,a2a,b2a)
theta_app1=c(0,0,0,0)
k=0
while(sum((theta_app1-theta_app)^2)>0.0005)
{
  if(k>0)
  {
    theta_app=theta_app1
  }
  theta=theta_app
f=matrix(0,nrow=25,ncol=4)
for (i in 1:25)
{
  for(j in 1:4)
  {
    if(j==1)
    {
      f[i,j]=exp(theta[2]*i)
    }
    if(j==2)
    {
      f[i,j]=theta[1]*i*exp(theta[2]*i)
    }
    if(j==3)
    {
      f[i,j]=exp(theta[4]*i)
    }
    if(j==4)
    {
      f[i,j]=theta[3]*i*exp(theta[4]*i)
    }
  }
}
func=function(theta)
{
  t=1:25
  return (theta[1]*exp(t*theta[2])+theta[3]*exp(t*theta[4]))
}
da=solve(t(f)%*%f)%*%t(f)%*%(data$y-func(theta_app))
theta_app1=da+theta_app
k=k+1
print(k)
}
t=1:25
yhat1=theta[1]*exp(theta[2]*t)+theta[3]*exp(theta[4]*t)
plot(yhat1,type="l",ylab="Predicted values",xlab="t",main ="Predicted values from Gauss-Newton Method")
lines(y,col="red")
yhat1
plot(y-yhat1,ylab="Error",xlab="t",main ="Error from Gauss Newton Method")

#Part-5

u1t=cbind(diag(23),matrix(0,nrow=23,ncol=2))
u2t=cbind(matrix(0,nrow=23,ncol=1),diag(23),matrix(0,nrow=23,ncol=1))
u3t=cbind(matrix(0,nrow=23,ncol=2),diag(23))
U=list(u1t,u2t,u3t)
ga=c(0.01707969,-0.36000500,0.93279402)
ga=ga/sqrt(sum(ga^2))
ga1=c(0,0,0)
l=0
while(sqrt(sum(ga-ga1)^2)>0.000001 ) 
{
  if(l>0)
  {
    ga=ga1
  }
Gt=ga[1]*u1t+ga[2]*u2t+ga[3]*u3t
b=matrix(0,nrow=3,ncol=3)
for (k in 1:3)
{
  for(j in 1:3)
  {
    b[k,j]=(t(y)%*%t(U[[k]])%*%solve(Gt%*%t(Gt))%*%U[[j]]%*%y)+(t(y)%*%t(U[[j]])%*%solve(Gt%*%t(Gt))%*%U[[k]]%*%y)+(t(y)%*%t(Gt)%*%solve(Gt%*%t(Gt))%*%(U[[k]]%*%t(U[[j]])+U[[j]]%*%t(U[[k]]))%*%solve(Gt%*%t(Gt))%*%Gt%*%y)                                         
  }
}
e=eigen(b)
ga1=e$vectors[,3]
ga1=ga1/sqrt(sum(ga1^2))
l=l+1
print(l)
}
library(polynom)
p=as.polynomial(c(ga1[1],ga1[2],ga1[3]))
pz=solve(p)
pz
bhat1=log(pz)
x=matrix(0,nrow=25,ncol=2)
for (i in 1:25)
{
  for (j in 1:2)
  {
    x[i,j]=exp(bhat1[j]*i)
  }
}
ahat1=solve(t(x)%*%x)%*%t(x)%*%y



#Part-8
error1=function(k1)
{
  t=1:25
S1=array(0)
for (p in 1:25)
{
  y_temp=y[-p]
  t_temp=t[-p]
  n=length(y_temp)
k=2
z=matrix(0,ncol=(k+1),nrow=(n-k))

for (i in 1:(n-k))
{
  for (j in 1:(k+1))
  {
    z[i,j]=y_temp[i+j-1]
  }
}
v=t(z)%*%(z)
e=eigen(v)
e$vectors
e$values
p1=as.polynomial(e$vectors[,k+1])
pz=solve(p1)
pz
bhat=log(abs(pz))
x=matrix(0,nrow=n,ncol=k)
for (i in 1:length(t_temp))
{
  for (j in 1:k)
  {
    x[i,j]=exp(bhat[j]*i)
  }
}
ahat=ginv(t(x)%*%x)%*%t(x)%*%y_temp
S1[p]=(y[p]-t(exp(bhat*p))%*%ahat)^2
}
return(S1)
}
error1(2)

e1=array(0)
sq=seq(2,20,1)
for (g in sq)
{
  e1[g-1]=mean(error1(g))
  print(g)
}
plot(e1,type="l",ylab="CVE(k)",xlab="k",main="Plot of Cross Validation Error against k")


#Confidence interval of parameters

t=1:25
Mse=sum((y-ahat[1]*exp(bhat[1]*t)-ahat[2]*exp(bhat[2]*t))^2)
Mse=Mse/21
Mse
theta=c(ahat[1],bhat[1],ahat[2],bhat[2])
f=matrix(0,nrow=25,ncol=4)
for (i in 1:25)
{
  for(j in 1:4)
  {
    if(j==1)
    {
      f[i,j]=exp(theta[2]*i)
    }
    if(j==2)
    {
      f[i,j]=theta[1]*i*exp(theta[2]*i)
    }
    if(j==3)
    {
      f[i,j]=exp(theta[4]*i)
    }
    if(j==4)
    {
      f[i,j]=theta[3]*i*exp(theta[4]*i)
    }
  }
}
f
var=solve(t(f)%*%f)
var=Mse*25*var
a1_conf=c(ahat[1]-sqrt(var[1,1])*qt(0.95,24,lower.tail = T)/5,ahat[1]+sqrt(var[1,1])*qt(0.95,24,lower.tail = T)/5)
a2_conf=c(ahat[2]-sqrt(var[3,3])*qt(0.95,24,lower.tail = T)/5,ahat[2]+sqrt(var[3,3])*qt(0.95,24,lower.tail = T)/5)
b1_conf=c(bhat[1]-sqrt(var[2,2])*qt(0.95,24,lower.tail = T)/5,bhat[1]+sqrt(var[2,2])*qt(0.95,24,lower.tail = T)/5)
b2_conf=c(bhat[2]-sqrt(var[4,4])*qt(0.95,24,lower.tail = T)/5,bhat[2]+sqrt(var[4,4])*qt(0.95,24,lower.tail = T)/5)


















