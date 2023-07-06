
#-----------------------------------------------------------------------------------
# Author: Irene Tubikanec
# Date:  2023-05-23
#-----------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
# Visualizaton of SMC-ABC results for the JRNMM
#-----------------------------------------------------------------------------------


N<-4 #number of populations in the JRNMM
M<-500 #number of kept samples

#--------------------------------------------------------
#Read results

filename<-"ABC_Results" 

#weights
weights_c<-as.vector(t(read.table(paste(filename,"/norm_weights_c.txt",sep=""),header=F)))
weights_c<-weights_c/sum(weights_c)

#A's
A_kept<-matrix(nrow=M,ncol=N)
for(i in 1:N){
  A_kept[,i]<-as.vector(t(read.table(paste(filename,"/A",i,"vec.txt",sep=""),header=F)))
}

#L,c
L_kept<-as.vector(t(read.table(paste(filename,"/Lvec.txt",sep=""),header=F)))
c_kept<-as.vector(t(read.table(paste(filename,"/cvec.txt",sep=""),header=F)))

#Rho's
p_kept<-matrix(0,nrow=M,ncol=N*(N-1))
counter<-0
for(j in 1:N){
  for(k in 1:N){
    if(j!=k){
      counter<-counter+1
      p_kept[,counter]<-as.vector(t(read.table(paste(filename,"/p",j,k,"vec.txt",sep = ""),header=F)))
    }
  }
}

#--------------------------------------------------------
#Prior distribution for continuous model parameters
Pr_cont<-matrix(0,nrow=N+2,ncol=2)
for(j in 1:N){
  Pr_cont[j,]<-c(2,4) #A's
}
Pr_cont[N+1,]<-c(100,2000) #L
Pr_cont[N+2,]<-c(0.5,1) #c

#--------------------------------------------------------
#Calculation of estimates: weighted ABC posterior means

#A's
hat_As_mean<-rep(0,N)
for(i in 1:N){
  hat_As_mean[i]<-weighted.mean(A_kept[,i],w=weights_c)
}

#L
hat_L_mean<-weighted.mean(L_kept,w=weights_c)

#c
hat_c_mean<-weighted.mean(c_kept,w=weights_c)

#Rho's
p_mode<-rep(0,N*(N-1))
counter<-0
for(j in 1:N){
  for(k in 1:N){
    if(j!=k){
      counter<-counter+1
      p_hat<-sum(p_kept[,counter])/M
      if(p_hat>(1/2)){
        p_mode[counter]<-1
      } 
    }
  }
}

#--------------------------------------------------------
# Plot of ABC results: continuous model parameters

par(mfrow=c(2,5),mai = c(0.3, 0.3, 0.3, 0.08))

#A1
i<-1
Al<-Pr_cont[i,1]
Ar<-Pr_cont[i,2]
plot(density(A_kept[,i],weights=weights_c,from=Al,to=Ar),col="blue",xlim=c(Al,Ar),xlab="",main="")
abline(v=3.6,col="green",lty=1,lwd=1)
curve(dunif(x,min=Al,max=Ar),add=TRUE,col="red",lwd=1)
leg<-c( expression(paste(pi(A[1]))), expression(paste(pi[ABC])(paste(A[1],"|",y))), expression(A[1]), expression(hat(A)[1]) )
legend("topleft",legend= leg  ,lty=c(1,1,1,3),col=c("red","blue","green","black"),cex=1.1,lwd=c(1,1,1,1.5),seg.len=1.0,bty="n")
abline(v=hat_As_mean[i],lty=3,col="black",lwd=1.5)

#A2
i<-2
Al<-Pr_cont[i,1]
Ar<-Pr_cont[i,2]
plot(density(A_kept[,i],weights=weights_c,from=Al,to=Ar),col="blue",xlim=c(Al,Ar),xlab="",main="")
abline(v=3.25,col="green",lty=1,lwd=1)
curve(dunif(x,min=Al,max=Ar),add=TRUE,col="red",lwd=1)
leg<-c( expression(paste(pi(A[2]))), expression(paste(pi[ABC])(paste(A[2],"|",y))), expression(A[2]), expression(hat(A)[2]) )
legend("topleft",legend= leg  ,lty=c(1,1,1,3),col=c("red","blue","green","black"),cex=1.1,lwd=c(1,1,1,1.5),seg.len=1.0,bty="n")
abline(v=hat_As_mean[i],lty=3,col="black",lwd=1.5)

#A3
i<-3
Al<-Pr_cont[i,1]
Ar<-Pr_cont[i,2]
plot(density(A_kept[,i],weights=weights_c,from=Al,to=Ar),col="blue",xlim=c(Al,Ar),xlab="",main="")
abline(v=3.25,col="green",lty=1,lwd=1)
curve(dunif(x,min=Al,max=Ar),add=TRUE,col="red",lwd=1)
leg<-c( expression(paste(pi(A[3]))), expression(paste(pi[ABC])(paste(A[3],"|",y))), expression(A[3]), expression(hat(A)[3]) )
legend("topleft",legend= leg  ,lty=c(1,1,1,3),col=c("red","blue","green","black"),cex=1.1,lwd=c(1,1,1,1.5),seg.len=1.0,bty="n")
abline(v=hat_As_mean[i],lty=3,col="black",lwd=1.5)

#A4
i<-4
Al<-Pr_cont[i,1]
Ar<-Pr_cont[i,2]
plot(density(A_kept[,i],weights=weights_c,from=Al,to=Ar),col="blue",xlim=c(Al,Ar),xlab="",main="")
abline(v=3.25,col="green",lty=1,lwd=1)
curve(dunif(x,min=Al,max=Ar),add=TRUE,col="red",lwd=1)
leg<-c( expression(paste(pi(A[4]))), expression(paste(pi[ABC])(paste(A[4],"|",y))), expression(A[4]), expression(hat(A)[4]) )
legend("topleft",legend= leg  ,lty=c(1,1,1,3),col=c("red","blue","green","black"),cex=1.1,lwd=c(1,1,1,1.5),seg.len=1.0,bty="n")
abline(v=hat_As_mean[i],lty=3,col="black",lwd=1.5)

#L
Ll<-Pr_cont[N+1,1]
Lr<-Pr_cont[N+1,2]
plot(density(L_kept,weights=weights_c,from=Ll,to=Lr),col="blue",xlim=c(Ll,Lr),xlab="",main="")
abline(v=700,col="green",lty=1,lwd=1)
curve(dunif(x,min=Ll,max=Lr),add=TRUE,col="red",lwd=1)
legend("topright",legend=c(expression(pi(L)),expression(paste(pi[ABC])(paste(L,"|",y))),expression(L),expression(hat(L))),
       lty=c(1,1,1,3),col=c("red","blue","green","black"),cex=1.1,lwd=c(1,1,1,1.5),seg.len=1.0,bty="n")
abline(v=hat_L_mean,lty=3,col="black",lwd=1.5)

#c
cl<-Pr_cont[N+2,1]
cr<-Pr_cont[N+2,2]
plot(density(c_kept,weights=weights_c,from=cl,to=cr),col="blue",xlim=c(cl,cr),xlab="",main="",ylim=c(0,4))
abline(v=0.8,col="green",lty=1,lwd=1)
curve(dunif(x,min=cl,max=cr),add=TRUE,col="red",lwd=1)
legend("topleft",legend=c(expression(pi(c)),expression(paste(pi[ABC])(paste(c,"|",y))),expression(hat(c))),
       lty=c(1,1,3),col=c("red","blue","black"),cex=1.1,lwd=c(1,1,1.5),seg.len=1.0,bty="n")
abline(v=hat_c_mean,lty=3,col="black",lwd=1.5)

#--------------------------------------------------------
# Plot of ABC results: Bernoulli network parameters

par(mfrow=c(3,6),mai=c(0.3,0.3,0.25,0.01)) 

counter<-0
for(j in 1:N){
  for(k in 1:N){
    if(j!=k){
      
      counter<-counter+1
      
      if(p_mode[counter]==0){
        hist(p_kept[,counter],c(-0.5,0.5,1.5),freq=FALSE,main="",ylab="",xlab="",col="grey",xaxt="n",yaxt="n",ylim=c(0,1))
        mtext(bquote(paste(pi[ABC])(paste(rho[paste(.(j),.(k))],"|",y))), line = 0, side = 3, outer = F,cex=0.9)
        axis(1,c(0,1))
        axis(2,c(0,0.5,1))
      }
      if(p_mode[counter]==1){
        hist(p_kept[,counter],c(-0.5,0.5,1.5),freq=FALSE,main="",ylab="",xlab="",col="grey",xaxt="n",yaxt="n",ylim=c(0,1))
        mtext(bquote(paste(pi[ABC])(paste(rho[paste(.(j),.(k))],"|",y))), line = 0, side = 3, outer = F,cex=0.9)
        axis(1,c(0,1))
        axis(2,c(0,0.5,1))
      }
      if(k==j+1||(j==3&&k==2)||(j==1&&k==3)){
        lines(c(1,1),c(0,1),col="green",lty=1,lwd=2) #line for 1
      } else {
        lines(c(0,0),c(0,1),col="green",lty=1,lwd=2) #line for 0
      }
      
    }
  }
}

