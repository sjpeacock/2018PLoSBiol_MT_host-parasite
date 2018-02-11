# rm(list=ls())	

library(deSolve)
library(dclone)
library(gplots)
source("models.R")
source('plotting_functions.R')

######################################################################################################
# Apr 4: TD model with LOGnormal prior on r, all temps, 15 clones
######################################################################################################

# load("***whatever you called the saved workspace from DT_model_fitting.R***")

#---------------------
# Plot chains (10 chains)
#---------------------
col.chains<-c(1:8, colors()[c(90,139)])
parnames<-c("mu", "beta1", "beta2", "r", "K", "alpha")
thin.x<-seq(10,2000,10)

leave.out<-list(
	T1 =list(K1=NA, K2=c(4,6,3,8,10), K3=c(1,3,5,8,9), K4=c(2,7), K5=c(2,4,7), K6=c(4,7), K7=c(1,8,9), K8=c(2,5,1,10), K9=c(2,4,6,7), K10=c(4,6,10), K11=c(4,5,6,9), K12=c(6,7), K13=c(3,7,9), K14=c(1,2,9,3,6), K15=c(1,3,8,9)),
	T2 =list(K1=NA, K2=NA, K3=c(2,9), K4=c(1,4,5), K5=NA, K6=NA, K7=c(2,5), K8=NA, K9=NA, K10=c(1,6), K11=c(4), K12=NA, K13=NA, K14=NA, K15=NA),
	T3 =list(K1=NA, K2=NA, K3=NA, K4=NA, K5=NA, K6=NA, K7=NA, K8=NA, K9=NA, K10=NA, K11=NA, K12=NA, K13=NA, K14=NA, K15=NA),
	T4 =list(K1=NA, K2=NA, K3=NA, K4=NA, K5=NA, K6=NA, K7=NA, K8=NA, K9=NA, K10=NA, K11=NA, K12=NA, K13=NA, K14=NA, K15=NA),
	T5 =list(K1=c(2,7), K2=c(4,5), K3=c(5,4,6,2,10), K4=c(5), K5=c(7), K6=c(4,10,7,8), K7=c(1,4,7), K8=c(6,7,9,10), K9=c(5), K10=c(10), K11=c(2,3,10), K12=c(1,7), K13=c(2,3,4), K14=c(10), K15=c(1,3,6,10)),
	T6 =list(K1=NA, K2=NA, K3=NA, K4=NA, K5=c(2,3), K6=NA, K7=NA, K8=NA, K9=NA, K10=c(9), K11=NA, K12=c(2,3,10), K13=NA, K14=NA, K15=NA),
	T7 =list(K1=NA, K2=c(2), K3=NA, K4=c(2,4,6,7,8), K5=c(5,6,7), K6=c(3,10), K7=c(4,5,2,3), K8=NA, K9=c(2,10), K10=c(6,8,9,10), K11=c(2,3,5), K12=c(3,9), K13=c(1,5,7,9), K14=c(9), K15=c(4,8,9,1,10)),
	T8 =list(K1=NA, K2=NA, K3=NA, K4=NA, K5=NA, K6=NA, K7=NA, K8=NA, K9=NA, K10=NA, K11=NA, K12=NA, K13=NA, K14=NA, K15=NA),
	T9 =list(K1=NA, K2=NA, K3=NA, K4=c(6), K5=NA, K6=NA, K7=NA, K8=NA, K9=NA, K10=NA, K11=NA, K12=NA, K13=NA, K14=c(7), K15=c(1))
	)
		
chains.in<-list(); length(chains.in)<-9*length(n.clones.all); dim(chains.in)<-c(9,length(n.clones.all))
for(j in 1:9){
	for(k in 1:length(n.clones.all)){
		chains.in[[j,k]]<-which(is.na(match(c(1:10), leave.out[[j]][[k]]))==TRUE) 
	}
}

plotJK()
 
#---------------------
# Plot mean and var for each clone
#---------------------
parnames.all<-paste(rep(parnames, each=9), "[", rep(1:9, 6), "]", sep='')

fit.chains.in<-list(); length(fit.chains.in)<-max(n.clones.all)*9; dim(fit.chains.in)<-c(9,max(n.clones.all))
for(k in 1:max(n.clones.all)){ #for each clone
	for(j in 1:9){ #for each temp
		fit.chains.in[[j,k]]<-list(); length(fit.chains.in[[j,k]])<-length(chains.in[[j,k]])
		for(i in 1:length(chains.in[[j,k]])) fit.chains.in[[j,k]][[i]]<-fit.1T[[j,k]][[chains.in[[j,k]][i]]]
	}}
	

parEst.all<-list(); length(parEst.all)<-max(n.clones.all)
for(k in 1:max(n.clones.all)){ #for each clone
	parEst.all[[k]]<-matrix(NA, nrow=length(parnames.all), ncol=2)
	rownames(parEst.all[[k]])<-parnames.all
	colnames(parEst.all[[k]])<-c("mean", "variance")
	for(j in 1:9){ #for each temp
		S<-summary(as.mcmc(fit.chains.in[[j,k]]))
		parEst.all[[k]][seq(j, 54, 9),1]<-S[[1]][match(parnames, rownames(S[[1]])),1]
		parEst.all[[k]][seq(j, 54, 9),2]<-S[[1]][match(parnames, rownames(S[[1]])),2]^2
		}
}

parEst.all<-list(); length(parEst.all)<-max(n.clones.all)
for(k in 1:max(n.clones.all)){ #for each clone
	parEst.all[[k]]<-matrix(NA, nrow=length(parnames.all), ncol=2)
	rownames(parEst.all[[k]])<-parnames.all
	colnames(parEst.all[[k]])<-c("mean", "variance")
	for(j in 1:9){ #for each temp
		S<-summary(as.mcmc(fit.chains.in[[j,k]]))
		parEst.all[[k]][seq(j, 54, 9),1]<-S[[1]][match(parnames, rownames(S[[1]])),1]
		parEst.all[[k]][seq(j, 54, 9),2]<-S[[1]][match(parnames, rownames(S[[1]])),2]^2
		}
}

n.iter<-dim(fit.chains.in[[j,k]][[1]])[1]
parEst.allLOG<-list(); length(parEst.allLOG)<-max(n.clones.all)
for(k in 1:max(n.clones.all)){ #for each clone
	parEst.allLOG[[k]]<-matrix(NA, nrow=length(parnames.all), ncol=2)
	rownames(parEst.allLOG[[k]])<-parnames.all
	colnames(parEst.allLOG[[k]])<-c("mean", "variance")
	for(j in 1:9){ #for each temp
		xijk<-matrix(nrow=length(fit.chains.in[[j,k]])*n.iter, ncol=dim(fit.chains.in[[j,k]][[1]])[2])
		for(i in 1:length(fit.chains.in[[j,k]])) xijk[((i-1)*n.iter+1):(i*n.iter),]<-fit.chains.in[[j,k]][[i]]
		
		parEst.allLOG[[k]][seq(j, 54, 9),1]<-apply(log(xijk), 2, mean)[match(parnames, rownames(S[[1]]))]
		parEst.allLOG[[k]][seq(j, 54, 9),2]<-apply(log(xijk), 2, var)[match(parnames, rownames(S[[1]]))]
		
		S<-summary(as.mcmc(fit.chains.in[[j,k]]))
		parEst.all[[k]][seq(j, 54, 9),1]<-S[[1]][match(parnames, rownames(S[[1]])),1]
		parEst.all[[k]][seq(j, 54, 9),2]<-S[[1]][match(parnames, rownames(S[[1]])),2]^2
		}
}



Rhat<-list(); length(Rhat)<-max(n.clones.all); dim(Rhat)<-max(n.clones.all)
for(k in 1:max(n.clones.all)){ #for each clone
	Rhat[[k]]<-numeric(length(parnames.all))
	names(Rhat[[k]])<-parnames.all
	for(j in 1:9){ #for each temp
		Rhat[[k]][seq(j, 54, 9)]<-gelman.diag(fit.chains.in[[j,k]])[[1]][c(5,3,4,6,1,2),1]
		}}		


parEst.byClone<-list(); length(parEst.byClone)<-length(parnames.all)
parEst.byCloneLOG<-list(); length(parEst.byCloneLOG)<-length(parnames.all)
for(i in 1:length(parnames.all)){
	parEst.byClone[[i]]<-data.frame(n.clones=n.clones.all, mean=rep(NA, max(n.clones.all)), var=rep(NA, max(n.clones.all)), r.hat=rep(NA, max(n.clones.all)))
	for(k in 1:max(n.clones.all)){
		parEst.byClone[[i]]$mean[k]<-parEst.all[[k]][i,1]
		parEst.byClone[[i]]$var[k]<-parEst.all[[k]][i,2]	
		parEst.byClone[[i]]$r.hat[k]<-Rhat[[k]][match(parnames.all[i], names(Rhat[[k]]))]
		
		parEst.byCloneLOG[[i]]$mean[k]<-parEst.allLOG[[k]][i,1]
		parEst.byCloneLOG[[i]]$var[k]<-parEst.allLOG[[k]][i,2]	
		parEst.byCloneLOG[[i]]$r.hat[k]<-Rhat[[k]][match(parnames.all[i], names(Rhat[[k]]))]

	}
}


# pdf("Figures/T1_paramsSE_4Apr2017_lognormal.pdf", onefile=TRUE, width=6, height=5)
# quartz(width=6, height=6)
par(mfrow=c(3,3), mar=c(2,2,1,1), oma=c(3,3,3,0))
for(m in 1:6){ #For each parameter
	for(j in 1:9){ #for each temperature
		i<-(m-1)*9+j
		plotCI(parEst.byClone[[i]]$n.clones, parEst.byClone[[i]]$mean, uiw=sqrt(parEst.byClone[[i]]$var), liw=sqrt(parEst.byClone[[i]]$var), bty="l", pt.bg=c(2,3)[as.numeric(parEst.byClone[[i]]$r.hat<1.1)+1], gap=0, pch=21, xlim=c(0.5,max(n.clones.all)+0.5), xaxt="n", ylim=range((parEst.byClone[[i]]$mean-sqrt(parEst.byClone[[i]]$var)), (parEst.byClone[[i]]$mean+sqrt(parEst.byClone[[i]]$var))))#, dat$prior.mean[i]-dat$prior.sd[i], dat$prior.mean[i]+dat$prior.sd[i]))
		if(j>6) axis(side=1, at=parEst.byClone[[i]]$n.clones) else axis(side=1, at=parEst.byClone[[i]]$n.clones, labels=FALSE)
		mtext(side=3, paste(parnames.all[i], "-", round(dat$temp[j]), "*C"), cex=par('cex'), font=2)
		
		# Plot prior overtop
		par(new=TRUE)
		y<-seq(min(parEst.byClone[[i]]$mean-sqrt(parEst.byClone[[i]]$var)), max(parEst.byClone[[i]]$mean+sqrt(parEst.byClone[[i]]$var)), length.out=100)
		x<-dlnorm(y, meanlog=dat$prior.mean[c(rep(1, 9), rep(2, 18), rep(3:5, each=9))[i]], sdlog=dat$prior.sd[c(rep(1, 9), rep(2, 18), rep(3:5, each=9))[i]])
		plot(x,y,"l", xaxt="n", yaxt="n", xlab="", ylab="", bty="n", col="#00000030", lwd=2)
		
		
		# abline(h=dat$prior.mean[i], col="#00000070", lwd=2)
		# abline(h=dat$prior.mean[i]-dat$prior.sd[i], col="#00000070", lty=2)
		# abline(h=dat$prior.mean[i]+dat$prior.sd[i], col="#00000070", lty=2)
	}}
mtext(side=1, outer=TRUE, "Number of clones (K)", line=1)
mtext(side=2, outer=TRUE, "Posterior estimate (mean +/- SE)", line=1)
	
dev.off()


# Plot all temps together

# #---------------------
# # Plot variance over number of clones
# #---------------------

# # pdf("Figures/T1_paramsScaledVar_4Apr2017_lognormal.pdf", onefile=TRUE, width=6, height=5)
par(mfrow=c(3,3), mar=c(2,2,1,1), oma=c(3,3,3,0))
for(m in 1:6){ #For each parameter
	for(j in 1:9){ #for each temperature
		i<-(m-1)*9+j
		plot(parEst.byClone[[i]]$n.clones, parEst.byClone[[i]]$var/parEst.byClone[[i]]$var[1], type="l", bty="l", xlim=c(0.5,max(n.clones.all)+0.5), xaxt="n", yaxt="n", ylim=c(0, max(parEst.byClone[[i]]$var/parEst.byClone[[i]]$var[1])))
		lines(parEst.byClone[[i]]$n.clones, 1/parEst.byClone[[i]]$n.clones, lty=3)
		points(parEst.byClone[[i]]$n.clones, parEst.byClone[[i]]$var/parEst.byClone[[i]]$var[1], bg=c(2,3)[as.numeric(parEst.byClone[[i]]$r.hat<1.1)+1], pch=21)
		if(i>18) axis(side=1, at=parEst.byClone[[i]]$n.clones) else axis(side=1, at=parEst.byClone[[i]]$n.clones, labels=FALSE)
		if(is.element(i, seq(1,21,3))) axis(side=2, at=c(0,1)) else axis(side=2, at=c(0,1), labels=FALSE)
		mtext(side=3, paste(parnames.all[i], "-", round(dat$temp[j]), "*C"), cex=par('cex'), font=2)
			
}}

mtext(side=1, outer=TRUE, "Number of clones (K)", line=1)
mtext(side=2, outer=TRUE, "Scaled variance in posterior", line=1)
# # dev.off()

#---------------------
# Plot TD estimates for k clones
#---------------------
k<-15

params.TD<-list(mean=matrix(parEst.all[[k]][,1], nrow=9, ncol=6), li=matrix(parEst.all[[k]][,1]-1.96*sqrt(k*parEst.all[[k]][,2]), nrow=9, ncol=6), ui=matrix(parEst.all[[k]][,1]+1.96*sqrt(k*parEst.all[[k]][,2]), nrow=9, ncol=6))

params.TDLOG<-list(mean=matrix(parEst.allLOG[[k]][,1], nrow=9, ncol=6), li=matrix(parEst.allLOG[[k]][,1]-1.96*sqrt(k*parEst.allLOG[[k]][,2]), nrow=9, ncol=6), ui=matrix(parEst.allLOG[[k]][,1]+1.96*sqrt(k*parEst.allLOG[[k]][,2]), nrow=9, ncol=6))

for(i in 1:3){
	colnames(params.TD[[i]])<-parnames
	colnames(params.TDLOG[[i]])<-parnames
	}

rhat.TD<-params.TD[[1]]
for(i in 1:6){for(j in 1:9){
		rhat.TD[j,i]<-parEst.byClone[[(i-1)*2+j]]$r.hat[k]
}}

# Matrix of point colors: red = no convergence, yellow = convergence, but not estimable, green = converged and estimable
point.col<-matrix(1, nrow=9, ncol=6)
colnames(point.col)<-parnames 
point.col[,'mu']<-rep(3, 9)
point.col[,'beta1']<-rep(3, 9)
point.col[,'beta2']<-rep(3, 9)
point.col[,'r']<-c(2,2,rep(3,6), "gold")
point.col[,'K']<-c(2,2,rep(3,5), rep("gold", 2))
point.col[,'alpha']<-c(2,2,rep(3,5), rep("gold", 2))

parnames.fig<-c(expression(mu), expression(beta[1]), expression(beta[2]), expression(italic(r)), expression(italic(K)), expression(alpha))

# MT.pred<-read.csv("~/Google Drive/Kirk/Code/Updated data/MT_pred_Mar22.csv")

ylims<-list(c(0, 0.35), c(0, 10), c(0,10), c(-1, 2), c(0, 200), c(0, 0.005))
# quartz(width=6, height=4)
par(mfrow=c(2,3), mar=c(3,3,2,0), oma=c(1,2,0,1))
for(i in 1:6){
	# plotCI(dat$temp, params.TD[[1]][,i], li=params.TD[[2]][,i], ui=params.TD[[3]][,i], xlab="", ylab="", las=1, bty="l", main=parnames.fig[i], pch=19, gap=0, ylim=ylims[[i]], xlim=c(3,38), col=point.col[,i])
	plotCI(dat$temp, exp(params.TDLOG[[1]][,i]), li=exp(params.TDLOG[[2]][,i]), ui=exp(params.TDLOG[[3]][,i]), xlab="", ylab="", las=1, bty="l", main=parnames.fig[i], pch=19, gap=0, ylim=ylims[[i]], xlim=c(3,38), col=point.col[,i])
	plotCI(dat$temp+0.5, params.TD[[1]][,i], li=params.TD[[2]][,i], ui=params.TD[[3]][,i], add=TRUE, cex=0.8, gap=0.3)
# lines(MT.pred[,1], MT.pred[,i+1])
	if(i==4) abline(h=0, lty=2)
	if(i==6) legend("topleft", pch=19, col=c(3,"gold", 2), c("Converged, estimable", "Converged, not estimable", "Not converged"), bty="n")
	}
mtext(side=1, outer=TRUE, expression(paste("Temperature (", degree, 'C)')))
mtext(side=2, outer=TRUE, "Parameter estimate", line=0.5)

#---------------------
# Table of parameter estimates
#---------------------

T1_est<-data.frame(parameter=rep(parnames, each=9), temp=rep(dat$temp, 6), mean=as.numeric(params.TD$mean), SE=sqrt(k*parEst.all[[k]][,2]), li=as.numeric(params.TD$li), ui=as.numeric(params.TD$ui), rhat=Rhat[[k]], estimable=as.numeric(point.col==3))

T1_estLOG<-data.frame(parameter=rep(parnames, each=9), temp=rep(dat$temp, 6), mean=as.numeric(params.TDLOG$mean), SE=sqrt(k*parEst.allLOG[[k]][,2]), meanEXP=exp(as.numeric(params.TDLOG$mean)), li=exp(as.numeric(params.TDLOG$li)), ui=exp(as.numeric(params.TDLOG$ui)), rhat=Rhat[[k]], estimable=as.numeric(point.col==3))

library(xlsx)

write.xlsx(T1_est, file="TD_ParameterEstimates_2017Apr4.xlsx")
