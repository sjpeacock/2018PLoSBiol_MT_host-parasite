setwd("~/Google Drive/Kirk/Code/Updated data")
rm(list=ls())	

library(deSolve)
library(dclone)
library(gplots)
library(xlsx)
source("models.R")
source('plotting_functions.R')


params.TD<-read.xlsx("DT/DT_ParameterEstimates_2017Apr4_lognormal.xlsx", sheetIndex=1)

######################################################################################################
######################################################################################################
# July 8, 2017: MT model with fixed EL_mu and E_K normal prior
######################################################################################################
######################################################################################################


# load('Workspaces/MTfit_fixedEL_mu_2017-07-08.RData')

plot(1:25, t1.all, "o", ylab="Time to fit (hrs)", xlab="Number of clones (K)", bty="l", las=1); abline(h=mean(t1.all), lty=2)
sum(t1.all)/24

#---------------------
# Plot chains (10 chains)
#---------------------
col.chains<-c(1:8, colors()[c(90,139)])

# Which chains to leave out because they did not converge, for clone K
leave.out<-list(
	K1=c(NA), 
	K2=c(1,4,10) , 
	K3=c(7), 
	K4=c(NA), 
	K5=c(1), 
	K6=c(7,8),
	K7=c(5,6),
	K8=c(3,9),
	K9=c(NA),
	K10=c(1,6),
	K11=c(10),
	K12=c(NA),
	K13=c(7,10),
	K14=c(NA),
	K15=c(1,8),
	K16=c(1,10),
	K17=c(2),
	K18=c(NA),
	K19=c(NA),
	K20=c(9),
	K21=c(6),
	K22=c(NA),
	K23=c(2),
	K24=c(1,2,8),
	K25=c(3,4,10)	
	) 
chains.in<-list(); length(chains.in)<-length(fit.MT)
for(k in 1:length(fit.MT)) chains.in[[k]]<-which(is.na(match(1:length(fit.MT[[k]]), leave.out[[k]]))==TRUE)
#for(k in 1:length(fit.MT)) chains.in[[k]]<-1:10

parnames<-as.character(priors$parameter)[c(1:2,4:21)]
thin.x<-seq(10,2000,10)

# pdf("MTmodel_all_clones_6Mar2017.pdf", onefile=TRUE, width=6, height=7)
k<-25# for(k in 1:length(fit.MT)){
	par(mfrow=c(5,4), mar=c(1,2,1,1), oma=c(2,2,3,0))
	for(i in 1:20){
		
		xi<-matrix(NA, nrow=length(thin.x), ncol=length(fit.MT[[k]]))
		for(j in 1:length(fit.MT[[k]])) xi[,j]<-fit.MT[[k]][[j]][thin.x,match(parnames[i], colnames(fit.MT[[k]][[1]]))]
		
		plot(thin.x, xi[,1], type='n', xaxt="n", bty="l", ylim=range(xi[,chains.in[[k]]]))
		for(j in chains.in[[k]]) lines(thin.x, xi[,j], col=col.chains[j])
		
		mtext(side=3, parnames[i])
	}
	legend("right", fill=col.chains[chains.in[[k]]], bg="white", legend=chains.in[[k]], ncol=2, border=NA)
	
	mtext(side=3, outer=TRUE, paste(k, "clones"), line=1)
# }
# dev.off()

#---------------------
# Plot mean and var for each clone
#---------------------

fit.chains.in<-list(); length(fit.chains.in)<-length(fit.MT)
for(k in 1:length(fit.MT)){
	fit.chains.in[[k]]<-list(); length(fit.chains.in[[k]])<-length(chains.in[[k]])
	for(j in 1:length(chains.in[[k]])) fit.chains.in[[k]][[j]]<-fit.MT[[k]][[chains.in[[k]][j]]]
	}

parEst.all<-list(); length(parEst.all)<-length(fit.MT)
for(k in 1:length(fit.MT)){
	parEst.all[[k]]<-MeanVarest(fit.MT[[k]], leave.out=leave.out[[k]], par.order= parnames)
}

Rhat<-list(); length(Rhat)<-length(fit.MT)
for(k in 1:length(fit.MT)) Rhat[[k]]<-gelman.diag(fit.chains.in[[k]])		

parEst.byClone<-list(); length(parEst.byClone)<-20
for(i in 1:20){
	parEst.byClone[[i]]<-data.frame(n.clones=c(1:length(fit.MT)), mean=rep(NA, length(fit.MT)), var=rep(NA, length(fit.MT)), r.hat=rep(NA, length(fit.MT)))
	for(k in 1:length(fit.MT)){
		parEst.byClone[[i]]$mean[k]<-parEst.all[[k]][i,1]
		parEst.byClone[[i]]$var[k]<-parEst.all[[k]][i,2]	
		parEst.byClone[[i]]$r.hat[k]<-Rhat[[k]][[1]][which(rownames(Rhat[[k]][[1]])==parnames[i]),1]
	}
}

col.title<-rep(1, 20)#c(3,3,2,2,2,2,3,3,3,2,3,3,2,2,3,3,2,3,2,3,3)
prior.type.MT<-c(0,0,0,1,1,0,0,0,1,0,0,0,1,0,1,0,0,1,1,0) # Changed for E_K ~ normal!!


parnamesMT.fig<-c(expression(mu[0]), expression(italic(E[mu])), expression(italic(E[H*mu])), expression(italic(T[L*mu])), expression(italic(T[H*mu])), expression(beta[0]), expression(italic(E[beta])), expression(italic(E[H*beta])), expression(italic(T[H*beta])), expression(italic(r[0])), expression(italic(E[r])), expression(italic(E[H*r])), expression(italic(T[H*r])), expression(italic(K[0])), expression(italic(E[K])), expression(italic(E[L*K])), expression(italic(E[H*K])), expression(italic(T[L*K])), expression(italic(T[H*K])), expression(alpha[0]))

pdf("Figures/MT_paramsSE_2017Jul8.pdf", width=8, height=9)
par(mfrow=c(5,4), mar=c(1.5,2,1,1), oma=c(2,3,1,0))
for(i in 1:20){
	plotCI(parEst.byClone[[i]]$n.clones, parEst.byClone[[i]]$mean, uiw=sqrt(parEst.byClone[[i]]$var), liw=sqrt(parEst.byClone[[i]]$var), bty="l", pt.bg=c(2,3)[as.numeric(parEst.byClone[[i]]$r.hat<1.1)+1], gap=0, pch=21, xlim=c(0.5,length(fit.MT)+0.5), xaxt="n", cex=0.8)
	
	if(i>16) axis(side=1, at=seq(0, 25, 5)) else axis(side=1, at=seq(0, 25, 5), labels=FALSE)
	mtext(side=3, parnamesMT.fig[i], cex=par('cex'), font=2)
	
	# Plot prior overtop
	par(new=TRUE)
	y<-seq(min(parEst.byClone[[i]]$mean-sqrt(parEst.byClone[[i]]$var)), max(parEst.byClone[[i]]$mean+sqrt(parEst.byClone[[i]]$var)), length.out=100)
	if(prior.type.MT[i]==1){ #if temp param, prior is normal
		x<-dnorm(y, mean=dat$prior.mean[i], sd=dat$prior.sd[i])		
	}else{# if other param, prior is lognormal
		x<-dlnorm(y, meanlog=dat$prior.mean[i], sdlog=dat$prior.sd[i])
		}
	plot(x,y,"l", xaxt="n", yaxt="n", xlab="", ylab="", bty="n", col="#00000030", lwd=2)
	}
	
	# abline(h=dat$prior.mean[i], col="#00000070", lwd=2)
	# abline(h=dat$prior.mean[i]-dat$prior.sd[i], col="#00000070", lty=2)
	# abline(h=dat$prior.mean[i]+dat$prior.sd[i], col="#00000070", lty=2)
mtext(side=1, outer=TRUE, "Number of clones (K)", line=1)
mtext(side=2, outer=TRUE, "Posterior estimate (mean +/- SE)", line=1)
dev.off()

#---------------------
# Plot variance over number of clones
#---------------------

pdf("Figures/MT_paramsScaledVar_2017Jul8.pdf", width=8, height=9)
par(mfrow=c(5,4), mar=c(1.5,1,1,0), oma=c(2,3,1,1))
for(i in 1:20){
	plot(parEst.byClone[[i]]$n.clones, parEst.byClone[[i]]$var/parEst.byClone[[i]]$var[1], type="l", bty="l", xlim=c(0.5,length(fit.MT)+0.5), xaxt="n", yaxt="n", ylim=c(0, max(parEst.byClone[[i]]$var/parEst.byClone[[i]]$var[1])))
	lines(parEst.byClone[[i]]$n.clones, 1/parEst.byClone[[i]]$n.clones, lty=3)
	points(parEst.byClone[[i]]$n.clones, parEst.byClone[[i]]$var/parEst.byClone[[i]]$var[1], bg=c(2,3)[as.numeric(parEst.byClone[[i]]$r.hat<1.1)+1], pch=21)
	if(i>16) axis(side=1, at=seq(0, 25, 5)) else axis(side=1, at=seq(0, 25, 5), labels=FALSE)
	if(is.element(i, seq(1,20,4))) axis(side=2, at=c(0,1), las=1) else axis(side=2, at=c(0,1), labels=FALSE)
	mtext(side=3, parnamesMT.fig[i], cex=par('cex'), font=2, col=col.title[i])	
}

mtext(side=1, outer=TRUE, "Number of clones (K)", line=1)
mtext(side=2, outer=TRUE, "Scaled variance in posterior", line=1)
dev.off()

#---------------------
# Plot MT predictions with TD estimates for 10 clones
#---------------------
k<-25
par.with.EL_mu<-c(parEst.all[[k]][1:2,1], 5*parEst.all[[k]][2,1], parEst.all[[k]][3:20,1])
names(par.with.EL_mu)<-as.character(priors$parameter)

MT.pred<-MT.parPred(par.with.EL_mu)
write.csv(MT.pred, file="MT_pred_Jul8.csv")

MT.pred_old<-read.csv("MT_pred_Apr14.csv")
# MT.pred_older<-read.csv("MT_pred_Mar22.csv")

TD.lognorm<-read.xlsx("~/Google Drive/Kirk/Code/Updated data/TD_ParameterEstimates_2017Apr4_lognormal.xlsx", sheetIndex=1)
TD.lognorm[is.na(TD.lognorm[,'col'])==TRUE,'col']<-colors()[90]
TD.lognorm[TD.lognorm[,'col']==3,'col']<-colors()[139]

parnames.fig<-c(expression(mu), expression(beta[1]), expression(beta[2]), expression(italic(r)), expression(italic(K)), expression(alpha))


ylims<-list(c(0, 0.35), c(0, 10), c(0,10), c(-0.5, 5), c(0, 200), c(0, 0.001))
# quartz(width=6, height=4)
par(mfrow=c(2,3), mar=c(3,3,2,0), oma=c(1,2,0,1))
for(i in 1:6){
	plotCI(dat$temp, TD.lognorm[((i-1)*9+1):(i*9),'mean'], li=TD.lognorm[((i-1)*9+1):(i*9),'li'], ui=TD.lognorm[((i-1)*9+1):(i*9),'ui'], xlab="", ylab="", las=1, bty="l", main=parnames.fig[i], pch=21, pt.bg="white", gap=0, ylim=ylims[[i]], xlim=c(3,38), col=TD.lognorm[((i-1)*9+1):(i*9),'col'])
	lines(rownames(MT.pred), MT.pred[,i])
	lines(rownames(MT.pred), MT.pred_old[,i+1], lty=3, lwd=2)
	#if(i==4) abline(h=0, lty=2)
	if(i==6) legend("topleft", lty=c(1,2,3), lwd=c(1,2), c(expression(paste(E[K], "~ Norm (Jul 8)", sep="")), expression(paste(E[K], "~ logNorm (Apr 14)", sep=""))), bg="white")
	}
mtext(side=1, outer=TRUE, expression(paste("Temperature (", degree, 'C)')))
mtext(side=2, outer=TRUE, "Parameter estimate", line=0.5)

# #---------------------
# # Table of parameter estimates
# #---------------------
for(i in 1:20) parEst.byClone[[i]]<-cbind(parEst.byClone[[i]], rescaled_var=parEst.byClone[[i]][,'var']*parEst.byClone[[i]][,'n.clones'])

library(xlsx)

write.xlsx(parEst.byClone[[1]], file="MT_ParameterEstimates_2017Jul8.xlsx", sheetName=parnames[1])
for(i in 2:20) write.xlsx(parEst.byClone[[i]], file="MT_ParameterEstimates_2017Jul8.xlsx", sheetName=parnames[i], append=TRUE)

k<-25
BestMTEst<-data.frame(mean=parEst.all[[k]][,'mean'], SE=sqrt(parEst.all[[k]][,'variance']*k), rhat=round(Rhat[[k]][[1]][match(rownames(parEst.all[[k]]), rownames(Rhat[[k]][[1]])),1], 2), estimable=c(rep(1,8), 0, rep(1,2), 0, rep(1,3), 0, 1, 0, rep(1,2)))


write.csv(BestMTEst, file="MT_ParameterEstimates_25clones_8Jul2017.csv")

