setwd("~/Google Drive/Kirk/Manuscript/PLoS Biol")
setwd("2018PLoSBiol_MT_host-parasite")
rm(list=ls())	

library(deSolve)
library(dclone)
library(gplots)

# See http://journals.plos.org/plosbiology/s/figures
# Width: 789 – 2250 pixels (at 300 dpi). Height maximum: 2625 pixels (at 300 dpi)
# equals 2.63 - 7.5 in by 8.75 in
#################################################################################################
#################################################################################################
# MAIN TEXT FIGURES
#################################################################################################
#################################################################################################

#------------------------------------------------------------------------------------------------
# Figure 1: Data (host survival + parasites) for each temperature
#------------------------------------------------------------------------------------------------

# Read in mortality data AVAILABLE AT ???? --------------------------------

data<-read.csv("mortality_data.csv")

# Uninfected hosts have no parasites, but should not contribute to the parasite growth parameters since they were unexposed.  Therefore, the parasites for unexposed individuals are NA
data$total_spvs[data$treatment==0]<-NA

# Create data list for JAGS model fitting
n<-dim(data)[1]
dat<-list(
	n=n, 						# Number of individuals for which we have time of death
	nt=286, 					# Number of timesteps to simulation (max number of days in data, 0:285)
	y=rep(1, n), 				# Dummy vector of 1s for Bernoulli likelihood
	ind.death=data$death+1, 	# Index of day found dead (day[1] = 0)
	ind.death.back=data$death, 	# Day last seen alive (surveyed daily)
	h=1, 						# Timestep for Euler (doesn't change)
	temp.ind=as.numeric(as.factor(data$temp)), # Index for temp treatment
	temp=c(5.979,9.463,11.801,16.179,20.107,24.286,27.387,29.668,33.258), # true temperatures of treatments
	infected=as.numeric(data$infection_status=="Infected"), # indicator of which individuals are infected
	p=data$total_spvs 			# Parasite burden at death
	)

# Create list of the proportion of individuals surviving from days 1 to 285 for each temp
# and infection status
S.obs<-list(); length(S.obs)<-9*2; dim(S.obs)<-c(9,2)
for(i in 1:9){ # for each temp
	for(j in 1:2){ #for uninfected and infected
		mij<-subset(data, data$temp==unique(data$temp)[i]&data$treatment==(j-1))
		if(dim(mij)[1]>0){
			S<-rep(dim(mij)[1], 285)
			for(z in 1:(dim(mij)[1])){
				S[mij$death[z]:285]<-S[mij$death[z]:285]-1
				}
			S.obs[[i,j]]<-S
			rm(S)
		}
	}
}

# Parameters under TD and MT models --------------------------------
T.ind<-match(round(dat$temp, 1), round(T.all, 1))

params.both<-list(
	mu=cbind(TD=params.TD[[1]][,'mu'], MT=MT.pred[T.ind,'mu']), 
	beta1=cbind(TD=params.TD[[1]][,'beta1'], MT=MT.pred[T.ind,'beta1']),
	beta2=cbind(TD=params.TD[[1]][,'beta2'], MT=MT.pred[T.ind,'beta2']),
	r=cbind(TD=params.TD[[1]][,'r'], MT=MT.pred[T.ind,'r']),
	K=cbind(TD=params.TD[[1]][,'K'], MT=MT.pred[T.ind,'K']),
	alpha=cbind(TD=params.TD[[1]][,'alpha'], MT=MT.pred[T.ind,'alpha']))

params.both$K[1,2]<-10^-10 # Very small predicted K of 10^-86 gives errors.

# Predictions -------------------------------------------------------
deriv<-function(t,y,parms){
	dy<-c(
		-parms['beta1']*parms['mu']^parms['beta1']*t^(parms['beta1']-1)*y[1],
		-parms['beta2']*(parms['mu']+parms['alpha']*y[3])^parms['beta2']*t^(parms['beta2']-1)*y[2],
		y[3]*parms['r']*(1-y[3]/parms['K'])
		)
	return(list(dy))
}


pred<-list(); length(pred)<-9*2; dim(pred)<-c(9,2)
for(i in 1:9){
	for(j in 1:2){
		parms.ij<-c(mu=params.both[['mu']][i,j], beta1=params.both[['beta1']][i,j], beta2=params.both[['beta2']][i,j], r=params.both[['r']][i,j], K=params.both[['K']][i,j], alpha=params.both[['alpha']][i,j])
		names(parms.ij)<-names(params.both)		
		# if(i==1&j==2) y0<-c(1,1,0) else 
		y0<-c(1,1,1)
		pred[[i,j]]<-lsoda(y=y0, times=seq(0,285,1), func=deriv, parms=parms.ij)
	}
}

#----------------------------------------------
# Plots
#----------------------------------------------

stretch<-TRUE

xmax<-rep(285, 9)
if(stretch==TRUE) for(i in 1:9) xmax[i]<-which((S.obs[[i,1]]+S.obs[[i,2]])==0)[1]

# Note: there is no FE model predictions for exposed hosts at 6*C, 9.5*C, 29.7*C, and 33.3*C due to estimability problems with parasite-related parameters.

# Plot HOST data and predictions -------------------------------------------------------
# quartz(width = 7.5, height=6.5)
par(mfrow=c(3,3), mar=c(2,2,2,1), oma=c(2,2.5,4.5,1), mgp=c(2.5, 0.5, 0), tck=-0.02)
for(i in 1:9){
	plot(0:284, S.obs[[i,1]]/S.obs[[i,1]][1], "l", col="#0000FF40", bty="l", las=1, ylim=c(0,1), xaxt="n", lwd=3, xlim=c(0, xmax[i]))
	
	if(stretch==TRUE) axis(side=1) else{
		if(i>=7) axis(side=1, at=seq(0, 285, 75), tck=-0.04) else axis(side=1, at=seq(0, 285, 75), tck=-0.04, labels=FALSE)}
		
	lines(0:284, S.obs[[i,2]]/S.obs[[i,2]][1], col="#FF000040", lwd=3)
	
	for(k in 1:2){ #for FE (aka TD) and MTE models
		lines(pred[[i,k]][,1], pred[[i,k]][,2], col=4, lty=k+1)
		if((i>2&i<8)|k==2) lines(pred[[i,k]][,1], pred[[i,k]][,3], col=2, lty=k+1)
		}
	
	mtext(side=3, substitute(paste(b, degree, "C"), list(b=sprintf("%.1f", round(dat$temp[i], 1)))))

	if(i==2){
		legend(-100, 1.6, fill=c(4,2), c("Unexposed hosts", "Exposed hosts"), bty="n", xpd=NA, border=NA)
		legend(200, 1.65, lty=c(1,2,3), lwd=c(3,1,1), col=c(grey(0.8), 1, 1), c("Data", "Fixed-temperature model*", "MTE model"), bty="n", xpd=NA, border=NA)
		}
}
mtext(side=1, outer=TRUE, "Days", line=1)
mtext(side=2, outer=TRUE, line=1, "Proportion surviving")

# Plot PARASITE data and predictions -------------------------------------------------------
# quartz(width = 7.5, height=6.5)
par(mfrow=c(3,3), mar=c(2,2,2,1), oma=c(2,2.5,4.5,1), mgp=c(2.5, 0.5, 0), tck=-0.02)
for(i in 1:9){
	plot(dat$ind.death[dat$temp.ind==i&dat$infected==1]-1, dat$p[dat$temp.ind==i&dat$infected==1], "n", xlim=c(0, xmax[i]), bty="l", xaxt="n", ylim=c(0, max(1, dat$p[dat$temp.ind==i&dat$infected==1], na.rm=TRUE)), las=1)
	
	if(stretch==TRUE) axis(side=1) else{
		if(i>=7) axis(side=1, at=seq(0, 285, 75), tck=-0.04) else axis(side=1, at=seq(0, 285, 75), tck=-0.04, labels=FALSE)}
		
	points(dat$ind.death[dat$temp.ind==i&dat$infected==1&dat$p==0]-1, dat$p[dat$temp.ind==i&dat$infected==1&dat$p==0], pch=4, col="#00000080")
	points(dat$ind.death[dat$temp.ind==i&dat$infected==1&dat$p!=0]-1, dat$p[dat$temp.ind==i&dat$infected==1&dat$p!=0], cex=0.8, col="#00000080")

	
	for(k in 1:2){ #for FE (aka TD) and MTE models
		if((i>2&i<8)|k==2) lines(pred[[i,k]][,1], pred[[i,k]][,4], col=1, lty=k+1) 
		}
	
	mtext(side=3, substitute(paste(b, degree, "C"), list(b=sprintf("%.1f", round(dat$temp[i], 1)))))
	
	if(i==2){
		legend(-200, 1.65, col="#00000080", pch=c(4,1), c("Zero parasites at TOD of exposed hosts", "Non-zero parasites at TOD of exposed hosts"), bty="n", xpd=NA, title="Data")
		legend(200, 1.65, lty=c(2,3), lwd=c(1,1), col=c(1, 1), c("Fixed-temperature model*", "MTE model"), bty="n", xpd=NA, border=NA, title="Model predictions")
		}

}
mtext(side=1, outer=TRUE, "Days", line=1)
mtext(side=2, outer=TRUE, line=1, "Parasite burden")

#------------------------------------------------------------------------------------------------
# Figure 2: TD parameters (estimable only) with MT model predictions
#------------------------------------------------------------------------------------------------

load('workspaces/DT_fits_4Apr2017.RData')
load("workspaces/MTE_fits_Jul8.RData")
source("CI_calc.R") # source R ode that produces MCMC Conf Int on the MT predicitions
T.all<-seq(0,40,0.1)

TD<-list(exp(params.TDLOG[[1]]), exp(params.TDLOG[[2]]), exp(params.TDLOG[[3]])) # These have 95% CI calculated on the log scale and then transformed back so that it's positive.  
# TD<-params.TD #use this if you want to plot the 95% CI calculated as mean +/- 1.96*SE on the scale of the parameter, but note that the lower bound will be negative in some cases.

# quartz(width=7.5, height=4.5, pointsize=12)
tiff(filename="Fig2.tiff", width=2250, height=1350, res=300)
ylims<-list(c(0, 0.35), c(0, 10), c(0,10), c(0, 0.8), c(0, 200), c(0, 0.0008))

par(mfrow=c(2,3), mar=c(1,3,2,0), oma=c(3,2,0,1))

for(i in 1:6){
	
	ind<-which(point.col[,i]!=3);ind2<-which(point.col[,i]==3)
	colT<-rep("white", 9)
	coli<-point.col[,i]; coli[ind]<-NA
	
	plotCI(c(5.979,9.463,11.801,16.179,20.107,24.286,27.387,29.668,33.258), TD[[1]][,i], type="n", lwd=NA,li=TD[[2]][,i], ui=TD[[3]][,i], xlab="", ylab="", las=1, bty="l", ylim=ylims[[i]], xlim=c(3,38), yaxt="n", xaxt="n")
	if(i>3) axis(side=1, at=seq(10,30,10), las=1) else axis(side=1, at=seq(10,30,10), las=1, labels=FALSE)
	axis(side=1, at=seq(5,35,5), tck=-0.02, labels=FALSE)
	
	if(i==1){
		axis(side=2, at=c(0, 0.1, 0.2, 0.3), las=1)
		axis(side=2, at=seq(0, 0.35, 0.05), tck=-0.02, labels=FALSE)
		mtext(side=3, adj=0, expression(paste("  a) ", mu, sep="")))
	}
	
	if(i==2) mtext(side=3, adj=0, expression(paste("  b) ", beta[1], sep="")))
	
	if(i==3) mtext(side=3, adj=0, expression(paste("  c) ", beta[2], sep="")))
	
	if(i==2|i==3){
		axis(side=2, at=seq(0, 10, 2), las=1)
		axis(side=2, at=seq(0, 10, 1), tck=-0.02, labels=FALSE)
		abline(h=1, lty=2)
		}
	
	if(i==4){
		axis(side=2, at=seq(0, 1, 0.2), las=1)
		axis(side=2, at=seq(0, 1, 0.1), tck=-0.02, labels=FALSE)
		mtext(side=3, adj=0, expression(paste("  d) ", italic(r), sep="")))
		arrows(27.387, ylims[[i]]*0.99, 27.387, 1.1*ylims[[i]], length=0.08, xpd=NA)
	}
	
	if(i==5){
		axis(side=2, at=seq(0, 200, 50), las=1)
		axis(side=2, at=seq(0, 200, 25), tck=-0.02, labels=FALSE)
		mtext(side=3, adj=0, expression(paste("  e) ", italic(theta), sep="")))
	}
	
	if(i==6){
		axis(side=2, at=seq(0, 0.0008, 0.0002), las=1, labels=seq(0,8,2))
		axis(side=2, at=seq(0, 0.0008, 0.0001), tck=-0.02, labels=FALSE)
		mtext(side=3, adj=0, expression(paste("  f) ", alpha%*%10^5), sep=""))
		arrows(27.387, ylims[[i]]*0.99, 27.387, 1.1*ylims[[i]], length=0.08, xpd=NA)
	}
	
	lines(rownames(MT.pred), MT.pred[,i])
	polygon(x=c(T.all, rev(T.all)), y=c(MCMT.CI[[i]][1,], rev(MCMT.CI[[i]][2,])), col="#00000030", border=NA)
	plotCI(c(5.979,9.463,11.801,16.179,20.107,24.286,27.387,29.668,33.258)[ind2], TD[[1]][ind2,i], li=TD[[2]][ind2,i], ui=TD[[3]][ind2,i], pch=21, gap=0, pt.bg=colT[ind2], add=TRUE, cex=1.3)
	
	if(i==5) mtext(side=1, expression(paste("Temperature (", degree, 'C)')), line=3)
	
}

mtext(side=2, outer=TRUE, "Parameter estimate", line=0.5)


dev.off()

#--------------------
# FIgure S3: Including non-estimable parasite parameters:
tiff(filename="FigS3.tiff",  width=2250, height=800, res=300)
ylims<-list(c(0, max(TD[[3]][,1])), c(0, max(TD[[3]][,2])), c(0,max(TD[[3]][,3])), c(0, 6), c(0, 200), c(0, 0.03))
par(mfrow=c(1,3), mar=c(1,3,2,0), oma=c(3,2,2,1))

for(i in 4:6){
	
	ind<-NA; ind2<-1:9
	colT<-point.col[,i]
	colT[which(colT=="gold")]<-"red"
	colT[which(colT==2)]<-"gold"
	
	coli<-point.col[,i]; coli[ind]<-NA
	
	plotCI(c(5.979,9.463,11.801,16.179,20.107,24.286,27.387,29.668,33.258), TD[[1]][,i], type="n", lwd=NA,li=TD[[2]][,i], ui=TD[[3]][,i], xlab="", ylab="", las=1, bty="l", ylim=ylims[[i]], xlim=c(3,38), yaxt="n", xaxt="n")
	if(i>3) axis(side=1, at=seq(10,30,10), las=1) else axis(side=1, at=seq(10,30,10), las=1, labels=FALSE)
	axis(side=1, at=seq(5,35,5), tck=-0.02, labels=FALSE)
	
	if(i==4){
		axis(side=2, at=seq(0, 6, 1), las=1)
		axis(side=2, at=seq(0, 6, 0.5), tck=-0.02, labels=FALSE)
		mtext(side=3, adj=0, expression(paste("  a) ", italic(r), sep="")), line=1)
		arrows(c(5.979, 33.258), rep(ylims[[i]][2]*0.99, 2), c(5.979, 33.258), rep(1.1*ylims[[i]][2], 2), length=0.08, xpd=NA)
	}
	
	if(i==5){
		axis(side=2, at=seq(0, 200, 50), las=1)
		axis(side=2, at=seq(0, 200, 25), tck=-0.02, labels=FALSE)
		mtext(side=3, adj=0, expression(paste("  b) ", italic(theta), sep="")), line=1)
		arrows(x0=c(9.463, 29.668), rep(ylims[[i]][2]*0.99, 2), c(9.463, 29.668), rep(1.1*ylims[[i]][2], 2), length=0.08, xpd=NA)
	}
	
	if(i==6){
		axis(side=2, at=seq(0, 0.03, 0.005), las=1, labels=seq(0,30,5))
		axis(side=2, at=seq(0, 0.03, 0.001), tck=-0.02, labels=FALSE)
		mtext(side=3, adj=0, expression(paste("  c) ", alpha%*%10^3), sep=""), line=1)
		arrows(c(9.463,29.668,33.258) , rep(ylims[[i]][2]*0.99, 3), c(9.463,29.668,33.258), rep(1.1*ylims[[i]][2], 3), length=0.08, xpd=NA)
	}
	
	lines(rownames(MT.pred), MT.pred[,i])
	polygon(x=c(T.all, rev(T.all)), y=c(MCMT.CI[[i]][1,], rev(MCMT.CI[[i]][2,])), col="#00000030", border=NA)
	plotCI(c(5.979,9.463,11.801,16.179,20.107,24.286,27.387,29.668,33.258)[ind2], TD[[1]][ind2,i], li=TD[[2]][ind2,i], ui=TD[[3]][ind2,i], pch=21, gap=0, pt.bg=colT[ind2], add=TRUE, cex=1.3)
	
	if(i==5) mtext(side=1, expression(paste("Temperature (", degree, 'C)')), line=3)
	
}

mtext(side=2, outer=TRUE, "Parameter estimate", line=0.5)
dev.off()

#------------------------------------------------------------------------------------------------
# Figure 3: Predictions of average lifespan over temperature from MT and TD estimates
#------------------------------------------------------------------------------------------------


#################################################################################################
#################################################################################################
# SUPPORTING INFORMATION
#################################################################################################
#################################################################################################

#------------------------------------------------------------------------------------------------
# Figure S1: example of modified SS
#------------------------------------------------------------------------------------------------
temp<-seq(0, 40, 0.1)
E<-0.65
EH<-5*0.65
EL<-5*0.65
TH<-30
TL<-12
x_0<-1
k<-8.62e-05

SSU.dag<-x_0*exp(-E/k*(1/(temp+273.15)-1/(15+273.15)))*(1+exp(EH/k*(-1/(temp+273.15)+1/(TH+273.15))))
SSU<-x_0*exp(-E/k*(1/(temp+273.15)-1/(15+273.15)))/(1+exp(EH/k*(-1/(temp+273.15)+1/(TH+273.15))))

SS.dag<-x_0*exp(-E/k*(1/(temp+273.15)-1/(15+273.15)))*(1+exp(EL/k*(1/(temp+273.15)-1/(TL+273.15)))+exp(EH/k*(-1/(temp+273.15)+1/(TH+273.15))))
SS<-x_0*exp(-E/k*(1/(temp+273.15)-1/(15+273.15)))/(1+exp(EL/k*(1/(temp+273.15)-1/(TL+273.15)))+exp(EH/k*(-1/(temp+273.15)+1/(TH+273.15))))

quartz(width=6.8, height=2.5, pointsize=10)
par(mfrow=c(1,2), mar=c(4,1,1,0), oma=c(0,2,1,1))
plot(temp, SS, "l", yaxt="n", xlab="", ylab="", bty="l", ylim=range(SS, SSU))
a<-par('usr'); arrows(a[1], a[3], a[1], a[4], length=0.08, xpd=NA)
# lines(temp, SS, lty=2)
mtext(side=3, adj=0, "a) Original Sharpe-Schoolfied", line=0.5)
# legend('topleft', lty=c(1,2), c("SS.U", "SS"), bty="n")
mtext(side=2,  "Rate or parameter value", line=1)

plot(temp, SS.dag, "l", yaxt="n", xlab="", ylab="", bty="l")
a<-par('usr'); arrows(a[1], a[3], a[1], a[4], length=0.08, xpd=NA)
# lines(temp, SS.dag, lty=2)
mtext(side=3, adj=0, "b) Modified Sharpe-Schoolfied (\u2020)", line=0.5)

# legend('topleft', lty=c(1,2), c("SS.U", "SS"), bty="n")

mtext(side=1, outer=TRUE, expression(paste("Temperature (", degree, "C)")), line=-1)


#------------------------------------------------------------------------------------------------
# Fig S2: TD estimates; estimability
#------------------------------------------------------------------------------------------------

load('workspaces/DT_fits_4Apr2017.RData')

parnames.fig<-c(expression(mu), expression(beta[1]), expression(beta[2]), expression(italic(r)), expression(theta), expression(alpha))


# Scaled variance over number of clones - all temps plotted together
colT<-colorRampPalette(c("skyblue", 4,  "gold", 2, 1))(n=9)

# pdf("Figures/T1_paramsScaledVar_4Apr2017_lognormal.pdf", onefile=TRUE, width=6, height=5)
# quartz(width=7.7, height=4.8)
par(mfrow=c(2,3), mar=c(2,2,1,1), oma=c(2,2,1,6))
for(m in 1:6){ #For each parameter
	plot(parEst.byClone[[1]]$n.clones, 1/parEst.byClone[[1]]$n.clones, type="n", bty="l", xlim=c(0.5,max(n.clones.all)+0.5), ylim=c(0, c(1, 1, 1, 3, 1.6, 4)[m]), las=1)
	lines(parEst.byClone[[1]]$n.clones, 1/parEst.byClone[[1]]$n.clones, lwd=2, col=grey(0.8))
	for(j in 1:9){ #for each temperature
		i<-(m-1)*9+j
		lines(parEst.byClone[[i]]$n.clones, parEst.byClone[[i]]$var/parEst.byClone[[i]]$var[1], lty=3, col=colT[j])
		points(parEst.byClone[[i]]$n.clones, parEst.byClone[[i]]$var/parEst.byClone[[i]]$var[1], pch=c(1,19)[as.numeric(parEst.byClone[[i]]$r.hat<1.1)+1], col=colT[j])
		}
		
		# if(i>18) axis(side=1, at=parEst.byClone[[i]]$n.clones) else axis(side=1, at=parEst.byClone[[i]]$n.clones, labels=FALSE)
		# if(is.element(i, seq(1,21,3))) axis(side=2, at=c(0,1)) else axis(side=2, at=c(0,1), labels=FALSE)
		# mtext(side=3, paste(parnames.all[i], "-", round(dat$temp[j]), "*C"), cex=par('cex'), font=2)

	if(m==3){
		legend(17.5,0.7, pch=19, col=colT, legend=round(dat$temp, 1), lty=3, bty="n", title=expression(paste("Temp (", degree, "C)")), xpd=NA)	
		legend(14.5,1, pch=c(19,1,NA), lty=c(3,3, 1), lwd=c(1,1, 2), col=c(1,1,grey(0.8)), legend=c("Converged", "Not converged", "Ideal"), bty="n", xpd=NA)	
		}
	
	mtext(side=3, line=-1, parnames.fig[m])
		
}

mtext(side=1, outer=TRUE, "Number of clones (K)", line=1)
mtext(side=2, outer=TRUE, "Scaled variance in posterior", line=1)
# # dev.off()

#------------------------------------------------------------------------------------------------
# Fig S3: MT estimates: mean +/- SE
#------------------------------------------------------------------------------------------------
load("workspaces/MTE_fits_Jul8.RData")

parnamesMT.fig<-expression(mu[0], italic(E[mu]), italic(E[H * mu]), italic(T[L * 
    mu]), italic(T[H * mu]), beta[0], italic(E[beta]), italic(E[H * 
    beta]), italic(T[H * beta]), italic(r[0]), italic(E[r]), 
    italic(E[H * r]), italic(T[H * r]), italic(theta[0]), italic(E[theta]), 
    italic(E[L * theta]), italic(E[H * theta]), italic(T[L * theta]), italic(T[H * 
        theta]), alpha[0])

# pdf("FigS3_MTparamsSE_2017Jul8.pdf", width=8, height=9)
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
# dev.off()

#------------------------------------------------------------------------------------------------
# Fig S4: MT estimates: scaled variance over number of clones
#------------------------------------------------------------------------------------------------

# pdf("FigS4_MT_ScaledVar_2017Jul8.pdf", width=8, height=9)
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
# dev.off()

