# rm(list=ls())
# load('~/Google Drive/Kirk/Code/Updated data/Final figs etc/Jul8-EKNormal/MT_Jul8.RData')

# fit25<-fit.chains.in[[25]]

# length(fit25)

# fit25R<-fit25[[1]]
# for(j in 2:length(fit25)) fit25R<-rbind(fit25R, fit25[[j]])
	
# # 20 parameters. How "lognormal" do they look?  How do we calculate 95% CI?
# par(mfrow=c(5,4), mar=c(3,1,2,1))
# for(i in 1:20){
	# plot(density(fit25R[,i]), main=colnames(fit25[[1]])[i], bty="n", yaxt="n", xlab="", ylab="")
	# abline(v=mean(fit25R[,i]), lty=2)
# }

# # Monte Carlo simulation for confidence interval on output
# BestMTEst


set.seed(2389)
iter<-1000
MCMT<-matrix(rep(NA, iter*6*401));dim(MCMT)<-c(iter, 401, 6)

par.draw<-function(){
	par.i<-rnorm(20, mean=BestMTEst[,'mean'], sd=BestMTEst[,'SE'])
	par.i<-c(par.i[1:2], 5*par.i[2], par.i[3:20])
	names(par.i)<-as.character(priors$parameter)
	return(par.i)
	}

for(i in 1:iter){
	
	par.i<-par.draw()
	while(length(which(par.i[c(1:4,7:9,11:13,15,17,18,21)]<0))>0){
		par.i<-par.draw()
	}

	MCMT[i,,]<-MT.parPred(par.i)
	}

MCMT.CI<-list(); length(MCMT.CI)<-2
for(i in 1:6) MCMT.CI[[i]]<-apply(MCMT[,,i], MARGIN=2, quantile, c(0.025, 0.975))

