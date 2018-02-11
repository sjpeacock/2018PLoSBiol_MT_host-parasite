# setwd("~/Google Drive/Kirk/Manuscript/PLoS Biol/2018PLoSBiol_MT_host-parasite")
# rm(list=ls())	
library(deSolve)
library(dclone)
source("models.R")

######################################################################################################
# Read in updated data
######################################################################################################

data<-read.csv("mortality_data.csv")

# Uninfected hosts have no parasites, but should not contribute to the parasite growth parameters since they were unexposed.  Therefore, the parasites for unexposed individuals are NA
data$total_spvs[data$treatment==0]<-NA

#----------------------------------------
# Create data list for JAGS model fitting
#----------------------------------------

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


#------------------------------------------------------
# Add priors to data list
priors<-read.csv("DT_priors.csv")
dat$prior.mean<-priors$mean; names(dat$prior.mean)<-priors$parameter
dat$prior.sd<-priors$sd; names(dat$prior.sd)<-priors$parameter


#------------------------------------------------------

n.chains<-10 # Up from 8 for 1-4 clones
parnames<-c("mu", "beta1", "beta2", "r", "K", "alpha")
n.clones.all<-c(1:15)


####################################################################################################
# Loop through each temperature
####################################################################################################

	
fit.1T<-list(); length(fit.1T)<-length(n.clones.all)*9; dim(fit.1T)<-c(9, length(n.clones.all))
t1.all1T<-matrix(NA, nrow=9, ncol=length(n.clones.all))

for(temp in c(1:9)){
	
	#------------------------------------------------------
	# Subset temperature
	dat.temp<-dat
	dat.temp$y<-dat$y[dat$temp.ind==temp]
	dat.temp$ind.death<-dat$ind.death[dat$temp.ind==temp]
	dat.temp$ind.death.back<-dat$ind.death.back[dat$temp.ind==temp]
	dat.temp$temp.ind<-dat$temp.ind[dat$temp.ind==temp]
	dat.temp$infected<-dat$infected[dat$temp.ind==temp]
	dat.temp$p<-dat$p[dat$temp.ind==temp]
	dat.temp$n<-length(dat.temp$y[dat$temp.ind==temp])

	# JAGS requires that initial conditions be in list format, so convert:
	inits.TD<-list(); length(inits.TD)<-n.chains
	set.seed(c(894, 4395, 198237, 382, 23, 19034, 2198, 54, 80)[temp])
	for(i in 1:n.chains){
		
		inits.TD[[i]]<-list(); length(inits.TD[[i]])<-length(parnames)+2
		# add names, including random number generator initial conditions
		names(inits.TD[[i]])<-c(parnames, '.RNG.name', '.RNG.seed') 
		
		for(j in 1:length(parnames)){
			inits.TD[[i]][[j]]<-rlnorm(1, mean=priors$mean[c(1,2,2,3,4,5)[j]], sd=priors$sd[c(1,2,2,3,4,5)[j]])
			}
		
		# Set random number generator and seed for each chain
		inits.TD[[i]][length(parnames)+1]<-"base::Super-Duper"
		inits.TD[[i]][length(parnames)+2]<-c(73,35,657,48,3897,1235,135,53,9023,645)[i]
		}
	
	#------------------------------------------------------
	# Clone data and model in loop
	for(k in 1:length(n.clones.all)){
		datk<-dat.temp
		if(n.clones.all[k]>1){
			# Multiply:
			datk$n<-dat.temp$n*n.clones.all[k]
			
			# Clone:
			datk$y <- rep(dat.temp$y, n.clones.all[k])
			datk$ind.death <- rep(dat.temp$ind.death, n.clones.all[k])
			datk$ind.death.back <- rep(dat.temp$ind.death.back, n.clones.all[k])
			datk$temp.ind <- rep(dat.temp$temp.ind, n.clones.all[k])
			datk$infected <- rep(dat.temp$infected, n.clones.all[k])
			datk$p <- rep(dat.temp$p, n.clones.all[k])
			
			# Unchanged: "nt", "h", "temp", "alpha.prior", "mu.prior", "r.prior", "beta.prior", "K.prior", "k"
			}
	
		#------------------------------------------------------
		# Fit model
		
		t0<-proc.time()[3] # Start timer to know how long things take
		
		# Fit Bayesian model ---------------------------
		cl <- makePSOCKcluster(n.chains) # Create cluster to fit chains in parallel
		fit.1T[[temp,k]]<-jags.parfit(
			cl=cl, 
			data=datk, 
			params=c("mu", "beta1", "beta2", "r", "K", "alpha"), 
			model=model.1T, 
			n.chains=n.chains, 
			n.adapt=30000, 
			n.update=30000, 
			n.iter=2000) 
		
		# Close cluster
		stopCluster(cl)
		t1.all1T[temp,k]<-(proc.time()[3]-t0)/60
		cat("\nTime to fit (mins):", t1.all1T[temp,k], " temp", temp, ", clone: ", k) 
	
	}# end k clones
} # end temp (9) temperatures
#------------------------------------------------------
# Save workspace
# save.image(file=paste("DTfit_", Sys.Date(), "_rLOGNormal.RData", sep=""))

