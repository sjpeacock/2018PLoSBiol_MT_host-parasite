# rm(list=ls())	
library(deSolve)
library(dclone)
source("models.R")

######################################################################################################
# Read in updated data
######################################################################################################

data<-read.csv("motality_data.csv")

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

######################################################################################################
######################################################################################################
######################################################################################################
# MT model
######################################################################################################
######################################################################################################
######################################################################################################

#------------------------------------------------------
# Add priors to data list
priors<-read.csv("MT_priors.csv")
dat$prior.mean<-priors$mean; names(dat$prior.mean)<-priors$parameter
dat$prior.sd<-priors$sd; names(dat$prior.sd)<-priors$parameter

#------------------------------------------------------
# Add Boltzman constant to data list
dat$k <- 8.62e-05

#------------------------------------------------------

n.chains<-10 # Up from 8 for 1-4 clones
parnames.MT<-priors$parameter

# JAGS requires that initial conditions be in list format, so convert:
inits.MT<-list(); length(inits.MT)<-n.chains
set.seed(894)
for(i in 1:n.chains){
	
	inits.MT[[i]]<-list(); length(inits.MT[[i]])<-23
	# add names, including random number generator initial conditions
	names(inits.MT[[i]])<-c(as.character(priors$parameter), '.RNG.name', '.RNG.seed') 
	
	for(j in 1:21){
		# For normal params (temperature)
		if(is.element(j, grep("T", priors$parameter, value=FALSE))==TRUE){ 
			inits.MT[[i]][[j]]<-rnorm(1, mean=priors$mean[j], sd=priors$sd[j])
		# For lognormal params (all others)
		}else if(priors$parameter[j]=="E_K"){
			inits.MT[[i]][[j]]<-rnorm(1, mean=priors$mean[j], sd=priors$sd[j])
		}else{
			inits.MT[[i]][[j]]<-rlnorm(1, mean=priors$mean[j], sd=priors$sd[j])
			}
	}
	
	# Set random number generator and seed for each chain
	inits.MT[[i]][22]<-"base::Super-Duper"
	inits.MT[[i]][23]<-c(73,38,37,48,34,1235,135,68,9023,645)[i]
	}

#------------------------------------------------------
# Clone data and model in loop
n.clones.all<-c(1:25)
fit.MT<-list(); length(fit.MT)<-length(n.clones.all)
t1.all<-numeric(length(n.clones.all))

for(k in 1:length(n.clones.all)){
	datk<-dat
	if(n.clones.all[k]>1){
		# Multiply:
		datk$n<-dat$n*n.clones.all[k]
		
		# Clone:
		datk$y <- rep(dat$y, n.clones.all[k])
		datk$ind.death <- rep(dat$ind.death, n.clones.all[k])
		datk$ind.death.back <- rep(dat$ind.death.back, n.clones.all[k])
		datk$temp.ind <- rep(dat$temp.ind, n.clones.all[k])
		datk$infected <- rep(dat$infected, n.clones.all[k])
		datk$p <- rep(dat$p, n.clones.all[k])
		
		# Unchanged: "nt", "h", "temp", "alpha.prior", "mu.prior", "r.prior", "beta.prior", "K.prior", "k"
		}
	
	#------------------------------------------------------
	# Fit model
	
	t0<-proc.time()[3] # Start timer to know how long things take
	
	# Fit Bayesian model ---------------------------
	cl <- makePSOCKcluster(n.chains) # Create cluster to fit 3 chains in parallel
	fit.MT[[k]]<-jags.parfit(
		cl=cl, 
		data=datk, 
		#params=as.character(priors$parameter),
		params=as.character(priors$parameter)[which(as.character(priors$parameter)!="EL_mu")], 
		model=model.MT_fixedEL_mu, 
		n.chains=n.chains, 
		n.adapt=20000, 
		n.update=20000, 
		n.iter=2000) 
	
	# Close cluster
	stopCluster(cl)
	t1.all[k]<-(proc.time()[3]-t0)/60/60
	cat("\nTime to fit (hrs):", t1.all[k], "Clone: ", k)
	# cat("E_K =", summary(fit.MT[[k]])[[1]]['E_K',1]) #Keep track of E_K which was a problem parameter initially

}# end k clones
#------------------------------------------------------
# Save workspace
# save.image(file=paste("MTE_fit_", Sys.Date(), ".RData", sep=""))

