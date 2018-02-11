######################################################################################################
######################################################################################################
######################################################################################################
# TD model
######################################################################################################
######################################################################################################
######################################################################################################

model.TD<-function(){
    
    #----------------------------------------------------------------------------------------
    # Priors and constraints 
    #----------------------------------------------------------------------------------------
   for(i in 1:9){ # For each temperature
		mu[i] ~ dlnorm(prior.mean[1], prior.sd[1]^(-2))
		beta1[i] ~ dlnorm(prior.mean[2], prior.sd[2]^(-2))
		beta2[i] ~ dlnorm(prior.mean[2], prior.sd[2]^(-2))
		r[i] ~ dlnorm(prior.mean[3], prior.sd[3]^(-2))
		K[i] ~ dlnorm(prior.mean[4], prior.sd[4]^(-2))
		alpha[i] ~ dlnorm(prior.mean[5], prior.sd[5]^(-2))
		} # end i
	
	#----------------------------------------------------------------------------------------
	# Survival probabilities
	# Simulate model using Euler method
	# Note: now we have infected data at all temps, which makes this simpler	
	#----------------------------------------------------------------------------------------
    
    for(i in 1:9){ # For each temperature
		
		# Intial conditions:
		U[i,1]<-1 	# Uninfected
		I[i,1]<-1	# Infected
		P[i,1]<-1	# Parasite burden		 	
		
		for(t in 2:nt){	# For each timestep
			
			# Changes from time t-1 to time t (timestep h = 1 day)
			dU[i,t-1] <- -beta1[i]*mu[i]^beta1[i]*(t-h)^(beta1[i]-1)*U[i,t-1]
			dI[i,t-1]<- -beta2[i]*(mu[i] + alpha[i]*P[i,t-1])^beta2[i]*(t-h)^(beta2[i]-1)*I[i,t-1]
			dP[i,t-1]<- r[i]*P[i,t-1]*(1-P[i,t-1]/K[i])		
			
			# Update variables (Euler, timestep h = 1)
			U[i,t] <- max(0, U[i,t-1] + dU[i,t-1]*h)
			I[i,t] <- max(0, I[i,t-1] + dI[i,t-1]*h)
			P[i,t] <- max(0, P[i,t-1] + dP[i,t-1]*h)	
			
			} #end time t
		}# end temp i
	
	

	#----------------------------------------------------------------------------------------
    # Probability 
	#----------------------------------------------------------------------------------------
    
    for(z in 1:n){ # For each individual
		prob[z]<-max(
			
			10^-10, # Avoid having a probability of zero, which give NAs in likelihood
			
			# if individual is not infected (if infected==1, this will = 0):
			(1-infected[z])*(-U[temp.ind[z], ind.death[z]] + U[temp.ind[z], ind.death.back[z]]), 
			
			# if individual is infected (if infected==0, this will = 0):
			(infected[z])*(-I[temp.ind[z], ind.death[z]] + I[temp.ind[z], ind.death.back[z]])) 
		} #end z
	
	#----------------------------------------------------------------------------------------
    # Likelihood 
	#----------------------------------------------------------------------------------------
    
    for (z in 1:n){ # For each individual
		y[z] ~ dbern(prob[z]) # Likelihood of surviving is prob. of Bernoulli==1 given modelled probability
		p[z] ~ dpois(P[temp.ind[z], ind.death[z]]) #Likelihood of parasites at time of death is given by P (unexposed individuals have p[z]==NA, which does not contribute to likelihood)		
		} # end z

} #end model



######################################################################################################
######################################################################################################
######################################################################################################
# MT model
######################################################################################################
######################################################################################################
######################################################################################################

model.MT<-function(){
    
    #----------------------------------------------------------------------------------------
    # Priors and constraints 
    #----------------------------------------------------------------------------------------
    # Note: JAGS takes precision instead of variance, where
    # precision = sd^(-2)
    
   # mu - SS
    mu_0 ~ dlnorm(prior.mean[1], prior.sd[1]^(-2))
	E_mu ~ dlnorm(prior.mean[2], prior.sd[2]^(-2))
	EL_mu ~ dlnorm(prior.mean[3], prior.sd[3]^(-2))
	EH_mu ~ dlnorm(prior.mean[4], prior.sd[4]^(-2)) 
	TL_mu ~ dnorm(prior.mean[5], prior.sd[5]^(-2))
	TH_mu ~ dnorm(prior.mean[6], prior.sd[6]^(-2))
	
   	# beta - SS.U
    beta_0 ~ dlnorm(prior.mean[7], prior.sd[7]^(-2))
	E_beta ~ dlnorm(prior.mean[8], prior.sd[8]^(-2))
	EH_beta ~ dlnorm(prior.mean[9], prior.sd[9]^(-2))
	TH_beta ~ dnorm(prior.mean[10], prior.sd[10]^(-2))
	
	# r = SS.U
    r_0 ~ dlnorm(prior.mean[11], prior.sd[11]^(-2))
	E_r ~ dlnorm(prior.mean[12], prior.sd[12]^(-2))
	EH_r ~ dlnorm(prior.mean[13], prior.sd[13]^(-2))
	TH_r ~ dnorm(prior.mean[14], prior.sd[14]^(-2))

    # K - SS
	K_0 ~ dlnorm(prior.mean[15], prior.sd[15]^(-2))
	E_K ~ dnorm(prior.mean[16], prior.sd[16]^(-2))
	EL_K ~ dlnorm(prior.mean[17], prior.sd[17]^(-2))
	EH_K ~ dlnorm(prior.mean[18], prior.sd[18]^(-2)) 
	TL_K ~ dnorm(prior.mean[19], prior.sd[19]^(-2))
	TH_K ~ dnorm(prior.mean[20], prior.sd[20]^(-2))
	
	# alpha - TD
    alpha_0 ~ dlnorm(prior.mean[21], prior.sd[21]^(-2))
	
    
	
	for(i in 1:9){ # For each temperature, calculate parameter according to MT relationship
		
		mu[i] <- mu_0*exp(-E_mu/k*(1/(temp[i]+273.15)-1/(15+273.15)))*(1+exp(EL_mu/k*(1/(temp[i]+273.15)-1/(TL_mu+273.15)))+exp(EH_mu/k*(-1/(temp[i]+273.15)+1/(TH_mu+273.15))))
		
		beta1[i] <- beta_0*exp(-E_beta/k*(1/(temp[i]+273.15)-1/(15+273.15)))/(1+exp(EH_beta/k*(-1/(temp[i]+273.15)+1/(TH_beta+273.15))))
		
		beta2[i] <- beta1[i] # Assume beta2 = beta1
		
		alpha[i]<-alpha_0 # Constant
		
		r[i]<-r_0*exp(-E_r/k*(1/(temp[i]+273.15)-1/(15+273.15)))/(1+exp(EH_r/k*(-1/(temp[i]+273.15)+1/(TH_r+273.15))))

		K[i] <- K_0*exp(-E_K/k*(1/(temp[i]+273.15)-1/(15+273.15)))/(1+exp(EL_K/k*(1/(temp[i]+273.15)-1/(TL_K+273.15)))+exp(EH_K/k*(-1/(temp[i]+273.15)+1/(TH_K+273.15))))
		}

	
 	#----------------------------------------------------------------------------------------
    # Survival probabilities	
	#----------------------------------------------------------------------------------------
    
     for(i in 1:9){ # For each temperature
		
		# Intial conditions:
		U[i,1]<-1 	# Uninfected
		I[i,1]<-1			# Infected
		P[i,1]<-1			# Parasite burden		 	
		
		for(t in 2:nt){	# For each timestep
			
			# Changes from time t-1 to time t (timestep h = 1 day)
			dU[i,t-1] <- -beta1[i]*mu[i]^beta1[i]*(t-h)^(beta1[i]-1)*U[i,t-1]
			dI[i,t-1]<- -beta2[i]*(mu[i] + alpha[i]*P[i,t-1])^beta2[i]*(t-h)^(beta2[i]-1)*I[i,t-1]
			dP[i,t-1]<- r[i]*P[i,t-1]*(1-P[i,t-1]/K[i])		
			
			# Update variables (Euler, timestep h = 1)
			U[i,t] <- max(0, U[i,t-1] + dU[i,t-1]*h)
			I[i,t] <- max(0, I[i,t-1] + dI[i,t-1]*h)
			P[i,t] <- max(0, P[i,t-1] + dP[i,t-1]*h)	
			
			} #end time t
		}# end temp i

	#----------------------------------------------------------------------------------------
    # Probability 
	#----------------------------------------------------------------------------------------
    
    for(z in 1:n){ # For each individual
		prob[z]<-max(
			
			10^-10, # Avoid having a probability of zero, which give NAs in likelihood
			
			# if individual is not infected (if infected==1, this will = 0):
			(1-infected[z])*(-U[temp.ind[z], ind.death[z]] + U[temp.ind[z], ind.death.back[z]]), 
			
			# if individual is infected (if infected==0, this will = 0):
			(infected[z])*(-I[temp.ind[z], ind.death[z]] + I[temp.ind[z], ind.death.back[z]])) 
		} #end z

	#----------------------------------------------------------------------------------------
    # Likelihood 
	#----------------------------------------------------------------------------------------
    
    for (z in 1:n){ # For each individual
		y[z] ~ dbern(prob[z]) # Likelihood of surviving is prob. of Bernoulli==1 given modelled probability
		p[z] ~ dpois(max(10^-10, P[temp.ind[z], ind.death[z]])) #Likelihood of parasites at time of death is given by P (unexposed individuals have p[z]==NA, which does not contribute to likelihood)	
		# Note that JAGS gives errors if dpois has an expectation of zero, hence the max() argument.	
		} # end z
	

} #end model


######################################################################################################
######################################################################################################
######################################################################################################
# Single temp model
######################################################################################################
######################################################################################################
######################################################################################################

model.1T<-function(){
    
    #----------------------------------------------------------------------------------------
    # Priors and constraints 
    #----------------------------------------------------------------------------------------
  	mu ~ dlnorm(prior.mean[1], prior.sd[1]^(-2))
	beta1 ~ dlnorm(prior.mean[2], prior.sd[2]^(-2))
	beta2 ~ dlnorm(prior.mean[2], prior.sd[2]^(-2))
	# r ~ dnorm(prior.mean[3], prior.sd[3]^(-2))
	r ~ dlnorm(prior.mean[3], prior.sd[3]^(-2))
	K ~ dlnorm(prior.mean[4], prior.sd[4]^(-2))
	alpha ~ dlnorm(prior.mean[5], prior.sd[5]^(-2))
	
	#----------------------------------------------------------------------------------------
	# Survival probabilities
	# Simulate model using Euler method
	# Note: now we have infected data at all temps, which makes this simpler	
	#----------------------------------------------------------------------------------------
    
   # Intial conditions:
	U[1]<-1 	# Uninfected
	I[1]<-1	# Infected
	P[1]<-1	# Parasite burden		 	
		
	for(t in 2:nt){	# For each timestep
		
		# Changes from time t-1 to time t (timestep h = 1 day)
		dU[t-1] <- -beta1*mu^beta1*(t-h)^(beta1-1)*U[t-1]
		dI[t-1]<- -beta2*(mu + alpha*P[t-1])^beta2*(t-h)^(beta2-1)*I[t-1]
		dP[t-1]<- r*P[t-1]*(1-P[t-1]/K)		
		
		# Update variables (Euler, timestep h = 1)
		U[t] <- max(0, U[t-1] + dU[t-1]*h)
		I[t] <- max(0, I[t-1] + dI[t-1]*h)
		P[t] <- max(0, P[t-1] + dP[t-1]*h)	
		
		} #end time t
	
	

	#----------------------------------------------------------------------------------------
    # Probability 
	#----------------------------------------------------------------------------------------
    
    for(z in 1:n){ # For each individual
		prob[z]<-max(
			
			10^-10, # Avoid having a probability of zero, which give NAs in likelihood
			
			# if individual is not infected (if infected==1, this will = 0):
			(1-infected[z])*(-U[ind.death[z]] + U[ind.death.back[z]]), 
			
			# if individual is infected (if infected==0, this will = 0):
			(infected[z])*(-I[ind.death[z]] + I[ind.death.back[z]])) 
		} #end z
	
	#----------------------------------------------------------------------------------------
    # Likelihood 
	#----------------------------------------------------------------------------------------
    
    for (z in 1:n){ # For each individual
		y[z] ~ dbern(prob[z]) # Likelihood of surviving is prob. of Bernoulli==1 given modelled probability
		p[z] ~ dpois(P[ind.death[z]]) #Likelihood of parasites at time of death is given by P (unexposed individuals have p[z]==NA, which does not contribute to likelihood)		
		} # end z

} #end model


######################################################################################################
######################################################################################################
######################################################################################################
# MT model - fixed mu lower activation energy
######################################################################################################
######################################################################################################
######################################################################################################

model.MT_fixedEL_mu<-function(){
    
    #----------------------------------------------------------------------------------------
    # Priors and constraints 
    #----------------------------------------------------------------------------------------
    # Note: JAGS takes precision instead of variance, where
    # precision = sd^(-2)
    
   # mu - SS
    mu_0 ~ dlnorm(prior.mean[1], prior.sd[1]^(-2))
	E_mu ~ dlnorm(prior.mean[2], prior.sd[2]^(-2))
	EL_mu <- 5*E_mu
	EH_mu ~ dlnorm(prior.mean[4], prior.sd[4]^(-2)) 
	TL_mu ~ dnorm(prior.mean[5], prior.sd[5]^(-2))
	TH_mu ~ dnorm(prior.mean[6], prior.sd[6]^(-2))
	
   	# beta - SS.U
    beta_0 ~ dlnorm(prior.mean[7], prior.sd[7]^(-2))
	E_beta ~ dlnorm(prior.mean[8], prior.sd[8]^(-2))
	EH_beta ~ dlnorm(prior.mean[9], prior.sd[9]^(-2))
	TH_beta ~ dnorm(prior.mean[10], prior.sd[10]^(-2))
	
	# r = SS.U
    r_0 ~ dlnorm(prior.mean[11], prior.sd[11]^(-2))
	E_r ~ dlnorm(prior.mean[12], prior.sd[12]^(-2))
	EH_r ~ dlnorm(prior.mean[13], prior.sd[13]^(-2))
	TH_r ~ dnorm(prior.mean[14], prior.sd[14]^(-2))

    # K - SS
	K_0 ~ dlnorm(prior.mean[15], prior.sd[15]^(-2))
	E_K ~ dnorm(prior.mean[16], prior.sd[16]^(-2))
	EL_K ~ dlnorm(prior.mean[17], prior.sd[17]^(-2))
	EH_K ~ dlnorm(prior.mean[18], prior.sd[18]^(-2)) 
	TL_K ~ dnorm(prior.mean[19], prior.sd[19]^(-2))
	TH_K ~ dnorm(prior.mean[20], prior.sd[20]^(-2))
	
	# alpha - TD
    alpha_0 ~ dlnorm(prior.mean[21], prior.sd[21]^(-2))
	
    
	
	for(i in 1:9){ # For each temperature, calculate parameter according to MT relationship
		
		mu[i] <- mu_0*exp(-E_mu/k*(1/(temp[i]+273.15)-1/(15+273.15)))*(1+exp(EL_mu/k*(1/(temp[i]+273.15)-1/(TL_mu+273.15)))+exp(EH_mu/k*(-1/(temp[i]+273.15)+1/(TH_mu+273.15))))
		
		beta1[i] <- beta_0*exp(-E_beta/k*(1/(temp[i]+273.15)-1/(15+273.15)))/(1+exp(EH_beta/k*(-1/(temp[i]+273.15)+1/(TH_beta+273.15))))
		
		beta2[i] <- beta1[i] # Assume beta2 = beta1
		
		alpha[i]<-alpha_0 # Constant
		
		r[i]<-r_0*exp(-E_r/k*(1/(temp[i]+273.15)-1/(15+273.15)))/(1+exp(EH_r/k*(-1/(temp[i]+273.15)+1/(TH_r+273.15))))

		K[i] <- K_0*exp(-E_K/k*(1/(temp[i]+273.15)-1/(15+273.15)))/(1+exp(EL_K/k*(1/(temp[i]+273.15)-1/(TL_K+273.15)))+exp(EH_K/k*(-1/(temp[i]+273.15)+1/(TH_K+273.15))))
		}

	
 	#----------------------------------------------------------------------------------------
    # Survival probabilities	
	#----------------------------------------------------------------------------------------
    
     for(i in 1:9){ # For each temperature
		
		# Intial conditions:
		U[i,1]<-1 	# Uninfected
		I[i,1]<-1			# Infected
		P[i,1]<-1			# Parasite burden		 	
		
		for(t in 2:nt){	# For each timestep
			
			# Changes from time t-1 to time t (timestep h = 1 day)
			dU[i,t-1] <- -beta1[i]*mu[i]^beta1[i]*(t-h)^(beta1[i]-1)*U[i,t-1]
			dI[i,t-1]<- -beta2[i]*(mu[i] + alpha[i]*P[i,t-1])^beta2[i]*(t-h)^(beta2[i]-1)*I[i,t-1]
			dP[i,t-1]<- r[i]*P[i,t-1]*(1-P[i,t-1]/K[i])		
			
			# Update variables (Euler, timestep h = 1)
			U[i,t] <- max(0, U[i,t-1] + dU[i,t-1]*h)
			I[i,t] <- max(0, I[i,t-1] + dI[i,t-1]*h)
			P[i,t] <- max(0, P[i,t-1] + dP[i,t-1]*h)	
			
			} #end time t
		}# end temp i

	#----------------------------------------------------------------------------------------
    # Probability 
	#----------------------------------------------------------------------------------------
    
    for(z in 1:n){ # For each individual
		prob[z]<-max(
			
			10^-10, # Avoid having a probability of zero, which give NAs in likelihood
			
			# if individual is not infected (if infected==1, this will = 0):
			(1-infected[z])*(-U[temp.ind[z], ind.death[z]] + U[temp.ind[z], ind.death.back[z]]), 
			
			# if individual is infected (if infected==0, this will = 0):
			(infected[z])*(-I[temp.ind[z], ind.death[z]] + I[temp.ind[z], ind.death.back[z]])) 
		} #end z

	#----------------------------------------------------------------------------------------
    # Likelihood 
	#----------------------------------------------------------------------------------------
    
    for (z in 1:n){ # For each individual
		y[z] ~ dbern(prob[z]) # Likelihood of surviving is prob. of Bernoulli==1 given modelled probability
		p[z] ~ dpois(max(10^-10, P[temp.ind[z], ind.death[z]])) #Likelihood of parasites at time of death is given by P (unexposed individuals have p[z]==NA, which does not contribute to likelihood)	
		# Note that JAGS gives errors if dpois has an expectation of zero, hence the max() argument.	
		} # end z
	

} #end model

