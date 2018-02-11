# Function to estimate mode (does R not have one??)
estimate_mode <- function(s) {
  d <- density(s)
  d$x[which.max(d$y)]
}

#---------------------------------------------------------------------------------------------
# Get estimates from MT model fit

HPDest<-function(fit, leave.out=NA, par.order=NA){
	
	if(is.na(leave.out[1])) chains.in<-1:length(fit) else chains.in<-which(is.na(match(1:length(fit), leave.out))==TRUE)
	
	n2<-dim(fit[[1]])[2] 				# number of parameters
	
	params<-matrix(NA, nrow=n2, ncol=3)
	for(i in 1:n2){
		x<-numeric()
		for(j in chains.in) x<-c(x, fit[[j]][,i])
		
		params[i,1]<-estimate_mode(x)
		params[i,2:3]<-HPDinterval(as.mcmc(x))[1,]
	}
	
	rownames(params)<-colnames(fit[[1]])
	colnames(params)<-c("HPD", "2.5%", "97.5%")
	
	if(is.na(par.order[1])==FALSE) params<-params[match(par.order, rownames(params)),]
	
	return(params)	
	}


#---------------------------------------------------------------------------------------------
# Get estimates from MT model fit

MeanVarest<-function(fit, leave.out=NA, par.order=NA){
	
	if(is.na(leave.out[1])) chains.in<-1:length(fit) else chains.in<-which(is.na(match(1:length(fit), leave.out))==TRUE)
	
	n2<-dim(fit[[1]])[2] 				# number of parameters
	
	params<-matrix(NA, nrow=n2, ncol=2)
	for(i in 1:n2){
		x<-numeric()
		for(j in chains.in) x<-c(x, fit[[j]][,i])
		
		params[i,1]<-mean(x)
		params[i,2]<-var(x)
	}
	
	rownames(params)<-colnames(fit[[1]])
	colnames(params)<-c("mean", "variance")
	
	if(is.na(par.order[1])==FALSE) params<-params[match(par.order, rownames(params)),]
	
	return(params)	
	}


#---------------------------------------------------------------------------------------------
# MT predictions -----------------------------

MT.parPred<-function(par, T.all=seq(0,40,0.1)){
	MT.pred<-cbind(
		
		par['mu_0']*exp(-par['E_mu']/dat$k*(1/(T.all+273.15)-1/(15+273.15)))*(1+exp(par['EL_mu']/dat$k*(1/(T.all+273.15)-1/(par['TL_mu']+273.15)))+exp(par['EH_mu']/dat$k*(-1/(T.all+273.15)+1/(par['TH_mu']+273.15)))),
		
		par['beta_0']*exp(-par['E_beta']/dat$k*(1/(T.all+273.15)-1/(15+273.15)))/(1+exp(par['EH_beta']/dat$k*(-1/(T.all+273.15)+1/(par['TH_beta']+273.15)))),
		
		par['beta_0']*exp(-par['E_beta']/dat$k*(1/(T.all+273.15)-1/(15+273.15)))/(1+exp(par['EH_beta']/dat$k*(-1/(T.all+273.15)+1/(par['TH_beta']+273.15)))),
	
		par['r_0']*exp(-par['E_r']/dat$k*(1/(T.all+273.15)-1/(15+273.15)))/(1+exp(par['EH_r']/dat$k*(-1/(T.all+273.15)+1/(par['TH_r']+273.15)))), 
	
		par['K_0']*exp(-par['E_K']/dat$k*(1/(T.all+273.15)-1/(15+273.15)))/(1+exp(par['EL_K']/dat$k*(1/(T.all+273.15)-1/(par['TL_K']+273.15)))+exp(par['EH_K']/dat$k*(-1/(T.all+273.15)+1/(par['TH_K']+273.15)))), 
		
		rep(par['alpha_0'], length(T.all))
		)
		
		rownames(MT.pred)<-T.all
		colnames(MT.pred)<-c("mu", 'beta1', 'beta2', 'r', 'K', 'alpha')
		
		return(MT.pred)		
	}

#---------------------------------------------------------------------------------------------
# plot MT predictions -----------------------------

plot.MTpred<-function(MT.parPred, params.TD, T.all=seq(0,40,0.1), ymax){
	
	par.names<-c('mu', 'beta1', 'beta2', 'r', 'K', 'alpha')
	
	if(missing(ymax)) ymax<-apply(MT.parPred, 2, max, na.rm=TRUE)
	quartz(width=6, height=4)
	par(mfrow=c(2,3), mar=c(3,3,2,0), oma=c(1,2,0,1))
	for(i in 1:6){
		
		plotCI(dat$temp, params.TD[[1]][,i], type="n", liw=params.TD[[2]][,i], uiw=params.TD[[2]][,i], xlab="", ylab="", las=1, bty="l", main=par.names[i], barcol=NA, ylim=c(0, ymax[i]))
		
		lines(T.all, MT.parPred[,i], lwd=3, col=grey(0.7))
		
		plotCI(dat$temp, params.TD[[1]][,i], li=params.TD[[2]][,i], ui=params.TD[[3]][,i], gap=0, pch=21, add=TRUE, pt.bg=c('white'), cex=1.2)
		
		
	}
	mtext(side=1, outer=TRUE, "Temperature (*C)", cex=par('cex')*1.2, line=-0.5)
	mtext(side=2, outer=TRUE, "Parameter estimate", cex=par('cex')*1.2, line=0.5)
}

#-----------------------------------------------------------------------------------------------------
# Plot the trace plots for temp j and clone k from DT fits

plotJK<-function(j, k){# temp = j, clone = k
	par(mfrow=c(3,2), mar=c(1,2,1,1), oma=c(2,2,3,0))		
	for(i in 1:6){
		xi<-matrix(NA, nrow=length(thin.x), ncol=length(fit.1T[[j,k]]))
		for(m in 1:length(fit.1T[[j,k]])) xi[,m]<-fit.1T[[j,k]][[m]][thin.x, match(parnames[i], colnames(fit.1T[[j,k]][[1]]))]
			
			plot(thin.x, xi[,1], type='n', xaxt="n", bty="l", ylim=range(xi[,chains.in[[j,k]]]))
			for(m in chains.in[[j,k]]) lines(thin.x, xi[,m], col=col.chains[m])
			
			# mtext(side=3, paste(paste(k, "clones, "), sep=""))
		
			mtext(side=3, paste(parnames[i], sep=""))
			mtext(side=3, adj=0, line=-1, paste("   ", round(mean(xi[,chains.in[[j,k]]]), 6)))
			abline(h=mean(xi[,chains.in[[j,k]]]), lty=2)
			} # end i		
		
		 mtext(side=3, outer=TRUE, paste(k, "clones, ", j, "temp"), line=1)
		legend("right", fill=col.chains[chains.in[[j,k]]], bg="white", legend=chains.in[[j,k]], ncol=2, border=NA)
}
 

	