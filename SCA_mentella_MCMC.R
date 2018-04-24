# SCA_mentella_MCMC
# Benjamin Planque April 2018

# devtools::install_github("kaskr/tmbstan/tmbstan") # may be required to install TMBSTAN
require(tmbstan)

load('SCA_mentella_model.Rdata')            # Load model results from TMB
model$fit <- tmbstan(model$obj,par = names(model$obj$par),init = list("last.par.best"),control = list(adapt_delta = 0.85,max_treedepth = 10),iter = 20000,warmup = 5000,thin = 10,chains = 1)
save(model,file='SCA_mentella_model_MCMC.Rdata') # save model object with the MCMC chains
save(model,file=paste('SCA_mentella_model_MCMC',model$date.flag,'.Rdata',sep="")) # save a copy with the date flag

# load(file='SCA_mentella_model_MCMC.Rdata') # load model object with the MCMC chains
# cnames = unique(names(model$fit))
cnames = unique(names(model$fit))[18:36]
cnames = unique(names(model$fit))[c(29,33,43,45,46)]
#pdf(paste('SCA_mentella_MCMC_',model$date.flag,'.pdf',sep=""))
plot(model$fit,plotfun = "hist",pars = cnames)
traceplot(model$fit,pars = cnames)
#dev.off()
