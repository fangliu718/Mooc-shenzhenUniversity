res=read.csv(paste0('resp.csv'))
res<-as.matrix(res[2:length(res)])
dim(res)
size = dim(res)
N = size[1]
J = size[2]


require(nimble)
joint <- nimbleCode({
  for(i in 1:N){
    for (j in 1:J){
      ind[i, j] <- a[j]*(theta[i] - b[j])
      p[i, j] <- exp(ind[i, j])/(1+exp(ind[i, j]))
      resp[i, j] ~ dbern(prob = p[i, j])
    }
  }
  for(i in 1:N){ theta[i] ~ dnorm(mean = 0, sd = 1) }
  for(j in 1:J){
    a[j] ~ dlnorm(meanlog = 0, sdlog = 1)
    b[j] ~ dnorm(mean = 0, var = 100)
  }
})
constants <- list(N = N, J = J)
jointModel <- nimbleModel(joint, constants = constants,check = TRUE)
data <- list(resp = res)
set.seed(1)
inits <- list(theta = rnorm(N),a = rep(1, J), b = rnorm(J))
jointModel$setData(data)
jointModel$setInits(inits)
samples <- nimbleMCMC(model = jointModel,
                        niter = 25000, nchains = 1, nburnin = 5000,thin=2,
                        monitors = c("a", "b","theta"))
 
require(coda)
label_a=paste0('a[',1:J,']')
label_b=paste0('b[',1:J,']')
label_theta=paste0('theta[',1:N,']')
bmcmc=mcmc(samples[,label_b])
amcmc=mcmc(samples[,label_a])
thetamcmc=mcmc(samples[,label_theta])
dev.new()
plot(bmcmc,smooth=TRUE)
plot(amcmc,smooth=TRUE)
plot(thetamcmc,smooth=TRUE)


eap = colMeans(samples)
sd = apply(samples, 2, sd)
require(boa)
hpd = apply(samples, 2, boa.hpd, alpha = 0.05)
hpd = t(hpd)
summ = matrix(, nrow = ncol(samples), ncol = 4)
rownames(summ) = colnames(samples)
colnames(summ) = c('eap', 'sd', 'Lower','Upper')
summ[, 'eap'] = eap
summ[, 'sd'] = sd
summ[, c('Lower','Upper')] = hpd
write.csv(summ, paste0('summ_aconly.csv'))
  
 



