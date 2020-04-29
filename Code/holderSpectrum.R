require('wmtsa')
require('plyr')
require('lattice')

#dev.new()

setwd("~/Dropbox/2014 - Memory")
ts1<-read.delim('ts1.txt',col.names="White")
ts2<-read.delim('ts2.txt',col.names="Pink")
ts3<-read.delim('ts3.txt',col.names="Brownian")
ts=list(cbind(ts1$White[1:1024]),cbind(ts2$Pink[1:1024]),cbind(ts3$Brownian[1:1024]))
y <- sapply(t(ts), function(t) (t-min(t))/max(t) )

s.n = 256
s.rng = c(1,s.n)


# layout(matrix(c(1,1,1,2,2,2,3,3,3), 1, 9, byrow = TRUE))
# layout.show()
# Brownian -----------------------------------------------------------------------------------------------
W <- wavCWT(y[,3], wavelet="gaussian1",scale.range=s.rng,n.scale=s.n)
W.tree <- wavCWTTree(W,type='extrema')

plot(W, grid.size=2^11, col = palette(gray(seq(0.1,.9,len=s.n))),ylim=c(0.1,9.5),xlim=c(-10,1030),yaxt="n",bty="n",ylab="Scale")#,main="Brownian noise")
lines(1:1024, log1p(y[,3]*256), lwd=3,type="l",col="gray",yaxt="n",bty="n")
plot(W.tree,pch=20,ylim=c(0.1,9.5),xlim=c(-10,1030),add=T,yaxt="n",bty="n")

# Pink ---------------------------------------------------------------------------------------------------
W <- wavCWT(y[,2], wavelet="gaussian1",scale.range=s.rng,n.scale=s.n)
W.tree <- wavCWTTree(W,type='extrema')

plot(W, grid.size=2^11, col = palette(gray(seq(0.1,.9,len=s.n))),ylim=c(0.1,9.5),xlim=c(-10,1030),yaxt="n",bty="n",ylab="Scale",main="A Pink Noise Analogon")
lines(1:1024, log1p(y[,2]*256), lwd=3,type="l",col="gray",yaxt="n",bty="n")
plot(W.tree,pch=20,ylim=c(0.1,9.5),xlim=c(-10,1030),add=T,yaxt="n",bty="n")

# White ---------------------------------------------------------------------------------------------------
W <- wavCWT(y[,1], wavelet="gaussian1",scale.range=s.rng,n.scale=s.n)
W.tree <- wavCWTTree(W,type='extrema')

plot(W, grid.size=2^11, col = palette(gray(seq(0.1,.9,len=s.n))),ylim=c(0.1,9.5),xlim=c(-10,1030),yaxt="n",bty="n",ylab="Scale",main="White noise")
lines(1:1024, log1p(y[,1]*256), lwd=3,type="l",col="gray",yaxt="n",bty="n")
plot(W.tree,pch=20,ylim=c(0.1,9.5),xlim=c(-10,1030),add=T,yaxt="n",bty="n")


# dat <- ldply(llply(W.tree,function(Wt) sapply(Wt,cbind)))
#
# xyplot(scale[iscale] ~ time[itime],data=dat)
#
#
# levelplot(dat$extrema~dat$time[dat$itime]*dat$scale[dat$iscale],region=T)
#
# ds<-dat[[1]]
# ds$scale
# ## estimate Holder exponents
# holder <- holderSpectrum(W.tree)
# print(holder)
#
#
#
#
# nscale <- attributes(W)$n.scale
# n      <- attributes(W)$n.sample
# scales <- attributes(W)$scale
#
# maxmap <- matrix(ncol=nscale,nrow=n)
# t      <- 1:n
# tplus  <- c(t[n],1:n-1)
# tminus <- c(t[2:n],t[1])
# cwt   <- abs(W);
#
# # WMMT
# for(k in 1:nscale){
#   localmax    <-  (cwt[t,k] > cwt[tplus,k])&(cwt[t,k] > cwt[tminus,k])
#   maxmap[t,k] <- localmax[t]
# }
#
# # Partition function
# nq <- 61
# q  <- seq(-20,20,length=nq)
# z  <- matrix(ncol=nscale,nrow=nq)
#
# for(k in 1:nscale){
#   j <- which(maxmap[t,k])
#   if(length(j)>0){
#     C = abs(W[j,k])
#     for(i in 1:nq){
#       z[i,k] <- sum(C^q[i])
#     }
#   }
#   else
#   {
#     z[ ,k] = (1e-15)^q
#   }
# }
#
#
# # Moment generating function
# loscale <- 2
# hiscale <- 512
# scaled  <- seq(min(scales),nscale,length=nq)
#
# tau    <- matrix(ncol=1,nrow=nq)
# window <- (loscale <= scaled)&(scaled <= hiscale)
# ix     <- which(window)
# if(length(ix)>1){
#   for(kq in 1:nq){
#     y  <- log(z[kq,ix])
#     y  <- y-mean(y)
#     x  <- log(scaled[ix])
#     x  <- x-mean(x)
#     tau[kq] <- sum(y*x) / sum(x*x)
#   }
# }
#
# # q(h)
# alpha  <- seq(.05,.95,length=21)
# nalpha <- length(alpha);
# f      <- matrix(ncol=1,nrow=nalpha)
# for(kalpha in 1:nalpha){f[kalpha] <- min(alpha[kalpha]*q - tau)}
#
#
# # D(h)
# 	if(length(scaled)==1&scaled[1]==1){
# 		d = log(abs(tau))/(q-1)
# 	}	else {
# 	d = log(abs(tau))/((q-1)*log(scaled))
# }
#
#
# par(mfrow=c(1,3))
# plot(q,tau,xlab="q",ylab=expression(tau(q)))
# plot(q,d,xlab="q",ylab="D(q,s)")
# plot(alpha,f,xlab="h",ylab="f(h)")
# par(mfrow=c(1,3))

