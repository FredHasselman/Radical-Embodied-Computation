library(plyr)
library(tidyverse)
library(invctr)
library(casnet)


# op<-par(xaxt="n",yaxt="n",bty = "n", mar = c(0,0,0,0))
# plot(1.0*sin(seq(0,4*pi,length.out = 1000)),xlab = "",ylab="",ylim = c(-1,1),type = "l",lwd=5)
# plot(.3*sin(seq(0,4*pi,length.out = 1000)),xlab = "",ylab="", ylim = c(-1,1),type = "l",lwd=5)
# plot(0.3*sin(seq(0,8*pi,length.out = 1000)),xlab = "",ylab="", ylim = c(-1,1),type = "l",lwd=5)
# plot(1.0*sin(seq(0,8*pi,length.out = 1000)),xlab = "",ylab="", ylim = c(-1,1),type = "l",lwd=5)
#
# plot(1.0*sin(seq(0,6*pi,length.out = 1000)),xlab = "",ylab="", ylim = c(-1,1),type = "l",lwd=5)
# par(op)

# s.n = 256
# s.rng = c(1,s.n)
#
# inDir <- "~/Documents/Projects/Self-Affine Resonator/2014 - Memory"
# ts1<-read.delim(file.path(inDir,'ts1.txt'),col.names="White")
# ts2<-read.delim(file.path(inDir,'ts2.txt'),col.names="Pink")
# ts3<-read.delim(file.path(inDir,'ts3.txt'),col.names="Brownian")
# ts=list(White=cbind(ts1$White[1:2047]),Pink=cbind(ts2$Pink[1:2047]),Brownian = cbind(ts3$Brownian[1:2047]))
# y <- sapply(t(ts), function(t) (t-min(t))/max(t) )


library(wmtsa)
library(lattice)
library(tuneR)
library(EMD)

inDir  <- "~/Documents/GitHub/tempscaleGraph/STIMULI"
normal <- file.path(inDir,paste0("BAKDAK", c(1,10),".WAV"))
# slowed <- file.path(inDir,paste0("BAKDAKV",c(1,10),".WAV"))
# amped  <- file.path(inDir,paste0("BAKDAKA",c(1,10),".WAV"))
# both   <- file.path(inDir,paste0("BAKDAKB",c(1,10),".WAV"))

ts <- list(bAk = normal[1], dAk = normal[2]) #,slowed,amped,both)

S1  <- llply(normal,function(s) readWave(s))
#S1d <- llply(S1,function(s) llply(s,downsample,samp.rate=4096))

specOut <- llply(S1,function(d) spectro(data = d, return_data = TRUE))

bAk <- specOut[[1]]$p




# s.n <- 512
# cnt <- 0
# tt  <- tsy <- yy <- list()
#
# for(w in 1:length(S1d)){
#   for(s in 1:length(S1d[[w]])){
#     cnt <- cnt+1
#     cat(cnt)
#     tt  <- cbind((1:length(S1d[[w]][[s]]@left))*(1/S1d[[w]][[s]]@samp.rate))
#      y  <- hilbertspec(xt=cbind(S1d[[w]][[s]]@left),tt=tt)
#    #plot(ts(y$instantfreq))
#      yy[[cnt]] <- cbind(y=S1d[[w]][[s]]@left, time=tt)
#      tsy[[cnt]] <- ts(data=abs(y$amplitude)/max(abs(y$amplitude)),start=0,frequency=S1d[[w]][[s]]@samp.rate)
#     # s.rng    <- deltat(tsy[[cnt]]) * c(1, length(tsy[[cnt]]))
#     # W[[cnt]]      <- wavCWT(tsy[[cnt]], wavelet="gaussian2",scale.range=s.rng,n.scale=s.n)
#     # W.tree[[cnt]] <- wavCWTTree(W[[cnt]],type='extrema')
#   }
# }


#y <- sapply(tt, function(t) (t-min(t))/max(t) )
y <- sapply(yy, function(t) (t[,1]-min(t[,1]))/max(t[,1]) )
N <- NROW(y)

# Number of frequencies estimated cannot be set! (defaults to Nyquist)
# Use Tukey window: cosine taper with r = 0.5

# fast = TRUE ensures padding with zeros to optimize FFT to highly composite number.
# However, we just pad to nextPow2, except if length already is a power of 2.
npad  <- (stats::nextn(N,factors=2)-N)
if(npad<2){
  npad <- stats::nextn(N,factors=2)
}

Tukey <- sapa::taper(type="raised cosine", flatness = 0.5, n.sample = N)
CAgridList <- gList <- list()




for(tn in 1:NCOL(y)){

  #psd  <- sapa::SDF(y[,tn], method = "multitaper", recenter = TRUE)
  psd  <- sapa::SDF(elascer(y[,tn]), taper. = Tukey)

  powspec <- cbind.data.frame(freq.norm = attr(psd, "frequency")[-1],  size = attr(psd, "frequency")[-1]*stats::frequency(y), bulk = as.matrix(psd)[-1])

  #unique(round(abs(sin(seq(0,2**pi,length.out = Nfreq))) * Nfreq)) #unique(round(abs(cos(seq(0,Nfreq*pi,length.out = Nfreq)+.5)) * Nfreq))

  Nfreq <- 512

  #freqVec <- rep(0,Nfreq)
  freqVec <- (powspec$size*N)[1:Nfreq]

  #ampVec <- rep(0,Nfreq)
  ampVec <- elascer(sqrt(powspec$bulk[1:Nfreq])) #abs(ifultools:: decibel(powspec$bulk[1:Nfreq]))

  #plot(log10(freqVec),log10(ampVec))

   # plot(ts(freqVec))
   # plot(ts(ampVec))

  # plot(log(freqVec),log(ampVec))

  Ntime <- Nfreq*4

  lb <- 1
  ub <- 1
  ampWinSize <- 3

  CAgrid      <- data.frame(matrix(data = 0, nrow = Nfreq, Ntime, dimnames = list(paste0("f",rev(1:Nfreq)), paste0("t",1:Ntime))))
  CAgrid$freq <- paste0(rownames(CAgrid))
  Input       <- CAgrid

  thresh <- .95
  ampWin  <- laply(1:ampWinSize,function(w) thresh * .5^w)
  ampWin <- c(rev(ampWin),thresh,ampWin)

  sinList <- list()

  for(f in 1:Nfreq){
    sinList[[f]] <- sin(seq(0,freqVec[f]*pi,length.out = Ntime))
    Input[(Nfreq-f+1),1:Ntime] <- ifelse(sinList[[f]]>min(ampWin)*(1-ampVec[f]),1,0)
  }

  Input_long         <- gather(Input, key = "time", value = "state",-freq)
  Input_long$freq_n  <- as.numeric(gsub("f","",Input_long$freq))
  Input_long$time_n  <- as.numeric(gsub("t","",Input_long$time))
  Input_long$state_f <- factor(Input_long$state)


  # ggplot(Input_long, aes(x = time_n, y = freq_n)) +
  #   geom_raster(aes(fill = state_f)) +
  #   scale_fill_manual(values = c("0"="white","1"="black")) +
  #   scale_x_continuous(expand = c(0,0)) +
  #   scale_y_continuous(expand = c(0,0)) +
  #   theme_bw()

   df_amps <- data.frame(t(ldply(sinList)))
  # plotTS_multi(df = df_amps)

  for(f in 1:Nfreq){
    for(t in 1:Ntime){
      if(Input[f,t]==1){
        if( f < lb){
          lbt <- f
        } else{
          lbt <- lb
        }
        if((f + ub) > Nfreq){
          ubt <- (Nfreq - f)
        } else {
          ubt <- ub
        }
        ifelse(sinList[[f]]>thresh*(1-ampVec[f]),1,0)
        CAgrid[c((f-lbt):(f+ubt)),t] <- 1 #rep(1,lb+ub+1) * c()
      }
    }
  }

  CAgridList[[tn]] <-  CAgrid

  CAgrid_long <- gather(CAgrid, key = "time", value = "state",-freq)
  CAgrid_long$freq_n <- as.numeric(gsub("f","",CAgrid_long$freq))
  CAgrid_long$time_n <- as.numeric(gsub("t","",CAgrid_long$time))
  CAgrid_long$state_f <- factor(CAgrid_long$state)

  gList[[tn]] <- ggplot(CAgrid_long, aes(x = time_n, y = freq_n)) +
    geom_raster(aes(fill = state_f)) +
    scale_fill_manual("State",values = c("0"="white","1"="black")) +
    scale_x_continuous("Time",expand = c(0,0)) +
    scale_y_reverse("Frequency / Scale",expand = c(0,0)) +
    ggtitle(label=names(ts)[[tn]]) +
    theme_bw()

  rm(CAgrid_long)

 ggsave(filename = paste0(names(ts)[[tn]],"_CA.tiff"), plot = gList[[tn]])
}

identical(CAgridList[[1]],CAgridList[[2]])
CAgridList[[1]]
CAgridList[[2]]

CAgrid1 <- gather(CAgridList[[1]], key = "time", value = "state",-freq)
CAgrid2 <- gather(CAgridList[[2]], key = "time", value = "state",-freq)

CAgrid <- CAgrid1
CAgrid$state <- CAgrid1$state-CAgrid2$state

CAgrid_long <-CAgrid  #gather(CAgrid, key = "time", value = "state",-freq)
CAgrid_long$freq_n <- as.numeric(gsub("f","",CAgrid_long$freq))
CAgrid_long$time_n <- as.numeric(gsub("t","",CAgrid_long$time))
CAgrid_long$state_f <- factor(CAgrid_long$state, levels = c(-1,0,1), labels = c("active for dAk","same", "active for bAk"))

ggplot(CAgrid_long, aes(x = time_n, y = freq_n)) +
  geom_raster(aes(fill = state_f)) +
  scale_fill_manual("State",values = c("active for dAk"="red3","same"="white","active for bAk"="blue")) +
  scale_x_continuous("Time",expand = c(0,0)) +
  scale_y_reverse("Frequency / Scale",expand = c(0,0)) +
  ggtitle(label=paste0(names(ts)[[1]],"-",names(ts)[[2]])) +
  theme_bw()
