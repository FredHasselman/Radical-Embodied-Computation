library(signal, warn.conflicts = F, quietly = T) # signal processing functions
library(oce, warn.conflicts = F, quietly = T) # image plotting functions and nice color maps
library(tuneR, warn.conflicts = F, quietly = T) # nice functions for reading and manipulating .wav files
#
# # # define path to audio file
# # fin = '~/Data/online_sound_examples/right_whale_upcall1.wav'
#
# inDir  <- "~/Documents/GitHub/tempscaleGraph/STIMULI"
# normal <- file.path(inDir,paste0("BAKDAK", c(1,10),".WAV"))
#
# fin <- normal[1]
#
# # read in audio file
# data = readWave(fin)
#
# # extract signal
# snd = data@left
#
# # determine duration
# dur = length(snd)/data@samp.rate
# dur # seconds
# ## [1] 3.588
#
# # determine sample rate
# fs = data@samp.rate
# fs # Hz
# ## [1] 2000
#
# # demean to remove DC offset
# snd = snd - mean(snd)
#
# # plot waveform
# plot(snd, type = 'l', xlab = 'Samples', ylab = 'Amplitude')
#
# # number of points to use for the fft
# nfft=1024
#
# # window size (in points)
# window=256
#
# # overlap (in points)
# overlap=128
#
# # create spectrogram
# spec = specgram(x = snd,
#                 n = nfft,
#                 Fs = fs,
#                 window = window,
#                 overlap = overlap
# )
#
# # discard phase information
# P = abs(spec$S)
#
# # normalize
# P = P/max(P)
#
# # convert to dB
# P = 10*log10(P)
#
# # config time axis
# t = spec$t
#
# # plot spectrogram
# imagep(x = t,
#        y = spec$f,
#        z = t(P),
#        col = oce.colorsViridis,
#        ylab = 'Frequency [Hz]',
#        xlab = 'Time [s]',
#        drawPalette = T,
#        decimate = F
# )
#
#
#
# specOut <- spectro(data = data, return_data = TRUE)
#

spectro = function(data, nfft=1024, window=256, overlap=128, t0=0, plot_spec = T, normalize = F, return_data = F,...){

  library(signal)
  library(oce)

  # extract signal
  snd = data@left

  # demean to remove DC offset
  snd = snd-mean(snd)

  # determine duration
  dur = length(snd)/data@samp.rate

  # create spectrogram
  spec = specgram(x = snd,
                  n = nfft,
                  Fs = data@samp.rate,
                  window = window,
                  overlap = overlap
  )

  # discard phase info
  P = abs(spec$S)

  # normalize
  if(normalize){
    P = P/max(P)
  }

  # convert to dB
  P = 10*log10(P)

  # config time axis
  if(t0==0){
    t = as.numeric(spec$t)
  } else {
    t = as.POSIXct(spec$t, origin = t0)
  }

  # rename freq
  f = spec$f

  if(plot_spec){

    # change plot colour defaults
    par(bg = "black")
    par(col.lab="white")
    par(col.axis="white")
    par(col.main="white")

    # plot spectrogram
    imagep(t,f, t(P), col = oce.colorsViridis, drawPalette = T,
           ylab = 'Frequency [Hz]', axes = F,...)

    box(col = 'white')
    axis(2, labels = T, col = 'white')

    # add x axis
    if(t0==0){

      axis(1, labels = T, col = 'white')

    }else{

      axis.POSIXct(seq.POSIXt(t0, t0+dur, 10), side = 1, format = '%H:%M:%S', col = 'white', las = 1)
      mtext(paste0(format(t0, '%B %d, %Y')), side = 1, adj = 0, line = 2, col = 'white')

    }
  }

  if(return_data){

    # prep output
    spec = list(
      t = t,
      f = f,
      p = t(P)
    )

    return(spec)
  }
}



AffineResonator <- function(inputGrid, amplitudes = rep(1,NROW(inputGrid)), frequencies = 1/(1:NCOL(inputGrid)),
                            timeLabels = 1:NCOL(inputGrid), freqLabels = 1:NROW(inputGrid), freqWindow = c(1,1),
                            ampThreshold = .99, ampWindow = 3){


  #Ntime   <- length(frequencies)*4
  Ntime  <- NROW(inputGrid)
  Nfreqs <- NCOL(inputGrid)

  lb <- freqWindow[1]
  ub <- freqWindow[2]

  freqVec <- frequencies
  ampVec  <- elascer(amplitudes) #abs(ifultools:: decibel(powspec$bulk[1:Nfreq]))

  CAgrid      <- data.frame(matrix(data = 0, nrow = Nfreq, ncol = Ntime, dimnames = list(freqLabels, timeLabels)))
  CAgrid$freq <- freqLabels
  Input       <- CAgrid

  ampWin  <- laply(1:ampWindow, function(w) ampThreshold*.5^w)
  ampWin <- c(rev(ampWin),ampThreshold, ampWin)

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


