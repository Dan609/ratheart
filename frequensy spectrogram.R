# load data obtained from ImageJ

calcium <- read.csv2('20x20.csv') # 24x24 automatic ROI grid
calcium <- read.csv2('fluo4d.csv') # 40 cells are manually selected
calcium <- read.csv2('40cells.csv')

# spectrogram

par(mfrow=c(1,1), mar=c(5,5,5,5))

xcal <- calcium$ROI94

x <- ts(xcal, start = calcium$Time[1])

spectrogram(x, 14.29745, windowlength = 96000, 0.01,
            timestep = 2000,
            maxfreq = 0.5,
            colors = TRUE,
            dynamicrange = 55, nlevels = dynamicrange, maintitle = "", 
            show = TRUE, window = 'gaussian', windowparameter = 0.3, 
            quality = TRUE)


# spectral density

xcal <- calcium$ROI377

calcium_ts <- ts(xcal, start = calcium$Time[1])

del <- 0.1
x.spec <- spectrum(calcium_ts, log = "no", span = 10, plot = FALSE)
spx <- x.spec$freq/del
spy <- 2*x.spec$spec

plot(spy ~ spx, xlab="frequency (Hz)", ylab="spectral density",
     type="l", lwd = 2, ylim = c(0, 10000))


###

??spectrogram
?spectrum
?ts
??fpeaks
?fpeaks

update.packages()

library(sm)
library(oce)
library(pdc)
library(xts)
library(Rwave)
library(tsDyn)
library(stats)
library(arfima)
library(signal)
library(seewave)
library(deSolve)
library(ggplot2)
library(TSclust)
library(devtools)
library(biwavelet)
library(phonTools)
library(tseriesChaos)
library(tseriesEntropy)
library(nonlinearTseries)
library(fractal)
library(fracdiff)
library(fractaldim) 
library(MSMVSampEn)
library(WaveletComp)
library(scatterplot3d)
library(XML)
library(stringi)
library(quantmod)
library(pracma)


## Peaks detection, lenght = 109.6 s

s.rate <- 14.29745

z <- ts(calcium$ROI37, start = calcium$Time[1])

fpeaks(spec(z, f=s.rate), nmax = 5,
       threshold = 10,
       plot = T,
       title = TRUE,
       xlab = "Frequency (Hz)", ylab = "Amplitude",
       labels = TRUE, legend = TRUE, digits = 2)*1000

z <- ts(calcium$ROI377, start = calcium$Time[1])

spec(z, f=s.rate, wl = 512, wn = "hanning", fftw = FALSE, norm = TRUE,
     scaled = FALSE, PSD = FALSE, PMF = FALSE, correction="none", dB = NULL, dBref = NULL,
     at = NULL, from = NULL, to = NULL,
     identify = FALSE, col = "black", cex = 1,
     plot = 1, flab = "Frequency (kHz)",
     alab = "Amplitude", flim = NULL,
     alim = NULL, type="l")


localpeaks(spec(z), f=s.rate, bands = 10, mel = FALSE, plot = TRUE,
           xlab = NULL, ylab = "Amplitude", labels = TRUE)


findPeaks(z, thresh=0)

?findpeaks

findpeaks(z, nups = 1, ndowns = nups, zero = "0", 
          peakpat = NULL, minpeakheight = -Inf, 
          minpeakdistance = 1, 
          threshold = 0, npeaks = 0, 
          sortstr = FALSE)

findpeaks(calcium$ROI377, npeaks=3, threshold=4, sortstr=TRUE)

?spec

spec(z,f=14.29745,col="red",plot=2,flab="",yaxt="n")

spec(z,f=14.29745)

fpeaks(spec(z,f=14.29745), nmax = 50,
       threshold = 1,
       plot = T,
       title = TRUE,
       xlab = "Frequency (Hz)", ylab = "Amplitude",
       labels = TRUE, legend = TRUE, digits = 2)*1000
