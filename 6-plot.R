# Analisys of calcium oscillations by Dan Bobkov, 2018

library(WaveletComp)
library(pdc)
library(tseriesEntropy)

# load data obtained from ImageJ

calcium <- read.csv2('20x20.csv')

# choose variable

xcal <- calcium$ROI416

# time series analysis

calcium_ts <- ts(xcal, start = calcium$Time[1])

del <- 0.1
x.spec <- spectrum(calcium_ts, log = "no", span = 10, plot = FALSE)
spx <- x.spec$freq/del
spy <- 2*x.spec$spec

# plot results

par(mfrow=c(2,3))

plot(xcal ~ calcium$Time, 
     xlab = "time (s)", ylab = "fluorescence intensity (a.u.)", t="l", 
     xlim = c(0, 110), ylim = c(0, 150), main = "Total signal")

plot(xcal ~ calcium$Time, 
     xlab = "time (s)", ylab = "fluorescence intensity (a.u.)", t="l", 
     xlim = c(0, 20), ylim = c(0, 150), main = "Signal segment")

boxplot(xcal, ylab = "fluorescence intensity (a.u.)",
        main = "Mean intensity", ylim = c(0, 150))

plot(spy ~ spx, xlab="frequency (Hz)", ylab="spectral density",
     type="l", lwd = 2)

acf(xcal, lag.max = 100, main = "Autocorrelation")

# Entropy measure

w <- as.integer(xcal)

Srho(w, lag.max = 10, stationary = TRUE, plot = TRUE,
     version = "FORTRAN")



