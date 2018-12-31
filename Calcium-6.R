# Analisys of calcium oscillations by Dan Bobkov, 2018

library(WaveletComp)
library(pdc)

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
     xlim = c(0, 20), ylim = c(0, 150), main = "Time-series")
plot(spy ~ spx, xlab="frequency (Hz)", ylab="spectral density", type="l", main = "Spectrum")
boxplot(xcal, ylab = "fluorescence intensity (a.u.)", main = "Mean intensity", ylim = c(0, 150))
acf(xcal, lag.max = 100, main = "Autocorrelation")
plot(entropyHeuristic(calcium_ts), ylim = c(0.6, 1), main = "Entropy")

# The wavelet transform of x is computed as follows:

x <- ts(xcal, start = calcium$Time[1])
my.data <- data.frame(x = x)
my.w <- analyze.wavelet(my.data, "x",
                        loess.span = 0,
                        dt = 1, 
                        dj = 1/250,
                        lowerPeriod = 2,
                        upperPeriod = 128,
                        make.pval = TRUE, n.sim = 10)

# A series with a variable period

wt.image(my.w, n.levels = 250,
         legend.params = list(lab = "wavelet power levels", label.digits = 2),
         color.palette = "gray((n.levels):1/n.levels)",
         col.ridge = "blue",
         main = "Wavelet power spectrum")



