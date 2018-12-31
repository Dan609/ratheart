# Analisys of calcium oscillations by Dan Bobkov, 2018
# load data obtained from ImageJ

calcium <- read.csv2('20x20.csv')
calcium <- read.csv2('fluo4d.csv') # 40 cells are manually selected

# choose ROI

xcal <- calcium$ROI481

# plot results

par(mfrow=c(2,3))

# time series analysis


calcium_ts <- ts(xcal, start = calcium$Time[1])

del <- 0.1
x.spec <- spectrum(calcium_ts, log = "no", span = 10, plot = FALSE)
spx <- x.spec$freq/del
spy <- 2*x.spec$spec

plot(spy ~ spx,
     xlab="frequency (Hz)",
     ylab="spectral density",
     type="l", main = "Spectrum")

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

# Entropy measure

w <- as.integer(xcal)  
Srho(w, lag.max = 10, stationary = TRUE, plot = TRUE,
     version = "FORTRAN")

## time series signal

plot(xcal ~ calcium$Time, 
     xlab = "time (s)", ylab = "fluorescence intensity (a.u.)", t="l", 
     xlim = c(0, 20), ylim = c(0, 100), main = "Time-series")

# spectrogram

x <- ts(xcal, start = calcium$Time[1])

spectrogram(x, 14.29745, windowlength = 96000, 0.01,
            timestep = 5000,
            maxfreq = 0.3,
            colors = TRUE,
            dynamicrange = 50, nlevels = dynamicrange, maintitle = "", 
            show = TRUE, window = 'gaussian', windowparameter = 0.3, 
            quality = TRUE)

# A series with a variable period

wt.image(my.w, n.levels = 250,
         legend.params = list(lab = "wavelet power levels", label.digits = 2),
         color.palette = "gray((n.levels):1/n.levels)",
         col.ridge = "blue",
         main = "Wavelet power spectrum")

# The reconstructed series
## my.rec <- reconstruct(my.w)
## x.rec <- my.rec$series$x.r

# Phase difference

x <- ts(calcium$ROI481, start = calcium$Time[1])
y <- ts(calcium$ROI555, start = calcium$Time[1])

my.data <- data.frame(x = x, y = y)

my.wc <- analyze.coherency(my.data, my.pair = c("x", "y"),
                           loess.span = 0,
                           dt = 1, 
                           dj = 1/100,
                           lowerPeriod = 4, upperPeriod = 64,
                           make.pval = TRUE, n.sim = 10)

wc.phasediff.image(my.wc, which.contour = "wp", use.sAngle = TRUE,
                   n.levels = 100, siglvl = 0.1,
                   legend.params = list(lab = "phase difference levels",
                                        lab.line = 3),
                   color.palette = "gray((n.levels):1/n.levels)",
                   timelab = "",
                   max.contour.segments = 2000)





