library(sm)
library(oce)
library(pdc)
library(Rwave)
library(tsDyn)
library(stats)
library(fftw)
library(deSolve)
library(biwavelet)
library(tseriesChaos)
library(tseriesEntropy)
library(nonlinearTseries)
library(fractal)
library(fracdiff)
library(fractaldim)
library(arfima)
library(seewave)
library(phonTools)
library(signal)
library(ggplot2)
library(TSclust)
library(devtools)
library(MSMVSampEn)
library(WaveletComp)
library(scatterplot3d)

# load data obtained from ImageJ

calcium <- read.csv2('20x20.csv')

head(calcium)

# insert ROI 57 cells

x <- ts(calcium$ROI409, start = calcium$Time[1])
y <- ts(calcium$ROI415, start = calcium$Time[1])
z <- ts(calcium$ROI333, start = calcium$Time[1])

# choose variable

xcal <- calcium$ROI416

# time series analysis

calcium_ts <- ts(xcal, start = calcium$Time[1])

del <- 0.1
x.spec <- spectrum(calcium_ts, log = "no", span = 10, plot = FALSE)
spx <- x.spec$freq/del
spy <- 2*x.spec$spec

?spectrum

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

# Entropy tests of serial and cross dependence for time series
## WARNING: computationally intensive, increase B with caution!

Srho.test.ts.p(z, B = 2, 
               lag.max = 10,
               ci.type = c("perm"),
               quant = 0.95, plot = TRUE)


Srho.ts(z, lag.max = 10, method = "integral", plot = TRUE,
        maxpts = 0, tol = 0.001)


w <- as.integer(calcium$ROI222)

ent <- Srho(w, lag.max = 10, stationary = TRUE,
            version = "FORTRAN")
ent


MSMVSampEn(mat = z, 
           M = 5, tau = 10, r = 1, eps = 1567, scaleMat = T)

print(which.max(spx))

print(which.max(spy))

print(which.max(x.spec))

print(x.spec)

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

# The wavelet transform of x is computed as follows:

par(mfrow=c(1,1))
my.data <- data.frame(x = x)
my.w <- analyze.wavelet(my.data, "x",
                        loess.span = 0,
                        dt = 1, 
                        dj = 1/250,
                        lowerPeriod = 2,
                        upperPeriod = 128,
                        make.pval = TRUE, n.sim = 10)

# Plot the wavelet power spectrum, a series with a constant period

wt.image(my.w, color.key = "quantile", n.levels = 250, 
         legend.params = list(lab = "wavelet power levels", mar = 4.7))

# A series with a variable period
par(mfrow=c(2,1))
wt.image(my.w, n.levels = 250,
         legend.params = list(lab = "wavelet power levels"))

# The reconstructed series

my.rec <- reconstruct(my.w)
x.rec <- my.rec$series$x.r   #x: name of original series

# A series with two periods

my.data <- data.frame(x = x)
my.w <- analyze.wavelet(my.data, "x",
                        loess.span = 0,
                        dt = 1, 
                        dj = 1/250,
                        lowerPeriod = 2,
                        upperPeriod = 128,
                        make.pval = TRUE, n.sim = 10)
wt.image(my.w, n.levels = 250,
         legend.params = list(lab = "wavelet power levels"))

# Reconstruction, using only period 8

reconstruct(my.w, sel.period = 8, plot.waves = TRUE, lwd = c(1,2), legend.coords = "bottomleft")

reconstruct(my.w, sel.period = 16, plot.waves = TRUE, lwd = c(1,2), legend.coords = "bottomleft")

# actually uses the period closest to 10

my.w$Period[(my.w$Period > 9) & (my.w$Period < 11)]

# Compare two series with average powers calling:

x <- ts(calcium$ROI415, start = calcium$Time[1])
y <- ts(calcium$ROI416, start = calcium$Time[1])

my.data <- data.frame(x = x, y = y)
my.wx <- analyze.wavelet(my.data, "x", loess.span = 0,
                         dt = 1, dj = 1/20,
                         lowerPeriod = 2, upperPeriod = 64,
                         make.pval = TRUE,n.sim = 10)
my.wy <- analyze.wavelet(my.data, "y", loess.span = 0,
                         dt = 1, dj = 1/20,
                         lowerPeriod = 2, upperPeriod = 64,
                         make.pval = TRUE,n.sim = 10)

maximum.level = 1.001*max(my.wx$Power.avg, my.wy$Power.avg)
par(mfrow=c(1,2))
wt.avg(my.wx, maximum.level = maximum.level)
wt.avg(my.wy, maximum.level = maximum.level)

# White noise method

par(mfrow=c(1,1))

my.data <- data.frame(x = x)
my.w <- analyze.wavelet(my.data, "x",
                        method = "white.noise",
                        loess.span = 0,
                        dt = 1, 
                        dj = 1/250,
                        lowerPeriod = 2,
                        upperPeriod = 256,
                        make.pval = TRUE, n.sim = 10)
wt.image(my.w, color.key = "interval", n.levels = 250,
         legend.params = list(lab = "wavelet power levels"))

# Fourier randomization method

my.data <- data.frame(x = x)
my.w <- analyze.wavelet(my.data, "x",
                        method = "Fourier.rand",
                        loess.span = 0,
                        dt = 1, 
                        dj = 1/250,
                        lowerPeriod = 2,
                        upperPeriod = 256,
                        make.pval = TRUE, n.sim = 10)
wt.image(my.w, color.key = "interval", n.levels = 250,
         legend.params = list(lab = "wavelet power levels"))

# Plotting the power spectrum

par(mfrow=c(1,1))

wt.image(my.w, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", label.digits = 2))


# Grayscale

wt.image(my.w, n.levels = 250,
         legend.params = list(lab = "wavelet power levels", label.digits = 2),
         color.palette = "gray((n.levels):1/n.levels)",
         col.ridge = "blue")

# Time axis

my.data <- data.frame(x = x)

my.w <- analyze.wavelet(my.data, "x",
                        method = "Fourier.rand",
                        loess.span = 0,
                        dt = 1, 
                        dj = 1/250,
                        lowerPeriod = 2,
                        upperPeriod = 256,
                        make.pval = TRUE, n.sim = 10)

wt.image(my.w, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", label.digits = 2),
         timelab = "time")

# Analysis of bivariate time series

x <- ts(calcium$ROI415, start = calcium$Time[1])
y <- ts(calcium$ROI416, start = calcium$Time[1])

par(mfrow=c(1,1))

my.data <- data.frame(x = x, y = y)

my.wc <- analyze.coherency(my.data, my.pair = c("x", "y"),
                           loess.span = 0,
                           dt = 1, 
                           dj = 1/100,
                           lowerPeriod = 4, upperPeriod = 32,
                           make.pval = TRUE, n.sim = 10)


wc.image(my.wc, n.levels = 250, color.key = "interval",
         siglvl.contour = 0.1, siglvl.arrow = 0.05, which.arrow.sig = "wt",
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = "")

######


wc.image(my.wc, n.levels = 250,
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = "", periodlab = "period")


# A plot of the time-averaged cross-wavelet power

wc.avg(my.wc, siglvl = 0.01, sigcol = "red", sigpch = 20,
       periodlab = "period")
?wc.avg
# Coherence

wc.image(my.wc, which.image = "wc", color.key = "interval", n.levels = 250,
         siglvl.contour = 0.1, siglvl.arrow = 0.05,
         legend.params = list(lab = "wavelet coherence levels"),
         timelab = "")


my.wc2 <- analyze.coherency(my.data, my.pair = c("x", "y"),
                            loess.span = 0, dt = 1, dj = 1/100,
                            lowerPeriod = 2, upperPeriod = 32,
                            window.type.t = 1, window.type.s = 1,
                            window.size.t = 5, window.size.s = 1,
                            make.pval = TRUE, n.sim = 10)

wc.image(my.wc2, which.image = "wc", color.key = "interval", n.levels = 250,
         siglvl.contour = 0.1, siglvl.arrow = 0.05,
         legend.params = list(lab = "wavelet coherence levels"),
         timelab = "")

my.wc3 <- analyze.coherency(my.data, my.pair = c("x", "y"),
                            loess.span = 0, dt = 1, dj = 1/100,
                            lowerPeriod = 2, upperPeriod = 64,
                            window.type.t = 3, window.type.s = 3,
                            window.size.t = 5, window.size.s = 1,
                            make.pval = TRUE, n.sim = 10)

wc.image(my.wc3, which.image = "wc", color.key = "interval", n.levels = 250,
         siglvl.contour = 0.1, siglvl.arrow = 0.05,
         legend.params = list(lab = "wavelet coherence levels"),
         timelab = "")

# Phase difference

wc.phasediff.image(my.wc, which.contour = "wc", use.sAngle = TRUE,
                   n.levels = 250, siglvl = 0.1,
                   legend.params = list(lab = "phase difference levels",
                                        lab.line = 3),
                   timelab = "")

## Peaks detection

par(mfrow=c(1,1))

s.rate <- 14.29745
# lenght 109.6 s

409
412
414
415
416
422

z <- ts(calcium$ROI415, start = calcium$Time[1])

# peaks detection

fpeaks(spec(z, f=s.rate), nmax = 5,
       threshold = 10,
       plot = T,
       title = TRUE,
       xlab = "Frequency (Hz)", ylab = "Amplitude",
       labels = TRUE, legend = TRUE, digits = 2)*1000


###
z <- ts(calcium$ROI222, start = calcium$Time[1])

specgram(z, 250, 14.29745, 10,
         overlap = ceiling(length(window)/2))

specgram(z, 500, 14.29745, 11,
         overlap = ceiling(length(window)/2))

specgram(z, 500, 14.29745, 12, 10)

specgram(z, 500, 14.29745, 11)

specgram(z, 500, 14.29745, 11)


?specgram
col = heat.colors(256)

###

spectro3D(z, 14.29745, wl = 512, wn = "hanning", zp = 0,
          ovlp = 0, norm = TRUE, correction = "none", fftw = FALSE,
          plot = TRUE,
          magt = 2000, magf = 1000, maga = 10000,
          palette = reverse.terrain.colors)

?spectro3D
### spectrogram for the signal

par(mfrow=c(2,2), mar=c(3,1,1,2))

z <- ts(calcium$ROI5, start = calcium$Time[1])

spectrogram(z, 14.29745, windowlength = 95000, 0.01,
            timestep = 2000,
            maxfreq = 3,
            colors = TRUE,
            dynamicrange = 50, nlevels = dynamicrange, maintitle = "", 
            show = TRUE, window = 'kaiser', windowparameter = 4, 
            quality = TRUE)

spectrogram(z, 14.29745, windowlength = 95000, 0.01,
            timestep = 2000,
            maxfreq = 3,
            colors = TRUE,
            dynamicrange = 50, nlevels = dynamicrange, maintitle = "", 
            show = TRUE, window = 'gaussian', windowparameter = 0.3, 
            quality = TRUE)

spectrogram(z, 14.29745, windowlength = 95000, 0.01,
            timestep = 2000,
            maxfreq = 3,
            colors = TRUE,
            dynamicrange = 50, nlevels = dynamicrange, maintitle = "", 
            show = TRUE, window = 'rectangular', 
            quality = TRUE)

spectrogram(z, 14.29745, windowlength = 95000, 0.01,
            timestep = 1000,
            maxfreq = 3,
            colors = TRUE,
            dynamicrange = 40, nlevels = dynamicrange, maintitle = "", 
            show = TRUE, window = 'hann', 
            quality = TRUE)

spectrogram(z, 14.29745, windowlength = 95000, 0.01,
            timestep = 2000,
            maxfreq = 3,
            colors = TRUE,
            dynamicrange = 40, nlevels = dynamicrange, maintitle = "", 
            show = TRUE, window = 'hamming', 
            quality = TRUE)

spectrogram(z, 14.29745, windowlength = 95000, 0.01,
            timestep = 2000,
            maxfreq = 3,
            colors = TRUE,
            dynamicrange = 40, nlevels = dynamicrange, maintitle = "", 
            show = TRUE, window = 'cosine', 
            quality = TRUE)

spectrogram(z, 14.29745, windowlength = 95000, 0.01,
            timestep = 2000,
            maxfreq = 3,
            colors = TRUE,
            dynamicrange = 50, nlevels = dynamicrange, maintitle = "", 
            show = TRUE, window = 'bartlett', 
            quality = TRUE)

?spectrogram

####

spectro(z, 14.29745, wl = 512, wn = "hanning", zp = 16,
        ovlp = 50, complex = FALSE, norm = TRUE, correction="none",
        fftw = TRUE, plot = TRUE,
        flog = FALSE, grid = TRUE, osc = TRUE, scale = TRUE, cont = TRUE,
        collevels = NULL, palette = spectro.colors,
        contlevels = NULL, colcont = "black",
        colbg = "white", colgrid = "black",
        colaxis = "black", collab="black",
        cexlab = 1, cexaxis = 1, 
        tlab = "Time (s)",
        flab = "Frequency (kHz)",
        alab = "Amplitude",
        scalelab = "Amplitude",
        main = NULL,
        scalefontlab = 1, scalecexlab =0.75,
        axisX = TRUE, axisY = TRUE, tlim = NULL, trel = TRUE,
        flim = NULL, flimd = NULL,
        widths = c(6,1), heights = c(3,1),
        oma = rep(0,4),
        listen=FALSE)

install.packages("fftw")

?spectro
###

spectro(z, f=14.29745, osc=TRUE, flim=c(0,0.0001))

spectro(z, f=14.29745, osc=TRUE)

spectro(z,f=14.29745,ovlp=85,zp=16,osc=TRUE,
        cont=TRUE,contlevels=seq(-30,0,20),colcont="red",
        lwd=1.5,lty=2,palette=reverse.terrain.colors)

### не похоже на првильный результат:

which(Mod(fft(z)) == max(abs(Mod(fft(z)))))*s.rate/length(z)

##

ggspectro(z, 14.29745)

###

z <- ts(calcium$ROI416, start = calcium$Time[1])
plot(z)
plot(fft(z), ylim = c(-2000, 2000))


#####



z <- ts(calcium$ROI543, start = calcium$Time[1])
z1 <- ts(calcium$ROI543, start = calcium$Time[1])
z2 <- ts(calcium$ROI333, start = calcium$Time[1])
RoverS(z, n.block.min = 2, scale.ratio = 2, scale.min = 8)
?RoverS
##

i <- c(1:10)

for (i in 1:10) {
  print(i)
}


## no lines for long series
z <- ts(calcium$ROI410, start = calcium$Time[1])
?ts
z3 <- ts(cbind(calcium$ROI147, calcium$ROI543, calcium$ROI333),
         start = calcium$Time[1])


par(mfrow=c(1,1))

plot(z)

plot(z3, plot.type = 's', col = 1:3)

z_mult <- ts(cbind(calcium[,410:420]),
             start = calcium$Time[1])

plot(z_mult, plot.type = 's', col = 1:10)


####

z3 <- ts(cbind(calcium$ROI147, calcium$ROI416, calcium$ROI409),
         start = calcium$Time[1])

plot(z3, plot.type = 's', col = 1:3)

plot(z3, plot.type = 'm', col = 1:3)

plot(z3, plot.type = 'm', mar=c(gap=0.3, 5.1, gap=0.3, 2.1))

####
wf(z, 14.29745)
####

oscillo(z, 14.29745)


#####

install.packages("biwavelet")

plot(ts(z), xlab = NA, ylab = NA)
dim_z <- dim(z)

plot(z, type = "power")

plot(z, type = "power.corr.norm", plot.cb=TRUE, plot.phase=FALSE)

## returns the cepstrum of a time wave

ceps(z, 14.29745)

##

coh(z1, z2, 14.29745)

# Time series creation

z <- ts(calcium$ROI147, start = calcium$Time[1])

# Explorative analysis
# First, we log transform the data:

x <- log10(z)

# Plot of the time series and time-inverted time series:

par(mfrow=c(2,1), mar=c(0,0,0,0))
plot(x, ax=F)
box()
plot(x[length(x):1], type="l", ax=F)
box()

# Nonparametric regression function of 
#VXt versus Xt−1 and of Xt versus Xt−3 (kernel estimation):

par(mfrow=c(2,1), mar=c(3,3,1,1))
autopairs(x, lag=1, type="regression")
autopairs(x, lag=3, type="regression")

# The marginal histogram of data shows unimodality:

hist(x, br=13)

# Global and partial autocorrelation:

par(mfrow=c(2,1), mar=c(2,4,0,0))
acf(x)
pacf(x)

# The Average Mutual Information

mutual(x)

# Recurrence plot

recurr(x, m=3, d=1, levels=c(0,0.2,1)) # computation time!

# Directed lines are a tipycal tool for time series explorations. 

lag.plot(x, lags=3, layout=c(1,3))

# conditional mutual independence and linearity for lags 2 and 3:

delta.test(x)
delta.lin.test(x)

# Model selection

mod.ar <- linear(x, m=2)
mod.ar

# a SETAR  model with threshold delay δ = 1

mod.setar <- setar(x, m=2, mL=2, mH=2, thDelay=1)
mod.setar

# fit different linear and nonlinear models
# and directly compare some measures of their fit:

mod <- list()
mod[["linear"]] <- linear(x, m=2)
mod[["setar"]] <- setar(x, m=2, thDelay=1)
mod[["lstar"]] <- lstar(x, m=2, thDelay=1)

mod[["nnetTs"]] <- nnetTs(x, m=2, size=3)
mod[["aar"]] <- aar(x, m=2)

sapply(mod, AIC)
sapply(mod, MAPE)

summary(mod[["setar"]])

plot(mod[["setar"]])

# Out-of-sample forecasting

set.seed(10)
mod.test <- list()
x.train <- window(x, end=1924)
x.test <- window(x, start=1925)
mod.test[["linear"]] <- linear(x.train, m=2)
mod.test[["setar"]] <- setar(x.train, m=2, thDelay=1)
mod.test[["lstar"]] <- lstar(x.train, m=2, thDelay=1, trace=FALSE, control=list(maxit=1e5))
mod.test[["nnet"]] <- nnetTs(x.train, m=2, size=3, control=list(maxit=1e5))

mod.test[["aar"]] <- aar(x.train, m=2)

# Inspecting model skeleton

x.new <- predict(mod[["linear"]], n.ahead=100)
lag.plot(x.new, 1)

x.new <- predict(mod[["setar"]], n.ahead=100)
lag.plot(x.new, 1)

# A stable periodic cycle

x.new <- predict(mod[["nnetTs"]], n.ahead=100)
lag.plot(x.new, 1)

# Possibly chaotic systems:

mod.chaos1 <- setar(x, m=5, mL=5, mH=3, thDelay=1, th=0)
lag.plot(predict(mod.chaos1, n.ahead=100))

mod.chaos2 <- setar(x, m=5, mL=5, mH=3, thDelay=1, th=0)
lag.plot(predict(mod.chaos2, n.ahead=100))

# estimating the maximal Lyapunov exponent with the Kantz algorithm

N <- 1000
x.new <- predict(mod[["setar"]], n.ahead=N)
x.new <- x.new + rnorm(N, sd=sd(x.new)/100)
ly <- lyap_k(x.new, m=2, d=1, t=1, k=2, ref=750, s=200, eps=sd(x.new)/10)
plot(ly)

# There is no scaling region, so the maximal Lyapunov exponent can assumed to be ≤ 0.

x.new <- predict(mod.chaos2, n.ahead=N)
x.new <- x.new + rnorm(N, sd=sd(x.new)/100)
ly <- lyap_k(x.new, m=5, d=1, t=1, k=2, ref=750, s=200, eps=sd(x.new)/10)
plot(ly)

# Here there is a scaling region. The final λ estimate for this time series is the
# slope of the plotted curve in that region:

lyap(ly,start=0,end=5)

### FRACTAL

z <- ts(calcium$ROI415, start = calcium$Time[1])
b <- ts(calcium$ROI1, start = calcium$Time[1])


FNN(z, dimension = 5, tlag = NULL, rtol=10, atol = 2, olag = 1)

FNS(z, dimension = 5, tlag = NULL, image.tol=1, atol = 1, olag = 1)

KDE(z, at=NULL, n.grid = 100)

#Estimate the Hurst coefficient 

RoverS(z, n.block.min = 2, scale.ratio = 2, scale.min = 8)

#Detrended fluctuation analysis

?DFA

DFA(z, detrend = "poly1", sum.order = 0,
    overlap = 0,
    scale.max = trunc(length(z)/2), scale.min = NULL,
    scale.ratio = 2, verbose = FALSE)

plot(DFA(z, detrend = "poly1", sum.order = 1))

#Hurst coefficient by Whittle method

FDWhittle(z, method = "continuous", dc = F, freq.max = 0.5,
          delta.min = -1, delta.max = 2.5, sdf.method = "direct")

# Estimete Hurst parameter H

hurstBlock(z, method = "aggAbs", scale.min = 8, scale.max = NULL,
           scale.ratio = 2, weight = function(z) rep(1, length(z)),
           fit = lm)

#Information Dimension

infoDim(z, dimension = 5,tlag = NULL,
        olag = 0, n.density = 100, metric = Inf,
        max.neighbors = as.integer(min(c(round(length(z)/3),100))),
        n.reference = as.integer(round(length(z)/20)))

# PSR test for nonstationarity

stationarity(z, n.taper = 5,
             n.block = max(c(2, floor(logb(length(z), 
                                           base = 2)))),
             significance = 0.05,
             center = TRUE, recenter = FALSE)

# the existence of determenistic structure

determinism(z, dimension = 6, tlag = NULL,
            olag = 1, scale.min = NULL, scale.max = NULL,
            resolution = NULL, method = "ce",
            n.realization = 10, attach.summary = TRUE,
            seed = 0)
# create a map using the extrema of scalar time series

f <- poincareMap(z, extrema = "min", denoise = FALSE)
f <- embedSeries(f$amplitude, tlag = 1, dimension = 2)
plot(f, pch=1, cex=1)

# 
?spaceTime

spaceTime(z, dimension = 2, tlag = timeLag(z, method = "acfdecor"),
          olag.max = as.integer(min(500,length(z)/20)),
          probability = 0.1)

plot(spaceTime(z, dimension = 2, tlag = timeLag(z, method = "acfdecor"),
               olag.max = as.integer(min(500,length(z)/20)),
               probability = 0.1))

# Для регулярных движений L<=0
# В хаотических режимах L>0

lyapunov(z, tlag = NULL, dimension = 5, local.dimension = 3,
         reference = NULL, n.reference = NULL, olag = 2,
         sampling.interval = NULL, polynomial.order = 3,
         metric = Inf, scale = NULL)

plot(lyapunov(z, tlag = NULL, dimension = 5, local.dimension = 3,
              reference = NULL, n.reference = NULL, olag = 2,
              sampling.interval = NULL, polynomial.order = 3,
              metric = Inf, scale = NULL))

# Для регулярных движений L<=0
# В хаотических режимах L>0
?lyapunov
lyapunov(z, dimension = 5, local.dimension = 3, olag = 2,
         polynomial.order = 3,
         metric = Inf)


plot(lyapunov(z, dimension = 5, local.dimension = 3, olag = 2,
              polynomial.order = 3,
              metric = Inf))

## evaluate the maximal Lyapunov exponent 

z <- ts(calcium$ROI147, start = calcium$Time[1])

output <- lyap_k(z, m=5, d=3, s=1500, t=10, ref=1000, k=3, eps = 10)
plot(output, ylim = c(1, 4))

lyap(output, 0.73, 2.47)

?lyap_k

## ML Estimates for Fractionally-Differenced ARIMA (p,d,q) models

fracdiff(z, nar = 0, nma = 0,
         ar = rep(NA, max(1)), ma = rep(NA, max(1)),
         dtol = NULL, drange = c(0, 0.5), M = 100, trace = 0)

##

estimateEmbeddingDim(z, number.points = length(z),
                     time.lag = 1, max.embedding.dim = 10, threshold = 0.95,
                     max.relative.change = 0.1, do.plot = TRUE,
                     main = "Computing the embedding dimension", xlab = "dimension (d)",
                     ylab = "E1(d) & E2(d)", ylim = NULL, xlim = NULL)

## Space-time separation plot
z <- ts(calcium$ROI222, start = calcium$Time[1])
stplot(z, m=3, d=10, idt=1, 500)

## Plot time series against lagged versions of themselves

lag.plot(z, 4, diag.col = "forest green", type = "l", main = "Lag plot")

## Multivariate (but non-stationary! ...)

lag.plot(calcium[,409:417], lags = 1, type = "l")

## no lines for long series
z <- ts(calcium$ROI410, start = calcium$Time[1])

lag.plot(sqrt(z), set = c(1:4, 9:12), pch = ".", col = "black")


###

##

fracdiff(z, nar = 0, nma = 0,
         ar = rep(NA, max(1)), ma = rep(NA, max(1)),
         dtol = NULL, drange = c(0, 0.5), M = 100, trace = 0)

##

estimateEmbeddingDim(z, number.points = length(z),
                     time.lag = 1, max.embedding.dim = 10, threshold = 0.95,
                     max.relative.change = 0.1, do.plot = TRUE,
                     main = "Computing the embedding dimension", xlab = "dimension (d)",
                     ylab = "E1(d) & E2(d)", ylim = NULL, xlim = NULL)

##

xyz <- embedd(z, m=3, d=8)

scatterplot3d(xyz, type = "p")

?scatterplot3d

### Для регулярных движений L<=0
# В хаотических режимах L>0

z <- ts(calcium$ROI543, start = calcium$Time[1])
output <- lyap_k(z, m=5, d=3, s=1500, t=10, ref=1000, k=3, eps = 10)
plot(output, ylim = c(1, 4))

lyap(output, 0.73, 2.47)

?lyap_k

## Space-time separation plot
z <- ts(calcium$ROI410, start = calcium$Time[1])
stplot(z, m=3, d=10, idt=1, 500)

##Plot time series against lagged versions of themselves
?lag.plot

lag.plot(z, 1, layout = NULL, set.lags = 1:lags,
         main = NULL, asp = 1,
         diag = TRUE, diag.col = "gray", type = "p", oma = NULL,
         ask = NULL, do.lines = (n <= 150), labels = do.lines)

lag.plot(z, 4, diag.col = "forest green", type = "l", main = "Lag plot")

lag.plot(z, 5)

lag.plot(z, 6, layout = c(2,1), asp = NA, col.main = "blue")

lag.plot(calcium[,409:417], lags = 1, type = "l")

## no lines for long series
z <- ts(calcium$ROI415, start = calcium$Time[1])

lag.plot(sqrt(z), set = c(1:4, 9:12), pch = ".", col = "black")

### PCA

calcium2 <- cbind(calcium$ROI553,
                  calcium$ROI461,
                  calcium$ROI418,
                  calcium$ROI333,
                  calcium$ROI222,
                  calcium$ROI148)

pairs(calcium2)

arc.pca1 <- princomp(calcium2, scores = TRUE, cor = TRUE)

arc.pca2 <- prcomp(calcium2)

summary(arc.pca2)

par(mfrow=c(1,1))

plot(arc.pca2)

biplot(arc.pca2)

arc.pca1$loadings

arc.pca1$scores

df  <- mtcars

df$vs  <- factor(df$vs  , labels = c("V", "S"))
df$am  <- factor(df$am  , labels = c("Auto", "Manual"))

hist(df$mpg, breaks = 10, xlab = "MPG")
boxplot(mpg ~ am, df, ylab = "MPG", xlab = "AM")
plot(df$mpg, df$hp)

hist(df$mpg, breaks = 20, xlab = "MPG", main ="Histogram of MPG", 
     col = "green", cex.lab = 1.3, cex.axis = 1.3)

plot(density(df$mpg), xlab = "MPG", main ="Density of MPG", 
     col = "green", cex.lab = 1.3, cex.axis = 1.3)

boxplot(mpg ~ am, df, ylab = "MPG", main ="MPG and AM", 
        col = "green", cex.lab = 1.3, cex.axis = 1.3)

boxplot(df$mpg[df$am == "Auto"], df$mpg[df$am == "Manual"], ylab = "MPG", main ="MPG and AM", 
        col = "green", cex.lab = 1.3, cex.axis = 1.3)


plot(df$mpg, df$hp, xlab = "MPG", ylab ="HP" , main ="MPG and HP", pch = 22)

plot(~ mpg + hp, df) 


#Step 2, 3: Library ggplot2

library(ggplot2)

ggplot(df, aes(x = mpg))+
  geom_histogram(fill = "white", col = "black", binwidth = 0.1)

ggplot(df, aes(x = mpg))+
  geom_histogram(fill = "white", col = "black", binwidth = 2)+
  xlab("Miles/(US) gallon")+
  ylab("Count")+
  ggtitle("MPG histogram")

ggplot(df, aes(x = mpg, fill = am))+
  geom_dotplot()+
  xlab("Miles/(US) gallon")+
  ylab("Count")+
  scale_fill_discrete(name="Transmission type")+
  ggtitle("MPG dotplot")


ggplot(df, aes(x = mpg))+
  geom_density(fill = "red")

ggplot(df, aes(x = mpg, fill = am))+
  geom_density(alpha = 0.5)+
  xlab("Miles/(US) gallon")+
  ylab("Count")+
  scale_fill_discrete(name="Transmission type")+
  ggtitle("MPG density plot")


ggplot(df, aes(x = am, y = hp, fill = vs))+
  geom_boxplot()+
  xlab("Transmission type")+
  ylab("Gross horsepower")+
  scale_fill_discrete(name="Engine type")+
  ggtitle("Gross horsepower and engine type")


ggplot(df, aes(x = mpg, y = hp, size = qsec))+
  geom_point()+
  xlab("Miles/(US) gallon")+
  ylab("Gross horsepower")+
  scale_size_continuous(name="1/4 mile time")+
  ggtitle("Miles/(US) gallon and Gross horsepower")


my_plot  <- ggplot(df, aes(x = mpg, y = hp, col = vs, size = qsec))+
  geom_point()

my_plot2  <- ggplot(df, aes(x = am, y = hp, fill = vs))

my_plot2 + geom_boxplot()

##

cell <- calcium[, c(392:400, 416:424)]
head(cell)
cell 
class(cell)

cell$mean <- NULL
mean(cell[, c(2:3)])
a <- data.frame(c(1:2), 0, 5)
a <- data.frame(c(1:2), 0, 5, 6)
a
mean(a$X5)
class(a)
mean(a[1,])
mean(a[,1])
mean(a[,2])
mean(a[,3])
mean(a[,4])
mean(a[,3], a[,2], a[,1])
print(a[,1:3])
mean(a[,1:3])

?data.frame
str(cell
    )

x <- ts(cell, start = calcium$Time[1])
y <- ts(calcium$ROI543, start = calcium$Time[1])

# choose variable

xcal <- x

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

# Entropy tests of serial and cross dependence for time series
## WARNING: computationally intensive, increase B with caution!

Srho.test.ts.p(z, B = 2, 
               lag.max = 10,
               ci.type = c("perm"),
               quant = 0.95, plot = TRUE)


Srho.ts(z, lag.max = 10, method = "integral", plot = TRUE,
        maxpts = 0, tol = 0.001)


w <- as.integer(calcium$ROI222)

ent <- Srho(w, lag.max = 10, stationary = TRUE,
            version = "FORTRAN")
ent


MSMVSampEn(mat = z, 
           M = 5, tau = 10, r = 1, eps = 1567, scaleMat = T)

print(which.max(spx))

print(which.max(spy))

print(which.max(x.spec))

print(x.spec)

#### Another six , delete some...

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

# The wavelet transform of x is computed as follows:

par(mfrow=c(1,1))
my.data <- data.frame(x = x)
my.w <- analyze.wavelet(my.data, "x",
                        loess.span = 0,
                        dt = 1, 
                        dj = 1/250,
                        lowerPeriod = 2,
                        upperPeriod = 128,
                        make.pval = TRUE, n.sim = 10)

# Plot the wavelet power spectrum, a series with a constant period

wt.image(my.w, color.key = "quantile", n.levels = 250, 
         legend.params = list(lab = "wavelet power levels", mar = 4.7))

# A series with a variable period
par(mfrow=c(2,1))
wt.image(my.w, n.levels = 250,
         legend.params = list(lab = "wavelet power levels"))

# The reconstructed series

my.rec <- reconstruct(my.w)
x.rec <- my.rec$series$x.r   #x: name of original series

# A series with two periods

my.data <- data.frame(x = x)
my.w <- analyze.wavelet(my.data, "x",
                        loess.span = 0,
                        dt = 1, 
                        dj = 1/250,
                        lowerPeriod = 2,
                        upperPeriod = 128,
                        make.pval = TRUE, n.sim = 10)
wt.image(my.w, n.levels = 250,
         legend.params = list(lab = "wavelet power levels"))

# Reconstruction, using only period 8

reconstruct(my.w, sel.period = 8, plot.waves = TRUE, lwd = c(1,2), legend.coords = "bottomleft")

reconstruct(my.w, sel.period = 16, plot.waves = TRUE, lwd = c(1,2), legend.coords = "bottomleft")

# actually uses the period closest to 10

my.w$Period[(my.w$Period > 9) & (my.w$Period < 11)]

# Compare two series with average powers calling:

x <- ts(calcium$ROI415, start = calcium$Time[1])
y <- ts(calcium$ROI416, start = calcium$Time[1])

my.data <- data.frame(x = x, y = y)
my.wx <- analyze.wavelet(my.data, "x", loess.span = 0,
                         dt = 1, dj = 1/20,
                         lowerPeriod = 2, upperPeriod = 64,
                         make.pval = TRUE,n.sim = 10)
my.wy <- analyze.wavelet(my.data, "y", loess.span = 0,
                         dt = 1, dj = 1/20,
                         lowerPeriod = 2, upperPeriod = 64,
                         make.pval = TRUE,n.sim = 10)

maximum.level = 1.001*max(my.wx$Power.avg, my.wy$Power.avg)
par(mfrow=c(1,2))
wt.avg(my.wx, maximum.level = maximum.level)
wt.avg(my.wy, maximum.level = maximum.level)

# White noise method

par(mfrow=c(1,1))

my.data <- data.frame(x = x)
my.w <- analyze.wavelet(my.data, "x",
                        method = "white.noise",
                        loess.span = 0,
                        dt = 1, 
                        dj = 1/250,
                        lowerPeriod = 2,
                        upperPeriod = 256,
                        make.pval = TRUE, n.sim = 10)
wt.image(my.w, color.key = "interval", n.levels = 250,
         legend.params = list(lab = "wavelet power levels"))

# Fourier randomization method

my.data <- data.frame(x = x)
my.w <- analyze.wavelet(my.data, "x",
                        method = "Fourier.rand",
                        loess.span = 0,
                        dt = 1, 
                        dj = 1/250,
                        lowerPeriod = 2,
                        upperPeriod = 256,
                        make.pval = TRUE, n.sim = 10)
wt.image(my.w, color.key = "interval", n.levels = 250,
         legend.params = list(lab = "wavelet power levels"))

# Plotting the power spectrum

par(mfrow=c(1,1))

wt.image(my.w, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", label.digits = 2))


# Grayscale

wt.image(my.w, n.levels = 250,
         legend.params = list(lab = "wavelet power levels", label.digits = 2),
         color.palette = "gray((n.levels):1/n.levels)",
         col.ridge = "blue")

# Time axis

my.data <- data.frame(x = x)

my.w <- analyze.wavelet(my.data, "x",
                        method = "Fourier.rand",
                        loess.span = 0,
                        dt = 1, 
                        dj = 1/250,
                        lowerPeriod = 2,
                        upperPeriod = 256,
                        make.pval = TRUE, n.sim = 10)

wt.image(my.w, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", label.digits = 2),
         timelab = "time")

# Analysis of bivariate time series

x <- ts(calcium$ROI415, start = calcium$Time[1])
y <- ts(calcium$ROI416, start = calcium$Time[1])

par(mfrow=c(1,1))

my.data <- data.frame(x = x, y = y)

my.wc <- analyze.coherency(my.data, my.pair = c("x", "y"),
                           loess.span = 0,
                           dt = 1, 
                           dj = 1/100,
                           lowerPeriod = 4, upperPeriod = 32,
                           make.pval = TRUE, n.sim = 10)


wc.image(my.wc, n.levels = 250, color.key = "interval",
         siglvl.contour = 0.1, siglvl.arrow = 0.05, which.arrow.sig = "wt",
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = "")

######


wc.image(my.wc, n.levels = 250,
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = "", periodlab = "period")


# A plot of the time-averaged cross-wavelet power

wc.avg(my.wc, siglvl = 0.01, sigcol = "red", sigpch = 20,
       periodlab = "period")

# Coherence

wc.image(my.wc, which.image = "wc", color.key = "interval", n.levels = 250,
         siglvl.contour = 0.1, siglvl.arrow = 0.05,
         legend.params = list(lab = "wavelet coherence levels"),
         timelab = "")


my.wc2 <- analyze.coherency(my.data, my.pair = c("x", "y"),
                            loess.span = 0, dt = 1, dj = 1/100,
                            lowerPeriod = 2, upperPeriod = 32,
                            window.type.t = 1, window.type.s = 1,
                            window.size.t = 5, window.size.s = 1,
                            make.pval = TRUE, n.sim = 10)

wc.image(my.wc2, which.image = "wc", color.key = "interval", n.levels = 250,
         siglvl.contour = 0.1, siglvl.arrow = 0.05,
         legend.params = list(lab = "wavelet coherence levels"),
         timelab = "")

my.wc3 <- analyze.coherency(my.data, my.pair = c("x", "y"),
                            loess.span = 0, dt = 1, dj = 1/100,
                            lowerPeriod = 2, upperPeriod = 64,
                            window.type.t = 3, window.type.s = 3,
                            window.size.t = 5, window.size.s = 1,
                            make.pval = TRUE, n.sim = 10)

wc.image(my.wc3, which.image = "wc", color.key = "interval", n.levels = 250,
         siglvl.contour = 0.1, siglvl.arrow = 0.05,
         legend.params = list(lab = "wavelet coherence levels"),
         timelab = "")

# Phase difference

wc.phasediff.image(my.wc, which.contour = "wc", use.sAngle = TRUE,
                   n.levels = 250, siglvl = 0.1,
                   legend.params = list(lab = "phase difference levels",
                                        lab.line = 3),
                   timelab = "")

## Peaks detection

par(mfrow=c(1,1))

s.rate <- 14.29745
lenght 109.6 s

409
412
414
415
416
422

z <- ts(calcium$ROI415, start = calcium$Time[1])

# peaks detection

fpeaks(spec(z, f=s.rate), nmax = 5,
       threshold = 10,
       plot = T,
       title = TRUE,
       xlab = "Frequency (Hz)", ylab = "Amplitude",
       labels = TRUE, legend = TRUE, digits = 2)*1000


###
z <- ts(calcium$ROI333, start = calcium$Time[1])

specgram(z, 250, 14.29745, 10,
         overlap = ceiling(length(window)/2))

specgram(z, 500, 14.29745, 11,
         overlap = ceiling(length(window)/2))

specgram(z, 500, 14.29745, 12, 10)

specgram(z, 500, 14.29745, 11)

specgram(z, 500, 14.29745, 11)


?specgram
col = heat.colors(256)

###

spectro3D(z, 14.29745, wl = 512, wn = "hanning", zp = 0,
          ovlp = 0, norm = TRUE, correction = "none", fftw = FALSE,
          plot = TRUE,
          magt = 100, magf = 100, maga = 20,
          palette = reverse.terrain.colors)

?spectro3D
###

spectrogram(z, 14.29745, windowlength = 5000, 0.1,
            timestep = 1000,
            maxfreq = 10,
            colors = TRUE)

spectrogram(z, 14.29745)

?spectrogram
###

spectro(z, 14.29745, wl = 512, wn = "hanning", zp = 0,
        ovlp = 0, complex = FALSE, norm = TRUE, correction="none",
        fftw = FALSE, plot = TRUE,
        flog = FALSE, grid = TRUE, osc = FALSE, scale = TRUE, cont = FALSE,
        collevels = NULL, palette = spectro.colors,
        contlevels = NULL, colcont = "black",
        colbg = "white", colgrid = "black",
        colaxis = "black", collab="black",
        cexlab = 1, cexaxis = 1, 
        tlab = "Time (s)",
        flab = "Frequency (kHz)",
        alab = "Amplitude",
        scalelab = "Amplitude",
        main = NULL,
        scalefontlab = 1, scalecexlab =0.75,
        axisX = TRUE, axisY = TRUE, tlim = NULL, trel = TRUE,
        flim = NULL, flimd = NULL,
        widths = c(6,1), heights = c(3,1),
        oma = rep(0,4),
        listen=FALSE)

?spectro
###

spectro(z, f=14.29745, osc=TRUE, flim=c(0,0.0001))

spectro(z, f=14.29745, osc=TRUE)

spectro(z,f=22050,ovlp=85,zp=16,osc=TRUE,
        cont=TRUE,contlevels=seq(-30,0,20),colcont="red",
        lwd=1.5,lty=2,palette=reverse.terrain.colors)

### не похоже на првильный результат:

which(Mod(fft(z)) == max(abs(Mod(fft(z)))))*s.rate/length(z)

##

ggspectro(z, 14.29745)

###

z <- ts(calcium$ROI416, start = calcium$Time[1])
plot(z)
plot(fft(z), ylim = c(-2000, 2000))


#####



z <- ts(calcium$ROI543, start = calcium$Time[1])
z1 <- ts(calcium$ROI543, start = calcium$Time[1])
z2 <- ts(calcium$ROI333, start = calcium$Time[1])
RoverS(z, n.block.min = 2, scale.ratio = 2, scale.min = 8)

##

i <- c(1:10)

for (i in 1:10) {
  print(i)
}


## no lines for long series
z <- ts(calcium$ROI410, start = calcium$Time[1])
?ts
z3 <- ts(cbind(calcium$ROI147, calcium$ROI543, calcium$ROI333),
         start = calcium$Time[1])


par(mfrow=c(1,1))

plot(z)

plot(z3, plot.type = 's', col = 1:3)

z_mult <- ts(cbind(calcium[,410:420]),
             start = calcium$Time[1])

plot(z_mult, plot.type = 's', col = 1:10)


####

z3 <- ts(cbind(calcium$ROI147, calcium$ROI416, calcium$ROI409),
         start = calcium$Time[1])

plot(z3, plot.type = 's', col = 1:3)

plot(z3, plot.type = 'm', col = 1:3)

plot(z3, plot.type = 'm', mar=c(gap=0.3, 5.1, gap=0.3, 2.1))

####
wf(z, 14.29745)
####

oscillo(z, 14.29745)


#####

install.packages("biwavelet")

plot(ts(z), xlab = NA, ylab = NA)
dim_z <- dim(z)

plot(z, type = "power")

plot(z, type = "power.corr.norm", plot.cb=TRUE, plot.phase=FALSE)

## returns the cepstrum of a time wave

ceps(z, 14.29745)

##

coh(z1, z2, 14.29745)

# Time series creation

z <- ts(calcium$ROI147, start = calcium$Time[1])

# Explorative analysis
# First, we log transform the data:

x <- log10(z)

# Plot of the time series and time-inverted time series:

par(mfrow=c(2,1), mar=c(0,0,0,0))
plot(x, ax=F)
box()
plot(x[length(x):1], type="l", ax=F)
box()

# Nonparametric regression function of 
#VXt versus Xt−1 and of Xt versus Xt−3 (kernel estimation):

par(mfrow=c(2,1), mar=c(3,3,1,1))
autopairs(x, lag=1, type="regression")
autopairs(x, lag=3, type="regression")

# The marginal histogram of data shows unimodality:

 ags 2 and 3:

delta.test(x)
delta.lin.test(x)

# Model selection

mod.ar <- linear(x, m=2)
mod.ar

# a SETAR  model with threshold delay δ = 1

mod.setar <- setar(x, m=2, mL=2, mH=2, thDelay=1)
mod.setar

# fit different linear and nonlinear models
# and directly compare some measures of their fit:

mod <- list()
mod[["linear"]] <- linear(x, m=2)
mod[["setar"]] <- setar(x, m=2, thDelay=1)
mod[["lstar"]] <- lstar(x, m=2, thDelay=1)

mod[["nnetTs"]] <- nnetTs(x, m=2, size=3)
mod[["aar"]] <- aar(x, m=2)

sapply(mod, AIC)
sapply(mod, MAPE)

summary(mod[["setar"]])

plot(mod[["setar"]])

# Out-of-sample forecasting

set.seed(10)
mod.test <- list()
x.train <- window(x, end=1924)
x.test <- window(x, start=1925)
mod.test[["linear"]] <- linear(x.train, m=2)
mod.test[["setar"]] <- setar(x.train, m=2, thDelay=1)
mod.test[["lstar"]] <- lstar(x.train, m=2, thDelay=1, trace=FALSE, control=list(maxit=1e5))
mod.test[["nnet"]] <- nnetTs(x.train, m=2, size=3, control=list(maxit=1e5))

mod.test[["aar"]] <- aar(x.train, m=2)

# Inspecting model skeleton

x.new <- predict(mod[["linear"]], n.ahead=100)
lag.plot(x.new, 1)

x.new <- predict(mod[["setar"]], n.ahead=100)
lag.plot(x.new, 1)

# A stable periodic cycle

x.new <- predict(mod[["nnetTs"]], n.ahead=100)
lag.plot(x.new, 1)

# Possibly chaotic systems:

mod.chaos1 <- setar(x, m=5, mL=5, mH=3, thDelay=1)

# mod.chaos1 <- setar(x, m=5, mL=5, mH=3, thDelay=1, th=-1)
lag.plot(predict(mod.chaos1, n.ahead=100))

mod.chaos2 <- setar(x, m=5, mL=5, mH=3, thDelay=1)

#mod.chaos2 <- setar(x, m=5, mL=5, mH=3, thDelay=1, th=0)
lag.plot(predict(mod.chaos2, n.ahead=100))

# estimating the maximal Lyapunov exponent with the Kantz algorithm

N <- 1000
x.new <- predict(mod[["setar"]], n.ahead=N)
x.new <- x.new + rnorm(N, sd=sd(x.new)/100)
ly <- lyap_k(x.new, m=2, d=1, t=1, k=2, ref=750, s=200, eps=sd(x.new)/10)
plot(ly)

# There is no scaling region, so the maximal Lyapunov exponent can assumed to be ≤ 0.

x.new <- predict(mod.chaos2, n.ahead=N)
x.new <- x.new + rnorm(N, sd=sd(x.new)/100)
ly <- lyap_k(x.new, m=5, d=1, t=1, k=2, ref=750, s=200, eps=sd(x.new)/10)
plot(ly)

# Here there is a scaling region. The final λ estimate for this time series is the
# slope of the plotted curve in that region:

lyap(ly,start=0,end=5)

### FRACTAL

z <- ts(calcium$ROI5, start = calcium$Time[1])
b <- ts(calcium$ROI1, start = calcium$Time[1])


FNN(z, dimension = 5, tlag = NULL, rtol=10, atol = 2, olag = 1)

FNS(z, dimension = 5, tlag = NULL, image.tol=1, atol = 1, olag = 1)

KDE(z, at=NULL, n.grid = 100)

#Estimate the Hurst coefficient 

RoverS(z, n.block.min = 2, scale.ratio = 2, scale.min = 8)

#Detrended fluctuation analysis

?DFA

DFA(z, detrend = "poly1", sum.order = 0,
    overlap = 0,
    scale.max = trunc(length(z)/2), scale.min = NULL,
    scale.ratio = 2, verbose = FALSE)

plot(DFA(z, detrend = "poly1", sum.order = 1))

#Hurst coefficient by Whittle method

FDWhittle(z, method = "continuous", dc = F, freq.max = 0.5,
          delta.min = -1, delta.max = 2.5, sdf.method = "direct")

# Estimete Hurst parameter H

hurstBlock(z, method = "aggAbs", scale.min = 8, scale.max = NULL,
           scale.ratio = 2, weight = function(z) rep(1, length(z)),
           fit = lm)

#Information Dimension

infoDim(z, dimension = 5,tlag = NULL,
        olag = 0, n.density = 100, metric = Inf,
        max.neighbors = as.integer(min(c(round(length(z)/3),100))),
        n.reference = as.integer(round(length(z)/20)))

# PSR test for nonstationarity

stationarity(z, n.taper = 5,
             n.block = max(c(2, floor(logb(length(z), 
                                           base = 2)))),
             significance = 0.05,
             center = TRUE, recenter = FALSE)

# the existence of determenistic structure

determinism(z, dimension = 6, tlag = NULL,
            olag = 1, scale.min = NULL, scale.max = NULL,
            resolution = NULL, method = "ce",
            n.realization = 10, attach.summary = TRUE,
            seed = 0)
# create a map using the extrema of scalar time series

f <- poincareMap(z, extrema = "min", denoise = FALSE)
f <- embedSeries(f$amplitude, tlag = 1, dimension = 2)
plot(f, pch=1, cex=1)

# 
?spaceTime

spaceTime(z, dimension = 2, tlag = timeLag(z, method = "acfdecor"),
          olag.max = as.integer(min(500,length(z)/20)),
          probability = 0.1)

z <- ts(calcium$ROI543, start = calcium$Time[1])

plot(spaceTime(z, dimension = 2, tlag = timeLag(z, method = "acfdecor"),
               olag.max = as.integer(min(500,length(z)/20)),
               probability = 0.1))

# Для регулярных движений L<=0
# В хаотических режимах L>0

lyapunov(z, tlag = NULL, dimension = 5, local.dimension = 3,
         reference = NULL, n.reference = NULL, olag = 2,
         sampling.interval = NULL, polynomial.order = 3,
         metric = Inf, scale = NULL)

z <- ts(calcium$ROI543, start = calcium$Time[1])

plot(lyapunov(z, tlag = NULL, dimension = 10, local.dimension = 10,
              reference = NULL, n.reference = NULL, olag = 2,
              sampling.interval = NULL, polynomial.order = 3,
              metric = Inf, scale = NULL))

# Для регулярных движений L<=0
# В хаотических режимах L>0
?lyapunov
lyapunov(z, dimension = 5, local.dimension = 3, olag = 2,
         polynomial.order = 3,
         metric = Inf)

z <- ts(calcium$ROI543, start = calcium$Time[1])

plot(lyapunov(z, dimension = 5, local.dimension = 3, olag = 2,
              polynomial.order = 3,
              metric = Inf))

## evaluate the maximal Lyapunov exponent 

z <- ts(calcium$ROI147, start = calcium$Time[1])

output <- lyap_k(z, m=5, d=3, s=1500, t=10, ref=1000, k=3, eps = 10)
plot(output, ylim = c(1, 4))

lyap(output, 0.73, 2.47)

?lyap_k

## ML Estimates for Fractionally-Differenced ARIMA (p,d,q) models

fracdiff(z, nar = 0, nma = 0,
         ar = rep(NA, max(1)), ma = rep(NA, max(1)),
         dtol = NULL, drange = c(0, 0.5), M = 100, trace = 0)

##

estimateEmbeddingDim(z, number.points = length(z),
                     time.lag = 1, max.embedding.dim = 10, threshold = 0.95,
                     max.relative.change = 0.1, do.plot = TRUE,
                     main = "Computing the embedding dimension", xlab = "dimension (d)",
                     ylab = "E1(d) & E2(d)", ylim = NULL, xlim = NULL)

?estimateEmbeddingDim
## Space-time separation plot
z <- ts(calcium$ROI409, start = calcium$Time[1])
stplot(z, m=3, d=10, idt=1, 500)

z <- ts(calcium$ROI543, start = calcium$Time[1])
stplot(z, m=33, d=11, idt=3, 100)
?stplot
## Plot time series against lagged versions of themselves

lag.plot(z, 4, diag.col = "forest green", type = "l", main = "Lag plot")
?lag.plot
## Multivariate (but non-stationary! ...)

lag.plot(calcium[,409:417], lags = 1, type = "l")
?lag.plot
## no lines for long series
z <- ts(calcium$ROI410, start = calcium$Time[1])

lag.plot(sqrt(z), set = c(1:4, 9:12), pch = ".", col = "black")


###

##

fracdiff(z, nar = 0, nma = 0,
         ar = rep(NA, max(1)), ma = rep(NA, max(1)),
         dtol = NULL, drange = c(0, 0.5), M = 100, trace = 0)

##

estimateEmbeddingDim(z, number.points = length(z),
                     time.lag = 1, max.embedding.dim = 10, threshold = 0.95,
                     max.relative.change = 0.1, do.plot = TRUE,
                     main = "Computing the embedding dimension", xlab = "dimension (d)",
                     ylab = "E1(d) & E2(d)", ylim = NULL, xlim = NULL)

##

xyz <- embedd(z, m=5, d=88)

scatterplot3d(xyz, type = "p")

?scatterplot3d

### Для регулярных движений L<=0
# В хаотических режимах L>0

z <- ts(calcium$ROI420, start = calcium$Time[1])
output <- lyap_k(z, m=5, d=3, s=1500, t=10, ref=1000, k=3, eps = 10)
plot(output, ylim = c(1, 4))

lyap(output, 0.73, 2.47)

?lyap_k

## Space-time separation plot
z <- ts(calcium$ROI410, start = calcium$Time[1])
stplot(z, m=3, d=10, idt=1, 500)

##Plot time series against lagged versions of themselves
?lag.plot

lag.plot(z, 1, layout = NULL, set.lags = 1:lags,
         main = NULL, asp = 1,
         diag = TRUE, diag.col = "gray", type = "p", oma = NULL,
         ask = NULL, do.lines = (n <= 150), labels = do.lines)

lag.plot(z, 4, diag.col = "forest green", type = "l", main = "Lag plot")

lag.plot(z, 5)

lag.plot(z, 6, layout = c(2,1), asp = NA, col.main = "blue")

lag.plot(calcium[,409:417], lags = 1, type = "l")

## no lines for long series
z <- ts(calcium$ROI415, start = calcium$Time[1])

lag.plot(sqrt(z), set = c(1:4, 9:12), pch = ".", col = "black")

  ### PCA

calcium2 <- cbind(calcium$ROI553,
                  calcium$ROI461,
                  calcium$ROI418,
                  calcium$ROI333,
                  calcium$ROI222,
                  calcium$ROI148)

pairs(calcium2)

arc.pca1 <- princomp(calcium2, scores = TRUE, cor = TRUE)

arc.pca2 <- prcomp(calcium2)

summary(arc.pca2)

par(mfrow=c(1,1))

plot(arc.pca2)

biplot(arc.pca2)

arc.pca1$loadings

arc.pca1$scores

df  <- mtcars

df$vs  <- factor(df$vs  , labels = c("V", "S"))
df$am  <- factor(df$am  , labels = c("Auto", "Manual"))

hist(df$mpg, breaks = 10, xlab = "MPG")
boxplot(mpg ~ am, df, ylab = "MPG", xlab = "AM")
plot(df$mpg, df$hp)

hist(df$mpg, breaks = 20, xlab = "MPG", main ="Histogram of MPG", 
     col = "green", cex.lab = 1.3, cex.axis = 1.3)

plot(density(df$mpg), xlab = "MPG", main ="Density of MPG", 
     col = "green", cex.lab = 1.3, cex.axis = 1.3)

boxplot(mpg ~ am, df, ylab = "MPG", main ="MPG and AM", 
        col = "green", cex.lab = 1.3, cex.axis = 1.3)

boxplot(df$mpg[df$am == "Auto"], df$mpg[df$am == "Manual"], ylab = "MPG", main ="MPG and AM", 
        col = "green", cex.lab = 1.3, cex.axis = 1.3)


plot(df$mpg, df$hp, xlab = "MPG", ylab ="HP" , main ="MPG and HP", pch = 22)

plot(~ mpg + hp, df) 


#Step 2, 3: Library ggplot2

library(ggplot2)

ggplot(df, aes(x = mpg))+
  geom_histogram(fill = "white", col = "black", binwidth = 0.1)

ggplot(df, aes(x = mpg))+
  geom_histogram(fill = "white", col = "black", binwidth = 2)+
  xlab("Miles/(US) gallon")+
  ylab("Count")+
  ggtitle("MPG histogram")

ggplot(df, aes(x = mpg, fill = am))+
  geom_dotplot()+
  xlab("Miles/(US) gallon")+
  ylab("Count")+
  scale_fill_discrete(name="Transmission type")+
  ggtitle("MPG dotplot")


ggplot(df, aes(x = mpg))+
  geom_density(fill = "red")

ggplot(df, aes(x = mpg, fill = am))+
  geom_density(alpha = 0.5)+
  xlab("Miles/(US) gallon")+
  ylab("Count")+
  scale_fill_discrete(name="Transmission type")+
  ggtitle("MPG density plot")


ggplot(df, aes(x = am, y = hp, fill = vs))+
  geom_boxplot()+
  xlab("Transmission type")+
  ylab("Gross horsepower")+
  scale_fill_discrete(name="Engine type")+
  ggtitle("Gross horsepower and engine type")


ggplot(df, aes(x = mpg, y = hp, size = qsec))+
  geom_point()+
  xlab("Miles/(US) gallon")+
  ylab("Gross horsepower")+
  scale_size_continuous(name="1/4 mile time")+
  ggtitle("Miles/(US) gallon and Gross horsepower")


my_plot  <- ggplot(df, aes(x = mpg, y = hp, col = vs, size = qsec))+
  geom_point()

my_plot2  <- ggplot(df, aes(x = am, y = hp, fill = vs))

my_plot2 + geom_boxplot()

##

cell <- calcium[, c(392:400, 416:424)]
head(cell)
cell

boxplot(a)
a[,1]
mean(a[1, ])
a[1, ]
b <- c(a[1,])
b
mean(b)
boxplot(b)



mean(data.frame(c(1:100)))


q <- data.frame(c(1:100))
q
mean(q$c.1.100.)
aggregate(q)

?aggregate

##  tabula rasa

a <- read.csv2('a.csv')

# Explorative analysis
# First, we log transform the data:

z <- ts(calcium$ROI9, start = calcium$Time[1])
x <- log10(z)

# Plot of the time series and time-inverted time series:

par(mfrow=c(2,1), mar=c(0,0,0,0))
plot(x, ax=F)
box()
plot(x[length(x):1], type="l", ax=F)
box()

# Nonparametric regression function of 
#VXt versus Xt−1 and of Xt versus Xt−3 (kernel estimation):

par(mfrow=c(2,1), mar=c(3,3,1,1))
autopairs(x, lag=1, type="regression")
autopairs(x, lag=3, type="regression")

# The marginal histogram of data shows unimodality:

hist(x, br=13)

# Global and partial autocorrelation:

par(mfrow=c(2,1), mar=c(2,4,0,0))
acf(x)
pacf(x)

# The Average Mutual Information

mutual(x)

# Recurrence plot

recurr(x, m=3, d=1, levels=c(0,0.2,1)) # computation time!

# Directed lines are a tipycal tool for time series explorations. 

lag.plot(x, lags=3, layout=c(1,3))

# conditional mutual independence and linearity for lags 2 and 3:

delta.test(x)
delta.lin.test(x)

# 
?spaceTime

spaceTime(z, dimension = 2, tlag = timeLag(z, method = "acfdecor"),
          olag.max = as.integer(min(500,length(z)/20)),
          probability = 0.1)

plot(spaceTime(z, dimension = 2, tlag = timeLag(z, method = "acfdecor"),
               olag.max = as.integer(min(500,length(z)/20)),
               probability = 0.1))

# Для регулярных движений L<=0
# В хаотических режимах L>0

lyapunov(z, tlag = NULL, dimension = 5, local.dimension = 3,
         reference = NULL, n.reference = NULL, olag = 2,
         sampling.interval = NULL, polynomial.order = 3,
         metric = Inf, scale = NULL)

plot(lyapunov(z, tlag = NULL, dimension = 5, local.dimension = 3,
              reference = NULL, n.reference = NULL, olag = 2,
              sampling.interval = NULL, polynomial.order = 3,
              metric = Inf, scale = NULL))

# Для регулярных движений L<=0
# В хаотических режимах L>0
?lyapunov
lyapunov(z, dimension = 5, local.dimension = 3, olag = 2,
         polynomial.order = 3,
         metric = Inf)


plot(lyapunov(z, dimension = 5, local.dimension = 3, olag = 2,
              polynomial.order = 3,
              metric = Inf))

## evaluate the maximal Lyapunov exponent 

z <- ts(calcium$ROI147, start = calcium$Time[1])

output <- lyap_k(z, m=5, d=3, s=1500, t=10, ref=1000, k=3, eps = 10)
plot(output, ylim = c(1, 4))

lyap(output, 0.73, 2.47)

?lyap_k

## ML Estimates for Fractionally-Differenced ARIMA (p,d,q) models

fracdiff(z, nar = 0, nma = 0,
         ar = rep(NA, max(1)), ma = rep(NA, max(1)),
         dtol = NULL, drange = c(0, 0.5), M = 100, trace = 0)

##

estimateEmbeddingDim(z, number.points = length(z),
                     time.lag = 1, max.embedding.dim = 10, threshold = 0.95,
                     max.relative.change = 0.1, do.plot = TRUE,
                     main = "Computing the embedding dimension", xlab = "dimension (d)",
                     ylab = "E1(d) & E2(d)", ylim = NULL, xlim = NULL)

## Space-time separation plot
z <- ts(calcium$ROI2, start = calcium$Time[1])
stplot(z, m=3, d=10, idt=1, 500)

## Plot time series against lagged versions of themselves

lag.plot(z, 4, diag.col = "forest green", type = "l", main = "Lag plot")

## Multivariate (but non-stationary! ...)

lag.plot(calcium[,2:41], lags = 1, type = "l")

## no lines for long series
z <- ts(calcium[,2:41], start = calcium$Time[1])

lag.plot(sqrt(z), set = c(1:4, 9:12), pch = ".", col = "black")


###

##
z <- ts(calcium$ROI2, start = calcium$Time[1])

fracdiff(z, nar = 0, nma = 0,
         ar = rep(NA, max(1)), ma = rep(NA, max(1)),
         dtol = NULL, drange = c(0, 0.5), M = 100, trace = 0)

##

estimateEmbeddingDim(z, number.points = length(z),
                     time.lag = 1, max.embedding.dim = 10, threshold = 0.95,
                     max.relative.change = 0.1, do.plot = TRUE,
                     main = "Computing the embedding dimension", xlab = "dimension (d)",
                     ylab = "E1(d) & E2(d)", ylim = NULL, xlim = NULL)

##

xyz <- embedd(z, m=3, d=8)

scatterplot3d(xyz, type = "p")

?scatterplot3d

### Для регулярных движений L<=0
# В хаотических режимах L>0

z <- ts(calcium$ROI543, start = calcium$Time[1])
output <- lyap_k(z, m=5, d=3, s=1500, t=10, ref=1000, k=3, eps = 10)
plot(output, ylim = c(1, 4))

lyap(output, 0.73, 2.47)

?lyap_k

## Space-time separation plot
z <- ts(calcium$ROI410, start = calcium$Time[1])
stplot(z, m=3, d=10, idt=1, 500)

##Plot time series against lagged versions of themselves
?lag.plot

lag.plot(z, 1, layout = NULL, set.lags = 1:lags,
         main = NULL, asp = 1,
         diag = TRUE, diag.col = "gray", type = "p", oma = NULL,
         ask = NULL, do.lines = (n <= 150), labels = do.lines)

lag.plot(z, 4, diag.col = "forest green", type = "l", main = "Lag plot")

lag.plot(z, 5)

lag.plot(z, 6, layout = c(2,1), asp = NA, col.main = "blue")

lag.plot(calcium[,409:417], lags = 1, type = "l")

## no lines for long series
z <- ts(calcium$ROI415, start = calcium$Time[1])

lag.plot(sqrt(z), set = c(1:4, 9:12), pch = ".", col = "black")

### PCA

calcium2 <- cbind(calcium$ROI553,
                  calcium$ROI461,
                  calcium$ROI418,
                  calcium$ROI333,
                  calcium$ROI222,
                  calcium$ROI148)

pairs(calcium2)
### PCA
arc.pca1 <- princomp(calcium, scores = TRUE, cor = TRUE)
biplot(arc.pca2)
### PCA
arc.pca1 <- princomp(z, scores = TRUE, cor = TRUE)
biplot(arc.pca2)

?princomp
arc.pca2 <- prcomp(calcium2)

summary(arc.pca1)

par(mfrow=c(1,1))

plot(arc.pca1)

biplot(arc.pca2)

arc.pca1$loadings

arc.pca1$scores


#####



z <- ts(calcium$ROI543, start = calcium$Time[1])
z1 <- ts(calcium$ROI543, start = calcium$Time[1])
z2 <- ts(calcium$ROI333, start = calcium$Time[1])
RoverS(z, n.block.min = 2, scale.ratio = 2, scale.min = 8)

?RoverS

##

### Multivariate time series 

m <- ts(calcium[,2:41], start = calcium$Time[1])

mm <- ts(calcium[,2:11], start = calcium$Time[1])

z3 <- ts(cbind(calcium$ROI7, calcium$ROI16, calcium$ROI40),
         start = calcium$Time[1])

head(m)

par(mfrow=c(1,1), mar=c(2,5,1,1))

plot(m, plot.type = "s")

# plot(mm, plot.type = "m")



x <- ts(calcium$ROI27, start = calcium$Time[1])

RoverS(x, n.block.min = 2, scale.ratio = 2, scale.min = 8)

#Detrended fluctuation analysis

?DFA


z <- ts(calcium$ROI21, start = calcium$Time[1])

DFA(z, detrend = "poly1", sum.order = 0,
    overlap = 0,
    scale.max = trunc(length(z)/2), scale.min = NULL,
    scale.ratio = 2, verbose = FALSE)

plot(DFA(z, detrend = "poly1", sum.order = 1))


#Hurst coefficient by Whittle method


z <- ts(calcium$ROI7, start = calcium$Time[1])

FDWhittle(z, method = "continuous", dc = F, freq.max = 0.5,
          delta.min = -1, delta.max = 2.5, sdf.method = "direct")



## no lines for long series
z <- ts(calcium[,2:41], start = calcium$Time[1])

lag.plot(sqrt(z), set = c(1:4, 9:12), pch = ".", col = "black")



