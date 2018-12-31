# Calcium oscillations analysis, Dan Bobkov, 2018

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

# load data obtained from ImageJ

calcium <- read.csv2('20x20.csv') # 24x24 automatic ROI grid
calcium <- read.csv2('fluo4d.csv') # 40 cells are manually selected
# calcium <- read.csv2('ser017dot.csv') # ��� �� ��������
calcium <- read.csv2('ser017comma.csv')
potential <- read.csv2('di-8-anepps.csv')

head(calcium)

boxplot(calcium[,2:41])

# time series analysis

xcal <- calcium$ROI94

calcium_ts <- ts(xcal, start = calcium$Time[1])

del <- 0.1
x.spec <- spectrum(calcium_ts, log = "no", span = 10, plot = FALSE)
spx <- x.spec$freq/del
spy <- 2*x.spec$spec

par(mfrow=c(2,3))  # plot results

plot(xcal ~ calcium$Time, 
     xlab = "time (s)", ylab = "fluorescence intensity (a.u.)", t="l", 
     xlim = c(0, 50), ylim = c(0.7, 1.4), main = "Total signal")

plot(xcal ~ calcium$Time, 
     xlab = "time (s)", ylab = "fluorescence intensity (a.u.)", t="l", 
     xlim = c(0, 5), ylim = c(0.7, 1.4), main = "Signal segment")

boxplot(xcal, ylab = "fluorescence intensity (a.u.)",
        main = "Mean intensity", ylim = c(0.7, 1.4))

plot(spy ~ spx, xlab="frequency (Hz)", ylab="spectral density",
     type="l", lwd = 2)

acf(xcal, lag.max = 100, main = "Autocorrelation")

w <- as.integer(xcal)  # Entropy measure

Srho(w, lag.max = 10, stationary = TRUE, plot = TRUE,
     version = "FORTRAN")
?Srho

## Peaks detection, lenght = 109.6 s

par(mfrow=c(5,5), mar=c(2,3,1,1))

s.rate <- 14.29745

z <- ts(calcium$ROI94, start = calcium$Time[1])

fpeaks(spec(z, f=s.rate), nmax = 5,
       threshold = 10,
       plot = T,
       title = TRUE,
       xlab = "Frequency (Hz)", ylab = "Amplitude",
       labels = TRUE, legend = TRUE, digits = 2)*1000

## Space-time separation plot

par(mfrow=c(4,1), mar=c(2,3,1,1))

z <- ts(calcium$ROI94, start = calcium$Time[1])

stplot(z, m=3, d=10, idt=1, 500)

## Plot time series against lagged versions of themselves

par(mfrow=c(2,2), mar=c(2,3,1,1))

z <- ts(calcium$ROI94, start = calcium$Time[1])

lag.plot(z, 4, diag.col = "forest green", type = "l", main = "Lag plot")

## Multivariate (but non-stationary! ...)

lag.plot(calcium[,16:25], lags = 1, type = "l")

## evaluate the maximal Lyapunov exponent 
### Для регулярных движений L<=0
# В хаотических режимах L>0

z <- ts(calcium$ROI500, start = calcium$Time[1])

output <- lyap_k(z, m=2, d=10, s=10, t=100, ref=700, k=3, eps = 10)

plot(output, ylim = c(0, 0.1))

lyap(output, 0.73, 2.47)

plot(lyap(output, -1, 1))

?lyap_k

## 

z <- ts(calcium$ROI94, start = calcium$Time[1])

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
# spectrogram and periodogram

par(mfrow=c(4,1), mar=c(2,3,1,1))

xcal <- calcium$ROI24

x <- ts(xcal, start = calcium$Time[1])

spectrogram(x, 14.29745, windowlength = 95000, 0.01,
            timestep = 2000,
            maxfreq = 0.3,
            colors = TRUE,
            dynamicrange = 55, nlevels = dynamicrange, maintitle = "", 
            show = TRUE, window = 'gaussian', windowparameter = 0.3, 
            quality = TRUE)

my.data <- data.frame(x = x)

my.w <- analyze.wavelet(my.data, "x",
                        loess.span = 0,
                        dt = 1, 
                        dj = 1/500,
                        lowerPeriod = 2,
                        upperPeriod = 64,
                        make.pval = TRUE, n.sim = 10)

wt.image(my.w, n.levels = 250)

# The reconstructed series

my.rec <- reconstruct(my.w)
x.rec <- my.rec$series$x.r   #x: name of original series

# Reconstruction, using only period 16

reconstruct(my.w, sel.period = 16, plot.waves = TRUE, 
            lwd = c(1,2), legend.coords = "bottomleft")

# Reconstruction, using only period 8

reconstruct(my.w, sel.period = 5, plot.waves = TRUE, 
            lwd = c(1,2), legend.coords = "bottomleft")

# Compare two series with average powers calling:

par(mfrow=c(1,1), mar=c(2,3,1,1))

x <- ts(calcium$ROI9, start = calcium$Time[1])
y <- ts(calcium$ROI20, start = calcium$Time[1])

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

# Coherence analysis of bivariate time series

par(mfrow=c(1,1))

x <- ts(calcium$ROI8, start = calcium$Time[1])
y <- ts(calcium$ROI9, start = calcium$Time[1])

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

wc.image(my.wc, n.levels = 250,
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = "", periodlab = "period")


#

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

#

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

par(mfrow=c(2,1))

x <- ts(calcium$ROI8, start = calcium$Time[1])
y <- ts(calcium$ROI24, start = calcium$Time[1])

my.data <- data.frame(x = x, y = y)

my.wc <- analyze.coherency(my.data, my.pair = c("x", "y"),
                           loess.span = 0,
                           dt = 1, 
                           dj = 1/100,
                           lowerPeriod = 4, upperPeriod = 32,
                           make.pval = TRUE, n.sim = 10)

wc.phasediff.image(my.wc, which.contour = "wc", use.sAngle = TRUE,
                   n.levels = 250, siglvl = 0.1,
                   legend.params = list(lab = "phase difference levels",
                                        lab.line = 3),
                   timelab = "")

#Estimate the Hurst coefficient 

RoverS(z, n.block.min = 2, scale.ratio = 2, scale.min = 8)



# Explorative analysis
# First, we log transform the data:

z <- ts(calcium$ROI9, start = calcium$Time[1])
x <- log10(z)

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

z <- ts(calcium$ROI24, start = calcium$Time[1])


FNN(z, dimension = 5, tlag = NULL, rtol=10, atol = 2, olag = 1)

plot(FNN(z, dimension = 5, tlag = NULL, rtol=10, atol = 2, olag = 1))

FNS(z, dimension = 5, tlag = NULL, image.tol=1, atol = 1, olag = 1)

plot(FNS(z, dimension = 5, tlag = NULL, image.tol=1, atol = 1, olag = 1))

KDE(z, at=NULL, n.grid = 100)



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

determinism(z, dimension = 5, tlag = NULL,
            olag = 1, scale.min = NULL, scale.max = NULL,
            resolution = NULL, method = "phase",
            n.realization = 10, attach.summary = TRUE,
            seed = 0)

?determinism

# create a map using the extrema of scalar time series

f <- poincareMap(z, extrema = "min", denoise = FALSE)
f <- embedSeries(f$amplitude, tlag = 1, dimension = 2)
plot(f, pch=1, cex=1)

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

par(mfrow=c(3,3), mar=c(2,2,1,1))

z <- ts(calcium$ROI20, start = calcium$Time[1])

output <- lyap_k(z, m=5, d=3, s=1500, t=10, ref=1000, k=3, eps = 10)
plot(output, ylim = c(0, 4))

lyap(output, 0.73, 2.47)

?lyap_k

## ML Estimates for Fractionally-Differenced ARIMA (p,d,q) models

z <- ts(calcium$ROI20, start = calcium$Time[1])

fracdiff(z, nar = 0, nma = 0,
         ar = rep(NA, max(0, 1)), ma = rep(NA, max(1)),
         dtol = NULL, drange = c(0, 1), 
         h = min(0.1, 1.05e-8), 
         M = 100, trace = 0)
?fracdiff

##

z <- ts(calcium$ROI24, start = calcium$Time[1])

estimateEmbeddingDim(z, number.points = length(z),
                     time.lag = 1, max.embedding.dim = 10,
                     threshold = 0.95,
                     max.relative.change = 0.1, 
                     do.plot = TRUE,
                     main = "Computing the embedding dimension", 
                     xlab = "dimension (d)",
                     ylab = "E1(d) & E2(d)", 
                     ylim = NULL, xlim = NULL)


## Space-time separation plot

z <- ts(calcium$ROI2, start = calcium$Time[1])

stplot(z, m=3, d=10, idt=1, 500)

## Plot time series against lagged versions of themselves

lag.plot(z, 4, diag.col = "forest green", type = "l", main = "Lag plot")

## Multivariate (but non-stationary! ...)

lag.plot(calcium[,2:41], lags = 1, type = "l")

?lag.plot



i <- c(1:40)

for (i in 1:40) {
  lyap(lyap_k(ts(calcium[,i+1],
                 start = calcium$Time[1]),
              m=5, d=3, s=1500, t=10, ref=1000, 
              k=3, eps = 10), 0.73, 2.47)
}

#

lyap(lyap_k(ts(calcium[,i+1],
               start = calcium$Time[1]),
            m=5, d=3, s=1500, t=10, ref=1000, 
            k=3, eps = 10), 0.73, 2.47)

##

par(mfrow=c(3,3), mar=c(2,2,1,1))

z <- ts(calcium$ROI24, start = calcium$Time[1])

output <- lyap_k(ts(calcium$ROI24,
                    start = calcium$Time[1]),
                 m=5, d=3, s=1500, t=10, ref=1000, k=3, eps = 10)

lyap(lyap_k(ts(calcium$ROI24,
               start = calcium$Time[1]),
            m=5, d=3, s=1500, t=10, ref=1000, 
            k=3, eps = 10), 0.73, 2.47)


####

for (i in 1:40) {
  FDWhittle(ts(calcium[,i+1], start = calcium$Time[1]),
            method = "continuous", dc = F, freq.max = 0.5,
            delta.min = -1, delta.max = 2.5, sdf.method = "direct")
}


###

for (i in 1:40) {
  z <- ts(calcium[,i+1], start = calcium$Time[1])
  
  FDWhittle(z, method = "continuous", dc = F, freq.max = 0.5,
            delta.min = -1, delta.max = 2.5, sdf.method = "direct")
}




