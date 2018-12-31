library(fractal)
library(fracdiff)
library(arfima)
library(nonlinearTseries)
library(tseriesChaos)

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











