calcium <- read.csv2(file.choose()) # select file

# Phase difference

par(mfrow=c(2,1))

x <- ts(calcium$ROI1, start = calcium$Time[1])
y <- ts(calcium$ROI2, start = calcium$Time[1])

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

?wc.phasediff.image
