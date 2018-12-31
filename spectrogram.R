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