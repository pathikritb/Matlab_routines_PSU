% Analyze frequency content of data, use stacked template
Fs = 25e6;
for i = 1
    i
    data = fulstack(847:909);
    data = padarray(data,nextpow2(length(data)));
    data = data/(max(data)-min(data));
    window = hann(length(data));
    windat = window.*data;
    [psdxwin,freq] = periodogram(data,window,0:0.1:Fs/2,Fs);
    psdx     = periodogram(data);
    %psdx        = psdx(1:end-1);
    %psdxwin     = psdxwin(1:end-1);
    if i==1
        h = subplot(1,2,1);
        semilogx(h,freq,psdx,'r-')
        hwin = subplot(1,2,2);
        semilogx(hwin,freq,psdxwin,'b--')
        set(h,'NextPlot','Add')
        set(hwin,'NextPlot','Add')
    else
        semilogx(h,freq,psdx,'r-')
        semilogx(hwin,freq,psdxwin,'b--')
    end
end
        
    