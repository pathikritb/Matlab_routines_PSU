% Power spectral density of signals
function [psdx,freq] = powerspecdens(x,Fs)
N = length(x);
xdft = fft(x);
xdft = xdft(1:round(N/2)+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(x):Fs/2;