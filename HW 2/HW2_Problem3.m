%% HW 2 problem 3
close all; clear variables;
data = load( 'hw2_3.mat' );
N = 2^12;
win = hamming(N);%Hamming window
wss = sqrt(sum(win.^2));
win = repmat(win/wss, 1,2^16/N);
d = reshape(data.d(:,2), N,[]);
% d = d-mean(d);
times = reshape((data.d(:,1)),N,[]);
windowedData = win.*d;%Columns of this matrix are the windowed data
%% Welch Periodogram
dhat = fft( windowedData );
dt = mean(diff( data.d(:,1) ));df = 1/(N*dt);
sHat = abs(dhat).^2;
s = mean(sHat.');
freqs = df*(0:N-1);freqs = freqs-mean(freqs);
figure;semilogx( freqs,10*log10(fftshift(s/sqrt(N))) );grid on;
xlabel('Frequency, Hz')
ylabel('Welch Periodogram, dB')
title('Welch Periodogram estimate')
saveas( gcf,'welch estimate.png' )
savefig(gcf, 'welch estimate.fig')
% ylim([0,50])
% BLackman Tukey Estimate
ac = [];
sBlackmanTukey = [];
for ii = 1:size(d,2)
    ac(ii,:) = conv(windowedData(:,ii),flip(windowedData(:,ii))) ;%Autocorrelation
    sBlackmanTukey(ii,:) = fft(ac(ii,:));
end
dfBT = df/2;
freqsBT = dfBT*(0:2*N-2);freqsBT = freqsBT-mean(freqsBT);
figure;semilogx(freqsBT, 10*log10(fftshift(abs(mean(sBlackmanTukey/sqrt(N))))));grid on
xlabel('Frequency, Hz')
ylabel('Blackman-Tukey PSD, dB')
title('Blackman-Tukey estimate')
saveas( gcf,'blackman tukey.png' )
savefig(gcf, 'blackman tukey.fig')
% ylim([0,50])
lags = (0:length(mean(ac))-1);
figure;plot(lags-mean(lags),mean(ac));
xlabel('Lag');ylabel('Autocorrelation')
title('Estimate of the signal autocorrelation function')
grid on
saveas( gcf,'autocorrelation function estimate.png' )
savefig(gcf, 'autocorrelation function esstimate.fig')
%% Autoregressive estimate
%Run on the entire data stream at once. 
%all-pole model
acYuleWalker = conv( data.d(:,2) ,flip(data.d(:,2) ) );
acYuleWalker = acYuleWalker/max(acYuleWalker);
m = 200;%Number of poles

lag0Ndx = 2^16;

nDataPoints = 2^16-m-1;
X = zeros(nDataPoints, m);
for ii = 1:m
    X(:,ii) = acYuleWalker( lag0Ndx-ii+1:lag0Ndx+nDataPoints-ii );
end
b = -flip( acYuleWalker( lag0Ndx-nDataPoints:lag0Ndx-1 ) );

a = (X'*X)\X'*b;
a0 = acYuleWalker( lag0Ndx)+sum(acYuleWalker( lag0Ndx+1:lag0Ndx+m).*a);

z = exp(pi*1j*freqsBT/(freqsBT(end)));
Z = repmat(z.', 1, m);
for ii = 1:m
    Z(:,ii) = Z(:,ii).^ii;
end

S = a0./abs(1+Z*a).^2;
figure;semilogx( freqsBT, 10*log10(S) )
xlabel('Frequency,Hz')
ylabel('Spectral estimate, dB')
title('All-pole spectral estimate')
grid on;
% xlim([-50,50])
saveas( gcf,'all pole estimate.png' )
savefig(gcf, 'all pole esstimate.fig')