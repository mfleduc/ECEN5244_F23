%% HW 3 Problem 4 part 3
clear variables; close all;
data = load( 'hw3_4.mat' );

Rxx = zeros(512);
ixm=257;iym=ixm;
for ii = 1:512
    for jj = 1:512
        ndx = 1+min( [255, round(sqrt((ii-ixm)^2+(jj-iym)^2))] );
        Rxx(ii,jj) = data.Rxx_est(ndx);
    end
end
% Rxx = 10*Rxx*round( Rxx(257,257)/10 )/Rxx(257,257);
Rxx = Rxx/Rxx(257,257);
Sxx = (abs(fft2(Rxx)));
Sxx = data.Rxx_est(1)*Sxx./sum(sum(Sxx));
%% Smoothing the image

Snn = ones( 512 )/512^2;
Snn = data.sigman^2*Snn;

W = psf2otf( data.pattern, size(data.dcpn) );

denom = abs(W).^2.*Sxx+Snn;
numer = conj(W).*Sxx;
D = numer./denom;

xhat = ifft2( D.*fft2(data.dcpn) );

figure;subplot(1,3,3)
imagesc(abs(xhat));colormap('gray');title('Deconvolved');colorbar;
subplot(1,3,1)
imagesc(data.d);colormap('gray');title('Original');colorbar;
subplot(1,3,2)
imagesc(data.dcpn);colormap('gray');title( 'Blurring and noise' );colorbar;
sgtitle( 'Image deconvolution: Noncausal Wiener Filter' )

noise_ncwf = var(data.d(:)-xhat(:));
signal_ncwf = var(data.d(:));
snr_ncwf = 10*log10(signal_ncwf/noise_ncwf);%4.76 dB: About a factor of 3/2 better
%%%%%

multipliers = 2.^linspace( -2,2, 25 );
snr_tmp = zeros(1,length(multipliers));
for ii = 1:length(multipliers)
    Snn_tmp = Snn*multipliers(ii);
    
    denom = abs(W).^2.*Sxx+Snn_tmp;
    numer = conj(W).*Sxx;
    D_tmp = numer./denom;
%     x_tmp = deconvwnr(data.dcpn, data.pattern, 1./SNR_tmp);
    x_tmp = ifft2( D_tmp.*fft2(data.dcpn) );
    snr_tmp(ii) = 10*log10(signal_ncwf/var(x_tmp(:)-data.d(:)));
end
figure;
plot( multipliers*data.sigman^2, snr_tmp,'b', 'linewidth', 2 );
grid on;
xline(data.sigman^2, 'k', 'linewidth', 2);
xlabel('Assumed process variance')
ylabel('SNR achieved')
title('SNR achieved vs assumed variance')
legend( 'SNR', 'True process variance' )



