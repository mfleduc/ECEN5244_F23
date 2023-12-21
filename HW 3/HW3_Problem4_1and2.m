%% HW 3 Problem 4
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
D = (Sxx./(Sxx+Snn));
d = ifftshift((ifft2(D)));
xhat = conv2( data.dpn,d,'same' );
figure; subplot(1,3,3)
imagesc(abs(xhat));colormap('gray');title('Filtered');colorbar;
subplot(1,3,1)
imagesc(data.d);colormap('gray');title('Uncorrupted');colorbar;
subplot(1,3,2)
imagesc(data.dpn);colormap('gray');title( 'Corrupted' );colorbar;
sgtitle( 'Image smoothing: Noncausal Wiener Filter' )
noise_ncwf = var(data.d(:)-xhat(:));
signal_ncwf = var(data.d(:));
snr_ncwf = 10*log10(signal_ncwf/noise_ncwf);

multipliers = 2.^linspace( -2, 2, 15 );
snr_tmp = zeros(1,length(multipliers));
for ii = 1:length(multipliers)
    Snn_tmp = Snn*multipliers(ii);
    D_tmp = (Sxx./(Sxx+Snn_tmp));
    d_tmp = ifftshift((ifft2(D_tmp)));
    x_tmp = conv2( data.dpn,d_tmp,'same' );
    snr_tmp(ii) = 10*log10(signal_ncwf/var(x_tmp(:)-data.d(:)));
end
figure;
plot( multipliers*data.sigman^2, snr_tmp,'b', 'linewidth', 2 );
grid on;
xline(data.sigman^2, 'k', 'linewidth', 2)
xlabel('Assumed process variance')
ylabel('SNR achieved')
title('SNR achieved vs assumed variance')
legend( 'SNR', 'True process variance' )

