%% Homework 3 Problem 3
clear variables; close all
data = load( 'hw3_3.mat' );
N = size(data.Y,1);
EY = mean(data.Y).';
dat = data.Y-EY.';
Sigmahat = cov(dat);
Rhat = corr(dat);
%Data appears to be highly correlated
% figure;subplot(1,2,1)
% imagesc( Sigmahat );colorbar;
% title('Covariance matrix')
% subplot(1,2,2)
% imagesc( Rhat );colorbar;
% title('Correlation matrix')
[V,D] = eig( Sigmahat );%eigenvalue decomposition
% [~,ndcs] = (sort(diag(D)));
% V = V(flip(ndcs),flip(ndcs));
% D = D(flip(ndcs),flip(ndcs));
figure;plot( fliplr(V(:,end-2:end) ));
xlabel('i')
ylabel('v_i')
title('Three most dominant eigenmodes')
legend('1st mode', '2nd mode', '3rd mode')
grid on

f3PCS = dat*V(:,end-2:end); %First 3 principal components of the data

figure;plot( fliplr(f3PCS)) ;
xlabel('Sample')
ylabel('Principal component')
title('Three most dominant principal components')
legend('1st PC', '2nd PC', '3rd PC')
grid on