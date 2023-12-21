clear variables; close all;
rng(5684587)
data = load( 'hw2_2.mat' );
%% Wish to minimize error formula as given in the lecture 10 notes

H = data.d(:,2)+1j*data.d(:,3);

nZeros = 1;nPoles = 2; %Inputs: Number of desired poles and zeros

[a,b] = getTransferFn(data.d(:,1) ,H,nPoles, nZeros);
s = ( 1j*data.d(:,1 ) );%%s = i\omega!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Sz = zeros(length(s), nZeros+1);
for ii = 0:nZeros
    Sz(:,ii+1) = s.^(ii);
end
Hnum = Sz*a;

Sp = zeros(length(s), 1+nPoles);
for ii = 0:nPoles
    Sp(:,ii+1) = s.^(ii);
end
Hdenom = Sp*b;

Hest = Hnum./Hdenom;

% figure;
% sgtitle('Real and imaginary parts of given transfer function data with fit')
% subplot(2,1,1)
% 
% semilogx(data.d(:,1),real(H),'b','linewidth',2);
% hold on;semilogx( data.d(1:5:end,1),real(Hest(1:5:end)),'r.' )
% grid on
% legend('Actual real part','Estimated real part')
% title('Real part')
% xlabel('\omega')
% ylabel('Re(H)')
% subplot(2,1,2);
% semilogx(data.d(:,1),imag(H),'b','linewidth',2);
% hold on;semilogx( data.d(1:5:end,1),imag(Hest(1:5:end)),'r.' )
% grid on;
% legend('Actual imaginary part','Estimated imaginary part')
% title('Imaginary part')
% xlabel('\omega')
% ylabel('Im(H)')
% ylim([-2 2])
fprintf('No added noise\n');
fprintf('The transfer function has zeros at %.3f\n',roots(flip(a)));
poles = roots(flip(b));
fprintf('The transfer function has poles at %.3f and %.3f\n',poles(1),poles(2));

sigmas = (linspace(0.000000001,0.0005,1001));
poleEsts = zeros( length(sigmas),nPoles );
zeroEsts = zeros( length(sigmas),nZeros );
rhos = zeros( length(sigmas)+1,1 );
rhos(1) = abs(sum( conj(Hest).*H)) /sqrt( sum( abs(H).^2 )*sum(abs(Hest).^2)); 
for ii = 1:length(sigmas)
    noise = sigmas(ii)/sqrt(2)*(randn(1,length(H))+1j*randn(1,length(H))).';
    [a,b] = getTransferFn(data.d(:,1) ,H+noise,nPoles, nZeros);
    zeroEsts(ii) = roots(flip(a));
    if any(~isreal(b))
        keyboard
    end
    poleEsts(ii,:) = 1/2*[-b(2)+sqrt(b(2)^2-4*b(3)),-b(2)-sqrt(b(2)^2-4*b(3))]/b(3);
    Hfit = (Sz*a)./Sp*b;
    rhos(ii+1) = abs(sum( conj(Hfit).*H)) /sqrt( sum( abs(H).^2 )*sum(abs(Hfit).^2)); 
end

figure;
subplot(3,1,1)
plot(sigmas, (real((poleEsts))))
ylabel('Re(z)')
xlabel('\sigma')
title('Real part of the estimated poles')
grid on
subplot(3,1,2)
plot(sigmas, (imag((poleEsts))))
ylabel('Im(z)')
xlabel('\sigma')
title('Imaginary part of the estimated poles')
grid on
subplot(3,1,3)
plot(sigmas, zeroEsts);
grid on;
title('Estimated zeroes')
xlabel('\sigma')
ylabel('Root estimate')


figure;
subplot(3,1,1)
plot(sigmas, (real((poleEsts))))
ylabel('Re(z)')
xlabel('\sigma')
title('Real part of the estimated poles')
grid on
xlim([0,3*10^-5])
subplot(3,1,2)
plot(sigmas, (imag((poleEsts))))
ylabel('Im(z)')
xlabel('\sigma')
title('Imaginary part of the estimated poles')
grid on
xlim([0,3*10^-5])
subplot(3,1,3)
plot(sigmas, zeroEsts);
grid on;
title('Estimated zeroes')
xlabel('\sigma')
ylabel('Root estimate')
xlim([0,3*10^-5])
