%HW4_6.M
%AJG 11/16/17

%The Matlab file hw4_6.m and associated SOHO (Solar and Heliospheric Observatory) .jpg image files can be used to read and 
%analyze the pdf of four-band imagery of the sun. These data were observed using the SOHO Extreme ultraviolet Imaging Telescope (EIT) 
%on November 20, 2011 at ~0100 UTC. Run this code to load the images and
%perform a preliminary rendering and scatter analysis of the data recorded
%by EIT at the following UV wavelengths:
%  Plane 1: 171A 
%  Plane 2: 195A 
%  Plane 3: 284A 
%  Plane 4: 495A 
%For this homework build on this Matlab code to study the entropy of the imagery,
%develop a principle components analysis, estimate magnetic field anomolies, and detect 
%sunspots using statistical LLMSE estimation applied to EIT data. Once developed, your 
%algorithm should be tested on real (to-date) SOHO data that you can download at 
%http://sohowww.nascom.nasa.gov/data/realtime-images.html .
%Note that the test (unkwown) data will be assumed to be data from another
%instrument on SOHO, the SDO Helioseismic and Magnetic Imager (SDO/HMI), which
%observes at a longer wavelength (6173A) than does EIT.

clear
close all

%Load SOHO Extreme ultraviolet Imaging Telescope (EIT) images for 11/20/2011
lambda=[]; y=[]; x=[];
yo = imread('20111120_0100_eit171_512.jpg');
y=cat(3,y,sum(yo,3)/3);
lambda=[lambda;171];
yo = imread('20111120_0113_eit195_512.jpg');
y=cat(3,y,sum(yo,3)/3);
lambda=[lambda;195];
yo = imread('20111120_0106_eit284_512.jpg');
y=cat(3,y,sum(yo,3)/3);
lambda=[lambda;284];
yo = imread('20111120_0119_eit304_512.jpg');
y=cat(3,y,sum(yo,3)/3);
lambda=[lambda;304];
tystring='SOHO 20111120 0100-0119 EIT';
y=y(16:495,16:495,:);
[N1,N2,Nc]=size(y);

%Load SOHO Helioseismic and Magnetic Imager (SDO/HMI) images for 11/20/2011
xo = imread('20111120_0130_hmiigr_512.jpg');
x=cat(3,x,sum(xo,3)/3);
xo = imread('20111120_0130_hmimag_512.jpg');
x=cat(3,x,sum(xo,3)/3);
txstring='SOHO 20111120 0130 SDO/HMI';
x=x(16:495,16:495,:);
tstring='SOHO 20111120 0100-0130 UTC';

%Convert data to standard real double precision (8 bytes) for processing
xx=reshape(x,N1*N2,2);
yy=reshape(y,N1*N2,4);
%Estimate pdf's for each image plane
p=hist(yy,256)/(N1*N2);
p(find(p==0))=eps;

%Display false color images
figure('Color',[1,1,1]);
subplot(2,2,1);
image(uint8(cat(3,y(:,:,1),zeros(size(y(:,:,[1:2])))))); 
axis('equal')
v=axis; axis([1,N2,1,N1]);
title(['\lambda=',num2str(lambda(1)),char(197)]);
subplot(2,2,2);
image(uint8(cat(3,zeros(size(y(:,:,1))),y(:,:,2),zeros(size(y(:,:,1)))))); 
axis('equal')
v=axis; axis([1,N2,1,N1]);
title(['\lambda=',num2str(lambda(2)),char(197)]);
subplot(2,2,3);
image(uint8(cat(3,zeros(size(y(:,:,[1:2]))),y(:,:,3)))); 
axis('equal')
v=axis; axis([1,N2,1,N1]);
title(['\lambda=',num2str(lambda(3)),char(197)]);
subplot(2,2,4);
image(uint8(cat(3,y(:,:,4),y(:,:,4),zeros(size(y(:,:,1)))))); 
axis('equal')
v=axis; axis([1,N2,1,N1]);
title(['\lambda=',num2str(lambda(4)),char(197)]);
suptitle(tystring);

%Display histograms
fc=[[1,0,0];[0,1,0];[0,0,1];[1,1,0]];
figure('Color',[1,1,1]);
for ii=1:4
    subplot(2,2,ii);
    h=bar(p(:,ii),1);
    shading flat
    set(h,'facecolor',fc(ii,:));
    grid on
    title(['p_k: \lambda=',num2str(lambda(ii)),char(197),' (R)']);
end
suptitle(tystring);

%Scatter plot of image planes
figure('Color',[1,1,1]);
for ii=1:3
    for jj=(ii+1):4
        subplot(3,3,3*(ii-1)+jj-1);
        plot(yy(:,ii),yy(:,jj),'.');
        axis('equal')
        v=axis; axis([1,256,1,256]);
        xlabel(['\lambda=',num2str(lambda(ii)),char(197)]);
        ylabel(['\lambda=',num2str(lambda(jj)),char(197)]);
    end
end
suptitle(tystring);

%From here, use x,y,xx,yy,p,lambda,N1,N2,Nc, etc... to solve the following
%problems:

%1) From our discussion of the entropy of an image, find the information content in each of the image planes by
%calculating the expected number of bits that would be required to code %each pixel. (Note that this would be the 
%entropy of each bit plane defined using the natural logarithm divided %by ln(2).) How does this compare 
%to the number of bits required to encode random images (i.e., with values uniformly distributed over 
%all levels of each image plane)?

%2) Using the multispectral pixels as a set of four-dimensional vectors,
%compute the covariance matrix for each band. How does the black background affect the 
%covariance matrix calculation? Can you avoid using pixels with intensity value 3 or less? 

%3) Using the covariance matrix from (2), compute the principal compnent
%imagery. Render this data in both image form and using scatter plots as
%was rendered the original data. Comment on the correlation between the
%principal component amplitudes and imagery. What are the standard
%deviations of the principal component amplitudes?

%4) Render a single false-colored (RGB) image of the first three (most
%dominant) principal components? Comment on the utility of this image.

%5) Render the HMI continuum and magnetic anomaly images (planes 1 and 2 of
%data x, respectively) and compare visually with the rendered data from EIT
%and its principal components. 

%6) Compute the D-matrix to estimate the HMI continuum and magnetic anomaly
%data from the EIT data.

%7) Render a scatter plot of the continuum and magnetic anomaly data. Is
%the strong orthogonality condition achieved? The weak condition? 

%% Question 1
entropy = zeros(1,4);
for kk = 1:size(yy,2)
    entropy(kk) = -1*sum(p(:,kk).*log(p(:,kk)));
end
entropy = entropy*size(yy,1)/log(2);
uniformEntropy = log(size(yy,1))/log(2);
%% Question 2
mn = mean(yy);
% mnZero = yy-mn;
Sigma = cov(yy);
yy(yy<=3)=0;
SigmaNoBg = cov(yy);
%% Question 3
[V,L] = eig(Sigma);
pcs = (yy-mn)*V;
pcs = reshape(pcs, size(y));
figure;
subplot(2,2,1)
pcolor(1:480,1:480,pcs(:,:,4));shading flat;colorbar;
title('First principal component')
subplot(2,2,2)
pcolor(1:480,1:480,pcs(:,:,3));shading flat;colorbar;
title('Second principal component')
subplot(2,2,3)
pcolor(1:480,1:480,pcs(:,:,2));shading flat;colorbar;
title('Third principal component')
subplot(2,2,4)
pcolor(1:480,1:480,pcs(:,:,1));shading flat;colorbar;
title('Fourth principal component')
%scatterplot
figure;
pcsvec = reshape(pcs,[],4);
sgtitle('Principal component scatterplot')
for ii=1:3
    for jj=(ii+1):4
        subplot(3,3,3*(ii-1)+jj-1);
        plot(pcsvec(:,ii),pcsvec(:,jj),'.');
        axis('equal')
%         v=axis; axis([1,256,1,256]);
        xlabel(['PC ', num2str(ii)]);
        ylabel(['PC ', num2str(jj)]);
        grid on
    end
end

%% Question 4
f3pcs = pcs(:,:,2:end);
im2plot = sum(f3pcs,3);
figure;
image(uint8(flipud(im2plot)))

%% Question 6
%Need to correlate magnetic anomaly and HMI continuum data with the images
x(428:end,1:218,:)=0;
x(428:end,1:240,2)=0;
xxnew = reshape(x,[],2);%Removing the watermark from the image
%  pcsnorm = pcsvec./vecnorm(pcsvec);
Rxymag = zeros(1,4);
Rxyhmi=Rxymag;
for ii=1:4
    Rxymag(ii) = xcov( yy(:,ii),xxnew(:,2),0 , 'unbiased');
    Rxyhmi(ii) = xcov( yy(:,ii),xxnew(:,1),0 , 'unbiased');
end
Sigmainv = inv(SigmaNoBg);
Dmag = Rxymag*Sigmainv;
Dhmi = Rxyhmi*Sigmainv;

xhatMag = (Dmag*yy')';%xhatMag = uint8(reshape(xhatMag, 480,480));
xhatHmi = (Dhmi*yy')';%xhatHmi = uint8(reshape(xhatHmi, 480,480));

errMag = uint8(xhatMag-xxnew(:,2));errMag = reshape(errMag,480,480);
errHmi = uint8(xhatHmi-xxnew(:,1));errHmi = reshape(errHmi,480,480);

figure;
subplot(1,2,1)
plot(errMag(:), uint8(xxnew(:,2)), 'b.')
xlabel('Error');ylabel( 'Data');title('Magnetic anomaly data')
grid on
subplot(1,2,2)
plot(errHmi(:),uint8( xxnew(:,1)),'b.')
grid on
xlabel('Error'); ylabel('Data');title('HMI continuum data')
sgtitle('False-color image data/error scatterplot')

figure;imagesc(1:480,1:480, abs(errMag));colorbar;title('Absolute error in reconstructing magnetic anomaly data')
figure;imagesc(1:480,1:480,abs(errHmi));colorbar;title('Absolute error in reconstructing HMI continuum data')

%Q7
ynew=[];
yo = imread('eit 171 now.jpg');
ynew=cat(3,ynew,sum(yo,3)/3);
% lambda=[lambda;171];
yo = imread('eit 195 now.jpg');
ynew=cat(3,ynew,sum(yo,3)/3);
% lambda=[lambda;195];
yo = imread('eit 284 now.jpg');
ynew=cat(3,ynew,sum(yo,3)/3);
% lambda=[lambda;284];
yo = imread('eit 304 now.jpg');
ynew=cat(3,ynew,sum(yo,3)/3);
% lambda=[lambda;304];
% tystring='SOHO 20111120 0100-0119 EIT';
ynew=ynew(16:495,16:495,:);
[N1,N2,Nc]=size(ynew);

yynew = reshape(ynew, [],4);
yynew(yynew<=3)=0;
xnew = [];
x0 = imread('hmi now.jpg');
xnew=cat(3, xnew,sum(x0,3)/3);
x0 = imread('mag now.jpg');
xnew = cat(3, xnew,sum(x0,3)/3);
xnew=xnew(16:495,16:495,:);
xxnew = reshape(xnew, [],2);

xnewMag = (Dmag*yynew')';%xhatMag = uint8(reshape(xhatMag, 480,480));
xnewHmi = (Dhmi*yynew')';%xhatHmi = uint8(reshape(xhatHmi, 480,480));

errNewMag = uint8(xnewMag-xxnew(:,2));errNewMag = reshape(errNewMag,480,480);
errNewHmi = uint8(xnewHmi-xxnew(:,1));errNewHmi = reshape(errNewHmi,480,480);

figure;imagesc(1:480,1:480, abs(errNewMag));colorbar;title('Absolute error in reconstructing magnetic anomaly data')
figure;imagesc(1:480,1:480,abs(errNewHmi));colorbar;title('Absolute error in reconstructing HMI continuum data')

figure;imagesc(1:480,1:480,reshape(xnewMag,480,480));title('Predicted magnetic anomaly data')
figure;imagesc(1:480,1:480,reshape(xnewHmi,480,480));title('Predicted HMI data')