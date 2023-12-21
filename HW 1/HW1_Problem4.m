%% ECEN 5244 HW 1 Problem 4
% dbstop if nan
clear variables;close all
rng(24572467)
data = load( 'hw1_4.mat' );
sigma = 0.5;
%%%%%Model is a1*cos(a2*x)*exp(a3*x)
%% First: Apply the steepest descent method
b = [1, 1, -1].';%[2.9, 2, -0.2].';
del = 1;
tol = 1e-3;
cnt=0;
while cnt<100 && del>tol 
    %First: Find the feasible region with vanilla gradient descent
    cnt=cnt+1;
    gradf = getGradient( b, data )/sigma^2;
    db = -0.0001*gradf;
    b = b+db;
    del = max(abs(db));
end
tol = tol/1000;
cnt = 0;
hessMat = eye(3);
while cnt<500 && del>tol && cond(hessMat)<10000 %Now use big boy gradient descent
    cnt=cnt+1;
    gradf = getGradient( b, data )/sigma^2;
    hessMat = getHessian(b, data)/sigma^2;
    db = -hessMat\gradf;
    b = b+db;
    del = max(abs(db));
end
chi2_steep = getChi2(b, data, sigma);
%% Second: Apply simulated annealing
a = linspace(1, 7, 2000);%Amplitude
w = linspace(0, 4, 2000);%Frequency
d = linspace(-2, 0, 2000);%Damping
ind = randi( 2000,[1,3]);
z = [a(ind(1)), w(ind(2)), d(ind(3))];%Parameter vector
maxIters = 1e5;
t = 1:maxIters;
cooling = 1./sqrt(t); %Cooling schedule
pHop = 0.9; %Transition probability
chiVals = zeros([1, maxIters+1]);
chiVals(1) = getChi2(z, data, 0.5);
for jj = t
    delind = (randi(3,[1,3])-2);%Choosing a random neighbor
    if all(delind==0)
        delind = (randi(3,[1,3])-2);
    end
    ind2 = ind+delind;
    mask = (ind2>2000);
    ind2(mask) = 2000;
    chiVals(jj+1) = getChi2( [a(ind2(1)), w(ind2(2)), d(ind2(3))], data, 0.5 );
    changeProb = pHop*exp(-1/cooling(jj)*max(0, chiVals(jj+1)-chiVals(jj)));
    changeBool = rand(1)<changeProb;%Update chi2 and switch based on changeProb
    if changeBool
        ind = ind2;
    end
end

saf = figure;semilogx( t, chiVals(2:end), 'k', 'Linewidth', 2 );
grid on;
xlabel('Time')
ylabel('$\chi^2$', 'interpreter', 'latex')
title('$\chi^2$ vs time step, simulated annealing','interpreter', 'latex')
saveas(saf, 'sim_anneal.png')
savefig( saf, 'sim_anneal.fig' )
b_sa = [a(ind(1)), w(ind(2)), d(ind(3))];%Final fitted parameters
%% Last: Apply the Prony Method
%Prefilter the data. Sensitive to the filter!
filtLength = 50;
tFilt = bartlett(filtLength*2+1);%abs(sawtooth( linspace(-pi, pi, 2*filtLength+1 )));
wss = (sum(tFilt.^2));
tFilt = tFilt/wss;
% if ~iscola(tFilt, 2*filtLength)
%     keyboard
% end
data.y = conv(tFilt, data.y);
data.y = data.y(2*filtLength+1:end-2*(filtLength));
%M=4
dx = mean(diff(data.x));
M = 2;
A = [data.y(1:end-2),data.y(2:end-1)];
beta = -(A'*A)\A'*data.y(3:end);
uj = roots([1,flip(beta.')]);
alpha = log( uj )/dx;
X = exp(data.x(filtLength+1:end-filtLength)*alpha.');

amp = (X'*X)\X'*data.y(1:end);
expPart = [];
for ii = 1:M
    expPart(ii,:) = exp(alpha(ii)*data.x(filtLength+1:end-filtLength));
end
sigma_prony = sqrt(1/wss)*sigma ;

chi2_prony = 1/sigma_prony^2*sum(abs(data.y-amp(1)*exp(alpha(1)*data.x(filtLength+1:end-filtLength))...
    -amp(2)*exp( alpha(2)*data.x(filtLength+1:end-filtLength) )).^2);

d2 = load('hw1_4.mat');
f = figure;plot(d2.x, d2.y,'bo','markersize',10);hold on;
plot(data.x(filtLength+1:end-filtLength), data.y, 'g', 'linewidth',2);
plot(data.x(filtLength+1:end-filtLength),(sum((amp.'.*expPart.').') ),'k', 'linewidth',2);
plot(data.x, b(1)*cos(b(2)*data.x).*exp(b(3)*data.x), 'r', 'linewidth',2);
plot(data.x, b_sa(1)*cos(b_sa(2)*data.x).*exp(b_sa(3)*data.x), 'm', 'linewidth',2);
grid on;
xlabel('x')
ylabel('y')
title('Fitted models vs actual data')
legend('Unfiltered data', 'Filtered data', 'Fit from Prony''s algorithm', 'Fit using steepest descent','Fit using simulated annealing')

%% Functions
function gradF = getGradient( b, data )%Taking advantage of the known model here
x = data.x;
y = data.y;
gF = zeros( length(b), length(x) );
F = y-b(1)*exp(b(3)*x).*cos(b(2)*x);
gF(1,:) = -cos(b(2)*x).*exp(b(3)*x).*F;
gF(2,:) = b(1)*x.*sin( b(2)*x).*exp(b(3)*x).*F;
gF(3,:) = -b(1)*x.*cos( b(2)*x).*exp(b(3)*x).*F;
gradF = 2*sum( gF.' ).';
end
function chi2 = getChi2(b, data, sigma)
chi2 = sum( abs(data.y-b(1)*exp(b(3)*data.x).*cos(b(2)*data.x)).^2)/sigma^2;
end
function H = getHessian( b, data )
x=data.x;y=data.y;
F = y-b(1)*cos(b(2)*x).*exp(b(3)*x);
H = zeros(length(b),length(b),length(data.x));
H(1,1,:) = 2*cos(b(2)*x).^2.*exp(2*b(3)*x);
H(1,2,:) = 2*x.*sin(b(2)*x).*exp(b(3)*x).*F-2*x*b(1).*exp(2*b(3)*x).*cos(b(2)*x).*sin(b(2)*x);
H(2,1,:) = H(1,2,:);
H(1,3,:) = -2*x.*cos( b(2)*x ).*exp(b(3)*x).*F+2*x*b(1).*cos(b(2)*x).^2.*exp(2*b(3)*x);
H(3,1,:) = H(1,3,:);
H(2,2,:) = 2*b(1)*x.^2.*cos(b(2)*x).*exp(b(3)*x).*F+2*b(1)^2*x.^2.*sin(b(2)*x).^2.*exp(2*b(3)*x);
H(2,3,:) = 2*b(1)*x.^2.*sin(b(2)*x).*exp(b(3)*x).*F-2*b(1)^2*x.^2.*sin(b(2)*x).*cos(b(2)*x).*exp(2*b(3)*x);
H(3,2,:) = H(2,3,:);
H(3,3,:) = -2*b(1)*x.^2.*cos(b(2)*x).*exp(b(3)*x).*F+2*b(1)^2*x.^2.*cos(b(2)*x).^2.*exp(2*b(3)*x);
H = sum(H,3);
end