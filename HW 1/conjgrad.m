function x0CG = conjGradLsq(W,obs, tol, maxIters)
%% Conjugate gradient time!
x0CG = ones(size(x(:)));%randn(numel(x), 1 );%figure;imagesc(reshape(x0, size(x)));
Wx = W*x0CG;
p = W.'*obs(:) - W.'*Wx;
% p = W*tmp;
r = p;
k = 0;
alpha = 0;
tol = 10^-6;
dx = 1+tol;
maxIters = max(size(W));
resid = zeros(1, maxIters+1);
resid(1) = sqrt(sum( r.^2 ));

while dx > tol && k<maxIters
    k=k+1;
    Wp = W*p;
    alpha = resid(k)^2/sum(Wp.^2);
    x0CG = x0CG+alpha*p;
    dx = mean( (alpha*p).^2 );
    r1 = r - alpha*( W.'*Wp );
    resid(k+1) = sqrt(sum(r1.^2));
    beta = resid(k+1)^2/resid(k)^2;
    p = r1+beta*p;
    r = r1;
    
end

if k ==maxIters
    fprintf('Desired tolerance not reached, loop reached max iterations\n');
else
    fprintf('Reached desired tolerance in %d iterations\n',k);
end