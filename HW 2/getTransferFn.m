function [a,b] = getTransferFn( w,H,nPoles, nZeros )
%Uses input w (radial frequency), H (values of the transfer function) data to approximate the given transfer function by a
%rational function with nPoles poles and nZeros zeros.\
%varargin: 1 if w is given in terms of radial frequency, 0 if w is given as
%values on the unit circle
%Returns a (vector of coefficients on the numerator) and b (vector of
%coefficients on the denominator)
N = length(w);
z = ( 1j*w );%Values on the unit circle
z = reshape(z,N,1);%Ensuring they are column vectors
H = reshape(H,N,1);

Az = zeros(N,nZeros+1);Ap = zeros(N,nPoles);

for ii = 0:nZeros
    Az(:,ii+1) = z.^(ii);
end
for ii = 1:nPoles
    Ap(:,ii) = H.*z.^(ii);
end
A = [Az,Ap];%The whole design matrix

A2 = [real(A);imag(A)];
H2 = [real(H);imag(H)];
x = (A2'*A2)\(A2'*H2);%MATLAB backslash to solve. Want real coefficients
a = x(1:nZeros+1);%Numerator coefficients
b = [1;-x(nZeros+2:end)];%Denominator coefficients
end