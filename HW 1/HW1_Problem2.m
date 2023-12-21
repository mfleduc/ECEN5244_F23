%% ECEN 5244 HW 1 Problem 2
clear variables;close all;
rng(452467365)
N = 500;%Number of points to simulate
R = diag([1,2,3]).^2;% Constructing covariance matrix
R(1,2) = 0.5*sqrt(R(1,1)*R(2,2));R(2,1)=R(1,2);
R(2,3) = -0.5*sqrt(R(2,2)*R(3,3)); R(3,2)=R(2,3);
%% Begin simulation
x = randn(3,N);
% L = chol(R, 'lower');
[E,Lambda] = eig(R);
z = E*sqrt(Lambda)*x;
figure;
subplotCounter = 1;
for ii=1:2
    for jj = ii+1:3
        subplot(1, 3, subplotCounter)
        scatter( z(ii,:), z(jj,:) )
        grid on
        xlabel(sprintf('x_{%d}', ii))
        ylabel(sprintf('x_{%d}', jj))
        subplotCounter = subplotCounter+1;
    end
end
sgtitle('Scatterplots of x_i and x_j for all pairs i,j')
% savefig( 'Q2_scatterplot.fig')
% saveas(gcf, 'Q2_scatterplot.png')

R_sample = z*z'/(N-1); %Sample covariance matrix

figure;imagesc( R );colorbar;caxis([-3,9])
figure;imagesc( R_sample );colorbar;caxis([-3, 9])