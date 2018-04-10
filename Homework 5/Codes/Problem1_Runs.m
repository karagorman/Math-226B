% Math 226B - Homework #5
% Problem 1
% Script to run Poisson2DGMRES for all 5 gamma values
clear all
gammaVec = [1,10,50,100,1000];
m=100;
h = 1/(m+1);


for i = 1:length(gammaVec)
    figure(1)
    hold on
    subplot(2,3,i)
    [u,relres,resvec] = Poisson2DGMRES(m,gammaVec(i));
    x=h:h:1-h;
    y=h:h:1-h; y=y';
    surf(x,y,reshape(u,[m,m]),'EdgeColor','none'); 
    %view(2); 
    %colormap hot
    colorbar
    title(strcat("\gamma = ", num2str(gammaVec(i))))
    set(gca,'FontSize',15)
end
hold off

for i = 1:length(gammaVec)
    
    figure(2)
    hold on
    [u,relres,resvec] = Poisson2DGMRES(m,gammaVec(i));
    plot(1:length(resvec),log(resvec),'lineWidth',2)
    
end

legend('\gamma=1','\gamma=10','\gamma=50','\gamma=100','\gamma=1000')
xlabel('Number of Iterations')
ylabel('Log of Relative Residual')
title('Relative Residual of GMRES, with Preconditioner A_0')


