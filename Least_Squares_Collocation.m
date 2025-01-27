
              %%% ***** Interpolation by inverse methode ***** %%%

% This code allows a user to evaluate the behavior of the solution based on the correlation length
% variation in the case of interpolation of a one variable problem.


% When the correlation length is small, the solution tends to be forced to pass
% through the observed points and elsewhere tends to resemble the a priori solution.
% The standard deviation of the solution is close to the standard deviation of the data error. 
% When the correlation length is small, the observed points cannot share their information
% with the neighboring points.
% However, if we start to increase the correlation length, the solution tends to become smoother and smoother.
% A smooth solution is a solution whose values do not change abruptly from one point to another. 
% The standard deviation a posteriori of the solution tends to be uniform around the region
% where the observed data are given. In this case, the information contained in precise observed data points
% can propagate to close interpolated points that are less precise.
% The standard deviation a posteriori of the solution at the point of lesser precision
% is lower than the standard deviation of the a priori error of the data.

% This script reproduce the solution of Olivier Francis :
% Introduction aux problemes inverses, 1991

% Programmed by -- Renaldo Sauveur

% Last Modified 06-August-2024


% Begin initialization code
delta = 3.1 ; % Correlation length
sigma_priori = 2 ; % A priori solution uncertainty
rest_point = 0:0.1:20 ; % Solution location
sigma = [0.2,0.4,0.4,0.1,0.2,0.2,0.3,0.25,0.3] ; % Standard deviation(Ecarts type) of the obseved data 
data = [1,1.5,1.25,0.85,0.25,0.5,1.35,0.4,-0.5] ; % observed data
dist_pt = [5,5.5,7.5,8.25,10,11,12.5,13,15] ; % Observed data location
  
% Variance matrix for observed data        
     ME = diag(sigma.*sigma) ; % error matrix
    
     
% Distance matrix for observed point
    D = zeros(length(data));
for i = 1:length(data)
    for j = 1:length(data)
           D(i,j)= abs(dist_pt(i) - dist_pt(j)) ;
    end
end
 
sig = sigma_priori*sigma_priori ;
del = delta*delta ;delo = -0.5/del ;
        
%covariance function 
    cov_dd = zeros(length(data));
for i = 1:length(data)     
    for j = 1:length(data)          
           cov_dd(i,j) = sig*exp(D(i,j)*D(i,j)*delo) ;   % Covariance matrix for observed data
    end    
end

S = ME + cov_dd ;
     
% Distance matrix between observed point and calculated point
    PD = zeros(length(data));
for i = 1:length(rest_point)
    for j = 1:length(data)
           PD(i,j)= abs(rest_point(i) - dist_pt(j)) ;
    end
end

   
% Covariance function between observed point and calculated point 
    cov_PD = zeros(length(data));
for i = 1:length(rest_point)   
    for j = 1:length(data)       
       cov_PD(i,j) = sig*exp(PD(i,j)*PD(i,j)*delo) ;           
    end    
end
   
   
% Symmetric definite positive function
opts.SYM = true ;
opts.POSDEF = true  ;
X = linsolve(S,data',opts) ;

% Final Solution 
solution = cov_PD*X ;

plot(rest_point,solution,'-b','Linewidth',5) ;
hold on ;
errorbar(dist_pt,data,sigma,'ok', 'MarkerFaceColor','k','MarkerSize',15,'LineWidth',3) ;
box on;
grid on;


% %%%%%%%%%%%%%%%%%%%%%%
% A priori uncertainty

     cov_DP = transpose(cov_PD) ;  
     
%function symmetric definite positive
opts.SYM = true ;
opts.POSDEF = true  ;
X_1 = linsolve(S,cov_DP,opts) ;   
    
% Matrix distance of calculated point 
    PP = zeros(length(rest_point));
for i = 1:length(rest_point)
    for j = 1:length(rest_point)
           PP(i,j)= abs(rest_point(i) - rest_point(j)) ;
    end
end
   
% covariance function of calculated point  
    cov_PP = zeros(length(rest_point));
for i = 1:length(rest_point)   
    for j = 1:length(rest_point)        
       cov_PP(i,j) = sig*exp(PP(i,j)*PP(i,j)*delo) ;           
    end    
end

% Covariance model
covariance_modele_1 = cov_PD*X_1 ;
covariance_modele_1_F = cov_PP - covariance_modele_1 ;     
hold on
 plot(rest_point,solution + diag(sqrt(covariance_modele_1_F)),':k','Linewidth',2)    ;    
hold on
 plot(rest_point,solution - diag(sqrt(covariance_modele_1_F)),':k','Linewidth',2)    ;  
 
 
title(strcat('Length of correlation =',num2str(delta))) ;
xlabel('Data location', 'FontWeight', 'bold','FontSize',10) ;
ylabel('Observed data', 'FontWeight', 'bold','FontSize',10) ; 
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',20,'FontWeight','bold','LineWidth',3) ;
set(gca,'XTickLabelMode','auto') ;
 
ylim([-3.25 2.75]) ; % set value limit into x-axis
xlim([0 20]) ;% set distance limit into x-axis
yline(0,'--g','A priori solution','LineWidth',3,'FontSize',16) ;
yline(-2,'--r','A priori Solution uncertainty','LineWidth',3,'FontSize',16) ;
yline(2,'.--r','A priori Solution uncertainty','LineWidth',3,'FontSize',16) ;
legend('Solution','Observed data','Solution +- Standard deviation','Location','south','FontSize',17) ;

% Plot margin cut
pos = get(gca, 'Position');
pos(1) = 0.045; %Left
pos(2) = 0.105;  %Bottom
pos(3) = 0.94;  %Right
pos(4) = 0.84;  %Top
set(gca, 'Position', pos)
 
 

