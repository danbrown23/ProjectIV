%% SUMMARY
% ------------------------------------------------------------------------%
% Plots the ecological functions featured in evolution of resource
% preference and niche width
%% PRELIMINARIES

clear

set(0,'defaultTextFontName', 'Arial')
set(0,'defaultaxesfontsize', 20); % 27 for 1X3, 20 for 1X2
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','latex');
set(0,'defaultAxesXGrid','off');
set(0,'defaultAxesYGrid','off');
set(0,'defaultAxesTickDir','out');
set(0,'defaultAxesLineWidth',1.5);

%% Evolution of Niche Width

function r = GrowthRate(x,s,R0,w,T)
% given a phenotype matrix, returns the growth rate vector of species x
    % r=exp(norm(x));
    r = (R0*w*s)/(sqrt(1+s^2))*exp(-1/2*(x.^2)./(1+s^2))-T;
end

function K = CarryingCapacity(x,var_K)
% Input: x is the (traits x species) phenotype configuration matrix
% Output: K is the carrying capacity vector of the n species
D2 = sum(x.^2,1);
K = exp(-D2./(2.*var_K));
K = K(:);
end

function alpha = Competition(x,var_a)
    % Input: x is the (traits x species) phenotype configuration matrix
    % Output: alpha is the (species x species) competition matrix. alpha_ij
    %  is the competition species i experiences from species j
    D2 = pdist2(x', x', 'euclidean').^2; % pairwise squared distances
    alpha = exp(-D2 ./ (2 * var_a)); % symmetric Gaussian competition
end


%% population dynamics

function dXdt = populationdynamics(X,r,K,alpha)
% given an instantaneous population vector and r, K, and alpha precalculated, returns the dynamics vector 
    dXdt = r.*X.*(1-alpha*X./K);
end

%% Parameter Regime

R0 =1;
w = 1;
T = 0.05;
s = 1;

%% Plots

x = -5:0.1:5;
r = GrowthRate(x,s,R0,w,T);

figure

plot(x,r,LineWidth=3)
xlabel('$x$')
ylabel('$r_X$')