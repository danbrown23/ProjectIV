%% Fitness plotter

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

%% phenotype parametrisations

function r = GrowthRate(x)
% given a phenotype matrix, returns the growth rate vector of species x
    % r=exp(norm(x));
    r = 1;
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

function X = EffectiveCarryingCapacity(x,X0,var_K,var_a,tspan) % (allowing extinction)
% given a phenotype matrix and an initial population vector, returns the effective carrying capacity vector 
    r = GrowthRate(x);
    K = CarryingCapacity(x,var_K);
    alpha = Competition(x,var_a);
    [~, X] = ode45(@(t,X) populationdynamics(X,r,K,alpha),tspan,X0);
    X = X(end,:);
    X = X(:);
    X(X<0.001) = 0;
end

%% fitness

function f = Fitness(x,y,X,var_K,var_a)
    % given a phenotype matrix and a mutant matrix, returns the fitness
    % vector of the mutant species in an environment created by community x
    n = size(x,2);
    phenotype = [x y];
    X = X(:);
    % define the bottom left corner of the competition matrix we need
    alpha = Competition(phenotype,var_a);
    beta = alpha(n+1:2*n,1:n);
    f = GrowthRate(y).*(1-beta*X./CarryingCapacity(y,var_K));
end

%% parameter configuration

% configuration space dimensions
n0 = 1; % initial number of species
d = 1; % dimension of phenotype space 

% system parameters
var_a = 0.25; % competition var
var_K = 1; % carrying capacity var
divlim = 1000; % maximal diversity

% initial conditions
figure
x = -2:0.01:2;
y = -2:0.01:2;
for i = 1:length(x)
    xi = x(i);
    K = CarryingCapacity(xi,var_K);
    for j = 1:length(y)
        yj = y(j);
        f(i,j) = Fitness(xi,yj,K,var_K,var_a);
    end
end
[X, Y] = meshgrid(x,y);
surf(X,Y,f,'EdgeColor','none','FaceAlpha',0.7)
xlabel('resident')
ylabel('mutant')
zlabel('fitness')
hold on
% Use fourth input for color scale.
patch([2 -2 -2 2], [2 2 -2 -2], [0 0 0 0], 'white')  
hold off

figure
x0 = 0;
K = CarryingCapacity(x0,var_K);
for j = 1:length(y)
    yj = y(j);
    F(j) = Fitness(x0,yj,K,var_K,var_a);
end
plot(y,F, linewidth = 3)