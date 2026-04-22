%% SUMMARY:
% ------------------------------------------------------------------------%
% Implements the canonical equation of adaptive dynamics to model evolution.
%
% Features:
%   -Extinction
%   -Population Tracking
%   -Branching via splitting a randomly chosen population in half every b
%   time steps
%   -Merging of sufficiently similar species
%   -The option for an artificially imposed maximal diversity 
% ------------------------------------------------------------------------%
%% PRELIMINARIES
% ------------------------------------------------------------------------%
clear

set(0,'defaultTextFontName', 'Arial')
set(0,'defaultaxesfontsize', 20); % 27 for 1X3, 20 for 1X2
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','latex');
set(0,'defaultAxesXGrid','off');
set(0,'defaultAxesYGrid','off');
set(0,'defaultAxesTickDir','out');
set(0,'defaultAxesLineWidth',1.5);

import Functions.viridis;
% ------------------------------------------------------------------------%
%% Phenotype Parametrisations
% ------------------------------------------------------------------------%
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
% ------------------------------------------------------------------------%
%% Parameter Configuration
% ------------------------------------------------------------------------%
% configuration space dimensions
n0 = 1; % initial number of species
d = 1; % dimension of phenotype space 

% system parameters
var_a = 0.25; % competition variance
var_K = 1; % carrying capacity variance
divlim = 1000; % maximal diversity

% initial conditions
x0 = randn(d,n0); % (trait x species) phenotype matrix
X0 = 0.1.*ones(n0,1); % initial population vector

% numerical parameters
frames = 5; % surveys of the population distribution taken
stepsize = 0.01; % variance of the random normal mutation step
steps = 1000; % evolutionary timespan
tspan = [0,1000]; % ecological integration timespan
branchtime = 1; % steps between each branch event
% ------------------------------------------------------------------------%
%% Adaptive Dynamics
% ------------------------------------------------------------------------%
% population dynamics
function dXdt = populationdynamics(X,r,K,alpha)
% given an instantaneous population vector and r, K, and alpha 
% precalculated, returns the dynamics vector 
    dXdt = r.*X.*(1-alpha*X./K);
end

function X = EffectiveCarryingCapacity(x,X0,var_K,var_a,tspan) % (allowing 
% extinction) given a phenotype matrix and an initial population vector, returns the 
% effective carrying capacity vector 
    r = GrowthRate(x);
    K = CarryingCapacity(x,var_K);
    alpha = Competition(x,var_a);
    [~, X] = ode45(@(t,X) populationdynamics(X,r,K,alpha),tspan,X0);
    X = X(end,:);
    X = X(:);
    X(X<0.001) = 0;
end

% fitness
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

% initialise

nmax = n0 + floor(steps/branchtime);
trajectory = NaN(d,nmax,steps);
PopulationData = zeros(nmax,steps);
nlist = NaN(steps);
x = NaN(d,nmax);
X = zeros(1,nmax);
branchcount = 0;

% set initial conditions
x(:,1:n0) = x0;
nbound = n0; % number of species in the system data
X(1:n0) = EffectiveCarryingCapacity(x0,X0,var_K,var_a,tspan); 
survivors = find(X~=0);
xeff = x(:,survivors); % surviving species currently in the system
Xeff = X(survivors);
trajectory(:,survivors,1) = xeff;
PopulationData(survivors,1) = Xeff;
nlist(1) = n0;
n = length(survivors); % number of species currently in the system

