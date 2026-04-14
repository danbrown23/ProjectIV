%% SUMMARY:
% ------------------------------------------------------------------------%
% Implements a stochastic trait substitution sequence to model evolution.
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
d = 10; % dimension of phenotype space 

% system parameters
var_a = 0.25; % competition variance
var_K = 1; % carrying capacity variance
divlim = 1000; % maximal diversity

% initial conditions
x0 = [0.4;0.4;0.4]; % (trait x species) phenotype matrix
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

tic
% trait substitution algorithm
for i = 2:steps
    % disp(["step =",i]); 
    % disp(["species =",n]); 
    if i-branchcount > branchtime - 1 && n < divlim
        % branch event: randomly selects a species to branch then generates
        % two mutants, splitting the population between them
        % disp("branched!")
        branchcount = i;
        nbound = nbound + 1;
        r = randi(n,1);
        branch_species = survivors(r);
        x(:,nbound) = x(:,branch_species) + 0.25*stepsize*randn(d,1);
        X(nbound) = 0.5*X(branch_species);
        % x(:,branch_species) = x(:,branch_species) + 0.5*stepsize*randn(d,1);
        X(branch_species) = 0.5*X(branch_species);
        survivors = [survivors,nbound]; % add the branched species
    end
    xeff = x(:,survivors);
    Xeff = X(survivors);
    n = length(survivors);
    % mutation
    y = xeff + stepsize*randn(d,n); % mutation step
    f = Fitness(xeff,y,Xeff,var_K,var_a);
    invaders = find(f>0); % find which species outcompete
    xeff(:,invaders) = y(:,invaders);
    % reincorporate to the full system
    x(:,survivors) = xeff;
    X(survivors) = Xeff;
    % check for ecological extinction
    X(survivors) = EffectiveCarryingCapacity(xeff,Xeff,var_K,var_a,tspan); 
    survivors = find(X~=0);
    % update data
    n = length(survivors); % current number of species
    xeff = x(:,survivors);
    Xeff = X(survivors);
    trajectory(:,survivors,i) = xeff; 
    PopulationData(survivors,i) = Xeff; 
    nlist(i) = n;
end
toc
% ------------------------------------------------------------------------%
%% plots
% ------------------------------------------------------------------------%
% number of species against time
figure
plot(1:steps,nlist)

% populations against time
figure
hold on
for i = 1:nbound
    y = PopulationData(i,:);
    plot(1:steps,y,'o');
end
hold off

% total population against time
figure
popsums = sum(PopulationData,1);
plot(1:steps,popsums);

% trajectory
% figure
% plot3(0,0,0,'o') 
% hold on
% scatter3(trajectory(1,:,1),trajectory(2,:,1),trajectory(3,:,1),'o') 
% scatter3(trajectory(1,:,end),trajectory(2,:,end),trajectory(3,:,end),'x')
% for j = 1:nbound
%     y = trajectory(:,j,1:steps);
%     plot3(y(1,:),y(2,:),y(3,:),linewidth=4);
% end    
% hold off


% trajectory movie
% M(steps) = struct('cdata',[],'colormap',[]);
% colors = viridis(n);
% f = figure;
% 
% for i = 1:steps
%     clf
%     plot3(0,0,0,'o',color='black')
%     hold on
%     plot3(x0(1),x0(2),x0(3),'o')
%     for j = 1:nbound
%         x = trajectory(1,j,1:i);
%         y = trajectory(2,j,1:i);
%         z = trajectory(3,j,1:i);
%         plot3(x,y,z,linewidth = 2,color=colors(j,:))
%         plot3(x(i),y(i),z(i),'o',MarkerFaceColor=colors(j,:),MarkerEdgeColor=colors(j,:),linewidth=3) 
%     end
%     axis equal
%     xlim([-2,2])
%     ylim([-2,2])
%     zlim([-2,2])
%     view(0.25*(45+i),20)
%     M(i) = getframe(f, [300 30 740 580]);
% end
% hold off
% 
% tss = VideoWriter('traitsubsequence3d', 'MPEG-4');
% tss.FrameRate=60;
% 
% open(tss)
% writeVideo(tss,M)
% close(tss)