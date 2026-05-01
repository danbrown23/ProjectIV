%% SUMMARY:
% ------------------------------------------------------------------------%
% Implements a stochastic trait substitution sequence to model evolution.
% ------------------------------------------------------------------------%
% Features:
%   -Extinction
%   -Population Tracking
%   -Branching via splitting a randomly chosen population in half every b
%   time steps
%   -Merging of sufficiently similar species
%   -The option for an artificially imposed maximal diversity 
% ------------------------------------------------------------------------%

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

%% parameter configuration

% configuration space dimensions
n0 = 1; % initial number of species
d = 1; % dimension of phenotype space 

% model parameters
var_a = 0.25; % competition var
var_K = 1; % carrying capacity var
A = 2;
omega = 0.0015;
tstart = 3000;
divlim = 600; % maximal diversity

% initial conditions

x0 = 0.1;
X0 = 0.1; % initial population vector

% numerical parameters

steps = 24000; % steps
stepsize = 0.01; % step size
tspan = [0,1000]; % ecological integration timespan
branchtime = 1; % steps between each branch event

% plotting parameters

frames = 30;
step0 = 5000;
stepend = 10000;
xwindow = [-4,4];
ywindow = [0,2.5];

%% phenotype parametrisations

function r = GrowthRate(x)
% given a phenotype matrix, returns the growth rate vector of species x
    % r=exp(norm(x));
    r = 1;
end

function K = CarryingCapacity(x,var_K,A,omega,t,tstart)
% Input: x is the (traits x species) phenotype configuration matrix
% Output: K is the carrying capacity vector of the n species
if t >= tstart
    t = t-tstart;
    n = size(x,2);
    v = A*sin(omega*t);
    v = v(:);
    v = repmat(v,1,n);
    D2 = sum((x-v).^2,1);
else
    D2 = sum(x.^2,1);
end
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

function X = EffectiveCarryingCapacity(x,X0,var_K,var_a,tspan,A,omega,t,tstart) % (allowing extinction)
% given a phenotype matrix and an initial population vector, returns the effective carrying capacity vector 
    r = GrowthRate(x);
    K = CarryingCapacity(x,var_K,A,omega,t,tstart);
    alpha = Competition(x,var_a);
    [~, X] = ode45(@(t,X) populationdynamics(X,r,K,alpha),tspan,X0);
    X = X(end,:);
    X = X(:);
    X(X<0.001) = 0;
end

%% fitness

function f = Fitness(x,y,X,var_K,var_a,A,omega,t,tstart)
    % given a phenotype matrix and a mutant matrix, returns the fitness
    % vector of the mutant species in an environment created by community x
    n = size(x,2);
    phenotype = [x y];
    X = X(:);
    % define the bottom left corner of the competition matrix we need
    alpha = Competition(phenotype,var_a);
    beta = alpha(n+1:2*n,1:n);
    f = GrowthRate(y).*(1-beta*X./CarryingCapacity(y,var_K,A,omega,t,tstart));
end

%% adaptive dynamics

% initialise


nmax = min(n0 + floor(steps/branchtime),divlim);
trajectory = NaN(d,nmax,steps);
PopulationData = zeros(nmax,steps);
nlist = NaN(steps);
CCClist = zeros(steps,1);
x = NaN(d,nmax);
X = zeros(1,nmax);
branchcount = 0;
mergecount = 0;

% set initial conditions
x(:,1:n0) = x0;
nbound = n0; % number of species in the system data
X(1:n0) = EffectiveCarryingCapacity(x0,X0,var_K,var_a,tspan,A,omega,1,tstart);
survivors = find(X~=0);
xeff = x(:,survivors); % surviving species currently in the system
Xeff = X(survivors);
trajectory(:,survivors,1) = xeff;
PopulationData(survivors,1) = Xeff;
nlist(1) = n0;
n = length(survivors); % number of species currently in the system
i=1;

tic
% trait substitution algorithm
while i < steps && n < divlim && n > 0
    i = i+1;
    disp(["step =",i]); 
    % disp(["species =",n]); 
    if i-branchcount > branchtime - 1 && n < divlim
        % branch event: randomly selects a species to branch then generates
        % two mutants, splitting the population between them
        % disp("branched!")
        branchcount = i;
        r = randi(n,1);
        branch_species1 = survivors(r);
        [~,extinct] = setdiff(x,x(survivors),'stable');
        if isempty(x(extinct)) == false
            branch_species2 = extinct(1);
            x(:,branch_species2) = x(:,branch_species1) + 0.1*stepsize*randn(d,1);
            X(branch_species2) = 0.5*X(branch_species1);
        else
            nbound = nbound+1;
            branch_species2 = nbound;
            x(:,branch_species2) = x(:,branch_species1) + 0.1*stepsize*randn(d,1);
            X(branch_species2) = 0.5*X(branch_species1);
        end
        % nbound = nbound+1;
        % branch_species2 = nbound;
        % x(:,branch_species2) = x(:,branch_species1) + 0.1*stepsize*randn(d,1);
        % X(branch_species2) = 0.5*X(branch_species1);
        x(:,branch_species1) = x(:,branch_species1) + 0.1*stepsize*randn(d,1);
        X(branch_species1) = 0.5*X(branch_species1);
        survivors = [survivors,branch_species2]; % add the branched species
    end
    xeff = x(:,survivors);
    Xeff = X(survivors);
    n = length(survivors);
    % mutation
    y = xeff + stepsize*randn(d,n); % mutation step
    f = Fitness(xeff,y,Xeff,var_K,var_a,A,omega,i,tstart);
    invaders = find(f>0); % find which species outcompete
    xeff(:,invaders) = y(:,invaders);
    % reincorporate to the full system
    x(:,survivors) = xeff;
    X(survivors) = Xeff;
    % check for extinction
    X(survivors) = EffectiveCarryingCapacity(xeff,Xeff,var_K,var_a,tspan,A,omega,i,tstart); % ecological extinction
    survivors = find(X~=0);
    % adjust data
    n = length(survivors); % current number of species
    xeff = x(:,survivors);
    Xeff = X(survivors);
    trajectory(:,survivors,i) = xeff; 
    PopulationData(survivors,i) = Xeff; 
    nlist(i) = n;
    if i > tstart
    CCClist(i) = A*sin(omega*(i-tstart));
    end
end
toc

%% plots

set(0,'defaultaxesfontsize', 27);
% number of species against time
figure
plot(1:steps,nlist)
xlabel("Time steps")
ylabel("Number of Distinct Species")

% carrying capacity against time
figure
plot(1:steps,CCClist)
xlabel("Time steps")
ylabel("Centre of the carrying capacity")
% % populations against time % DONT DO THIS FOR LONG TIMESPANS IT CRASHES
% MATLAB AND MY LAPTOP
% figure
% hold on
% for i = 1:nmax
%     y = PopulationData(i,:);
%     plot(1:steps,y);
% end
% hold off
% xlabel("Time steps")
% ylabel("Species Population Densities")

% total population against time
figure
%PopulationData(:,1)
popsums = sum(PopulationData,1);
plot(1:steps,popsums);
xlabel("Time steps")
ylabel("Total Population Density")

% population profiles

tic
% step0 = 5000;
% stepend = 10000;
frame = step0;
figure
x = trajectory(:,:,frame);
CCC = CCClist(frame);
h = histogram(x,-2+CCC:0.05:2+CCC);
e = h.BinEdges;
popavg = zeros(1,h.NumBins-1);
for i = 1:h.NumBins
    popavg(i) = 0;
    for j = 1:n
        if e(i) < x(:,j) && x(:,j) < e(i+1)
            popavg(i) = popavg(i) + PopulationData(j,frame);
        end
    end
end
scaledcounts = h.BinCounts.*popavg;
bar(e(1:end-1), scaledcounts, 'hist');
hold on
plot(CCC,0,'-o',LineWidth=4)
hold off
xlim(xwindow)
ylim(ywindow)
xlabel("Phenotype")
ylabel("Population Density")
title(frame)

for k = 2:frames
    steprange = stepend-step0;
    frame = k*floor(steprange/frames)+step0;
    figure
    x = trajectory(:,:,frame);
    CCC = CCClist(frame);
    h = histogram(x,-2:0.05:2);
    e = h.BinEdges;
    % weighting process
    popavg = zeros(1,h.NumBins-1);
    for i = 1:h.NumBins
        popavg(i) = 0;
        for j = 1:n
            if e(i) < x(:,j) && x(:,j) < e(i+1)
                popavg(i) = popavg(i) + PopulationData(j,frame);
            end
        end
    end
    scaledcounts = h.BinCounts.*popavg;
    bar(e(1:end-1), scaledcounts, 'hist');
    hold on
    plot(CCC,0,'-o',LineWidth=4)
    hold off
    xlim(xwindow)
    ylim(ywindow)
    xlabel("Phenotype")
    ylabel("Population Density")
    title(frame)
end
toc

%% animations

% % population profile movie
% M(steps) = struct('cdata',[],'colormap',[]);
% f = figure;
% 
% for frame = 1:steps
%     clf
%     x = trajectory(:,:,frame);
%     h = histogram(x,-2:0.05:2);
%     e = h.BinEdges;
%     popavg = zeros(1,h.NumBins-1);
%     for i = 1:h.NumBins
%         popavg(i) = 0;
%         for j = 1:nmax
%             if e(i) < x(:,j) && x(:,j) < e(i+1)
%                 popavg(i) = popavg(i) + PopulationData(j,frame);
%             end
%         end
%     end
%     scaledcounts = h.BinCounts.*popavg;
%     bar(e(1:end-1), scaledcounts, 'hist');
% 
%     xlim([-2,2])
%     ylim([0,4])
% 
%     M(frame) = getframe(f);
% end
% hold off
% tss = VideoWriter('Symmetric_Population_Distribution_Evolution10000', 'MPEG-4');
% tss.FrameRate=60;
% 
% open(tss)
% writeVideo(tss,M)
% close(tss)


% % 1d phenotype space
% 
% trait substitution sequence
% figure
% hold on
% for i = 1:nmax
%     y = trajectory(1,i,:);
%     plot(1:steps,y(:))
% end
% hold off

% Initialise:

% M(steps) = struct('cdata',[],'colormap',[]);
% colors = viridis(n);
% f = figure;
% 
% for i = 1:steps
%     clf
%     plot3(0,0,0,'o',color='black')
%     hold on
%     plot3(x0(1),x0(2),x0(3),'o')
%     for j = 1:n
%         x = squeeze(trajectory(1,j,1:i));
%         y = squeeze(trajectory(2,j,1:i));
%         z = squeeze(trajectory(3,j,1:i));
%         plot3(x,y,z,linewidth = 2,color=colors(j,:))
%         plot3(x(i),y(i),z(i),'o',MarkerFaceColor=colors(j,:),MarkerEdgeColor=colors(j,:),linewidth=3)% plot trajectory up to i
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
% tss = VideoWriter('traitsubsequence5', 'MPEG-4');
% tss.FrameRate=60;
% 
% open(tss)
% writeVideo(tss,M)
% close(tss)