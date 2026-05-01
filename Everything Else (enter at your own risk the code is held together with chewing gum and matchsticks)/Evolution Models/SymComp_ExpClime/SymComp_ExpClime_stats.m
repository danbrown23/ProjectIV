%% PRELIMINARIES

clear

set(0,'defaultTextFontName', 'Arial')
set(0,'defaultaxesfontsize', 27); % 27 for 1X3, 20 for 1X2
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

% system parameters
var_a = 0.25; % competition var
var_K = 1; % carrying capacity var
% omegalist = [0.001,0.005,0.1];
% A = 3;
omega = 0.003;
Alist = [2, 3, 4];
tstart = 3000;
divlim = 600; % maximal diversity

% initial conditions

x0 = 0;
X0 = 0.1.*ones(n0,1); % initial population vector

% numerical parameters

stepsize = 0.01; % step size
steps = 10000; % steps
tspan = [0,1000]; % ecological integration timespan
branchtime = 1; % steps between each branch event

% plotting parameters

stepend = 4000;
step0 = 2000;
fstart = 500;
fend = 1000;
res = 0.001;
w = 3;

function beta = velocity(omega,A,t,tstart)
if t > tstart
    beta = A*sin(omega*(t-tstart));
else
    beta = 0;
end
end

%% phenotype parametrisations

function r = GrowthRate(~)
% given a phenotype matrix, returns the growth rate vector of species x
    % r=exp(norm(x));
    r = 1;
end

function K = CarryingCapacity(x,var_K,omega,A,t,tstart)
% Input: x is the (traits x species) phenotype configuration matrix
% Output: K is the carrying capacity vector of the n species
n = size(x,2);
v = velocity(omega,A,t,tstart)/A;
amp = velocity(omega,A,t-pi/2,tstart)+2*A;
v = v(:);
v = repmat(v,1,n);
if t > tstart
    D2 = sum((x-v).^2,1);
else
    D2 = sum(x.^2,1);
end
K = amp*exp(-D2./(2.*var_K))-A;
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

function X = EffectiveCarryingCapacity(x,X0,var_K,var_a,tspan,omega,A,t,tstart) % (allowing extinction)
% given a phenotype matrix and an initial population vector, returns the effective carrying capacity vector 
    r = GrowthRate(x);
    K = CarryingCapacity(x,var_K,omega,A,t,tstart);
    alpha = Competition(x,var_a);
    [~, X] = ode45(@(t,X) populationdynamics(X,r,K,alpha),tspan,X0);
    X = X(end,:);
    X = X(:);
    X(X<0.005) = 0;
end

%% fitness

function f = Fitness(x,y,X,var_K,var_a,omega,A,t,tstart)
    % given a phenotype matrix and a mutant matrix, returns the fitness
    % vector of the mutant species in an environment created by community x
    n = size(x,2);
    phenotype = [x y];
    X = X(:);
    % define the bottom left corner of the competition matrix we need
    alpha = Competition(phenotype,var_a);
    beta = alpha(n+1:2*n,1:n);
    f = GrowthRate(y).*(1-beta*X./CarryingCapacity(y,var_K,omega,A,t,tstart));
end

%% adaptive dynamics

popdata = zeros(steps,length(Alist));
ndata = zeros(steps,length(Alist));
mdata = NaN(d,steps,length(Alist));
tdata = zeros(steps,length(Alist));
tic

for sim = 1:length(Alist)

    A = Alist(sim);

    % initialise

    nmax = min(n0 + floor(steps/branchtime),divlim);
    trajectory = NaN(d,nmax,steps);
    PopulationData = zeros(nmax,steps);
    nlist = NaN(steps,1);
    mlist = NaN(d,steps);
    msh = length(-w:res:w);
    flist = NaN(msh,steps);
    CCClist = zeros(steps,1);
    x = NaN(d,nmax);
    X = zeros(1,nmax);
    branchcount = 0;
    mergecount = 0;
    xpts = NaN(msh,steps);

    % set initial conditions
    x(:,1:n0) = x0;
    nbound = n0; % number of species in the system data
    X(1:n0) = EffectiveCarryingCapacity(x0,X0,var_K,var_a,tspan,omega,A,1,tstart);
    survivors = find(X~=0);
    xeff = x(:,survivors); % surviving species currently in the system
    Xeff = X(survivors);
    trajectory(:,survivors,1) = xeff;
    PopulationData(survivors,1) = Xeff;
    nlist(1) = n0;
    mlist(:,1) = mean(x,2);
    n = length(survivors); % number of species currently in the system
    i=1;
    xpts(:,1) = -w:res:w;

    tic
    % trait substitution algorithm
    while i < steps && n < divlim && n > 0
        i = i+1;
        disp(["step =",i]);
        disp(["species =",n]);
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
        f = Fitness(xeff,y,Xeff,var_K,var_a,omega,A,i,tstart);
        invaders = find(f>0); % find which species outcompete
        xeff(:,invaders) = y(:,invaders);
        % reincorporate to the full system
        x(:,survivors) = xeff;
        X(survivors) = Xeff;
        % check for extinction
        X(survivors) = EffectiveCarryingCapacity(xeff,Xeff,var_K,var_a,tspan,omega,A,i,tstart); % ecological extinction
        survivors = find(X~=0);
        % adjust data
        n = length(survivors); % current number of species
        xeff = x(:,survivors);
        Xeff = X(survivors);
        trajectory(:,survivors,i) = xeff;
        PopulationData(survivors,i) = Xeff;
        nlist(i) = n;
        % CCC = v.*i;
        % xpts(:,i) = -w+CCC:res:w+CCC;
        % size(xpts(:,i))
        % size(xeff)
        % flist(:,i) = Fitness(xeff,xpts(:,i)',Xeff,var_K,var_a,omega,A,i,tstart);
        if i > tstart
            CCClist(i) = velocity(omega,A,i,tstart);
        end
        mlist(:,i) = mean(xeff,2);
    end
    toc
    popsums = sum(PopulationData,1);
    popdata(:,sim) = popsums;
    ndata(:,sim) = nlist;
    mdata(:,:,sim) = mlist;
    tdata(:,sim) = CCClist;
    figure
    hold on
    for i = 1:divlim
        y = trajectory(1,i,:);
        plot(1:steps,y(:),linewidth = 4)
        plot(1:steps,CCClist(1:steps))
    end
    hold off
end
%% plots

% number of species against time
figure
for i = 1:length(Alist)   
    % total population
    hold on
    plot(1:steps,popdata(:,i));
    xlabel("Time Steps")
    ylabel("Total Population")
end
hold off

figure
for i = 1:length(Alist)
    % number of species
    hold on
    plot(1:steps,ndata(:,i))
    xlabel("Time Steps")
    ylabel("Number of Distinct Species")
end
hold off
for i = 1:length(Alist)
    figure
    hold on
    plot(1:steps,mdata(1,:,i))
    plot(1:steps,tdata(:,i),color = 'black',linewidth = 3)
    hold off
end

