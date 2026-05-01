clear
%% Summary

% Implements the Gillespie SSA algorithm for the logistic competition model

%% The Model

function r = GrowthRate(x)
% given a phenotype matrix, returns the growth rate vector of species x
    % r=exp(norm(x));
    N = size(x,2);
    r = ones(N,1);
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
    alpha = 1-D2/var_a; % symmetric Gaussian competition
    % remove interspecies competition
    alpha = alpha - eye(size(alpha));
end

%% Parameter Configuration

% configuration space dimension
dim = 1;

% initial conditions
N = 400;
X0 = ones(N,1);
x0 = 0.1*randn(dim,N)+1;

% ecological parameters
var_a = 0.5;
var_K = 0.25;

% numerical parameters
var_m = 0.5;
sd_m = sqrt(var_m);
Nmax = 2000;
stepsmax = 2000;
threshhold = 10;

%% Gillespie Algorithm

% initialise
timedata = NaN(stepsmax,1);
popudata = zeros(Nmax,stepsmax);
phenodata = NaN(dim,Nmax,stepsmax);
X = popudata(:,1);
x = phenodata(:,:,1);
Nlist = NaN(stepsmax,1);

% set initial conditions
survivors = 1:N;
t = 0;
X(survivors) = X0;
x(:,survivors) = x0;
timedata(1) = t;
popudata(:,1) = X;
phenodata(:,:,1) = x;
step = 1;

% algorithm
while step < stepsmax & N < Nmax & N > threshhold
    step
    % 1. Generate random numbers r1,r2 uniformly:
    r1 = rand;
    r2 = rand;
    % 2. Compute the propensity function of the system:
    xeff = x(:,survivors);
    Xeff = X(survivors);
    r = GrowthRate(xeff);
    % d = r./CarryingCapacity(xeff,var_K).*sum(Competition(xeff,var_a),2);
    d = 1/N.*r./CarryingCapacity(xeff,var_K).*sum(Competition(xeff,var_a),2);
    R = sum(r);
    D = sum(d);
    propensity = R + D;
    % 3. Compute the next reaction time
    tau = 1/propensity*log(1/r1);
    t = t + tau;
    % 4. Find which of the n^2 reactions took place
    if propensity*r2 < R %birth
        disp("birth")
        % find the reproducing individual
        reproducer = survivors(randi(N));
        % find a population that has gone extinct for memory
        offspring = find(~X,1);
        % put in the baby
        X(offspring) = 1;
        x(:,offspring) = x(:,reproducer) + sd_m*randn(dim,1);
        survivors = [survivors, offspring];
    else % death
        disp("death")
        % find the individual destined to die
        i=1;
        dsum = R + d(i);
        while r2*propensity > dsum
            i = i+1;
            dsum = dsum + d(i);
        end
        % kill him
        deadeff = i;
        dead = survivors(i);
        X(dead) = 0;
        x(:,dead) = NaN(dim,1);
        survivors(i) = [];
    end
    % update system history
    step = step + 1;
    timedata(step) = t;
    popudata(:,step) = X;
    phenodata(:,:,step) = x;
    N = length(survivors);
    Nlist(step) = N;
end

tfinal = timedata(step);
frames = 10;
tframe = tfinal/frames;
stepframe = floor(step/frames);
% for f = 1:frames
%     figure
%     t = tframe*f;
%     i=1;
%     while t > timedata(i)
%         i = i+1;
%     end
%     xdata = phenodata(1,:,i);
%     ydata = phenodata(2,:,i);
%     scatter(xdata,ydata,'filled')
%     xlim([-5,5])
%     ylim([-5,5])
% end

figure
plot(timedata,Nlist)

% for f = 1:frames
%     figure
%     i = f*stepframe;
%     xdata = phenodata(1,:,i);
%     ydata = phenodata(2,:,i);
%     scatter(xdata,ydata,'filled')
%     xlim([-5,5])
%     ylim([-5,5])
% end

for f = 1:frames
    figure
    t = tframe*f;
    i=1;
    while t > timedata(i)
        i = i+1;
    end
    x = phenodata(:,:,i);
    histogram(x,-3:0.2:3)
end