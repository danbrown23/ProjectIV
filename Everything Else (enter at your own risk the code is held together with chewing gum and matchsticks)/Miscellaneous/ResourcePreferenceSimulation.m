%% PLOT STYLE

set(0,'defaultTextFontName', 'Arial')
set(0,'defaultaxesfontsize', 20); % 27 for 1X3, 20 for 1X2
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','latex');
set(0,'defaultAxesXGrid','off');
set(0,'defaultAxesYGrid','off');
set(0,'defaultAxesTickDir','out');
set(0,'defaultAxesLineWidth',1.5);

%% phenotype parametrisations

function r = growthrate(x,T)
% given a phenotype matrix, returns the growth rate vector of species x
    r = -T+1/sqrt(2)*exp(-x.^2./(4));
end

function K = carryingcapacity(x,T)
% given a phenotype matrix, returns the carrying capacity vector
    K = sqrt(2*pi)*sqrt(3)*exp(x.^2./3).*(-T+sqrt(2)^(-1)*exp(-x.^2./4));
    K = K(:);
end

function alpha = competition(x,i,j)
% given a phenotype matrix, returns the competition species i feels from species j
    alpha = exp((2.*x(:,i).^2-2.*x(:,i).^2+2.*x(:,i).*x(:,j)-(2).*x(:,j).^2)./(2.*3));
end

function alpha = competitionM(x,n)
% returns the competition matrix
alpha = zeros(n,n);
    for i = 1:n
        for j = 1:n
            alpha(i,j) = competition(x,i,j);  
        end
    end
end

%% population dynamics

function dXdt = populationdynamics(X,x,T,n)
% given a population vector and phenotype matrix, returns the dynamics vector 
    dXdt = growthrate(x,T).*X.*(1-competitionM(x,sd_a,n)*X./carryingcapacity(x,T));
end

% (assume no species go extinct)
function X = effectivecarryingcapacityM(x,T,n)
% given a phenotype matrix, returns the effective carrying capacity vector 
    carryingcapacity(x,T)
    X = inv(competitionM(x,n))*carryingcapacity(x,T);
end

%% fitness

function f = fitness(x,y,T,n)
    % given a phenotype matrix and a mutant matrix, returns the fitness
    % vector of the mutant species in an environment created by community x
    phenotype = [x y];
    beta = zeros(n,n);
    % define the bottom left corner of the competition matrix we need
    for j = n+1:2*n
        for k = 1:n
            beta(j-n,k) = competition(phenotype,j,k);
        end
    end
    % compute the fitness
    f = growthrate(y,T).*(1-beta*effectivecarryingcapacityM(x,T,n)./carryingcapacity(y,T));
end

%% generate n species with random phenotypes

% configuration space dimensions
n = 2; % number of species
d = 1; % dimension of phenotype space 

% parameter regime
w = 1; % mass per unit resource
R0 = 1; % resource carrying capacity mean 
s = 1; % resource carrying capacity standard deviation
T = 0.05; % threshhold mass of resource required to sustain consumer
r_R = 1; % resource growth rate in absence of consumer

% species initial conditions
x0 = 0.2*randn(d,n); % phenotype matrix. each column represents a species' phenotype
X0 = effectivecarryingcapacityM(x0,T,n); % population vector

%% adaptive dynamics

stepsize = 0.001; % step size
steps = 1000; % steps

% trait substitution sequence
x = x0;
trajectory = zeros(d,n,steps);
for i = 1:steps
    trajectory(:,:,i) = x;
    y = x+stepsize*randn(d,n);
    f = fitness(x,y,T,n);
    for j = 1:n
        if f(j) > 0
            x(:,j) = y(:,j);
        end
    end
end

figure

i = 1:steps;


hold on

plot(0,'o')
for j = 1:n
    plot(x0(j),'o')
end
for j = 1:n
    y = squeeze(trajectory(:,j,i)); 
    plot(i,y) 
end

hold off

