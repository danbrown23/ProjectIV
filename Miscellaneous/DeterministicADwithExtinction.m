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

function r = growthrate(x) 
% given a phenotype matrix, returns the growth rate vector of species x
    r = 1; % no
end

function K = carryingcapacity(x,i,sd_K)
% given a phenotype matrix, returns the carrying capacity of species i
    K = exp(-norm(x(:,i)).^2./(2.*sd_K^2));
end

function alpha = competition(x,i,j,sd_a)
% given a phenotype matrix, returns the competition species i feels from species j
    alpha = exp((-norm(x(:,i)-x(:,j)).^2)./(2.*sd_a^2)); % symmetric competition
end

%% ecological functional parameters

function alpha = competitionM(x,sd_a,n)
% returns the competition matrix
alpha = zeros(n,n);
    for i = 1:n
        for j = 1:n
            alpha(i,j) = competition(x,i,j,sd_a);  
        end
    end
end

function K = carryingcapacityM(x,sd_K,n)
% returns the carrying capacity vector
    K = zeros(n,1);
    for i = 1:n
        K(i) = carryingcapacity(x,i,sd_K);
    end
end

%% population dynamics

function dXdt = populationdynamics(X,x,sd_K,sd_a,n)
% given an instantaneous population vector and phenotype matrix, returns the dynamics vector 
    dXdt = growthrate(x).*X.*(1-competitionM(x,sd_a,n)*X./carryingcapacityM(x,sd_K,n));
end

function X = effectivecarryingcapacity(x,X0,sd_K,sd_a,n,tspan) % (allowing extinction)
% given a phenotype matrix and an initial population vector, returns the effective carrying capacity vector 
    [~, X] = ode45(@(t,X) populationdynamics(X,x,sd_K,sd_a,n),tspan,X0);
    X = X(end,:);
    X = X(:);
    X(X<0.001) = 0;
end

%% fitness

function f = fitness(x,y,X,sd_K,sd_a,n)
    % given a resident phenotype matrix, mutant phenotype matrix, and a
    % resident carrying capacity, returns the fitness vector of the mutant 
    % species in an environment created by community x
    phenotype = [x y];
    beta = zeros(n,n);
    % define the bottom left corner of the competition matrix we need
    for j = n+1:2*n
        for k = 1:n
            beta(j-n,k) = competition(phenotype,j,k,sd_a);
        end
    end
    % compute the fitness
    f = growthrate(y).*(1-beta*X./carryingcapacityM(y,sd_K,n));
end

%% parameter configuration

% configuration space dimensions
n0 = 1; % initial number of species
d = 2; % dimension of phenotype space 

% standard deviations
sd_a = 0.5; % competition sd
sd_K = 1; % carrying capacity sd

% species initial conditions
x0 = 0.2*randn(d,n0); % phenotype matrix. each column represents a species' phenotype
X0 = 0.1.*ones(n0,1); % initial population vector

% numerical parameters
stepsize = 0.01; % step size
steps = 200; % steps
tspan = [0,1000]; % ecological integration timespan

%% adaptive dynamics

% initialise
x = x0;
X = X0;
n = n0;
trajectory = zeros(d,n,steps);
trajectory(:,:,1) = x0;

% trait substitution algorithm
for i = 2:steps
    i
    % calculate resident ecological equilibrium
    X = effectivecarryingcapacity(x,X,sd_K,sd_a,n,tspan);
    % add mutations
    y = x+stepsize.*(randn(d,n));
    f = fitness(x,y,X,sd_K,sd_a,n);
    mutant = find(f>0 & X~=0);
    x(:,mutant) = y(:,mutant);
    X(mutant) = 0.1;
    trajectory(:,:,i) = x;
end

%% plots

for k = 1:n
    nonzerofirst = find(trajectory(1,k,:),1,'first');
    nonzerolast = find(trajectory(1,k,:),1,'last');
    xifirst = trajectory(:,k,nonzerofirst);
    xilast = trajectory(:,k,nonzerolast);
    for j = 1:nonzerofirst
        trajectory(:,k,j) = xifirst;
    end
    for j = nonzerolast:steps
        trajectory(:,k,j) = xilast;
    end
end

figure

i = 1:steps;

% 1d phenotype space
% hold on
% for j = 1:n
%     y = squeeze(trajectory(:,j,i)); 
%     plot(i,y) 
% end
% hold off

axis equal
zero = zeros(d,1);

% 2d phenotype space
hold on
% color = '#008080'
plot(zero(1),zero(2),'o') 
%
for j = 1:n0
    plot(x0(1,j),x0(2,j))
end
for j = 1:n
    y = squeeze(trajectory(:,j,i)); 
    plot(y(1,:),y(2,:),color = '#008080',linewidth = 4) 
end
scatter(trajectory(1,:,1),trajectory(2,:,1))
scatter(trajectory(1,:,end),trajectory(2,:,end),'x')

% 3d phenotype space
% plot3(zero(1),zero(2),zero(3),'o') 
% plot3(trajectory(1,:,1),trajectory(2,:,1),trajectory(3,:,1),'o')
% plot3(trajectory(1,:,end),trajectory(2,:,end),trajectory(3,:,end),'x')
% hold on
% for j = 1:n0
%     plot3(x0(1,j),x0(2,j),x0(3,j),'o')
% end
% for j = 1:n
%     y = squeeze(trajectory(:,j,i));
%     plot3(y(1,:),y(2,:),y(3,:));
% end    
% hold off