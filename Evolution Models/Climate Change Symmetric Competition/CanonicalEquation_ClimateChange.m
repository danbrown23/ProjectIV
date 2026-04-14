%% PLOT STYLE

set(0,'defaultTextFontName', 'Arial')
set(0,'defaultaxesfontsize', 20); % 27 for 1X3, 20 for 1X2
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','latex');
set(0,'defaultAxesXGrid','off');
set(0,'defaultAxesYGrid','off');
set(0,'defaultAxesTickDir','out');
set(0,'defaultAxesLineWidth',1.5);

import Functions.rainbow_cmap.*
%% phenotype parametrisations

function r = growthrate(x) 
% given a phenotype matrix, returns the growth rate vector of species x
    r = 1; % no
end

function K = carryingcapacity(x,i,sd_K,t)
% given a phenotype matrix, returns the carrying capacity of species i
    K = exp(-norm(x(:,i)-0.01*[1;0]*t).^2./(2.*sd_K^2));
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

function K = carryingcapacityM(x,sd_K,n,t)
% returns the carrying capacity vector
    K = zeros(n,1);
    for i = 1:n
        K(i) = carryingcapacity(x,i,sd_K,t);
    end
end

%% population dynamics

function X = effectivecarryingcapacity(x,sd_K,sd_a,n,t) % (assume no species go extinct)
% given a phenotype matrix, returns the effective carrying capacity vector 
    X = inv(competitionM(x,sd_a,n))*carryingcapacityM(x,sd_K,n,t);
end

%% adaptive dynamics

function f = fitness(x,y,X,sd_K,sd_a,n,t)
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
    f = growthrate(y).*(1-beta*X./carryingcapacityM(y,sd_K,n,t));
end

function s = selectiongradient(x,sd_K,sd_a,n,t)
    % given a resident phenotype matrix, finds the selection gradient
    % matrix
    d = length(x)/n;
    s = zeros(d,n);
    x = reshape(x,[d n]); % changes x from the vector ode45 uses to the matrix we use for intuition
    alpha = competitionM(x,sd_a,n);
    K = carryingcapacityM(x,sd_K,n,t);
    X = effectivecarryingcapacity(x,sd_K,sd_a,n,t);
    term2summand = zeros(d,n);
    for i = 1:n
        for j = 1:n
            term2summand(:,j) = alpha(i,j).*X(j).*x(:,j);
        end
        term2 = sum(term2summand(:,:),2);
        s(:,i) = (1/(sd_a)^2-1/(sd_K)^2).*x(:,i)-1/(sd_a^2.*K(i)).*term2;
    end   
    s = reshape(s,[n*d,1]);
end

%% parameter configuration

% configuration space dimensions
n = 10; % number of species
d = 2; % dimension of phenotype space 

% standard deviations
sd_a = 0.5; % competition sd
sd_K = 1; % carrying capacity sd

% species initial conditions
% x0 = [-0.1661    0.1459    0.1247   -0.0106 ;
%       0.2090    0.0457    0.1331   -0.2294]; % phenotype matrix. each column represents a species' phenotype
x0 = 2*randn(d,n);
X0 = 0.1.*ones(n,1); % initial population vector

% numerical parameters
stepsize = 0.01; % step size
steps = 20000; % steps
ecotspan = [0,1000]; % ecological integration timespan
adtspan = [0,steps/100]; % adaptive dynamics timespan

%% trait substitution sequence

% initialise
x = x0;
X = X0;
trajectory = zeros(d,n,steps);
trajectory(:,:,1) = x0;

% trait substitution algorithm
for i = 2:steps
    % calculate resident ecological equilibrium
    X = effectivecarryingcapacity(x,sd_K,sd_a,n,i);
    % add mutations
    y = x+stepsize.*(randn(d,n));
    f = fitness(x,y,X,sd_K,sd_a,n,i);
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

%% 2d phenotype space

% hold on
% plot(zero(1),zero(2),'o') 
% plot(trajectory(1,:,1),trajectory(2,:,1),'o')
% plot(trajectory(1,:,end),trajectory(2,:,end),'x')
% 
% % trait substitution sequence
% for j = 1:n
%     y = squeeze(trajectory(:,j,i)); 
%     plot(y(1,:),y(2,:),color = '#FFC107',linewidth = 2) 
% end
% % canonical equation of adaptive dynamics
% [t, x] = ode45(@(t,x) selectiongradient(x,sd_K,sd_a,n),adtspan,x0);
% for j = 1:n
%     y = x(:,(j-1)*d+1:j*d);
%     plot(y(:,1),y(:,2),color = '#008080',linewidth = 4)
% end
% xlabel('$x_1$')
% ylabel('$x_2$')
% lgd = legend("Stochastic","Deterministic");
% lgd.Location = 'southeast';
% lgd.FontSize = 30;
hold off

%% 3d phenotype space
% clf
% hold off
% plot3(0,0,0,'o') 
% hold on
% plot3(trajectory(1,:,1),trajectory(2,:,1),trajectory(3,:,1),'o')
% plot3(trajectory(1,:,end),trajectory(2,:,end),trajectory(3,:,end),'x')
% % trait substitution sequence
% for j = 1:n
%     plot3(x0(1,j),x0(2,j),x0(3,j),'o')
%     y = squeeze(trajectory(:,j,i));
%     plot3(y(1,:),y(2,:),y(3,:),color = '#FFC107', linewidth = 4);
% end  
% % canonical equation of adaptive dynamics
% [t, x] = ode45(@(t,x) selectiongradient(x,sd_K,sd_a,n),adtspan,x0);
% for j = 1:n
%     y = x(:,(j-1)*d+1:j*d);
%     plot3(y(:,1),y(:,2),y(:,3),color = '#008080', linewidth=4)
% end
% hold off

hold on
sims = 1;
% theta = 2*pi/sims*(1:sims);
colour = prism(sims);
scatter(0,0,'black',"filled")
for i = 1:sims
    x0 = 2*rand(d,n)-ones(2,n);
    % x0 = 0.2*[cos(theta(i));sin(theta(i));cos(theta(i)+pi/9);sin(theta(i)-pi)];
    [t, x] = ode45(@(t,x) selectiongradient(x,sd_K,sd_a,n,t),adtspan,x0);
    for j = 1:n
        y = x(:,(j-1)*d+1:j*d);
        plot(y(end,1),y(end,2),'o',linewidth = 4,color=colour(i,:))
        plot(y(:,1),y(:,2),linewidth = 2, color = colour(i,:))
    end
end
% xlim([-1.2,1.2])
% ylim([-1.2,1.2])
xlabel('$x_1$')
ylabel('$x_2$')

hold off
