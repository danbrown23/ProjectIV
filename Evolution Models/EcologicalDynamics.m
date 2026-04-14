clear

%% PLOT STYLE

set(0,'defaultTextFontName', 'Arial')
set(0,'defaultaxesfontsize', 30); % 27 for 1X3, 20 for 1X2
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','latex');
set(0,'defaultAxesXGrid','off');
set(0,'defaultAxesYGrid','off');
set(0,'defaultAxesTickDir','out');
set(0,'defaultAxesLineWidth',1.5);

%% Population Dynamics

function alpha = competition(x,i,j,sd_a)
% given a phenotype matrix, returns the competition species i feels from species j
    alpha = exp((-norm(x(i)-x(j)).^2)./(2.*sd_a^2)); % symmetric competition
    % alpha = eq(i,j); % no competition
end

function K = carryingcapacity(x,i,sd_K)
% given a phenotype matrix, returns the carrying capacity of species i
    K = exp(-norm(x(i)).^2./(2.*sd_K^2));
end

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

function r = growthrate(x)
% given a phenotype matrix, returns the growth rate vector of species x
    % r=exp(norm(x));
    r = 1;
end

function dXdt = populationdynamics(X,r,K,alpha)
% given an instantaneous population vector and phenotype matrix, returns the dynamics vector 
    dXdt = r.*X.*(1-alpha*X./K);
    % dXdt = growthrate(x).*X.*(1-competitionM(x,sd_a,n)*X./carryingcapacityM(x,sd_K,n));
end

%% parameter configuration

% configuration space dimensions
n = 500; % number of species

% standard deviations
sd_a = 0.5; % competition sd
sd_K = 1; % carrying capacity sd

% species initial conditions
x = rand(1,n)-0.5*ones(1,n); % phenotype matrix. each column represents a species' phenotype
X0 = 0.1*ones(1,n); % initial population vector
% x = [1 0.99];
% X0 = [0.6 0.01];

%numerical parameters
tspan = [0,10000]; % ecological integration timespan

figure(1)
% disp("vectorisation")
tic
r = growthrate(x);
K = carryingcapacityM(x,sd_K,n);
alpha = competitionM(x,sd_a,n);
[t, X] = ode45(@(t,X) populationdynamics(X,r,K,alpha),tspan,X0);
toc

% [t, X] = ode45(@(t,X) populationdynamics(X,x,sd_K,sd_a,n),tspan,X0);
clf
hold on
xlabel("Time")
ylabel("Population Density")

% 2 species
% plot(t,X(:,1), linewidth = 4,color='#008080')
% plot(t,X(:,2), linewidth = 4,color='#FFC107')
% 
% lgd = legend("x_1","x_2");
% lgd.Location = 'East';
% lgd.FontSize = 30;

% 5 species
plot(t,X, linewidth = 4)
% lgd = legend("x_1","x_2","x_3","x_4","x_5");
% lgd.Location = 'East';
% lgd.FontSize = 30;
hold off

figure(2)
clf
Xfinal = X(end,:);
feas = find(Xfinal >= 0 & Xfinal < 1 & -0.2 < x & x < 0.2);
x = x(feas);
Xfinal = Xfinal(feas);
scatter(x,Xfinal,'filled')