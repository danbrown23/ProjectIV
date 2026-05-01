%% Algorithm for simulating the evolution of niche width and niche position
clear

set(0,'defaultTextFontName', 'Arial')
set(0,'defaultaxesfontsize', 27); % 27 for 1X3, 20 for 1X2
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','latex');
set(0,'defaultAxesXGrid','off');
set(0,'defaultAxesYGrid','off');
set(0,'defaultAxesTickDir','out');
set(0,'defaultAxesLineWidth',1.5);

% ------------------------------------------------------------------------%
%% Phenotype Parametrisations
% ------------------------------------------------------------------------%

function beta = velocity(b,t,tstart)
if t > tstart
    beta = b*(t-tstart);
else
    beta = 0;
end
end

function r = GrowthRate(x,y,c,R_0,w,s,T,b,t,tstart)
r = R_0*w*s./(sqrt(s^2+y.^2)).*exp(-(c.*y)-(x-velocity(b,t,tstart)).^2./(2.*(s.^2+y.^2)))-T;
r = r(:);
end

function K = CarryingCapacity(x,y,c,R_0,w,s,T,r_R,b,t,tstart)
K = r_R.*y.*sqrt(2.*pi*(2.*s^2+y.^2)).*exp(c.*y).*exp((x-velocity(b,t,tstart)).^2./(2*s^2+y.^2)).*(exp(-(x-velocity(b,t,tstart)).^2./(2*(s^2+y.^2)))./(sqrt(s^2+y.^2))-(T.*exp(c.*y)/(R_0.*w.*s)));
K = K(:);
end

function alpha = Competition(x,y,c,s,b,t,tstart) % competition i experiences from j
n = length(x);
alpha = zeros(n,n);
for i =1:n
    for j = 1:n
        alpha(i,j) = y(i).*sqrt((2.*s.^2+y(i).^2)./(y(i).^2.*y(j).^2+s^2.*(y(i).^2+y(j).^2))).*exp(c.*(y(i)-y(j))+(x(i)-velocity(b,t,tstart)).^2./(2.*s^2+y(i).^2)-((x(j)-velocity(b,t,tstart)).^2.*y(i).^2+(x(i)-velocity(b,t,tstart)).^2.*y(j).^2+s^2.*(x(j)-x(i)).^2)./(2.*(y(j).^2.*y(i).^2+s^2.*(y(j).^2+y(i).^2))));
    end
end
end

% ------------------------------------------------------------------------%
%% Parameter Configuration
% ------------------------------------------------------------------------%
% configuration space dimensions
n0 = 1; % initial number of species
d = 2; % dimension of phenotype space

% system parameters
c = 1;
R_0 = 1.5;
w = 1;
s = 0.5;
T = 0.05;
r_R = 5;
% Alist = [1,3,5];
% omega = 0.0015;
A = 3;
blist = 0:0.0005:0.01;
tstart = 3000;
divlim = 25; % maximal diversity

% initial conditions
x0 = [0.4;0.4]; % (trait x species) phenotype matrix
X0 = 0.1.*ones(n0,1); % initial population vector

% numerical parameters
stepsize = 0.01; % variance of the random normal mutation step
steps = 6000; % evolutionary timespan
tspan = [0,1000]; % ecological integration timespan
branchtime = 5; % steps between each branch event

% plotting parameters
width = 3;
l = 3;
frames = 10;
fstart = tstart;
fend = 3000;

% ------------------------------------------------------------------------%
%% Adaptive Dynamics
% ------------------------------------------------------------------------%
% population dynamics
function dXdt = populationdynamics(X,r,K,alpha)
% given an instantaneous population vector and r, K, and alpha
% precalculated, returns the dynamics vector
dXdt = r.*X.*(1-alpha*X./K);
end

function X = EffectiveCarryingCapacity(x,y,c,R_0,w,s,T,r_R,X0,tspan,b,t,tstart) % (allowing
% extinction) given a phenotype matrix and an initial population vector, returns the
% effective carrying capacity vector
r = GrowthRate(x,y,c,R_0,w,s,T,b,t,tstart);
K = CarryingCapacity(x,y,c,R_0,w,s,T,r_R,b,t,tstart);
alpha = Competition(x,y,c,s,b,t,tstart);
[~, X] = ode45(@(t,X) populationdynamics(X,r,K,alpha),tspan,X0);
X = X(end,:);
X = X(:);
X(X<0.05) = 0;
end

% fitness
function f = Fitness(u,v,x,y,X,c,R_0,w,s,T,r_R,b,t,tstart)
n = length(x);
ux = [x,u];
vy = [y,v];
X = X(:);
r = GrowthRate(u,v,c,R_0,w,s,T,b,t,tstart);
K = CarryingCapacity(u,v,c,R_0,w,s,T,r_R,b,t,tstart);
alpha = Competition(ux,vy,c,s,b,t,tstart);
beta = alpha(n+1:2*n,1:n);
f = r.*(1-beta*X./K);
end

popdata = zeros(steps,length(blist));
ndata = zeros(steps,length(blist));
mdata = NaN(2,steps,length(blist));
vardata = NaN(2,steps,length(blist));
tdata = zeros(steps,length(blist));
for bindex = 1:length(blist)

    b = blist(bindex);

    % initialise

    nmax = min(n0 + floor(steps/branchtime),divlim);
    trajectory = NaN(d,nmax,steps);
    PopulationData = zeros(nmax,steps);
    nlist = NaN(steps,1);
    mlist = NaN(2,steps);
    varlist = NaN(2,steps);
    CCClist = zeros(steps,1);
    x = NaN(d,nmax);
    X = zeros(1,nmax);
    branchcount = 0;

    % set initial conditions
    x(:,1:n0) = x0;
    nbound = n0; % number of species in the system data
    X(1:n0) = EffectiveCarryingCapacity(x0(1,:),x0(2,:),c,R_0,w,s,T,r_R,X0,tspan,b,1,tstart);

    survivors = find(X~=0);
    xeff = x(:,survivors); % surviving species currently in the system
    Xeff = X(survivors);
    trajectory(:,survivors,1) = xeff;
    PopulationData(survivors,1) = Xeff;
    nlist(1) = n0;
    Xmean = mean(Xeff);
    npweightxeff = 1/Xmean.*xeff(1,:).*Xeff;
    nwweightxeff = 1/Xmean.*xeff(2,:).*Xeff;
    mlist(1,1) = mean(npweightxeff);
    mlist(2,1) = mean(nwweightxeff);
    varlist(1,1) = var(nwweightxeff);
    varlist(2,1) = var(nwweightxeff);
    n = length(survivors); % number of species currently in the system
    i=1;

    tic
    % trait substitution algorithm
    while i < steps && n > 0
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
            x(:,branch_species1) = x(:,branch_species1) + 0.1*stepsize*randn(d,1);
            X(branch_species1) = 0.5*X(branch_species1);
            survivors = [survivors,branch_species2]; % add the branched species
        end
        xeff = x(:,survivors);
        Xeff = X(survivors);
        n = length(survivors);
        % mutation
        u = xeff + stepsize*randn(d,n); % mutation step

        fxeff = xeff(1,:); % turn them into the right form
        fyeff = xeff(2,:);
        fueff = u(1,:);
        fveff = u(2,:);

        f = Fitness(fueff,fveff,fxeff,fyeff,Xeff,c,R_0,w,s,T,r_R,b,i,tstart);
        invaders = find(f>0); % find which species outcompete
        xeff(:,invaders) = u(:,invaders);
        % reincorporate to the full system
        x(:,survivors) = xeff;
        X(survivors) = Xeff;
        % check for ecological extinction
        X(survivors) = EffectiveCarryingCapacity(fxeff,fyeff,c,R_0,w,s,T,r_R,Xeff,tspan,b,i,tstart);
        survivors = find(X~=0);
        % check that niche width is bigger than 0
        infeasible = find(x(2,survivors)<0);
        X(survivors(infeasible)) = 0;
        survivors(infeasible) = []; 
        % update data
        n = length(survivors); % current number of species
        xeff = x(:,survivors);
        Xeff = X(survivors);
        trajectory(:,survivors,i) = xeff;
        PopulationData(survivors,i) = Xeff;
        nlist(i) = n;
        Xmean = mean(Xeff);
        npweightxeff = 1/Xmean.*xeff(1,:).*Xeff;
        nwweightxeff = 1/Xmean.*xeff(2,:).*Xeff;
        mlist(1,i) = mean(npweightxeff);
        mlist(2,i) = mean(nwweightxeff);
        varlist(1,i) = var(npweightxeff);
        varlist(2,i) = var(nwweightxeff);
        CCClist(i) = velocity(b,i,tstart);
    end
    toc

    % total population against time
    hold on
    popsums = sum(PopulationData,1);
    popdata(:,bindex) = popsums;
    ndata(:,bindex) = nlist;
    mdata(:,:,bindex) = mlist;
    vardata(:,:,bindex) = varlist;
    tdata(:,bindex) = CCClist;
    disp("yay")
    % for i = 1:divlim
    %     y = trajectory(:,i,:);
    %     plot3(y(1,:),y(2,:),1:steps)
    %     hold on
    %     plot3(CCClist,zeros(steps,1),1:steps)
    % end
    % view(0,5)
    % xlabel("Niche Position")
    % ylabel("Niche Width")
    % zlabel("Time Steps")
    % figure
    % for i = 1:divlim
    %     y = trajectory(:,i,:);
    %     plot3(y(1,:),y(2,:),1:steps)
    %     hold on
    %     plot3(CCClist,zeros(steps,1),1:steps)
    % end
    % view(-45,5)
    % xlabel("Niche Position")
    % ylabel("Niche Width")
    % zlabel("Time Steps")
    % hold off
    % figure
    % hold on
    % for i = 1:divlim
    %     y = trajectory(:,i,:);
    %     plot(1:steps,y(1,:))
    %     plot(1:steps,CCClist,color = 'blue')
    % end
    % hold off
    % figure
    % hold on
    % for i = 1:divlim
    %     y = trajectory(:,i,:);
    %     plot(1:steps,y(2,:))
    %     plot(1:steps,CCClist,color = 'blue')
    % end
    % hold off
    % xlabel("Time Steps")
    % ylabel("Phenotype")
end
%%
% figure
% for i = 1:length(blist)   
%     % total population
%     hold on
%     plot(1:steps,popdata(:,i));
%     xlabel("Time Steps")
%     ylabel("Total Population")
% end
% hold off
% 
% figure
% for i = 1:length(blist)
%     % number of species
%     hold on
%     plot(1:steps,ndata(:,i))
%     xlabel("Time Steps")
%     ylabel("Number of Distinct Species")
% end
% hold off
% for i = 1:length(blist)
%     figure
%     hold on
%     plot(1:steps,mdata(1,:,i))
%     plot(1:steps,mdata(2,:,i))
%     plot(1:steps,tdata(:,i),color = 'blue')
%     hold off
% end
% for i = 1:length(blist)
%     figure
%     hold on
%     plot(1:steps,vardata(1,:,i))
%     plot(1:steps,vardata(2,:,i))
%     plot(1:steps,tdata(:,i),color = 'blue')
%     hold off
% end

%%
figure
for i = 1:length(blist)
    hold on
    plot(1:steps,mdata(2,:,i))
end
hold off
xlabel("Time Steps")
ylabel("Mean Niche Width")
hold off

m = squeeze(mean(mdata(2,4000:6000,:),2));

figure
scatter(blist,m,80,'filled')
xlabel("CCC Speed")
ylabel("Mean Niche Width")

figure

nspec = mean(ndata(4000:6000,:),1);

scatter(blist,nspec,80,'filled')
xlabel("CCC Speed")
ylabel("Mean Number of Distinct Species")