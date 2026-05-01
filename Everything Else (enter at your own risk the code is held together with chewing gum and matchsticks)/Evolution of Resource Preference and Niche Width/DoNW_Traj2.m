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

function r = GrowthRate(x,y,c,R_0,w,s,T)
n = length(x);
r = zeros(n,1);
for i = 1:n
    r(i) = R_0*w*s./(sqrt(s^2+y(i).^2)).*exp(-(c.*y(i))-x(i).^2./(2.*(s.^2+y(i).^2)))-T;
end
end

function K = CarryingCapacity(x,y,c,R_0,w,s,T,r_R)
n = length(x);
K = zeros(n,1);
for i = 1:n
    K(i) = r_R.*y(i).*sqrt(2.*pi*(2.*s^2+y(i).^2)).*exp(c.*y(i)).*exp(x(i).^2./(2*s^2+y(i).^2)).*(exp(-x(i).^2./(2*(s^2+y(i).^2)))./(sqrt(s^2+y(i).^2))-(T.*exp(c.*y(i))/(R_0.*w.*s)));
end
end

function alpha = Competition(x,y,c,s) % competition x,y experiences from u,v
n = length(x);
alpha = zeros(n,n);
for i =1:n
    for j = 1:n
        alpha(i,j) = y(i).*sqrt((2.*s.^2+y(i).^2)./(y(i).^2.*y(j).^2+s^2.*(y(i).^2+y(j).^2))).*exp(c.*(y(i)-y(j))+x(i).^2./(2.*s^2+y(i).^2)-(x(j).^2.*y(i).^2+x(i).^2.*y(j).^2+s^2.*(x(j)-x(i)).^2)./(2.*(y(j).^2.*y(i).^2+s^2.*(y(j).^2+y(i).^2))));
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
c = 0.2;
R_0 = 1.5;
w = 1;
s = 0.5;
T = 0.05;
r_R = 1;
divlim = 1; % maximal diversity

% initial conditions
x0 = [0.4;0.4]; % 
X0 = 0.1.*ones(n0,1); % initial population vector

% numerical parameters
stepsize = 0.01; % variance of the random normal mutation step
steps = 400; % evolutionary timespan
tspan = [0,1000]; % ecological integration timespan
branchtime = 5; % steps between each branch event

% plotting parameters
res = 0.01;
edgeop = 0;
faceop = 0.5;
width = 1;
l = 3;
frames = 5;

% ------------------------------------------------------------------------%
%% Adaptive Dynamics
% ------------------------------------------------------------------------%
% population dynamics
function dXdt = populationdynamics(X,r,K,alpha)
% given an instantaneous population vector and r, K, and alpha 
% precalculated, returns the dynamics vector 
    dXdt = r.*X.*(1-alpha*X./K);
end

function X = EffectiveCarryingCapacity(x,y,c,R_0,w,s,T,r_R,X0,tspan) % (allowing 
% extinction) given a phenotype matrix and an initial population vector, returns the 
% effective carrying capacity vector 
    r = GrowthRate(x,y,c,R_0,w,s,T);
    K = CarryingCapacity(x,y,c,R_0,w,s,T,r_R);
    alpha = Competition(x,y,c,s);
    [~, X] = ode45(@(t,X) populationdynamics(X,r,K,alpha),tspan,X0);
    X = X(end,:);
    X = X(:);
    X(X<0.001) = 0;
end

% fitness
function f = Fitness(u,v,x,y,X,c,R_0,w,s,T,r_R)
n = length(x);
ux = [x,u];
vy = [y,v];
X = X(:);
r = GrowthRate(u,v,c,R_0,w,s,T);
K = CarryingCapacity(u,v,c,R_0,w,s,T,r_R);
alpha = Competition(ux,vy,c,s);
beta = alpha(n+1:2*n,1:n);
f = r.*(1-beta*X./K);
end

% initialise

nmax = min(n0 + floor(steps/branchtime),divlim);
trajectory = NaN(d,nmax,steps);
PopulationData = zeros(nmax,steps);
nlist = NaN(steps);
x = NaN(d,nmax);
X = zeros(1,nmax);
branchcount = 0;

% set initial conditions
x(:,1:n0) = x0;
nbound = n0; % number of species in the system data
X(1:n0) = EffectiveCarryingCapacity(x0(1,:),x0(2,:),c,R_0,w,s,T,r_R,X0,tspan);

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
    u = xeff + stepsize*randn(d,n); % mutation step

    fxeff = xeff(1,:); % turn them into the right form
    fyeff = xeff(2,:);
    fueff = u(1,:);
    fveff = u(2,:);

    f = Fitness(fueff,fveff,fxeff,fyeff,Xeff,c,R_0,w,s,T,r_R);
    invaders = find(f>0); % find which species outcompete
    xeff(:,invaders) = u(:,invaders);
    % reincorporate to the full system
    x(:,survivors) = xeff;
    X(survivors) = Xeff;
    % check for ecological extinction
    X(survivors) = EffectiveCarryingCapacity(fxeff,fyeff,c,R_0,w,s,T,r_R,Xeff,tspan);
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
% trait substitution sequence
i = 1:steps;
figure
hold on

for j = 1:n
    y = trajectory(:,j,i);
    plot(y(1,:),y(2,:),linewidth = 2) 
end
scatter(0,0.6093859637938,'filled','red')
hold off
xlim([-width,width])
ylim([0,l])
xlabel('Niche Position');
ylabel('Niche Width');

function r = GrowthRate1(x,y,c,R_0,w,s,T)
    r = R_0*w*s./(sqrt(s^2+y.^2)).*exp(-(c.*y)-x.^2./(2.*(s.^2+y.^2)))-T;
end
function K = CarryingCapacity1(x,y,c,R_0,w,s,T,r_R)
    K = r_R.*y.*sqrt(2.*pi*(2.*s^2+y.^2)).*exp(c.*y).*exp(x.^2./(2*s^2+y.^2)).*(exp(-x.^2./(2*(s^2+y.^2)))./(sqrt(s^2+y.^2))-(T.*exp(c.*y)/(R_0.*w.*s)));
end
function alpha = Competition1(x,y,u,v,c,s)
    alpha = y.*sqrt((2.*s.^2+y.^2)./(y.^2.*v.^2+s^2.*(y.^2+v.^2))).*exp(c.*(y-v)+x.^2./(2.*s^2+y.^2)-(u.^2.*y.^2+x.^2.*v.^2+s^2.*(u-x).^2)./(2.*(v.^2.*y.^2+s^2.*(v.^2+y.^2))));
end
function f = Fitness1(x,y,u,v,c,R_0,w,s,T,r_R)
f = GrowthRate1(u,v,c,R_0,w,s,T).*(1-Competition1(u,v,x,y,c,s).*CarryingCapacity1(x,y,c,R_0,w,s,T,r_R)./CarryingCapacity1(u,v,c,R_0,w,s,T,r_R));
end


[x,y] = meshgrid(-width:res:width,0:res:l);
frame = 1;
    % if i == 0
    %     frame = 1;
    % end
    xmin = x(1,1);
    xmax = x(end,end);
    ymin = y(1,1);
    ymax = y(end,end);
    % u = 0;
    % v = sqrt((s^5*R_0/(s*T))^(2/5)-s^2);
    u = trajectory(1,1,frame);
    v = trajectory(2,1,frame);

    figure
    f = Fitness1(u,v,x,y,c,R_0,w,s,T,r_R);
    r = GrowthRate1(x,y,c,R_0,w,s,T);
    surf(x,y,f,FaceAlpha = faceop, EdgeColor = 'interp',EdgeAlpha = edgeop);
    hold on
    scatter3(u,v,Fitness1(u,v,u,v,c,R_0,w,s,T,r_R),'filled','LineWidth',20,MarkerFaceColor='red',MarkerEdgeColor='red')
    contour3(x,y,f,[0,0],'white',LineWidth = 3)
    % contour3(x,y,r,[0,0],'white',LineWidth = 3)
    hold off
    % view(100,20);
    view(90,90) % heatmap
    xlim([xmin,xmax])
    ylim([ymin,ymax])
    xlabel('Niche Position');
    ylabel('Niche Width');
    zlabel('Fitness');
    title("Time ="+frame)
    colorbar;
for i = 1:frames
    frame = i*floor(steps/frames);
    % if i == 0
    %     frame = 1;
    % end
    xmin = x(1,1);
    xmax = x(end,end);
    ymin = y(1,1);
    ymax = y(end,end);
    % u = 0;
    % v = sqrt((s^5*R_0/(s*T))^(2/5)-s^2);
    u = trajectory(1,1,frame);
    v = trajectory(2,1,frame);

    figure
    f = Fitness1(u,v,x,y,c,R_0,w,s,T,r_R);
    r = GrowthRate1(x,y,c,R_0,w,s,T);
    surf(x,y,f,FaceAlpha = faceop, EdgeColor = 'interp',EdgeAlpha = edgeop);
    hold on
    scatter3(u,v,Fitness1(u,v,u,v,c,R_0,w,s,T,r_R),'filled','LineWidth',20,MarkerFaceColor='red',MarkerEdgeColor='red')
    contour3(x,y,f,[0,0],'white',LineWidth = 3)
    % contour3(x,y,r,[0,0],'white',LineWidth = 3)
    hold off
    % view(100,20);
    view(90,90) % heatmap
    xlim([xmin,xmax])
    ylim([ymin,ymax])
    xlabel('Niche Position');
    ylabel('Niche Width');
    zlabel('Fitness');
    title("Time ="+frame)
    colorbar;
end

disp(trajectory(:,:,end))
% figure
% [x,y] = meshgrid(-width:res:width,0:res:length);
% xmin = x(1,1);
% xmax = x(end,end);
% ymin = y(1,1);
% ymax = y(end,end);
% % u = 0;
% % v = sqrt((s^5*R_0/(s*T))^(2/5)-s^2);
% u = trajectory(1,:,end);
% v = trajectory(2,:,end);
% f = Fitness1(u,v,x,y,c,R_0,w,s,T,r_R);
% surf(x,y,f,FaceAlpha = faceop, EdgeColor = 'interp',EdgeAlpha = edgeop);
% hold on
% scatter3(u,v,Fitness1(u,v,u,v,c,R_0,w,s,T,r_R),'filled','LineWidth',5,MarkerFaceColor='red',MarkerEdgeColor='red')
% contour3(x,y,f,[0,0],'white',LineWidth = 3)
% hold off
% view(100,20);
% view(90,90) % heatmap
% xlim([xmin,xmax])
% ylim([ymin,ymax])
% xlabel('Niche Position');
% ylabel('Niche Width');
% zlabel('Fitness');
% colorbar;

%%
M(steps) = struct('cdata',[],'colormap',[]);
fig = figure;

for frame = 1:steps
    clf
    u = trajectory(1,1,frame);
    v = trajectory(2,1,frame);
    f = Fitness1(u,v,x,y,c,R_0,w,s,T,r_R);
    surf(x,y,f,FaceAlpha = faceop, EdgeColor = 'interp',EdgeAlpha = edgeop);
    hold on
    scatter3(u,v,Fitness1(u,v,u,v,c,R_0,w,s,T,r_R),'filled','LineWidth',5,MarkerFaceColor='red',MarkerEdgeColor='red')
    contour3(x,y,f,[0,0],'white',LineWidth = 3)
    hold off
    % view(100,20);
    view(90,90) % heatmap
    xlim([xmin,xmax])
    ylim([ymin,ymax])
    xlabel('Niche Position');
    ylabel('Niche Width');
    zlabel('Fitness');
    title("Time ="+frame)
    M(frame) = getframe(fig);
end
hold off

tss = VideoWriter('HeatmapMovie_c-02', 'MPEG-4');
% tss.FrameRate=60;

open(tss)
writeVideo(tss,M)
close(tss)
