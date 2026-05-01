%% Simulates the trajectory for 1 species under niche width evolution, not modelling extinction

%% parameter configuration

% configuration space dimensions
n = 1; % number of species
d = 2; % dimension of phenotype space 

% system parameters
c = 0.5;
R_0 = 1.5;
w = 1;
s = 0.5;
T = 0.05;
r_R = 1;

% species initial conditions
x0 = [0.4;0.4];
X0 = 0.1.*ones(n,1); % initial population vector

% numerical parameters
stepsize = 0.01; % step size
steps = 1000; % steps
tspan = [0,1000]; % ecological integration timespan
adtspan = [0,steps/100]; % adaptive dynamics timespan

%% Functions

function r = GrowthRate(x,y,c,R_0,w,s,T)
n = length(x);
r = zeros(n,1);
for i = 1:n
    r(i) = R_0*w*s./(sqrt(s^2+y(i).^2)).*exp(-(c.*y(i))-x(i).^2./(2.*(s.^2+y(i).^2)))-T;
end
end

function K = CarryingCapacity(x,y,c,R_0,w,s,T,r_R)
n = length(x);
K = zeros(1,n);
for i = 1:n
    K = r_R.*y(i).*sqrt(2.*pi*(2.*s^2+y(i).^2)).*exp(c.*y(i)).*exp(x(i).^2./(2*s^2+y(i).^2)).*(exp(-x(i).^2./(2*(s^2+y(i).^2)))./(sqrt(s^2+y(i).^2))-(T.*exp(c.*y(i))/(R_0.*w.*s)));
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
    [t, X] = ode45(@(t,X) populationdynamics(X,r,K,alpha),tspan,X0);
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

%% Algorithm


% initialise
x = x0;
X = X0;
trajectory = zeros(d,n,steps);
trajectory(:,:,1) = x0;
PopulationData = zeros(n,steps);
PopulationData(1:n,1) = X0;

% check feasible
r = GrowthRate(x0(1),x0(2),c,R_0,w,s,T);
if r < 0
    disp("UNFEASIBLE!")
end

% trait substitution algorithm
tic
for i = 2:steps
    % calculate resident ecological equilibrium
    X = EffectiveCarryingCapacity(x(1,:),x(2,:),c,R_0,w,s,T,r_R,X,tspan);
    if X > 0
        PopulationData(1:n,i) = X;
        % add mutations
        mut = x+stepsize.*(randn(d,n));
        f = Fitness(mut(1,:),mut(2,:),x(1,:),x(2,:),X,c,R_0,w,s,T,r_R);
        mutant = find(f>0);
        x(:,mutant) = mut(:,mutant);
        X(mutant) = 0.1;
        trajectory(:,:,i) = x;
    end
end
toc

%% Plots

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
xlim([-3,3])
ylim([0,3])
xlabel('Niche Position');
ylabel('Niche Width');

res = 0.1;
edgeop = 0.2;
faceop = 0.7;

function r = GrowthRate1(x,y,c,R_0,w,s,T)
    r = R_0*w*s./(sqrt(s^2+y.^2)).*exp(-(c.*y)-x.^2./(2.*(s.^2+y.^2)))-T;
end
function K = CarryingCapacity1(x,y,c,R_0,w,s,T,r_R)
    K = r_R.*y.*sqrt(2.*pi*(2.*s^2+y.^2)).*exp(c.*y).*exp(x.^2./(2*s^2+y.^2)).*(exp(-x.^2./(2*(s^2+y.^2)))./(sqrt(s^2+y.^2))-(T.*exp(c.*y)/(R_0.*w.*s)));
end
function alpha = Competition1(x,y,u,v,c,s)
    alpha = y.*sqrt((2.*s.^2+y.^2)./(y.^2.*v.^2+s^2.*(y.^2+v.^2))).*exp(c.*(y-v)+x.^2./(2.*s^2+y.^2)-(u.^2.*y.^2+x.^2.*v.^2+s^2.*(u-x).^2)./(2.*(v.^2.*y.^2+s^2.*(v.^2+y.^2))));
end
function f = Fitness1(u,v,x,y,c,R_0,w,s,T,r_R)
f = GrowthRate1(u,v,c,R_0,w,s,T).*(1-Competition1(u,v,x,y,c,s).*CarryingCapacity1(x,y,c,R_0,w,s,T,r_R)./CarryingCapacity1(u,v,c,R_0,w,s,T,r_R));
end

figure
[x,y] = meshgrid(-2:res:2,0:res:3);
xmin = x(1,1);
xmax = x(end,end);
ymin = y(1,1);
ymax = y(end,end);
% u = 0;
% v = sqrt((s^5*R_0/(s*T))^(2/5)-s^2);
u = x0(1,1);
v = x0(2,1);
f = Fitness1(u,v,x,y,c,R_0,w,s,T,r_R);
surf(x,y,f,FaceAlpha = faceop, EdgeColor = 'interp',EdgeAlpha = edgeop);
hold on
% plot3(u,v,Competition(u,v,u,v,c,s),'LineWidth',5)
scatter3(u,v,Fitness1(u,v,u,v,c,R_0,w,s,T,r_R),'filled','LineWidth',5,MarkerFaceColor='red',MarkerEdgeColor='red')
hold off
view(80,20);
% view(90,90) % heatmap
xlim([xmin,xmax])
ylim([ymin,ymax])
xlabel('Niche Position');
ylabel('Niche Width');
zlabel('Fitness');
colorbar;

figure
[x,y] = meshgrid(-2:res:2,0:res:3);
xmin = x(1,1);
xmax = x(end,end);
ymin = y(1,1);
ymax = y(end,end);
% u = 0;
% v = sqrt((s^5*R_0/(s*T))^(2/5)-s^2);
u = trajectory(1,:,end);
v = trajectory(2,:,end);
f = Fitness1(u,v,x,y,c,R_0,w,s,T,r_R);
surf(x,y,f,FaceAlpha = faceop, EdgeColor = 'interp',EdgeAlpha = edgeop);
hold on
% plot3(u,v,Competition(u,v,u,v,c,s),'LineWidth',5)
scatter3(u,v,Fitness1(u,v,u,v,c,R_0,w,s,T,r_R),'filled','LineWidth',5,MarkerFaceColor='red',MarkerEdgeColor='red')
hold off
view(80,20);
% view(90,90) % heatmap
xlim([xmin,xmax])
ylim([ymin,ymax])
xlabel('Niche Position');
ylabel('Niche Width');
zlabel('Fitness');
colorbar;
