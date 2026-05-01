clear
%% Summary

% Implements the Gillespie SSA algorithm for the logistic competition model

%% The Model

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
% remove interspecies competition
alpha = alpha - eye(size(alpha));
end

%% Parameter Configuration

% configuration space dimension
dim = 2;

% initial conditions
N = 200;
X0 = ones(N,1);
x0 = 0.01*randn(dim,N) + [0.4;0.4];

% ecological parameters
c = -0.1;
R_0 = 150;
w = 1;
s = 0.5;
T = 5;
r_R = 1;

% numerical parameters
var_m = 0.02;
sd_m = sqrt(var_m);
% get an estimate of the maximal population by integrating 
xpts = [-5:0.1:5,0:0.1:10];
ypts = CarryingCapacity(xpts(:,1),xpts(:,2),c,R_0,w,s,T,r_R);
NCC = trapz(ypts);
% Nmax = round(NCC + 50);
Nmax = 1000;
stepsmax = 1000;
threshhold = 5;

%% Gillespie Algorithm

tic

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
    r = GrowthRate(xeff(1,:),xeff(2,:),c,R_0,w,s,T);
    % d = r./CarryingCapacity(xeff,var_K).*sum(Competition(xeff,var_a),2);
    d = CarryingCapacity(xeff(1,:),xeff(2,:),c,R_0,w,s,T,r_R).*sum(Competition(xeff(1,:),xeff(2,:),c,s),2);
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

toc

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

for f = 0:frames
    figure
    t = tframe*f;
    i=1;
    while t > timedata(i)
        i = i+1;
    end
    x = phenodata(:,:,i);
    histogram2(x(1,:),x(2,:),-5:0.1:5,0:0.1:10)
end

% heatmap (?)
% xdata = reshape(phenodata(1,:,:),[Nmax*step,1]);
% ydata = repmat(timedata,1,Nmax);
% ydata = ydata';
% ydata = reshape(ydata,[Nmax*step,1]);
% % histogram2(xdata,ydata,-3:0.1:3,0:tfinal/step:tfinal)
% [counts, edges] = histcounts2(xdata,ydata,-2:0.1:1,timedata);
% figure
% imagesc(edges(1), edges(2), counts)
% colormap(jet)  % choose a colormap
% colorbar  % add a colorbar to the plot
% axis xy  % set the axis orientation to x-y
% xlabel('x')  % label the x-axis
% ylabel('y')  % label the y-axis