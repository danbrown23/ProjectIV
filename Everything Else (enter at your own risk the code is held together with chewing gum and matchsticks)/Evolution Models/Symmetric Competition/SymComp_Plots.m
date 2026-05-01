%% Plots for the evolution of niche width chapter in my dissertation

clear

set(0,'defaultTextFontName', 'Arial')
set(0,'defaultaxesfontsize', 27); % 27 for 1X3, 20 for 1X2
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','latex');
set(0,'defaultAxesXGrid','off');
set(0,'defaultAxesYGrid','off');
set(0,'defaultAxesTickDir','out');
set(0,'defaultAxesLineWidth',1.5);

%% Parameters

var_a = 0.25;
var_K = 1;


%% The Model

function r = GrowthRate(x)
% given a phenotype matrix, returns the growth rate vector of species x
    % r=exp(norm(x));
    N = size(x,2);
    r = ones(N,1);
end

function K = CarryingCapacity(x,var_K,K0)
% Input: x is the (traits x species) phenotype configuration matrix
% Output: K is the carrying capacity vector of the n species
D2 = sum(x.^2,1);
K = K0*exp(-D2./(2.*var_K));
K = K(:);
end

function alpha = Competition(x,var_a)
    % Input: x is the (traits x species) phenotype configuration matrix
    % Output: alpha is the (species x species) competition matrix. alpha_ij
    %  is the competition species i experiences from species j
    D2 = pdist2(x', x', 'euclidean').^2; % pairwise squared distances
    alpha = exp(-D2 ./ (2 * var_a)); % symmetric Gaussian competition
    % remove interspecies competition
    alpha = alpha - eye(size(alpha));
end

function f = Fitness(x,y,var_K,var_a)
    % given a phenotype matrix and a mutant matrix, returns the fitness
    % vector of the mutant species in an environment created by community x
    % define the bottom left corner of the competition matrix we need
    alpha = Competition(phenotype,var_a);
    beta = alpha(n+1:2*n,1:n);
    f = GrowthRate(y).*(1-beta*CarryingCapacity(x,var_K)./CarryingCapacity(y,var_K));
end

%% Plots

res = 0.1;
edgeop = 0.2;
faceop = 0.7;

figure
[x,y] = meshgrid(-2:res:2,0:res:4);
xmin = x(1,1);
xmax = x(end,end);
ymin = y(1,1);
ymax = y(end,end);
% u = 0;
% v = sqrt((s^5*R_0/(s*T))^(2/5)-s^2);
u = 0;
v = 0.1;
u = [u, v];
f = Fitness(x,u,var_K,var_a);
surf(x,y,f,FaceAlpha = faceop, EdgeColor = 'interp',EdgeAlpha = edgeop);
hold on
% plot3(u,v,Competition(u,v,u,v,c,s),'LineWidth',5)
scatter3(u,v,Fitness(u,u,var_K,var_a),'filled','LineWidth',5,MarkerFaceColor='red',MarkerEdgeColor='red')
hold on
contour3(x,y,f,[0,0],'white',LineWidth = 3)
hold off
hold off
view(145,20);
% view(90,90) % heatmap
xlim([xmin,xmax])
ylim([ymin,ymax])
xlabel('Niche Position');
ylabel('Niche Width');
zlabel('Fitness');
colorbar;

%% Plots 1d

res = 0.01;
edgeop = 0;
faceop = 0.7;

s = 1;
R_0 = 1;
steps = 100;

x = -5:0.01:5;

figure
r = GrowthRate(x,1,0,R_0,w,s,T);
plot(x,r)
xlabel('Niche Position');
ylabel('Growth Rate');
ylim([0,1])

figure
K = CarryingCapacity(x,1,0,R_0,w,s,T,r_R);
plot(x,K)
xlabel('Niche Position');
ylabel('Carrying Capacity');
ylim([0,5])

figure
X = 2;
z = -5:0.01:5;
[zpts, xpts] = meshgrid(z,x);
R = ResourceCarryingCapacity(zpts,xpts,X,R_0,r_R,s);
hold on
surf(z,x,R,FaceAlpha = faceop,EdgeAlpha = edgeop)
xlabel('Resource Phenotype');
ylabel('Consumer Phenotype');
title('Resource Carrying Capacity')
colorbar
view(0,90) % heatmap

figure
x = -2:res:2;
u = -2:res:2;
[xpts, upts] = meshgrid(x,u);
alpha = Competition(upts,1,xpts,1,0,s);
surf(xpts,upts,alpha,FaceAlpha = faceop,EdgeAlpha = edgeop)
xlabel('Resident Niche Position');
ylabel('Mutant Niche Position');
colorbar
view(0,90) % heatmap

figure
f = Fitness(upts,1,xpts,1,0,R_0,w,s,T,r_R);
res = -1.5;
for i = 1:2:steps
    reslist(i) = res;
    mutlist(i) = res;
    reslist(i+1) = res;
    mut = res + sqrt(0.02)*randn;
    if Fitness(mut,1,res,1,0,R_0,w,s,T,r_R) > 0
        res = mut;
    end
    mutlist(i+1) = mut;
end
hold on
surf(xpts,upts,f,FaceAlpha = faceop,EdgeAlpha = edgeop);
plot3(reslist,mutlist,ones(steps,1), '-o',color = 'white' ,linewidth = 3)
plot3(reslist(1),mutlist(1),1,'o',color = 'black', linewidth = 5)
contour3(xpts,upts,f,[0,0],'white',LineWidth = 3)
hold off
xlabel('Resident Niche Position');
ylabel('Mutant Niche Position');
colorbar
view(0,90) % heatmap

figure

% fitness profile at equilibrium
upts = -5:0.01:5;
f = Fitness(upts,1,0,1,0,R_0,w,s,T,r_R);
plot(upts,f)
xlabel('Mutant Niche Position');
ylabel('Fitness');
title("$x=0$")

figure

% fitness profile at equilibrium
upts = -5:0.01:5;
f = Fitness1(upts,1,0,1,0,R_0,w,s,T);
plot(upts,f)
xlabel('Mutant Niche Position');
ylabel('Fitness');
title("$x=0$")


