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

c = 0;
R_0 = 1.5;
w = 1;
s = 0.5;
T = 0.05;
r_R = 1;


%% Functions
function R = ResourceCarryingCapacity(z,x,X,R_0,r_R,s)
R = R_0/r_R*exp(-z.^2./(2*s^2)).*(r_R-X/sqrt(2*pi)*exp(-(z-x).^2/2));
end

function r = GrowthRate(x,y,c,R_0,w,s,T)
    r = R_0*w*s./(sqrt(s^2+y.^2)).*exp(-(c.*y)-x.^2./(2.*(s.^2+y.^2)))-T;
end

function K = CarryingCapacity(x,y,c,R_0,w,s,T,r_R)
    K = r_R.*y.*sqrt(2.*pi*(2.*s^2+y.^2)).*exp(c.*y).*exp(x.^2./(2*s^2+y.^2)).*(exp(-x.^2./(2*(s^2+y.^2)))./(sqrt(s^2+y.^2))-(T.*exp(c.*y)/(R_0.*w.*s)));
end

function alpha = Competition(x,y,u,v,c,s)
    alpha = y.*sqrt((2.*s.^2+y.^2)./(y.^2.*v.^2+s^2.*(y.^2+v.^2))).*exp(c.*(y-v)+x.^2./(2.*s^2+y.^2)-(u.^2.*y.^2+x.^2.*v.^2+s^2.*(u-x).^2)./(2.*(v.^2.*y.^2+s^2.*(v.^2+y.^2))));
end

function f = Fitness(u,v,x,y,c,R_0,w,s,T,r_R)
f = GrowthRate(u,v,c,R_0,w,s,T).*(1-Competition(u,v,x,y,c,s).*CarryingCapacity(x,y,c,R_0,w,s,T,r_R)./CarryingCapacity(u,v,c,R_0,w,s,T,r_R));
end

function f = Fitness1(u,v,x,y,c,R_0,w,s,T)
f = GrowthRate(u,v,c,R_0,w,s,T)-GrowthRate(x,y,c,R_0,w,s,T).*Competition(x,y,u,v,c,s);
end
%% Plots

res = 0.1;
edgeop = 0.2;
faceop = 0.7;

figure
[x,y] = meshgrid(-4:res:4,0:res:4);
xmin = x(1,1);
xmax = x(end,end);
ymin = y(1,1);
ymax = y(end,end);
r = GrowthRate(x,y,c,R_0,w,s,T);
r(r<-1)=-1;
surf(x,y,r,FaceAlpha = faceop, EdgeColor = 'interp',EdgeAlpha = edgeop);
hold on
contour3(x,y,r,[0,0],'white',LineWidth = 3)
hold off
view(125,20)
% view(90,90) % heatmap
xlim([xmin,xmax])
ylim([ymin,ymax])
xlabel('Niche Position');
ylabel('Niche Width');
zlabel('Growth Rate');
colorbar;

figure
[x,y] = meshgrid(-10:res:10,0:res:16);
xmin = x(1,1);
xmax = x(end,end);
ymin = y(1,1);
ymax = y(end,end);
K = CarryingCapacity(x,y,c,R_0,w,s,T,r_R);
K(K<-1) = -1;
surf(x,y,K,FaceAlpha = faceop, EdgeColor = 'interp',EdgeAlpha = edgeop);
hold on
contour3(x,y,K,[0,0],'white',LineWidth = 3)
hold off
view(80,20);
% view(90,90) % heatmap
xlim([xmin,xmax])
ylim([ymin,ymax])
zlim([-1,15])
xlabel('Niche Position');
ylabel('Niche Width');
zlabel('Carrying Capacity');
colorbar;

figure
[x,y] = meshgrid(-2:res:2,0:res:4);
xmin = x(1,1);
xmax = x(end,end);
ymin = y(1,1);
ymax = y(end,end);
u = 0;
v = sqrt((s^5*R_0/(s*T))^(2/5)-s^2);
alpha = Competition(x,y,u,v,c,s);
surf(x,y,alpha,FaceAlpha = faceop, EdgeColor = 'interp',EdgeAlpha = edgeop);
hold on
% plot3(u,v,Competition(u,v,u,v,c,s),'LineWidth',5)
scatter3(u,v,Competition(u,v,u,v,c,s),'filled','LineWidth',5,MarkerFaceColor='red',MarkerEdgeColor='red')
hold off
view(35,20);
% view(90,90) % heatmap
xlim([xmin,xmax])
ylim([ymin,ymax])
xlabel('Niche Position');
ylabel('Niche Width');
zlabel('Competition');
colorbar;

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
f = Fitness(x,y,u,v,c,R_0,w,s,T,r_R);
surf(x,y,f,FaceAlpha = faceop, EdgeColor = 'interp',EdgeAlpha = edgeop);
hold on
% plot3(u,v,Competition(u,v,u,v,c,s),'LineWidth',5)
scatter3(u,v,Fitness(u,v,u,v,c,R_0,w,s,T,r_R),'filled','LineWidth',5,MarkerFaceColor='red',MarkerEdgeColor='red')
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


