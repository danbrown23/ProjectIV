%% Plots for the evolution of niche position chapter in my dissertation

clear

set(0,'defaultTextFontName', 'Arial')
set(0,'defaultaxesfontsize', 20); % 27 for 1X3, 20 for 1X2
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','latex');
set(0,'defaultAxesXGrid','off');
set(0,'defaultAxesYGrid','off');
set(0,'defaultAxesTickDir','out');
set(0,'defaultAxesLineWidth',1.5);

%% Parameters

c = 0;
R_0 = 1;
w = 1;
s = 1;
T = 0.05;
r_R = 1;
y = 1;

%% Functions

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
%% Plots

res = 0.01;
edgeop = 0.2;
faceop = 0.7;

figure
x = -4:res:4;
xmin = x(1,1);
xmax = x(end,end);
r = GrowthRate(x,y,c,R_0,w,s,T);
plot(x,r);
xlim([xmin,xmax])
xlabel('Niche Position');
ylabel('Growth Rate');

figure
x = -4:res:4;
xmin = x(1,1);
xmax = x(end,end);
K = CarryingCapacity(x,y,c,R_0,w,s,T,r_R);
Kmax = max(K);
plot(x,K);
xlim([xmin,xmax])
ylim([0,Kmax+2])
xlabel('Niche Position');
ylabel('Consumer Carrying Capacity');

figure
x = -4:res:4;
xmin = x(1,1);
xmax = x(end,end);
hold on
for u = 0:2
v = 1;
alpha = Competition(x,y,u,v,c,s);
plot(x,alpha);
scatter(u,Competition(u,v,u,v,c,s),'filled','LineWidth',5,MarkerFaceColor='red',MarkerEdgeColor='red')
end
hold off

colormap(parula(5))
xlim([xmin,xmax])
xlabel('Niche Position');
ylabel('Competition');

figure
x = -2:res:2;
xmin = x(1,1);
xmax = x(end,end);
hold on
for u = 0:2
v = 1;
f = Fitness(x,y,u,v,c,R_0,w,s,T,r_R);
plot(x,f)
scatter(u,0,'filled','LineWidth',5,MarkerFaceColor='red',MarkerEdgeColor='red')
end
hold off
xlim([xmin,xmax])
xlabel('Niche Position');
ylabel('Fitness');