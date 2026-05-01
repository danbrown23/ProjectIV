% system parameters
c = 0.5;
R_0 = 1.5;
w = 1;
s = 0.5;
T = 0.05;
r_R = 1;

x0 = [1;2];

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
u = x0(1);
v = x0(2);
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
u = x0(1);
v = x0(2);
f = Fitness(u,v,x,y,c,R_0,w,s,T,r_R);
surf(x,y,f,FaceAlpha = faceop, EdgeColor = 'interp',EdgeAlpha = edgeop);
hold on
% plot3(u,v,Competition(u,v,u,v,c,s),'LineWidth',5)
scatter3(u,v,Fitness(u,v,u,v,c,R_0,w,s,T,r_R),'filled','LineWidth',5,MarkerFaceColor='red',MarkerEdgeColor='red')
hold off
view(145,20);
% view(90,90) % heatmap
xlim([xmin,xmax])
ylim([ymin,ymax])
xlabel('Niche Position');
ylabel('Niche Width');
zlabel('Fitness');
colorbar;