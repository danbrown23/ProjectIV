% The System:

function dxdt = onespeciescrm(x,c,a,w,T,r,K)
X = x(1);
R = x(2);
dxdt = [c.*X.*(a.*w.*R-T); 
    r.*R.*(1-R./K)-a.*X];
end

% R carrying capacity distribution

function K = cc(z,R_0,sigma_R)
K = R_0.*exp(-z.^2./(2.*sigma_R^2));
end

% a Interaction term

function a = inter(z,x)
a = 1./(sqrt(2*pi)).*exp(-(z-x).^2./2);
end

% R when at equilibrium given X

function R = resource(z,R_0,sigma_R,x,r,X)
R = cc(z,R_0,sigma_R)./r.*(r-inter(z,x).*X);
end

% Define Constants:

c = 0.5; % conversion efficiency of resource to biomass
a = 0.1; % probability of a fatal interaction
w = 4; % mass of resource eaten per fatal interaction
T = 0.2; % mass of resource required to sustain population
r = 0.1; % growth rate of resource
K = 1; % carrying capacity of the resource


figure(1);
[X,R] = meshgrid(0:0.1:1);
Xdot = c.*X.*(a.*w.*R-T);
Rdot = r.*R.*(1-R./K)-a.*R.*X;
quiver(R,X,Rdot,Xdot);
hold on;
plot([0,K],[0,0], "ro");
hold on;
plot(T/(a*w),r/a*(1-T/(a*w*K)), "r.")
axis equal;

% Define Constants:

c = 0.5; % conversion efficiency of resource to biomass
w = 1; % mass of resource eaten per fatal interaction
T = 0.2; % mass of resource required to sustain population
r = 1; % growth rate of resource

R_0 = 1; % amplitude of resource
sigma_R = 1; % sd of resource

x = 0.5; % preferred resource
X = 2.5; % current population size

z = -5:0.01:5;
K = cc(z,R_0,sigma_R); % carrying capacity of resource in absence of X

a = inter(z,x);
R = resource(z,R_0,sigma_R,x,r,X);

figure(2)

function r = r_X(x,c,w,R_0,sigma_R,T) % logistic growth rate of X
r = c.*((sqrt(2).*w.*R_0.*sigma_R)./(sqrt(1+sigma_R.^2)).*exp(-1/2.*x.^2./(1+sigma_R^2))-T);
end

function K = K_X(x,w,sigma_R,T,r_R) % logistic carrying capacity of X
% I_2 = w.*R_0.*sigma_R./(sqrt(2.*pi.*(2.*sigma_R.^2+1).*r_R).*exp(-x.^2./(1+2.*sigma_R.^2)));
% K = (r_X(x,c,w,R_0,sigma_R,T)./c)./I_2;
K = r_R.*sqrt(2*pi*(2.*sigma_R.^2+1)).*exp(x.^2./(1+2.*sigma_R.^2)).*(1/sqrt(sigma_R.^2+1).*exp(-1/2.*x.^2./(sigma_R.^2+1))-T/w);
end

function dxdt = consumer(X,x,c,w,R_0,sigma_R,T,r_R) % time evolution of X
dxdt = r_X(x,c,w,R_0,sigma_R,T).*X.*(1-X./K_X(x,w,sigma_R,T,r_R));
end

% Define Constants:

c = 1; % conversion efficiency of resource to biomass
w = 1; % mass of resource eaten per fatal interaction
T = 0.2; % mass of resource required to sustain population
r_R = 1; % growth rate of resource

R_0 = 1; % amplitude of resource
sigma_R = 1; % sd of resource

x = 0.5; % preferred resource
X_0 = 2.5; % initial population size

tspan = [0,40];
[t, X] = ode45(@(t,X) consumer(X,x,c,w,R_0,sigma_R,T,r_R), tspan, X_0);

plot(t,X)

figure(3)

function alpha = competition(x,y,sigma_R,r_R,w) % competition function
alpha = r_R./w.*exp(1/2.*x.^2./(1+2.*sigma_R.^2)).*exp(-1/2.*((1+sigma_R.^2).*x.^2-2.*sigma_R.^2.*x.*y+(1+sigma_R.^2).*y.^2)/(1+2.*sigma_R.^2));
end

y = -5:0.01:5;
hold on
n = 10; % number of curves to plot
for i = -5:10/n:5
    colour = [0 (i+5)/n 1-(i+5)/n];
    alpha = competition(i,y,sigma_R,r_R,w);
    plot(y,alpha,color = colour)
end

figure(4)

[x,y] = meshgrid(-5:0.2:5);
alpha = competition(x,y,sigma_R,r_R,w);
surf(x,y,alpha)

figure(5) % the plots start to look weird here, i assume i've got my fitness function wrong

function f = fitness(x,y,c,w,R_0,sigma_R,T,r_R)
f = r_X(y,c,w,R_0,sigma_R,T) - r_X(x,c,w,R_0,sigma_R,T).*competition(x,y,sigma_R,r_R,w);
end

x = 5;
y=-10:0.01:10;
f = fitness(y,x,c,w,R_0,sigma_R,T,r_R);
plot(y,f)

figure(6)

y = -5:0.01:5;
hold on
n = 10; % number of curves to plot
for i = -5:10/n:5
    colour = [0 (i+5)/n 1-(i+5)/n];
    f = fitness(y,i,c,w,R_0,sigma_R,T,r_R);
    plot(y,f,color = colour)
end
hold off

figure(7)

% parameter regime:
w = 1;
sigma_R = 1;
T = 0.01;
r_R = 1;
x = -5:0.01:5;
K = K_X(x,w,sigma_R,T,r_R);
plot(x,K)
axis([-5 5 0 inf])

figure(8)

[x,y] = meshgrid(-5:0.1:5);
f = fitness(y,x,c,w,R_0,sigma_R,T,r_R);
surf(x,y,f)

figure(9)

contour(x,y,f,[0 0])

figure(10)

function s = selection(x,r_R,sigma_R,w,c,R_0)
    s = -c.*w.*R_0.*sigma_R./sqrt(1+sigma_R.^2).*x./(1+sigma_R.^2).*exp(-1/2.*x.^2./(1+sigma_R.^2))+r_R./w.*x./(1+2.*sigma_R.^2).*exp(-1/2.*x.^2/(1+2.*sigma_R.^2));
end

% Parameter Regime:

c = 1; % conversion efficiency of resource to biomass
w = 1; % mass of resource eaten per fatal interaction
T = 0.2; % mass of resource required to sustain population
r_R = 1; % growth rate of resource
R_0 = 1; % amplitude of resource
sigma_R = 1; % sd of resource

x = -5:0.1:5;
s = selection(x,r_R,sigma_R,w,c,R_0);
plot(x,s,x,zeros(length(x)))