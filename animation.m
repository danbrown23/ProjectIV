function r = r_X(x,c,w,R_0,sigma_R,T) % logistic growth rate of X
r = c.*((w.*R_0.*sigma_R)./(sqrt(1+sigma_R.^2)).*exp(-1/2.*x.^2./(1+sigma_R^2))-T);
end

function K = K_X(x,c,w,R_0,sigma_R,T,r_R) % logistic carrying capacity of X
I_2 = w.*R_0.*sigma_R./(sqrt(2.*pi.*(2.*sigma_R.^2+1).*r_R).*exp(-x.^2./(1+2.*sigma_R.^2)));
K = (r_X(x,c,w,R_0,sigma_R,T)./c)./I_2;
end

function dxdt = Xdot(X,x,c,w,R_0,sigma_R,T,r_R) % time evolution of X
dxdt = r_X(x,c,w,R_0,sigma_R,T).*X.*(1-X./K_X(x,c,w,R_0,sigma_R,T,r_R));
end

function K = cc(z,R_0,sigma_R)
K = R_0.*exp(-z.^2./(2.*sigma_R^2));
end

function a = inter(z,x)
a = 1./(sqrt(2*pi)).*exp(-(z-x).^2./2);
end

function R = resource(z,R_0,sigma_R,x,r,X)
R = cc(z,R_0,sigma_R)./r.*(r-inter(z,x).*X);
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

tspan = [0,30];
[t, X] = ode45(@(t,X) Xdot(X,x,c,w,R_0,sigma_R,T,r_R), tspan, X_0);

obj = VideoWriter('My Animation');
open(obj);
z = -5:0.1:5;
axis([-5 5 0 0.6])
for i = 1:length(X)
    plot(z,resource(z,R_0,sigma_R,x,r_R,X(i)))
    M(i) = getframe;
    pause(0.1)
end
writeVideo(obj, M)
close(obj);