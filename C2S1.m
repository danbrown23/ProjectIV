% PARAMETER REGIME

w=1;
R0=1;
s=1;
T=0.05;
r_R=1;

% PLOT NEATING

set(0,'defaultTextFontName', 'Arial')
set(0,'defaultaxesfontsize', 20); % 27 for 1X3, 20 for 1X2
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','latex');
set(0,'defaultAxesXGrid','off');
set(0,'defaultAxesYGrid','off');
set(0,'defaultAxesTickDir','out');
set(0,'defaultAxesLineWidth',1.5);

%ylabel, rotation 0 

figure(1) % Growth Rate

function y = r_X(x,w,R0,s,T)
    y = -T+w*R0/sqrt(1+s^(-2))*exp(-x.^2./(2+2*s^2));
end

x = -5:0.01:5;
y = r_X(x,w,R0,s,T);
plot(x,y)

figure(2) % Carrying Capacity

function y = K_X(x,w,R0,s,T,r_R)
    y = sqrt(2*pi)*r_R*sqrt(2+s^(-2))*exp(x.^2./(1+2*s^2)).*(-T/(R0*w)+sqrt(1+s^(-2))^(-1)*exp(-x.^2./(2+2*s^2)));
end

y = K_X(x,w,R0,s,T,r_R);
plot(x,y)
axis([-5 5 0 inf])

figure(3) % Resource Curve

function y = R_X(z,x,s,R0,r_R,w,T)
    y = R0*exp(-1/2.*(z/s).^2).*(1-K_X(x,w,R0,s,T,r_R)./(sqrt(2*pi)*r_R)*exp(-1/2*(z-x).^2));
end

x=2;
z=-5:0.01:5;
y=R_X(z,x,s,R0,r_R,w,T);
plot(z,y)

figure(4) % Competition Kernel

function y = alpha(x,y,s) % competition experienced by species x from species y
    y = exp((2.*x.^2-(1+s^2).*x.^2+2.*s^2.*x.*y-(1+s^2).*y.^2)./(2.*(1+2*s^2)));
end

[x,y] = meshgrid(-5:0.1:5);
z = alpha(x,y,s);
surf(x,y,z)

figure(5) %% Competition Kernel 2D

y = -5:0.01:5;
hold on
n = 5; % number of curves to plot
for i = 0:5/n:5
    colour = [0 (i)/n 1-(i)/n];
    z = alpha(i,y,s);
    plot(y,z,color = colour)
    xline(i,color = colour)
end
hold off

% the competition profiles for species with 10 different phenotypes
