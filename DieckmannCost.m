function h = h(eta,alpha,beta) % cost function
h = alpha + beta.*eta.^2;
end

function g = g(S,a,eta,alpha,beta) % consumption rate function
g = (a.*S)./(1+a.*h(eta,alpha,beta).*S);
end

function dxdt = odefcn2(y,etares,D,S0,a,h) % uninvaded dynamical system
S = y(1);
X = y(2);
dxdt =  [D.*S0 - D.*S - g(S,a,h,alpha,beta).*X; -D.*X + etares.*g(S,a,h,alpha,beta).*X];
end

function dydt = odefcn3(y,etares,etainv,D,S0,a,h) % invaded dynamical system
S = y(1);
X = y(2);
Y = y(3);
dydt = [D.*S0 - D.*S - g(S,a,h).*(X+Y); -D.*X + etares.*g(S,a,h).*X; -D.*Y + etainv.*g(S,a,h).*Y];
end

D = 1; % death rate
S0 = 1; % constant flux of S
a = 1; % consumption parameter

etares = 0.2;
etainv = 0.3;

% for the nontrivial case of the uninvaded system, we have equilibria:

function Sbar = Sbar(eta,D,a,alpha,beta)
Sbar = D./(a.*(eta - D.*(alpha + beta.*eta.^2)));
end

function Xbar = Xbar(eta,S0,D,a,alpha,beta)
Xbar = eta.*(S0-Sbar(eta,D,a,beta,alpha));
end

% Then we require that Sbar < S0, which gives us conditions on the domain
% of h. We can choose alpha so that the minimal value of eta is 0, and beta
% such that the maximal value is 1

alpha = 0.2; 
beta = 1; 

% Now we can plot h:

figure(1)

eta = 0:0.001:1;
plot(eta,h(eta,alpha,beta))

figure(2)

plot(eta,Sbar(eta,D,a,alpha,beta))