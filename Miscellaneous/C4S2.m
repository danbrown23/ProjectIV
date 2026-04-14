%% PLOT STYLE

set(0,'defaultTextFontName', 'Arial')
set(0,'defaultaxesfontsize', 20); % 27 for 1X3, 20 for 1X2
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','latex');
set(0,'defaultAxesXGrid','off');
set(0,'defaultAxesYGrid','off');
set(0,'defaultAxesTickDir','out');
set(0,'defaultAxesLineWidth',1.5);

%% phenotype parametrisations

function alpha = competition(x,i,j,sd_a)
% given a phenotype matrix, returns the competition species i feels from species j
    alpha = exp((-norm(x(:,i)-x(:,j)).^2)./(2.*sd_a^2)); % symmetric competition
end

function K = carryingcapacity(x,i,sd_K)
% given a phenotype matrix, returns the carrying capacity of species i
    K = exp(-norm(x(:,i)).^2./(2.*sd_K^2));
end

function alpha = competitionM(x,sd_a,n)
% returns the competition matrix
alpha = zeros(n,n);
    for i = 1:n
        for j = 1:n
            alpha(i,j) = competition(x,i,j,sd_a);  
        end
    end
end

function K = carryingcapacityM(x,sd_K,n)
% returns the carrying capacity vector
    K = zeros(n,1);
    for i = 1:n
        K(i) = carryingcapacity(x,i,sd_K);
    end
end

function r = growthrate(x)
% given a phenotype matrix, returns the growth rate vector of species x
    % r=exp(norm(x));
    r = 1;
end

%% population dynamics

function dXdt = populationdynamics(X,x,sd_K,sd_a,n)
% given a population vector and phenotype matrix, returns the dynamics vector 
    dXdt = growthrate(x).*X.*(1-competitionM(x,sd_a,n)*X./carryingcapacityM(x,sd_K,n));
end

% (assume no species go extinct)
function X = effectivecarryingcapacityM(x,sd_K,sd_a,n)
% given a phenotype matrix, returns the effective carrying capacity vector 
    X = inv(competitionM(x,sd_a,n))*carryingcapacityM(x,sd_K,n);
end

%% fitness

function f = fitness(x,y,sd_K,sd_a,n)
    % given a phenotype matrix and a mutant matrix, returns the fitness
    % vector of the mutant species in an environment created by community x
    phenotype = [x y];
    beta = zeros(n,n);
    % define the bottom left corner of the competition matrix we need
    for j = n+1:2*n
        for k = 1:n
            beta(j-n,k) = competition(phenotype,j,k,sd_a);
        end
    end
    % compute the fitness
    f = growthrate(y).*(1-beta*effectivecarryingcapacityM(x,sd_K,sd_a,n)./carryingcapacityM(y,sd_K,n));
end

%% generate n species with random phenotypes

% configuration space dimensions
n = 2; % number of species
d = 2; % dimension of phenotype space 

% standard deviations
sd_a = 0.5; % competition sd
sd_K = 1; % carrying capacity sd

% species initial conditions
x0 = 2*randn(d,n); % phenotype matrix. each column represents a species' phenotype
% x0 = [0 0; -0.2 0.2];
% X0 = [0.1, 0.1]; % population vector
X0 = 0.1*ones(n,1);

% carrying capacity profile test
% d=1:
% figure
% x = -5:0.01:5;
% y = carryingcapacityM(x,sd_K,length(x));
% plot(x,y)
% d=2:
% figure
% x = meshgrid(-5:0.1:5)
% y = carryingcapacity

figure % note: this breaks if d+n > 10 for some reason
tspan = [0,100];
[t, X] = ode45(@(t,X) populationdynamics(X,x0,sd_K,sd_a,n),tspan,X0);

for i = 1:n
    plot(t,X(:,i),'-o')
    hold on
end
hold off


%% adaptive dynamics

stepsize = 0.02; % step size
steps = 10000; % steps

% trait substitution sequence
x = x0;
trajectory = zeros(d,n,steps);
for i = 1:steps
    trajectory(:,:,i) = x;
    y = x+stepsize*randn(d,n);
    f = fitness(x,y,sd_K,sd_a,n);
    for j = 1:n
        if f(j) > 0
            x(:,j) = y(:,j);
        end
    end
end

figure

i = 1:steps;
axis equal
zero = zeros(d,1);

% 1d phenotype space
% for j = 1:n
%     y = squeeze(trajectory(:,j,i)); 
%     plot(i,y) 
% end

% 2d phenotype space
hold on
for j = 1:n
    plot(x0(1,j),x0(2,j),'o')
end
for j = 1:n
    y = squeeze(trajectory(:,j,i)); 
    plot(y(1,:),y(2,:),linewidth=2) 
end
hold off

figure
clf
plot3(0,0,1:steps)
hold on
for j = 1:n
    y = squeeze(trajectory(:,j,i)); 
    plot3(y(1,:),y(2,:),1:steps,linewidth=2) 
end
hold off

% y1 = squeeze(trajectory(:,1,i));
% y2 = squeeze(trajectory(:,2,i));

% figure
% for j=1:steps
%     clf % clears plot
%     y1_j = y1(:,j);
%     y2_j = y2(:,j);
%     hold on % keeps all art
%     %% PLOT STYLE
% 
%     set(0,'defaultTextFontName', 'Arial')
%     set(0,'defaultaxesfontsize', 20); % 27 for 1X3, 20 for 1X2
%     set(0,'defaultAxesTickLabelInterpreter','none');
%     set(0,'defaulttextinterpreter','latex');
%     set(0,'defaultAxesXGrid','off');
%     set(0,'defaultAxesYGrid','off');
%     set(0,'defaultAxesTickDir','out');
%     set(0,'defaultAxesLineWidth',1.5);
% 
%     plot(zero(1),zero(2),'o') %origin
%     plot(x0(1,1),x0(2,1),'o',MarkerFaceColor='#008080') %starting pos of species
%     plot(x0(1,2),x0(2,2),'o',MarkerFaceColor='#ffc107') 
%     plot(y1(1,1:j),y1(2,1:j),color='#008080',linewidth = 2) % trajectory
%     plot(y2(1,1:j),y2(2,1:j),color='#ffc107',linewidth = 2)
%     plot(y1_j(1),y1_j(2),'go',color='#003336',linewidth=2, markersize = 15) % marks end
%     plot(y2_j(1),y2_j(2),'go',color='#ff6500',linewidth=2, markersize = 15) 
% 
%     axis equal
%     xlim([-0.6,0.6])
%     ylim([-0.6,0.6])
%     xlabel('x')
%     ylabel('y')
%     movie(j) = getframe;
% end
% hold off
% 
% tss = VideoWriter('traitsubsequence');
% tss.FrameRate=20;
% 
% open(tss)
% writeVideo(tss,movie)
% close(tss)

% 3d phenotype space
% plot3(zero(1),zero(2),zero(3),'o') 
% hold on
% for j = 1:n
%     plot3(x0(1,j),x0(2,j),x0(3,j),'o')
% end
% for j = 1:n
%     y = squeeze(trajectory(:,j,i));
%     plot3(y(1,:),y(2,:),y(3,:));
% end    
% hold off

