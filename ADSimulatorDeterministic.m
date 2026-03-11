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
% given an instantaneous population vector and phenotype matrix, returns the dynamics vector 
    dXdt = growthrate(x).*X.*(1-competitionM(x,sd_a,n)*X./carryingcapacityM(x,sd_K,n));
end

% function X = effectivecarryingcapacityM(x,sd_K,sd_a,n) % (assume no species go extinct)
% % given a phenotype matrix, returns the effective carrying capacity vector 
%     X = inv(competitionM(x,sd_a,n))*carryingcapacityM(x,sd_K,n);
% end

function X = effectivecarryingcapacity(x,X0,sd_K,sd_a,n) % (allowing extinction)
% given a phenotype matrix and an initial population vector, returns the effective carrying capacity vector 
    tspan = [0,1000];
    [t, X] = ode45(@(t,X) populationdynamics(X,x,sd_K,sd_a,n),tspan,X0);
    X = X(end,:);
    X = X(:);
    for i = 1:n
        if X(i) <0.001
            X(i) = 0;
        end
    end   
end

%% fitness

function f = fitness(x,y,X0,sd_K,sd_a,n)
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
    f = growthrate(y).*(1-beta*effectivecarryingcapacity(x,X0,sd_K,sd_a,n)./carryingcapacityM(y,sd_K,n));
end

%% generate n species with random phenotypes

% configuration space dimensions
n = 1; % number of species
d = 3; % dimension of phenotype space 

% standard deviations
sd_a = 0.5; % competition sd
sd_K = 1; % carrying capacity sd

% species initial conditions
x0 = 0.2*randn(d,n); % phenotype matrix. each column represents a species' phenotype
X0 = 0.1.*ones(n,1); % initial population vector

%% adaptive dynamics

stepsize = 0.01; % step size
steps = 400; % steps

% trait substitution sequence

x = x0;
X = X0;
m = n;
trajectory = NaN(d,n,steps);
count = 0;
for i = 1:steps
    i
    if m ~= size(trajectory,2)
        for j=m/2+1:m
            trajectory(:,j,1:i) = NaN(d,1,i);
        end
    end
    trajectory(:,:,i) = NaN(d,m,1);


    % check for extinction
    X = effectivecarryingcapacity(x,X,sd_K,sd_a,m);
    % survivors = find(X);
    % X = X(survivors);
    % x = x(:,survivors);
    % trajectory(:,survivors,i) = x;
    % check for branching
    for j = 1:m
        trajectory(:,j,i) = x(:,j);
    end
    xsize = zeros(1,20);
    if i-count > 100
        % for j = 1:50
        %     for k = 1:m
        %         xsize(k,j) = norm(trajectory(:,k,i-j));
        %     end
        % end
        % for k = 1:m
        %     for l =1:m
        %         if abs(mean(xsize(l,:))-mean(xsize(k,:))) < 0.05
        %             if x(1,m) ~= 0
        %                 x = [x,x(:,l)+stepsize*randn(d,1)];
        %                 X = [X;0.1];
        %             end
        %         end
        %     end
        %     if max(xsize(l,:))-min(xsize(l,:)) < 0.05
        %         if x(1,m) ~= 0
        %             x(:,m+1) = x(:,l)+stepsize*randn(d,1);
        %             X(m+1) = 0.1;
        %         end
        %     end
        % end
        % m = length(X);
        for j=m+1:2*m
            x(:,j) = x(:,j-m)+stepsize*randn(d,1);
            X(j) = 0.1;
        end
        m = 2*m;
        count = i;
    end
    % perform a mutation
    y = x+stepsize*randn(d,m);
    f = fitness(x,y,X,sd_K,sd_a,m);
    mut = find(f>0);
    x(:,mut) = y(:,mut);
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
% for k = 1:m
%     nonzero1 = find(trajectory(1,k,:),1,'first');
%     nonzerolast = find(trajectory(1,k,:),1,'last');
%     xi0 = trajectory(:,k,nonzero1);
%     xilast = trajectory(:,k,nonzerolast);
%     for j = 1:nonzero1
%         trajectory(:,k,j) = xi0;
%     end
%     for j = nonzerolast:steps
%         trajectory(:,k,j) = xilast;
%     end
% end
% 2d phenotype space
% hold on
% % color = '#008080'
% plot(zero(1),zero(2),'o') 
% for j = 1:n
%     plot(x0(1,j),x0(2,j),'o')
% end
% for j = 1:m
%     y = squeeze(trajectory(:,j,i)); 
%     plot(y(1,:),y(2,:), LineWidth = 4) 
% end
% scatter(trajectory(1,:,1),trajectory(2,:,1))
% scatter(trajectory(1,:,end),trajectory(2,:,end),'x')
% 3d phenotype space

trajectory
figure 
plot3(zero(1),zero(2),zero(3),'o') 
plot3(trajectory(1,:,1),trajectory(2,:,1),trajectory(3,:,1),'o')
plot3(trajectory(1,:,end),trajectory(2,:,end),trajectory(3,:,end),'x')
hold on
for j = 1:n
    plot3(x0(1,j),x0(2,j),x0(3,j),'o')
end
for j = 1:m
    y = squeeze(trajectory(:,j,i));
    plot3(y(1,:),y(2,:),y(3,:),linewidth = 2);
end    
hold off

figure

for j=1:steps
    clf
    plot3(zero(1),zero(2),zero(3),'o') 
    hold on
    % n = len(trajectory(:,:,i))
    for k = 1:n
        yk = squeeze(trajectory(:,k,i)); % all dims, kth species, ith timestep
        %% PLOT STYLE

        set(0,'defaultTextFontName', 'Arial')
        set(0,'defaultaxesfontsize', 27); % 27 for 1X3, 20 for 1X2
        set(0,'defaultAxesTickLabelInterpreter','none');
        set(0,'defaulttextinterpreter','latex');
        set(0,'defaultAxesXGrid','off');
        set(0,'defaultAxesYGrid','off');
        set(0,'defaultAxesTickDir','out');
        set(0,'defaultAxesLineWidth',1.5);

        plot3(squeeze(trajectory(1,k,1:j)),squeeze(trajectory(2,k,1:j)),squeeze(trajectory(3,k,1:j)),linewidth = 2)
        plot3(yk(1,k),yk(2,k),yk(3,k),'go',linewidth=2, markersize = 15) 
    end
    axis equal
    xlim([-0.5,0.5])
    ylim([-0.5,0.5])
    zlim([-0.5,0.5])
    xlabel('x')
    ylabel('y')    
    movie(j) = getframe;
end
hold off

tss = VideoWriter('traitsubsequence');
tss.FrameRate=20;

open(tss)
writeVideo(tss,movie)
close(tss)

