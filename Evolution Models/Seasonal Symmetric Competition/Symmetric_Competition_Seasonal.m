clear 
import 
%% parameter configuration

% configuration space dimensions
n0 = 1; % initial number of species
d = 2; % dimension of phenotype space 

% standard deviations
sd_a = 0.5; % competition sd
sd_K = 1; % carrying capacity sd

% species initial conditions
x0 = [0.4;0.4]; % phenotype matrix. each column represents a species' phenotype
X0 = 0.1.*ones(n0,1); % initial population vector

%numerical parameters
stepsize = 0.01; % step size
steps = 1200; % steps
tspan = [0,1000]; % ecological integration timespan
branchcd = 100; % steps between each branch event
mergedist = 0.5; % the threshhold average distance between two trajectories over 
mergetime = 10; % steps before they are merged for computational efficiency

%% adaptive dynamics

% initialise
nmax = 2^floor(steps/branchcd); % highest possible value of n
trajectory = NaN(d,nmax,steps);
x = NaN(d,nmax);
X = zeros(1,nmax);

% set initial conditions
x(:,1:n0) = x0;
nbound = n0; % upper bound for n
X(1:n0) = effectivecarryingcapacity(x0,X0,sd_K,sd_a,n0,tspan,1); 
survivors = find(X~=0);
xeff = x(:,survivors); 
Xeff = X(survivors);
trajectory(:,survivors,1) = x(:,survivors);
n = length(survivors); % number of species currently in the system
count = 0;



% trait substitution algorithm
for i = 2:steps
    disp(["step =",i]); 
    disp(["species =",n]); 
    % check for extinction/branching
    X(survivors) = effectivecarryingcapacity(xeff,Xeff,sd_K,sd_a,n,tspan,i); % ecological extinction
    survivors = find(X~=0);
    n = length(survivors);
    % if i > mergetime & nbound > 1 % merge branches
    %     for j = nbound/2:nbound
    %         if norm(trajectory(:,j,i-mergetime) - trajectory(:,2*j,i-mergetime)) < mergedist
    %             X(j) = 0; % kills off a branched species if it stays too close to its cousin
    %             disp("merged!");
    %         end
    %     end
    % end
    if i - count > branchcd-1 % branch event
        disp("branched!");
        x(:,nbound+survivors) = x(:,survivors) + stepsize.*(rand(d,n)); % duplicate each species
        X(nbound+survivors) = 0.01.*ones(1,n);
        X([survivors nbound+survivors]) = effectivecarryingcapacity(x(:,[survivors nbound+survivors]),...
                                                             X([survivors nbound+survivors]),sd_K,sd_a,...
                                                             length([survivors nbound+survivors]),tspan,i); % calculate ecological equilibrium
        nbound = 2*nbound;
        count = i;
    end
    survivors = find(X~=0);
    n = length(survivors);
    % mutate  
    % define the effective state of the system
    xeff = x(:,survivors); 
    Xeff = X(survivors); 
    % perform a trait substitution for the effective system
    y = xeff + stepsize*randn(d,n); % mutation step
    f = fitness(xeff,y,Xeff,sd_K,sd_a,n,i);
    mutant = find(f>0);
    xeff(:,mutant) = y(:,mutant);
    Xeff(mutant) = 0.01;
    % reincorporate to the full system
    x(:,survivors) = xeff;
    X(survivors) = Xeff;
    % update trajectory
    trajectory(:,survivors,i) = x(:,survivors); 
end

%% plots

figure(1)

axis equal
zero = zeros(d,1);

% 2d phenotype space
hold on
axis equal
plot(zero(1),zero(2),'o') 
for j = 1:n0
    plot(x0(1,j),x0(2,j),'o')
end
for j = 1:n
    y = squeeze(trajectory(:,j,:)); 
    plot(y(1,:),y(2,:), linewidth = 2) 
end

% 3d phenotype space
% i = 1:steps;
% plot3(zero(1),zero(2),zero(3),'o') 
% hold on
% plot3(trajectory(1,:,1),trajectory(2,:,1),trajectory(3,:,1),'o') % initial position
% plot3(trajectory(1,:,end),trajectory(2,:,end),trajectory(3,:,end),'x')
% % for j = 1:n0
% %     plot3(x0(1,j),x0(2,j),x0(3,j),'o')
% % end
% for j = 1:n
%     y = squeeze(trajectory(:,j,i));
%     plot3(y(1,end),y(2,end),y(3,end),'x')
%     plot3(y(1,:),y(2,:),y(3,:),linewidth=4);
% end    
% hold off

%% animations


% Initialise:

% M(steps) = struct('cdata',[],'colormap',[]);
% colors = viridis(n);
% f = figure;
% 
% for i = 1:steps
%     clf
%     plot3(0,0,0,'o',color='black')
%     hold on
%     plot3(x0(1),x0(2),x0(3),'o')
%     for j = 1:n
%         x = squeeze(trajectory(1,j,1:i));
%         y = squeeze(trajectory(2,j,1:i));
%         z = squeeze(trajectory(3,j,1:i));
%         plot3(x,y,z,linewidth = 2,color=colors(j,:))
%         plot3(x(i),y(i),z(i),'o',MarkerFaceColor=colors(j,:),MarkerEdgeColor=colors(j,:),linewidth=3)% plot trajectory up to i
%     end
%     axis equal
%     xlim([-2,2])
%     ylim([-2,2])
%     zlim([-2,2])
%     view(0.25*(45+i),20)
%     M(i) = getframe(f, [300 30 740 580]);
% end
% hold off
% 
% tss = VideoWriter('traitsubsequenceseasonal1', 'MPEG-4');
% tss.FrameRate=60;
% 
% open(tss)
% writeVideo(tss,M)
% close(tss)