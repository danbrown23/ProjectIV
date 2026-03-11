%% PLOT STYLE
clear

set(0,'defaultTextFontName', 'Arial')
set(0,'defaultaxesfontsize', 20); % 27 for 1X3, 20 for 1X2
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','latex');
set(0,'defaultAxesXGrid','off');
set(0,'defaultAxesYGrid','off');
set(0,'defaultAxesTickDir','out');
set(0,'defaultAxesLineWidth',1.5);

%% Evolution Dynamics

% function s = selectiongradient(x,b,a,d)
% % given a phenotype matrix, outputs the selection gradient matrix
%     s = zeros(d,1);
%     a2 = zeros(4,1);
%     for i = 1:d
%         a1 = zeros(d,d);
%         for j = 1:d
%             a1(i,j) = squeeze(a(i,j,:))'*x(:);
%         end
%         a2(i) = squeeze(a1(i,:))*x(:);
%         s(i) = b(i,:)*x(:) + a2(i) - x(i).^3;
%     end
% end

function s = selectiongradient(x,b,a) % chatgpt faster version

    x = x(:);
    d = length(x);

    a1 = reshape(a, d*d, d) * x;
    a1 = reshape(a1, d, d);

    a2 = a1 * x;

    s = b * x + a2 - x;
end

function m = MVCM(x,d)
% mutation variance covariance matrix
    m = eye(d,d);
end

function dxdt = canon(x,b,a,d)
% canonical equation of adaptive dynamics
    dxdt = MVCM(x,d)*selectiongradient(x,b,a);
end

%% parameter configuration

% configuration space dimensions
d = 2; % dimension of phenotype space 
n = 1; % number of species

% model parameters

% random
% b = randn(d,d);
% a = randn(d,d,d);

% lorenz
sigma = 10;
rho = 28;
beta = 8/3;
b = [1-sigma sigma 0; rho 0 0; 0 0 1-beta];
a = zeros(d,d,d);
a(2,1,3) = -1;
a(3,1,2) = 1;

% chen
% b = [6 0 0; 0 -9 0; 0 0 0.62];
% a = zeros(d,d,d);
% a(1,2,3) = -1;
% a(2,1,3) = 1;
% a(3,1,2) = 1/3;

% species initial conditions
x0 = d^2*randn(d,n); % random
% x0 = [1;2;7]; % lorenz

% numerical parameters
tspan = [0,50]; 


%% plot
figure(1)
clf
for i = 1:n
    [t, x] = ode45(@(t,x) canon(x,b,a,d),tspan,squeeze(x0(:,i)));
    % 2d
    plot(x0(1,i),x0(2,i),'o')
    hold on
    plot(x(:,1),x(:,2))
    % 3d
    % plot3(x0(1,i),x0(2,i),x0(3,i),'o')
    % hold on
    % plot3(x(:,1),x(:,2),x(:,3), color = '#008080',linewidth = 2)
end
hold off

%% animate

% [t, x] = ode45(@(t,x) canon(x,b,a,d),tspan,squeeze(x0(:,1)));
% 
% figh = figure;
% tic
% for j = 1:length(t)
%     clf
%     plot3(0,0,0,'o',color='#ffc107') 
%     hold on
%     % n = len(trajectory(:,:,i))
%     for k = 1:n
%         %% PLOT STYLE
% 
%         set(0,'defaultTextFontName', 'Arial')
%         set(0,'defaultaxesfontsize', 27); % 27 for 1X3, 20 for 1X2
%         set(0,'defaultAxesTickLabelInterpreter','none');
%         set(0,'defaulttextinterpreter','latex');
%         set(0,'defaultAxesXGrid','off');
%         set(0,'defaultAxesYGrid','off');
%         set(0,'defaultAxesTickDir','out');
%         set(0,'defaultAxesLineWidth',1.5);
% 
%         plot3(x(1:j,1),x(1:j,2),x(1:j,3),linewidth=2,color = '#008080')
%         plot3(x(j,1),x(j,2),x(j,3),linewidth=2, markersize = 8,color = '#008080') 
%     end
%     xlim([-25,25])
%     ylim([-25,25])
%     zlim([0,50])
%     view(45+j,20)
%     movie(j) = getframe(figh, [0 0 1000 620]);
% end
% toc
% hold off
% 
% chen = VideoWriter('lorenzsys','MPEG-4');
% chen.FrameRate=20;
% 
% open(chen)
% writeVideo(chen,movie)
% close(chen)