clear
%% Summary

% Implements the Gillespie SSA algorithm for the logistic competition model
% with symmetric competition

%% The Model

function r = GrowthRate(x)
% given a phenotype matrix, returns the growth rate vector of species x
    % r=exp(norm(x));
    N = size(x,2);
    r = ones(N,1);
end

function K = CarryingCapacity(x,var_K,K0)
% Input: x is the (traits x species) phenotype configuration matrix
% Output: K is the carrying capacity vector of the n species
D2 = sum(x.^2,1);
K = K0*exp(-D2./(2.*var_K));
K = K(:);
end

function alpha = Competition(x,var_a)
    % Input: x is the (traits x species) phenotype configuration matrix
    % Output: alpha is the (species x species) competition matrix. alpha_ij
    %  is the competition species i experiences from species j
    D2 = pdist2(x', x', 'euclidean').^2; % pairwise squared distances
    alpha = exp(-D2 ./ (2 * var_a)); % symmetric Gaussian competition
    % remove interspecies competition
    alpha = alpha - eye(size(alpha));
end

%% Parameter Configuration

% configuration space dimension
dim = 1;

% ecological parameters
var_a = 0.5;
var_K = 4;
K0 = 100;

% initial conditions
N = K0;
X0 = ones(N,1);
x0 = 1+0.01*randn(dim,N);

% numerical parameters
var_m = 0.01;
sd_m = sqrt(var_m);
% get an estimate of the maximal population by integrating 
% xpts = -10:0.1:10;
% ypts = CarryingCapacity(xpts,var_K,K0);
% NCC = trapz(ypts);
% Nmax = round(NCC + 50);
Nmax = 400;
stepsmax = 30000;
% extinction
threshhold = 10;

% plotting parameters
frames = 600;
framelength = floor(stepsmax/frames);

%% Gillespie Algorithm

tic

% initialise
timedata = NaN(frames,1);
popudata = zeros(Nmax,frames);
phenodata = NaN(dim,Nmax,frames);
X = popudata(:,frames);
x = phenodata(:,:,frames);
Nlist = NaN(frames,1);

% set initial conditions
survivors = 1:N;
t = 0;
X(survivors) = X0;
x(:,survivors) = x0;
timedata(1) = t;
popudata(:,1) = X;
phenodata(:,:,1) = x;
step = 1;
frame = 0;
frameno = 1;

% algorithm
while step < stepsmax && N < Nmax && N > threshhold
    % step
    % N
    % 1. Generate random numbers r1,r2 uniformly:
    r1 = rand;
    r2 = rand;
    % 2. Compute the propensity function of the system:
    xeff = x(:,survivors);
    Xeff = X(survivors);
    r = GrowthRate(xeff);
    % d = r./CarryingCapacity(xeff,var_K).*sum(Competition(xeff,var_a),2);
    d = r./CarryingCapacity(xeff,var_K,K0).*sum(Competition(xeff,var_a),2);
    R = sum(r);
    D = sum(d);
    propensity = R + D;
    % 3. Compute the next reaction time
    tau = 1/propensity*log(1/r1);
    t = t + tau;
    % 4. Find which of the n^2 reactions took place
    if propensity*r2 < R %birth
        disp("birth")
        % find the reproducing individual
        reproducer = survivors(randi(N));
        % find a population that has gone extinct for memory
        offspring = find(~X,1);
        % put in the baby
        X(offspring) = 1;
        x(:,offspring) = x(:,reproducer) + sd_m*randn(dim,1);
        survivors = [survivors, offspring];
    else % death
        disp("death")
        % find the individual destined to die
        i=1;
        dsum = R + d(i);
        while r2*propensity > dsum
            i = i+1;
            dsum = dsum + d(i);
        end
        % kill him
        deadeff = i;
        dead = survivors(i);
        X(dead) = 0;
        x(:,dead) = NaN(dim,1);
        survivors(i) = [];
    end
    % update system history
    step = step + 1;
    N = length(survivors);
    if step > frame
        timedata(frameno) = t;
        popudata(:,frameno) = X;
        phenodata(:,:,frameno) = x;
        Nlist(frameno) = N;
        frame = frame + framelength;
        frameno = frameno + 1;
    end
end

toc

tic
set(0,'defaultTextFontName', 'Arial')
set(0,'defaultaxesfontsize', 20); % 27 for 1X3, 20 for 1X2
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','latex');
set(0,'defaultAxesXGrid','off');
set(0,'defaultAxesYGrid','off');
set(0,'defaultAxesTickDir','out');
set(0,'defaultAxesLineWidth',1.5);

% tfinal = timedata(end);
% frames = 10;
% tframe = tfinal/frames;
% stepframe = floor(step/frames);

% Species against time

figure
plot(timedata,Nlist)
xlabel("Time")  % label the x-axis
ylabel("Number of Individuals")  % label the y-axis


% Snapshot Profiles
% 
% for f = 0:frames
%     figure
%     t = tframe*f;
%     i=1;
%     while t > timedata(i)
%         i = i+1;
%     end
%     x = phenodata(:,:,i);
%     histogram(x,-10:0.1:10)
%     xlabel("Phenotype")
%     ylabel("Population")
%     title(t)
%     xlim([-10,10])
%     ylim([0,45])
% end

%% Heatmap

figure
% Flatten phenotype data

Ypredata = phenodata(1,:,20:end);
Ydata = reshape(Ypredata(1,:,:), [], 1);

% Repeat time data to match

Xpredata = timedata(20:end);
Xdata = repelem(Xpredata(:), Nmax);

% Plot
h = histogram2(Xdata, Ydata, Xpredata, -10:0.1:10, ...
    'DisplayStyle','tile','ShowEmptyBins','off');

colorbar
xlabel('Time')
ylabel('Phenotype')
ylim([-5,5])
xlim([5,timedata(end)])
toc

%% MOVIE

% Initialise:

% M(frames-2) = struct('cdata',[],'colormap',[]);
% f = figure;
% 
% for i = 1:frames-2
%     clf
%     x = phenodata(:,:,i);
%     t = timedata(i);
%     histogram(x,-5:0.05:5)
%     xlabel("Phenotype")
%     ylabel("Population")
%     title("Time = "+t)
%     xlim([-3,3])
%     ylim([0,20])
%     M(i) = getframe(f);
% end
% 
% tss = VideoWriter('SymGillespie_vara025_varK1_200', 'MPEG-4');
% tss.FrameRate=60;
% 
% open(tss)
% writeVideo(tss,M)
% close(tss)