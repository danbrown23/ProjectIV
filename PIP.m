%% Standard PIP plotting

function y = r_X(x,w,R0,s,T) % Growth Rate
    y = -T+w*R0/sqrt(1+s^(-2))*exp(-x.^2./(2+2*s^2));
end

function y = K_X(x,w,R0,s,T,r_R) % Carrying Capacity
    y = sqrt(2*pi)*r_R*sqrt(2+s^(-2))*exp(x.^2./(1+2*s^2)).*(-T/(R0*w)+sqrt(1+s^(-2))^(-1)*exp(-x.^2./(2+2*s^2)));
end

function y = alpha(x,y,s) % competition experienced by species x from species y
    y = exp((2.*x.^2-(1+s^2).*x.^2+2.*s^2.*x.*y-(1+s^2).*y.^2)./(2.*(1+2*s^2)));
end

w=0.5;
R0=10;
s=1;
T=;
r_R=1;

% Define Phi(r, m) numerically as an anonymous function
Phi = @(r, m) r_X(m,w,R0,s,T) - r_X(r,w,R0,s,T).*alpha(r,m,s);

% Generate a mesh grid of r and m values
r_vals = linspace(-5, 5, 500);
m_vals = linspace(-5, 5, 500);
[R, M] = meshgrid(r_vals, m_vals);

% Compute Phi on the grid
Z = Phi(R, M);

% Produce the contour plot for Phi(r, m) = 0
figure;
contour(R, M, Z, [0 0], 'LineWidth', 2);
xlabel('r (resident trait value)');
ylabel('m','Rotation',0);
title('Contour Plot: $\Phi(r, m) = 0$','Interpreter','latex');
axis([-5 5 -5 5]);
grid on;
set(gca,'FontSize',24);

%% PIP with trait substitution sequence added

% Parameters
n_steps = 70; % number of mutations
step_size = 0.35;
initial_trait = -4;

% Trait values and meshgrid
traits = linspace(-5,5,200);
[X,Y] = meshgrid(traits, traits);

invasion_fitness = @(m, r) r_X(m,w,R0,s,T) - r_X(r,w,R0,s,T).*alpha(r,m,s);

% Simulate TSS
[residents, mutants] = simulate_tss_single_mutant(initial_trait, n_steps, step_size, traits, invasion_fitness);

% Compute PIP
pip = invasion_fitness(Y, X);

% Plot PIP and TSS
figure;
hold on;
scatter(residents(1), residents(1), 250, 'k','filled',...
    'MarkerEdgeColor','k');
contour(X, Y, pip, [0 0], 'LineWidth', 2);

%contour(X, Y, pip, [-10 0 10], 'LineColor', 'none', 'FaceAlpha',0.7);
scatter(residents(1:end-1), mutants, 110,'filled'); % Mutant steps
scatter(residents, residents, 110,'filled'); % Resident steps
scatter(residents(end), residents(end), 180, 'red','filled',...
    'MarkerEdgeColor','k');
xlabel('Resident trait');
ylabel('Mutant trait');
title('Trait Substitution Sequence on Pairwise Invasibility Plot');
legend('Starting trait value','PIP (zero invasion fitness)','Mutant trait','Resident trait',...
    'final trait value','Location','northeast');
hold off;

% Compute the PIP
function [residents, mutants] = simulate_tss_single_mutant(initial_trait, n_steps, step_size, traits, invasion_fitness)
    residents = zeros(n_steps+1, 1);
    mutants = zeros(n_steps, 1);
    residents(1) = initial_trait;
    for i = 1:n_steps
        res = residents(i);
        % Generate a single candidate mutant nearby
        mut = res + (2*rand-1) * step_size;
        mut = min(max(mut, traits(1)), traits(end));
        mutants(i) = mut;
        % If mutant is successful, update resident; else, resident doesn't change
        if invasion_fitness(mut, res) > 0
            residents(i+1) = mut;
        else
            residents(i+1) = res;
        end
    end
end