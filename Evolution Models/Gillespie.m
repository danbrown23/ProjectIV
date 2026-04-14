%% Summary
% Implements the Gillespie SSA algorithm for the logistic competition model

n = 2;
X0 = [1;1];
r = 1;
K = [10;10];
alpha = eye(n);
tfinal = 100;

X = X0;
t = 0;

while t < tfinal
    % 1. Generate random numbers r1,r2 uniformly:
    r1 = rand;
    r2 = rand;
    % 2. Compute the propensity function of the system:
    propensity = sum(r.*X + r.*X./K.*alpha*X);
    % 3. Compute the next reaction time
    tau = 1/propensity*log(1/r1);
    t = t + tau;
    % 4. Find which of the n(n+1) reactions took place
    propensities = zeros(n*(n+1),1);
    start = 1;
    for j = 1:n
        propensities(start:j*n) = r.*X./K.*alpha(j).*X(j); % death due to competition
        propensities(start+j-1)
        start = start + n;
    end
    propensities(n^2:n^2+n) = r.*X; % birth
    propensities(n:n:n^2) = r./K.*X.*(X-1); % death due to inter-species competition adjustment
    propensities
end