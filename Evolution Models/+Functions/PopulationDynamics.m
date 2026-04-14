function dXdt = PopulationDynamics(X,r,K,alpha)
% Input: population vector X, precalculated growth rate vector r, carrying
% capacity vector K, and competition matrix alpha
% Output: instantaneous rate of change of population vector
    dXdt = r.*X.*(1-alpha*X./K);
end