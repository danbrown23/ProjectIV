function X = EffectiveCarryingCapacity(x,X0,var_K,var_a,tspan,t) % (allowing extinction)
% Input: x is the (traits x species) phenotype configuration matrix, X0 is
% the initial population vector, tspan is the integration timespan
% Output: X is the final population vector 
import Functions.Growthrate
import Functions.CarryingCapacity
import Functions.Competition
import Functions.PopulationDynamics
    r = Growthrate(x);
    K = CarryingCapacity(x,var_K,t);
    alpha = Competition(x,var_a);
    [~, X] = ode45(@(t,X) PopulationDynamics(X,r,K,alpha),tspan,X0);
    X = X(end,:);
    X = X(:);
    X(X<0.001) = 0;
end