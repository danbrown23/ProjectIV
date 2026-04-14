function K = CarryingCapacity(x,var_K,t)
import Functions.CCC_Trajectory
% Input: x is the (traits x species) phenotype configuration matrix
% Output: K is the carrying capacity vector of the n species
d = size(x,1);
n = size(x,2);
v = repmat(CCC_Trajectory(t,d),1,n);
D2 = sum((x-v).^2,1);
K = exp(-D2./(2.*var_K));
end

