function f = Fitness(x,y,X,sd_K,sd_a,t)
import Functions.Growthrate
import Functions.CarryingCapacity
import Functions.Competition
% Input: resident phenotype matrix and mutant phenotype matrix
% vector of the mutant species in an environment created by community x
    n = size(x,2);
    phenotype = [x y];
    X = X(:);
% define the bottom left corner of the competition matrix we need
    beta = Competition(phenotype,sd_a);
    beta = beta(n+1:2*n,1:n);
% compute the fitness
    f = Growthrate(y).*(1-beta*X./CarryingCapacity(y,sd_K,t));
end