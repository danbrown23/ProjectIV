classdef SymmetricCompetition_Functions
    
    properties
        Value
    end

    methods(Static)
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

        function dXdt = populationdynamics(X,r,K,alpha)
            % given an instantaneous population vector and r, K, and alpha precalculated, returns the dynamics vector
            dXdt = r.*X.*(1-alpha*X./K);
        end

        function X = effectivecarryingcapacity(x,X0,sd_K,sd_a,n,tspan) % (allowing extinction)
            % given a phenotype matrix and an initial population vector, returns the effective carrying capacity vector
            r = self.growthrate(x);
            K = carryingcapacityM(x,sd_K,n);
            alpha = competitionM(x,sd_a,n);
            [~, X] = ode45(@(t,X) populationdynamics(X,r,K,alpha),tspan,X0);
            X = X(end,:);
            X = X(:);
            X(X<0.001) = 0;
        end

        %% fitness

        function f = fitness(x,y,X,sd_K,sd_a,n)
            % given a phenotype matrix and a mutant matrix, returns the fitness
            % vector of the mutant species in an environment created by community x
            phenotype = [x y];
            X = X(:);
            % define the bottom left corner of the competition matrix we need
            beta = zeros(n,n);
            for j = n+1:2*n
                for k = 1:n
                    beta(j-n,k) = competition(phenotype,j,k,sd_a);
                end
            end
            % compute the fitness
            f = growthrate(y).*(1-beta*X./carryingcapacityM(y,sd_K,n));
        end
    end
end
