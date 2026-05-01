function X = effectivecarryingcapacity(x,X0,sd_K,sd_a,n,tspan) % (allowing extinction)
% given a phenotype matrix and an initial population vector, returns the effective carrying capacity vector
r = 1;
D2 = sum((x).^2,1);
K = exp(-D2./(2.*sd_K^2));
K = K(:);
        alpha = zeros(n,n);
        for i = 1:n
            for j = 1:n
                alpha(i,j) = exp((-norm(x(:,i)-x(:,j)).^2)./(2.*sd_a^2));
            end
        end
[~, X] = ode45(@(t,X) r.*X.*(1-alpha*X./K),tspan,X0);
X = X(end,:);
X = X(:);
X(X<0.001) = 0;
end