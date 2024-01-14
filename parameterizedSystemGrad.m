function dydt = parameterizedSystemGrad(~, y, alpha, beta, gamma, delta, epsilon, lambda, N)

dydt = zeros(2*N,1);

E = delta/(8*sqrt(pi)*epsilon^(5/3));

for i=1:N
    xTerm = 0;
    yTerm = 0;
    for j=1:N
        if j ~= i
            % r^2 = (x1^2 + x2^2)
            % x1^2 = (y(i) - y(j)).^2
            % x2^2 = (y(i+N) - y(j+N)).^2
            
            x1 = (y(i) - y(j));
            y1 = (y(i+N) - y(j+N));
            r  = (x1.^2 + y1.^2).^(1/2);

            A = E*exp(-r^2/(4*epsilon));

            if alpha == 0
                B = lambda*r.^(-2);
            else
                B = alpha*lambda*r.^(-alpha-2);
            end

            if beta == 0
                C = gamma*r.^(-2);
            else
                C = beta*gamma*r.^(beta-2);
            end

            coeff = A - B + C;
            
            xTerm = xTerm + (coeff .* x1);
            yTerm = yTerm + (coeff .* y1);
        end
    end
    dydt(i) = (-1/N) .* xTerm;
    dydt(i+N) = (-1/N) .* yTerm;
end