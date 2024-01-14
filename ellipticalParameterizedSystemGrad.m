function dydt = ellipticalParameterizedSystemGrad(~, y, alpha, beta, gamma, delta, epsilon, lambda, a, b, N)

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
            
            % r is sqrt(x^2+y^2)
            r  = (x1.^2 + y1.^2).^(1/2);

            % relip = r_e = sqrt((ax)^2+(by)^2)
            relip = (a^2*(x1.^2) + b^2*(y1.^2)).^(1/2);

            A = E*exp(-relip^2/(4*epsilon));
%             A = E*exp(-r^2/(4*epsilon));
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

            xCoeff = (a^2*A) - B + C;
            yCoeff = (b^2*A) - B + C;
            
            xTerm = xTerm + (xCoeff .* x1);
            yTerm = yTerm + (yCoeff .* y1);
        end
    end
    dydt(i) = (-1/N) .* xTerm;
    dydt(i+N) = (-1/N) .* yTerm;
end