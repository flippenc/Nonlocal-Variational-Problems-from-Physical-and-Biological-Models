function dydt = normParameterizedSystemGrad(~, y, alpha, beta, gamma, delta, epsilon, lambda, a, b, c, N)

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
            
            % r = (x^c + y^c)^1/c
            r = (x1.^c + y1.^c).^(1/c);

            % relip = r_e = ((ax)^c+(by)^c)^1/c
            relip = ((a*x1).^c + (b*y1).^c).^(1/c);

            A = E*(c*a^c*x1.^(c-1)/2)*exp(-relip^2/(4*epsilon));

            if alpha == 0
                B = lambda*r.^(-c);
            else
                B = alpha*lambda*r.^(-alpha-c);
            end

            if beta == 0
                C = gamma*r.^(-c);
            else
                C = beta*gamma*r.^(beta-c);
            end

            xCoeff = (a^c*A) - B + C;
            yCoeff = (b^c*A) - B + C;
            
            xTerm = xTerm + (xCoeff .* (x1.^(c-1)));
            yTerm = yTerm + (yCoeff .* (y1.^(c-1)));
        end
    end
    dydt(i) = (-1/N) .* xTerm;
    dydt(i+N) = (-1/N) .* yTerm;
end