function genModel = parameterizedGenModelElliptical(alpha,beta,gamma,delta,epsilon,lambda,a,b,N,tMax)

alpha = str2double(alpha);
beta = str2double(beta);
gamma = str2double(gamma);
delta = str2double(delta);
epsilon = str2double(epsilon);
lambda = str2double(lambda);
a = str2double(a);
b = str2double(b);
N = str2double(N);
tMax = str2double(tMax);

% alpha = 0;
% beta = 0;
% gamma = 0;
% delta = 1;
% epsilon = 0.01;
% lambda = 0.5;
% a = 1;
% b = 0.5;
% N = 50;
% tMax = 100;

kernelParamList = strcat('alpha=',num2str(alpha),',beta=', ...
    num2str(beta),',gamma=',num2str(gamma),',delta=',num2str(delta), ...
    ',epsilon=',num2str(epsilon),',lambda=',num2str(lambda),',a=', ...
    num2str(a),',b=',num2str(b));

fullParamList = strcat(kernelParamList,',N=',num2str(N), ...
    ',tMax=',num2str(tMax));

% start the timer
tic

% start with the particles in a random circular distribution
r=zeros(N,2);
radius = 0.5;
xc = 0;
yc = 0;
theta = rand(1,N)*(2*pi);
r = sqrt(rand(1,N))*radius;
x = xc + r.*cos(theta);
y = yc + r.*sin(theta);
r(1:N,1)=x;
r(1:N,2)=y;
y0=zeros(2*N,1);  
y0(1:N,1)=r(:,1);
y0(N+1:2*N,1)=r(:,2);

tspan = [0 tMax];
sol = ode45(@(t,y) ellipticalParameterizedSystemGrad(t, y, alpha, beta, gamma, delta, epsilon, lambda, a, b, N), tspan, y0);

xpoints = sol.y(1:N,:);
ypoints = sol.y(N+1:2*N,:);
m = size(xpoints,2);

for i=0:m
    xcol = xpoints(:,m-i);
    ycol = ypoints(:,m-i);
    % scan the table from right to left. if there aren't any nan values,
    % then we have a valid set of points to plot
    if ~(any(isnan(xcol)) || any(isnan(ycol)))        
        f = figure(1);

        % determine max(d(x,y)) and min(d(x,y))
        maxD = -1;
        minD = Inf;
        for j=1:N
            for k=1:j-1
                p1x = xcol(j);
                p1y = ycol(j);
                p2x = xcol(k);
                p2y = ycol(k);
                dist = sqrt((p1x-p2x)^2+(p1y-p2y)^2);
                if dist > maxD
                    maxD = dist;
                end
                if dist < minD
                    minD = dist;
                end
            end
        end
        
        % determine x and y axis scale
        minX = min(xcol);
        maxX = max(xcol);
        minY = min(ycol);
        maxY = max(ycol);
        xLow = minX - abs(maxX - minX)/10;
        xHigh = maxX + abs(maxX - minX)/10;
        yLow = minY - abs(maxY - minY)/10;
        yHigh = maxY + abs(maxY - minY)/10;

        % determine runTime
        runTime = toc;

        % format plot title
        paramText = strcat('$\alpha = ',num2str(alpha),'\quad\beta = ', ...
            num2str(beta),'\quad\gamma = ',num2str(gamma),'\quad\delta = ', ...
            num2str(delta),'\quad\varepsilon = ',num2str(epsilon), ...
            '\quad\lambda = ',num2str(lambda),'\quad N = ',num2str(N), ...
            '\quad t_{max} = ',num2str(tMax),'$');

        ellipticalText = strcat('$(ax)^2 + (by)^2=(', ...
            num2str(a),'x)^2 + (',num2str(b),'y)^2$');

        distText = strcat('$\min(d(x_i,x_j)) = ',num2str(minD), ...
            '\quad\max(d(x_i,x_j)) = ',num2str(maxD),'$');
        
        distText2 = strcat('$\min_x = ',num2str(minX),'\quad\max_x = ',num2str(maxX), ...
            '\quad\min_y =',num2str(minY),'\quad\max_y = ',num2str(maxY),'$');

        timeText = strcat('$\mathrm{Iteration} = ',num2str(m-i), ...
            '\quad\mathrm{Running}\;\mathrm{time} =',num2str(runTime),'\;\mathrm{seconds}$'); 
        
        plot(xcol, ycol,'k.','MarkerSize',10);
        ax = gca;
        ax.TitleFontSizeMultiplier = 1.5;
        axis([xLow xHigh yLow yHigh])
        f.Position = [10 10 1000 735]; 
        axis equal
        axis on
        title({paramText ellipticalText distText distText2 timeText},'Interpreter','latex')
        print('-dpng','-r150', ...
            strcat(fullParamList,',iteration=',num2str(m-i),'.png'))
        close(1)
        break
    end
end
