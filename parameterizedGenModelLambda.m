function genModel = parameterizedGenModelLambda(alpha,beta,gamma,delta,epsilon,lambda,N,tMax)

alpha = str2double(alpha);
beta = str2double(beta);
gamma = str2double(gamma);
delta = str2double(delta);
epsilon = str2double(epsilon);
lambda = str2double(lambda);
% a = str2double(a);
% b = str2double(b);
N = str2double(N);
tMax = str2double(tMax);

% alpha = 0.1;
% beta = 0;
% gamma = 0;
% delta = 1;
% epsilon = 0.01;
% lambda = 0.3;
% % a = 4;
% % b = 1/4;
% N = 30;
% tMax = 30;

kernelParamList1 = strcat('alpha=',num2str(alpha),',beta=', ...
    num2str(beta),',gamma=',num2str(gamma),',delta=',num2str(delta), ...
    ',epsilon=',num2str(epsilon),',lambda=',num2str(lambda));

kernelParamList2 = strcat('alpha=',num2str(alpha),',beta=', ...
    num2str(beta),',gamma=',num2str(gamma),',delta=',num2str(delta), ...
    ',epsilon=',num2str(epsilon),',lambda=0');

fullParamList1 = strcat(kernelParamList1,',N=',num2str(N), ...
    ',tMax=',num2str(tMax));

fullParamList2 = strcat(kernelParamList2,',N=',num2str(N), ...
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
solLambda =   ode45(@(t,y) parameterizedSystemGrad(t, y, alpha, beta, gamma, delta, epsilon, lambda, N), tspan, y0);

xPoints1 = solLambda.y(1:N,:);
yPoints1 = solLambda.y(N+1:2*N,:);
mLambda = size(xPoints1,2);

xPoints2 = solNoLambda.y(1:N,:);
yPoints2 = solNoLambda.y(N+1:2*N,:);
mNoLambda = size(xPoints2,2);

% determine runTime
runTime = toc;

% code for making a graph of the kernel
r = linspace(0,10,200);

% for lambda =/= 0
lambdaStr1 = num2str(lambda);

gammaStr = strcat('(',num2str(gamma,'%G'),')');

if delta == 0
    A = '';
    Afunc = 0;
else
    A = strcat('-\frac{',num2str(delta),'}{4\sqrt{\pi}\left(', ...
        num2str(epsilon),'^{3/2}\right)}e^{-\frac{r^2}{', ...
        num2str(4*epsilon),'}}');
    Afunc = (-delta/(4*sqrt(pi)*epsilon))*exp(-r.^2/(4*epsilon));
end

% for lambda =/= 0
if alpha ~= 0
    B1 = strcat('+\frac{',lambdaStr1,'}{r^',num2str(alpha),'}');
    Bfunc1 = lambda*r.^(-alpha);
else
    B1 = strcat('+',lambdaStr,'\log\left(\frac{1}{r}\right)');
    Bfunc1 = lambda*log(r.^(-1));
end

% for lambda = 0
B2 = '';
Bfunc2 = 0;

if gamma == 0
    C = '$';
    Cfunc = 0;
elseif beta ~= 0
    C = strcat('+',gammaStr,'r^',num2str(beta),'$');
    Cfunc = gamma*r.^beta;
else
    C = strcat('+',gammaStr,'\log(r)$');
    Cfunc = gamma*log(r);
end

% lambda =/= 0 case
figure(1);
plot(r,Afunc+Bfunc1+Cfunc)
ax = gca;
ax.TitleFontSizeMultiplier = 1.5;
title(strcat('$K(r) = ',A,B1,C),'Interpreter','latex');
xlabel('r');
ylabel('K(r)');
print('-dpng','-r150', ...
    strcat('kernel-',kernelParamList1,'.png'))
close(1)

% lambda = 0 case
figure(1);
plot(r,Afunc+Bfunc2+Cfunc)
ax = gca;
ax.TitleFontSizeMultiplier = 1.5;
title(strcat('$K(r) = ',A,B2,C),'Interpreter','latex');
xlabel('r');
ylabel('K(r)');
print('-dpng','-r150', ...
    strcat('kernel-',kernelParamList2,'.png'))
close(1)

% both case
figure(1);
p1 = plot(r,Afunc+Bfunc1+Cfunc);
p1.Color = 'red';
hold on
p2 = plot(r,Afunc+Bfunc2+Cfunc);
p2.Color = 'blue';
hold off
ax = gca;
ax.TitleFontSizeMultiplier = 1.5;
title({strcat('In red: $K(r) = ',A,B1,C), ...
    strcat('In blue: $K(r) = ',A,B2,C)},'Interpreter','latex');
xlabel('r');
ylabel('K(r)');
print('-dpng','-r150', ...
    strcat('kernel2-',kernelParamList1,'.png'))
close(1)

% results for lambda =/= 0
for i=0:mLambda
    xcol1 = xPoints1(:,mLambda-i);
    ycol1 = yPoints1(:,mLambda-i);
    if ~(any(isnan(xcol1)) || any(isnan(ycol1)))        
        % determine max(d(x,y)) and min(d(x,y))
        maxD1 = -1;
        minD1 = Inf;
        for j=1:N
            for k=1:j-1
                p1x = xcol1(j);
                p1y = ycol1(j);
                p2x = xcol1(k);
                p2y = ycol1(k);
                dist = sqrt((p1x-p2x)^2+(p1y-p2y)^2);
                if dist > maxD1
                    maxD1 = dist;
                end
                if dist < minD1
                    minD1 = dist;
                end
            end
        end

        % determine x and y axis scale
        minX1 = min(xcol1);
        maxX1 = max(xcol1);
        minY1 = min(ycol1);
        maxY1 = max(ycol1);
        xLow1 = minX1 - abs(maxX1 - minX1)/10;
        xHigh1 = maxX1 + abs(maxX1 - minX1)/10;
        yLow1 = minY1 - abs(maxY1 - minY1)/10;
        yHigh1 = maxY1 + abs(maxY1 - minY1)/10;

        lambdaIter = mLambda-i;
        break
    end
end

% results for lambda = 0
for i=0:mNoLambda
    xcol2 = xPoints2(:,mNoLambda-i);
    ycol2 = yPoints2(:,mNoLambda-i);
    if ~(any(isnan(xcol2)) || any(isnan(ycol2)))
        % determine max(d(x,y)) and min(d(x,y))
        maxD2 = -1;
        minD2 = Inf;
        for j=1:N
            for k=1:j-1
                p1x = xcol2(j);
                p1y = ycol2(j);
                p2x = xcol2(k);
                p2y = ycol2(k);
                dist = sqrt((p1x-p2x)^2+(p1y-p2y)^2);
                if dist > maxD2
                    maxD2 = dist;
                end
                if dist < minD2
                    minD2 = dist;
                end
            end
        end

        % determine x and y axis scale
        minX2 = min(xcol2);
        maxX2 = max(xcol2);
        minY2 = min(ycol2);
        maxY2 = max(ycol2);
        xLow2 = minX2 - abs(maxX2 - minX2)/10;
        xHigh2 = maxX2 + abs(maxX2 - minX2)/10;
        yLow2 = minY2 - abs(maxY2 - minY2)/10;
        yHigh2 = maxY2 + abs(maxY2 - minY2)/10;

        noLambdaIter = mNoLambda-i;
        break
    end
end

% format plot title
lambdaParamText = strcat('$\alpha = ',num2str(alpha),'\quad\beta = ', ...
    num2str(beta),'\quad\gamma = ',num2str(gamma),'\quad\delta = ', ...
    num2str(delta),'\quad\varepsilon = ',num2str(epsilon), ...
    '\quad\lambda = ',num2str(lambda),'\quad N = ',num2str(N), ...
    '\quad t_{max} = ',num2str(tMax),'$');

noLambdaParamText = strcat('$\alpha = ',num2str(alpha),'\quad\beta = ', ...
    num2str(beta),'\quad\gamma = ',num2str(gamma),'\quad\delta = ', ...
    num2str(delta),'\quad\varepsilon = ',num2str(epsilon), ...
    '\quad\lambda = 0\quad N = ',num2str(N), ...
    '\quad t_{max} = ',num2str(tMax),'$');

% ellipticalText = strcat('$(ax)^2 + (by)^2=(', ...
%     num2str(a),'x)^2 + (',num2str(b),'y)^2$');

lambdaDistText = strcat('$\min(d(x_i,x_j)) = ',num2str(minD1), ...
    '\quad\max(d(x_i,x_j)) = ',num2str(maxD1),'$');

noLambdaDistText = strcat('$\min(d(x_i,x_j)) = ',num2str(minD2), ...
    '\quad\max(d(x_i,x_j)) = ',num2str(maxD2),'$');

lambdaTimeText = strcat('$\mathrm{Iteration} = ',num2str(lambdaIter), ...
    '\quad\mathrm{Running}\;\mathrm{time} =',num2str(runTime),'\;\mathrm{seconds}$'); 

noLambdaTimeText = strcat('$\mathrm{Iteration} = ',num2str(noLambdaIter), ...
    '\quad\mathrm{Running}\;\mathrm{time} =',num2str(runTime),'\;\mathrm{seconds}$'); 

timeText = strcat('$\mathrm{Iteration} = ',num2str(max(lambdaIter, noLambdaIter)), ...
    '\quad\mathrm{Running}\;\mathrm{time} =',num2str(runTime),'\;\mathrm{seconds}$'); 

% plot the lambda =/= 0 results
f = figure(1);
plot(xcol1, ycol1,'k.','MarkerSize',10);
ax = gca;
ax.TitleFontSizeMultiplier = 1.5;
axis([xLow1 xHigh1 yLow1 yHigh1])
f.Position = [10 10 1000 735]; 
axis equal
axis on
title({lambdaParamText lambdaDistText lambdaTimeText},'Interpreter','latex')
print('-dpng','-r150', ...
    strcat(fullParamList1,',iteration=',num2str(lambdaIter),'.png'))
close(1)

% plot the lambda = 0 results
f = figure(1);
plot(xcol2, ycol2,'k.','MarkerSize',10);
ax = gca;
ax.TitleFontSizeMultiplier = 1.5;
axis([xLow2 xHigh2 yLow2 yHigh2])
f.Position = [10 10 1000 735]; 
axis equal
axis on
title({noLambdaParamText noLambdaDistText noLambdaTimeText},'Interpreter','latex')
print('-dpng','-r150', ...
    strcat(fullParamList2,',iteration=',num2str(noLambdaIter),'.png'))
close(1)

% plot both results on the same plot
f = figure(1);
p1 = plot(xcol1, ycol1,'k.','MarkerSize',10);
p1.Color = 'red';
hold on
p2 = plot(xcol2, ycol2,'k.','MarkerSize',10);
p2.Color = 'blue';
hold off
ax = gca;
ax.TitleFontSizeMultiplier = 1.5;
axis([min(xLow1, xLow2) max(xHigh1, xHigh2) ...
    min(yLow1, yLow2) max(yHigh1, yHigh2)])
f.Position = [10 10 1000 835]; 
axis equal
axis on
title({strcat('In red:',lambdaDistText) lambdaParamText ...
    strcat('In blue:',noLambdaDistText) noLambdaParamText ...
    timeText},'Interpreter','latex')
print('-dpng','-r150', ...
            strcat(fullParamList1,'(andLambda=0),iteration=',num2str(max(lambdaIter, noLambdaIter)),'.png'))
close(1)
