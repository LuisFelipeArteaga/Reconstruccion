clc; clear
t = linspace(0,2*pi,100);
t2 = linspace(0,1,100);
c1x = 1*cos(t(1:end-1)); c1y = 1*sin(t(1:end-1));

pointsC2 =[-0.5716    0.1771;-0.3528   -0.0458;0.1214   -0.1742;0.1795   -0.2175;0.3009   -0.5162;0.6160   -0.3615;0.2928   -0.0554;0.2527    0.0020;0.0996    0.5063;-0.2879    0.5342;-0.5716    0.1771];
%Xc1 = .25*[3*cos(t)+cos(3*t); 3*sin(t)-sin(3*t)]'; Xc1 = Xc1(1:end-1,:);
Xc1 = .5*[c1x;c1y]';
%Xc1 = interparc(t2,pointsC2(:,1),pointsC2(:,2)); Xc1 = Xc1(1:end-1,:);
Xc2  = [c1x;c1y]';

%%
[LXc1,AXc1,pendenza1,n1]=preprocessing(Xc1);
[LXc2,AXc2,pendenza2,n2]=preprocessing(Xc2);

LCc1 = cumsum(LXc1)/sum(LXc1);
LCc2 = cumsum(LXc2)/sum(LXc2);
%%
dx1 = gradient (Xc1(:,1)');
dy1 = gradient (Xc1(:,2)');

dx2 = gradient (Xc2(:,1)');
dy2 = gradient (Xc2(:,2)');

%quiver(c1x,c1y,-dx,dy)
dXY1 = cross([dx1; dy1;zeros(1,numel(dx1))],repmat([0 0 1]',1,numel(dx1))); dXY1 = dXY1(1:2,:)';
dXY2 = cross([dx2; dy2;zeros(1,numel(dx2))],repmat([0 0 1]',1,numel(dx2)));  dXY2 = -dXY2(1:2,:)';
% figure(1); hold on
% quiver(Xc1(:,1),Xc1(:,2),dXY1(:,1),dXY1(:,2))
% quiver(Xc2(:,1),Xc2(:,2),dXY2(:,1),dXY2(:,2))
% hold off
%%
% X = [Xc1 dXY1];
% Y = [Xc2 dXY2];
X = Xc1;
Y = Xc2;
clearvars -except X Y
%%

opt.method = 'nonrigid';
opt.normalize = 10;
opt.max_it = 1;
opt.tol = 1e-5;
opt.viz = 1;
opt.corresp = 0;
opt.outliers = 0.1;
opt.fgt = 0;
opt.sigma2 = 0;

% reflejan la cantidad de regularización de la suavidad.
opt.beta = 2;
opt.lambda = 3;

%%
W0 = ones(numel(Y)/2,2);
outliers  = 0.1;
maxIter = 1000;

gauss = @(X,Y,s)(exp(-pdist2(X,Y).^2/(2*s^2)));

V = @(X,Y,s)mean(mean(gauss(X,Y,s)));

s = median(pdist(Y));

N = size(X,1);
M = size(Y,1);
D = size(X,2);

Xi = X;

b = 50;

cmap = parula(N);

for i = 1:maxIter
    
    Vxy = V(Xi,Y,s);
    Vx = V(Xi,Xi,s);
    c = Vxy*M/Vx/N;
    
    Kxy = gauss(Xi,Y,s);
    Kx = gauss(Xi,Xi,s);
    
    X1 = c*(1-b)/b*Kx*Xi;
    X2 = Kxy*Y;
    X3 = bsxfun(@times,c*(1-b)/b*sum(Kx,2),Xi);
    Xi = bsxfun(@plus,X1,X2-X3);
    Xi = bsxfun(@times,Xi,1./sum(Kxy,2));
    
    sx = median(pdist(Xi));
    ksig = -2*sx*sx;
    c2 = (outliers*M*((2*pi*sx)^(D/2))) / ((1- outliers)*N);
    

%     figure(4); clf;
    %imagesc(Kxy)
    %drawnow
    
    
    K =  pdist2(Xi,Y).^2;
    K = K/ksig;
    K = exp(K);
%     
%     
     a =1./(K'*ones(N,1)+ c2*ones(M,1));
     da = diag(a);
     P = K*da;
     figure(3); clf;
     imagesc(P)
     drawnow
    
    
    i
    if 0
        figure(2);
        clf
        hold on; axis equal
        scatter(X(:,1),X(:,2),[],cmap,'filled');
        scatter(Y(:,1),Y(:,2),'filled');
        scatter(Xi(:,1),Xi(:,2),[],cmap,'filled');
        
        for j=1:N
            plot([X(j,1) Xi(j,1)],[X(j,2) Xi(j,2)],'Color',cmap(j,:))
        end
        
        hold off;
        drawnow
        
    end
end

clf
hold on;
scatter(X(:,1),X(:,2),[],cmap,'filled');
scatter(Y(:,1),Y(:,2),'filled');
scatter(Xi(:,1),Xi(:,2),[],cmap,'filled');

for j=1:N
    plot([X(j,1) Xi(j,1)],[X(j,2) Xi(j,2)],'Color',cmap(j,:))
end

hold off;

%% CPD









