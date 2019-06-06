clc; clear; clear

%% varios contornos
% t = linspace(0,2*pi,100);
% x = cos(t(1:end-1))' + rand(99,1)*.2;
% y = sin(t(1:end-1))' + rand(99,1)*.2;
% scatter(x,y)
% X = [x,y];  Y = X;


%% Toroide
% theta = linspace(0,2*pi,30);
% gamma = theta;
% R = 1; r = .5;
% [Theta, Gamma ] = meshgrid(theta,gamma);
% x = cos(Theta).*(R+r*cos(Gamma)); x = x(:);
% y = sin(Theta).*(R+r*cos(Gamma)); y = y(:);
% z = r*sin(Gamma);  z = z(:);
% 
% scatter3(x,y,z); axis equal
% X = [x,y,z];
% Y = X;

%% Load point cloud;
% xx = pcread('/home/felipe/Desktop/Imagen_medic/PointCloud/Heart2.ply')
% X = xx.Location; X = downsample(X,3);
% Y = X;

%%


%X = [rand(50,2);rand(40,2)+[2,2]; rand(40,2)+[1,2];rand(100,2)+[4,4] ];
% 
% xx = pcread('/home/felipe/Desktop/Imagen_medic/PointCloud/HeartSlice.ply')
% X = xx.Location; Y = X;

load('/home/felipe/Desktop/Imagen_medic/PointCloud/imgSLICE_CONT.mat');
X = cell2mat(imgSLICE_CONT{1});
Y = cell2mat(imgSLICE_CONT{end});

figure; hold on
scatter(X(:,1),X(:,2),'.k');
scatter(Y(:,1),Y(:,2),'.r');
%% Normalizar and scalar datos
% n=size(X,1);
% m=size(Y,1);
% 
% normal.xd=mean(X);
% normal.yd=mean(Y);
% 
% X=X-repmat(normal.xd,n,1);
% Y=Y-repmat(normal.yd,m,1);
% 
% normal.xscale=sqrt(sum(sum(X.^2,2))/n);
% normal.yscale=sqrt(sum(sum(Y.^2,2))/m);
% 
% X=X/normal.xscale;
% Y=Y/normal.yscale;

%%
clearvars -except X Y Yc Xc  normal



%% Clustering


% %-------------normalization
n=size(X,1);
xd=mean(X);

X=X-repmat(xd,n,1);
xscale=sqrt(sum(sum(X.^2,2))/n);
X=X/xscale;


%------------------------

gauss = @(X,Y,s)(exp(-pdist2(X,Y).^2/(2*s^2)));
s = .5;

figure(3)

Xi = X;
pf = 1e6;
for ii = 1:200
    clf
    hold on
    scatter(X(:,1),X(:,2),'b')
    Kx = gauss(Xi,Xi,s);
    x_thao = (Kx*Xi);
    Xi = bsxfun(@times,x_thao,1./sum(Kx,2));
    scatter(Xi(:,1),Xi(:,2),'r')
    pf-mean(diag(pdist2(X,Xi)))
    if  1e-5 >  abs( pf-mean(diag(pdist2(X,Xi))))
        break;
    end
    pf = mean(diag(pdist2(X,Xi)));
    hold off
    pause(.5)
end
xd = mean(X);
pp = regionprops(X,'Centroid');
hold on;plot(xd(1),xd(2),'*g')
hold on;plot(pp.Centroid(1),pp.Centroid(2),'*m')

%%
%W0 = ones(numel(Y)/2,2);
maxIter = 1000;

gauss = @(X,Y,s)(exp(-pdist2(X,Y).^2/(2*s^2)));
V = @(X,Y,s)mean(mean(gauss(X,Y,s)));

s = .5;

N = size(X,1);
M = size(Y,1);
D = size(X,2);

Xi = X;

b = 10;

cmap = parula(N);
figure(2);
for i = 1:maxIter
    i
    Vxy = V(Xi,Y,s);
    Vx = V(Xi,Xi,s);
    c = Vxy*M/Vx/N;
    %c = 1.2;
    
    Kxy = gauss(Xi,Y,s);
    Kx = gauss(Xi,Xi,s);
    
    X1 = c*(1-b)/b*Kx*Xi;
    X2 = Kxy*Y;
    X3 = bsxfun(@times,c*(1-b)/b*sum(Kx,2),Xi);
    Xi = bsxfun(@plus,X1,X2-X3);
    Xi = bsxfun(@times,Xi,1./sum(Kxy,2));
    
    if  mod(i,1) == 0
        if size(X,2) == 3
            clf
            hold on; axis equal
            scatter3(X(:,1),X(:,2),X(:,3),'.k');
            scatter3(Xi(:,1),Xi(:,2),Xi(:,3),'.m');
            view([0 1 1])
            drawnow
            hold off;
            
        else
            clf
            hold on; axis equal
            scatter(X(:,1),X(:,2),'.k');
            scatter(Y(:,1),Y(:,2),'.r');
            scatter(Xi(:,1),Xi(:,2),'m','filled');
            drawnow
            hold off;
            
        end
    end
end









