 clear; clc

%load data/cpd_data2D_fish.mat
% t = linspace(0,2*pi,200);
% t2 = linspace(0,1,200);
% c1x = 1*cos(t(1:end-1)); c1y = 1*sin(t(1:end-1));
% 
% pointsC2 =[-0.5716    0.1771;-0.3528   -0.0458;0.1214   -0.1742;0.1795   -0.2175;0.3009   -0.5162;0.6160   -0.3615;0.2928   -0.0554;0.2527    0.0020;0.0996    0.5063;-0.2879    0.5342;-0.5716    0.1771];
% %Y = .2*[3*cos(t)+cos(3*t); 3*sin(t)-sin(3*t)]'; Y = Y(1:end-1,:);
% X = interparc(t2,pointsC2(:,1),pointsC2(:,2)); X = X(1:end-1,:);
% %X = .5*[c1x;c1y]';
% Y  = [c1x;c1y]';
load('/home/felipe/Desktop/Imagen_medic/PointCloud/imgSLICE_CONT.mat');
X = cell2mat(imgSLICE_CONT{1});
Y = cell2mat(imgSLICE_CONT{end});

clearvars -except X Y


%% LLE
%%
[M,D]=size(Y); [N, D2]=size(X);
opt.method = 'nonrigid';
opt.normalize = 0;
opt.max_it = 100;
opt.tol = 1e-5;
opt.viz = 1;
opt.corresp = 1;
opt.outliers = .01;
opt.fgt = 0;
opt.sigma2 = 1;

% reflejan la cantidad de regularización de la suavidad.
opt.beta = .1;
opt.lambda = 20;

addpath('mex/')
[Transform, C]=cpd_register(Y, X, opt);
cmap = parula(numel(C));
T = Transform.Y;
figure('Name','CPD'); hold on; axis equal
scatter(Y(:,1),Y(:,2),10,'b')
scatter(X(:,1),X(:,2),[], 'r')
for   ii =1: numel(C)
    plot([X(ii,1) T(ii,1)],[X(ii,2) T(ii,2)],'Color',cmap(ii,:))
end


