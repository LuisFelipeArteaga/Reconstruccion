%% Non-rigid Point Set Registration with Global-Local Topology Preservation
 clear; clc
%load data/cpd_data2D_fish.mat
t = linspace(0,2*pi,100);
t2 = linspace(0,1,100);
c1x = 1*cos(t(1:end-1)); c1y = 1*sin(t(1:end-1));

pointsC2 =[-0.5716    0.1771;-0.3528   -0.0458;0.1214   -0.1742;0.1795   -0.2175;0.3009   -0.5162;0.6160   -0.3615;0.2928   -0.0554;0.2527    0.0020;0.0996    0.5063;-0.2879    0.5342;-0.5716    0.1771];
%Xc1 = .25*[3*cos(t)+cos(3*t); 3*sin(t)-sin(3*t)]'; Xc1 = Xc1(1:end-1,:);
%X = interparc(t2,pointsC2(:,1),pointsC2(:,2)); X = X(1:end-1,:);
X  = .5*[c1x;c1y]';
Y  = [c1x;c1y]';

clearvars -except X Y

%%
[M,D]=size(Y); [N, D2]=size(X);
opt.method = 'nonrigid';
opt.normalize = 1;
opt.max_it = 250;
opt.tol = 1e-5;
opt.viz = 1;
opt.corresp = 1;
opt.outliers = .01;
opt.fgt = 0;
opt.sigma2 = 1;

% reflejan la cantidad de regularización de la suavidad.
opt.beta = 2;
opt.lambda = 1;
opt.K = 5;
opt.alpha = 3;

[Transform, C, T,Y]=GLTP_Algorithm(Y,X, opt);