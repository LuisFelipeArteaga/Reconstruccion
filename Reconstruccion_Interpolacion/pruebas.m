
clc; clear; close all
load XX.mat

t = linspace(0,2*pi,100);
t2 = linspace(0,1,100);
c1x = cos(t(1:end-1)); c1y = sin(t(1:end-1));

pointsC2 =[-0.5716    0.1771;-0.3528   -0.0458;0.1214   -0.1742;0.1795   -0.2175;0.3009   -0.5162;0.6160   -0.3615;0.2928   -0.0554;0.2527    0.0020;0.0996    0.5063;-0.2879    0.5342;-0.5716    0.1771];
%X = .25*[3*cos(t(1:end-1))+cos(3*t(1:end-1)); 3*sin(t(1:end-1))-sin(3*t(1:end-1))]';
X = [.6*c1x;.3*c1y]';
Y  = [.8*c1x;1.6*c1y]'; 
figure; hold on; axis equal
scatter(X(:,1),X(:,2))
scatter(Y(:,1),Y(:,2))


%%
mx = max([X;Y]);
mn = min([X;Y]);
dx = linspace(mx(1)+.2,mn(1)-.2,5)
dy = linspace(mx(2)+.2,mn(2)-.2,5)
[XX, YY] = meshgrid(dx,dy);
ctrl_pts =[XX(:) YY(:)];
scatter(XX(:),YY(:))


[param, model] = cpd(config,X,Y,ctrl_pts)
% DisplayPoints(model,Y,2);


