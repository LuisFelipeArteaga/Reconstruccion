%% 
clear; close all; clc

%%

t = linspace(0,2*pi,100);
t2 = linspace(0,1,100);
c1x = 1.2*cos(t(1:end-1)); c1y = .9*sin(t(1:end-1));
Xc1  = [c1x;c1y]';
pointsC2 =[-0.5716    0.1771;-0.3528   -0.0458;0.1214   -0.1742;0.1795   -0.2175;0.3009   -0.5162;0.6160   -0.3615;0.2928   -0.0554;0.2527    0.0020;0.0996    0.5063;-0.2879    0.5342;-0.5716    0.1771];
Xc2 = interparc(t2,pointsC2(:,1),pointsC2(:,2));

%%  

figure(2); hold on
plot(c1x,c1y,'.b')
plot(Xc2(:,1),Xc2(:,2),'.b')
plot(pointsC2(:,1),pointsC2(:,2),'+k')

%%  kernel 

Ds      = pdist(Xc1);
sm      = median(Ds);
%m      = kScaleOptimization_info_alpha(Ds,0.5);   
Kc1     = exp(-squareform(Ds).^2/(2*sm^2));



Ds      = pdist(Xc2);
sm      = median(Ds);
%m      = kScaleOptimization_info_alpha(Ds,0.5);   
Kc2    = exp(-squareform(Ds).^2/(2*sm^2));


