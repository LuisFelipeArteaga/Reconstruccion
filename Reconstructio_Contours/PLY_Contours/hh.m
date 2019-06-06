close all; clear all;clc
load Cont.mat;
pc1= pcread('c4.ply');
pc2= pcread('c2.ply');
pc3= pcread('c3.ply');
pc4= pcread('c4.ply');



aa = pointCloud(pc3.Location+ [.05 0  -.1]);
pcwrite(aa,'1c3c.ply')

% hold on;
% scatter3(cc{1}(:,1),zeros(100,1),cc{1}(:,2),'*')
% scatter3(cc{2}(:,1),zeros(100,1)+.5,cc{2}(:,2),'*')