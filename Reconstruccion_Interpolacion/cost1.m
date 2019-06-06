clc; clear; close all;

%X = [2 1.5; 2.5 2; 3.1 1];
%Y = [0.5 1; 2 4.1; 5.3 1];

t = linspace(0,2*pi,30);
t2 = linspace(0,1,30);
c1x = 1.2*cos(t(1:end-1)); c1y = .9*sin(t(1:end-1));
Xc1  = [c1x;c1y]';
pointsC2 =[-0.5716    0.1771;-0.3528   -0.0458;0.1214   -0.1742;0.1795   -0.2175;0.3009   -0.5162;0.6160   -0.3615;0.2928   -0.0554;0.2527    0.0020;0.0996    0.5063;-0.2879    0.5342;-0.5716    0.1771];
Xc2 = interparc(t2,pointsC2(:,1),pointsC2(:,2));

X = (Xc2(1:end-1,:)+ones(1,2)*5);
Y = (Xc1)+ones(1,2)*5;

X2 = X;
% Y2 = Y;

W0 = ones(29,2);

maxIter = 200;
lr      = 0.1;
figure(1),

loop = 1;
flag = 1;
wchange = 1;

for i = 1:maxIter
    i
    clf
    hold on;
    scatter(X2(:,1),X2(:,2),'b','filled');
    scatter(Y(:,1),Y(:,2),'filled');
     
    xlim([3 7])
    ylim([3 7])
    
    W1 = W0-lr*((1/size(X,1))*((W0.*X)-Y).*X);
    Xnew = X.*W1;
    
     scatter(Xnew(:,1),Xnew(:,2),'b','filled');
    hold off;
    
    W0 = W1;
    X = Xnew;
    pause(0.01);
end
hold on
for ii = 1:numel(c1x)
    line([X2(ii,1) X(ii,1)],[X2(ii,2) X(ii,2)])
end