clear; close all; clc
%%
addpath('Functions')
load Ro.mat

q0 = [1 0 0 0]';
q1 =rotm2quat(Roo{1})';
q0 = Quaternion.Unit(q0);
q1 = Quaternion.Unit(q1);

t = linspace(0,1,100);
qt = zeros(numel(t),4);
qt2 = zeros(numel(t),4);
 %% Slerp (Spline between  2 points)
for ii = 1:100
    qt(ii,:) = slerp(q0,q1,t(ii),1);
    q01 = Quaternion.Mult(Quaternion.Inv(q0),q1);
    q01 = Quaternion.Power(q01,t(ii));
    qt2(ii,:) = Quaternion.Mult(q0,q01);
end



%% Squad (Bezier Curve 4 points)

p = [1 0 0 0];
a = rotm2quat(Roo{1});
b =rotm2quat(Roo{2});
q =rotm2quat(Roo{3});

qtt = spline_Quaternion([p;a;b;q],[30 30  30])
