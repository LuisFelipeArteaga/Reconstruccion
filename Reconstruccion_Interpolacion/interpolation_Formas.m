close all; clear,clc
addpath('Functions','Files','ClosedCurves2D3D');
%%
% 
% t = linspace(0,2*pi,100);
% 
% c1x = cos(t); 
% c1y = sin(t);
% 
% c2x = .05*cos(t); 
% c2y = .1*sin(t);
% c3x = .2*cos(t)+.5; 
% c3y = .4*sin(t)+.5;
% 
% C1 = [c1x;c1y];
% C2 = [c2x;c2y];%[[c2x(1:2:end) c3x(1:2:end) ];[c2y(1:2:end) c3y(1:2:end)]];
%%
figure(3); hold on
s = linspace(0,pi,99);
ro1 = [(2+cos(s-pi/2));-1+sin(s-pi/2)];

s = linspace(pi,pi+pi/2,50); s = s(2:end);
ro2 = [2+cos(pi/2 - s);1+sin(pi/2 -s)];

s = linspace(pi+pi/2,2*pi+pi/2,100); s =s(2:end);
ro3  = [cos(s-3*pi/2);1+sin(s-3*pi/2)];

s = linspace(2*pi+pi/2,3*pi,50); s = s(2:end);
ro4 = [-2+cos(5*pi/2-s);1+sin(5*pi/2 -s)];

s = linspace(3*pi,4*pi,100); s = s(2:end-1);
ro5 = [-2+cos(s-5*pi/2);-1+sin(s-5*pi/2)];

s = linspace(-2,2,100); s =s(2:end);
ro6 = [s;-2*ones(1,numel(s))];

ro = [ro1'; ro2'; ro3';ro4';ro5';ro6']';
%ro = ro2;



s = linspace(0,pi,100);
ri1 = [1+cos(-s-pi/2);sin(-s-pi/2)];

s = linspace(pi,pi+pi/2,49); s = s(2:end);
ri2 = [1+.5*cos(2*s-5*pi/2);3/2+.5*sin(2*s-5*pi/2)];

s = linspace(1,-1,50); s = s(2:end-1);
ri3 = [s;2*ones(1,numel(s))];

s = linspace(3*pi/2+2,7*pi/2+2,200); s = s(1:end);
ri4 = [-1+2*cos(.5*s-1-pi/4);2*sin(.5*s-1-pi/4)];

s = linspace(-1,1,49); s = s(2:end-1);
ri5= [s;-2*ones(1,numel(s))];

s = linspace(7*pi/2+4,4*pi+4,50);s(2:end);
ri6 = [1+.5*cos(2*s-8-15*pi/2);-3/2+.5*sin(2*s-8-15*pi/2)];

%ri = [ri1'; ri2';ri3'; ri4';ri5'; ri6']';
ri = [ri6';ri1'; ri2';ri3'; ri4';ri5']';

C1 = ro(:,1:1:end);
C2 = ri(:,1:1:end);

%%
[ftc,~,~,Ro] = GeodesicElasticClosed2(C1,C2,5);

%%
clf
figure(1); hold on, axis equal
for  ii = 1:6
    aa = (ftc{ii}*Ro);
     plot(aa(:,1),aa(:,2),'.')
end
