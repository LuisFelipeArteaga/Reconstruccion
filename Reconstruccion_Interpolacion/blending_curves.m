%% Curvature-based blending of closed planar curves

close all; clear; clc;
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
ri = [ri1'; ri2';ri3'; ri4';ri5';ri6']';

p1 = ro(:,1:1:end);
p2 = ri(:,1:1:end);
plot(p1(1,:),p1(2,:),'.',p2(1,:),p2(2,:),'.')
plot(p1(1,1),p1(2,1),'+',p2(1,1),p2(2,1),'+')
%plot(ro5(1,:),ro5(2,:),'+',ro6(1,:),ro6(2,:),'.')
%% Parametric curves
% t = linspace(0,2*pi,1000);
% c1x = cos(t); p1x = c1x(1:2:end);
% c1y = sin(t);  p1y = c1y(1:2:end);
%
% c2x = .5*cos(t);
% c2y = sin(t);
%
% p1 =[[p1x p1x(1) ];[p1y p1y(1) ] ];
%
% % Polinomio
% ap2 = [c2x(1:end-1);c2y(1:end-1)];
% cx = sum(sqrt(sum((ap2 - [ap2(:,2:end) ap2(:,1)]).^2)));
% ui = cx./numel(p1x);
% [~,pos] = unique(ceil(cumsum(sqrt(sum((ap2 - [ap2(:,2:end) ap2(:,1)]).^2)))/ui));
% p2x = c2x(pos);
% p2y = c2y(pos);
% p2 =[[p2x p2x(1) ];[p2y p2y(1) ] ];
%%


% figure(1); axis equal
% plot(c1x,c1y,p1x,p1y,'+r')
%  hold on; plot(c2x,c2y,p2x,p2y,'+g')

%% lenght curve
L1 = sqrt(sum((p1 - [p1(:,2:end) p1(:,1)]).^2));
L2 = sqrt(sum((p2 - [p2(:,2:end) p2(:,1)]).^2));

%%
for ii = 1:numel(p1)/2
    sn1(ii) = sum(L1(1:ii))/ sum(L1) ;
    sn2(ii) = sum(L2(1:ii))/ sum(L2) ;
    t1(ii) = sum(L1(1:ii));
    t2(ii) = sum(L2(1:ii));
end
sn1 = [0,sn1];
sn2 = [0,sn2];
%%  Parametrizacion

for ii =2:numel(p1)/2-1
    ro1(:,ii-1) = p1(:,ii-1) + (t1(end)-t1(ii-1))*((p1(:,ii)-p1(:,ii-1))./(t1(ii)-t1(ii-1)));
    ro2(:,ii-1) = p2(:,ii-1) + (t2(end)-t2(ii-1))*((p2(:,ii)-p2(:,ii-1))./(t2(ii)-t2(ii-1)));
    
end
ro1 = ro1./max(ro1,[],2);
ro2 = ro2./max(ro2,[],2);

for ii =2:numel(p1)/2-1
    ri1(:,ii-1) = p1(:,ii-1) + (sn1(end)-sn1(ii-1))*((p1(:,ii)-p1(:,ii-1))./(sn1(ii)-sn1(ii-1)));
    ri2(:,ii-1) = p2(:,ii-1) + (sn2(end)-sn2(ii-1))*((p2(:,ii)-p2(:,ii-1))./(sn2(ii)-sn2(ii-1)));
end

alpha = geodesic_sphere_Full(ro1,ro2,3);
clf
figure(1); hold on, axis equal
for  ii = 1:4
    aa = alpha(:,:,ii);
     hold on, plot(aa(1,:),aa(2,:),'.')
end

%% curvature

for ii =2:numel(p1)/2-1
    pin = (p1(:,ii) - p1(:,ii-1));
    po = (p1(:,ii+1) - p1(:,ii));
    alpha1(ii-1) = (atan2(pin(2),pin(1)) -atan2(po(2),po(1)));
    if alpha1(ii-1) > pi || alpha1(ii-1) < -pi
        alpha1(ii-1) = alpha1(ii-1) -2*pi;
    end
    k1(ii-1) = 2*alpha1(ii-1)/(norm(pin)+norm(po));
end


for ii = 2:numel(p2)/2-1
    pin = (p2(:,ii) - p2(:,ii-1));
    po = (p2(:,ii+1) - p2(:,ii));
    alpha2(ii-1) = (atan2(pin(2),pin(1)) -atan2(po(2),po(1)));
    if alpha2(ii-1) > pi || alpha2(ii-1) < -pi
        alpha2(ii-1) = alpha2(ii-1) -2*pi;
    end
    k2(ii-1) = -2*alpha2(ii-1)/(norm(pin)+norm(po));
end
%% Interpolation curvature

tt = linspace(0,1,6)

figure(2);hold on, ylim([-3 -1.5 ])
for ii = 1:numel(tt)
    kt = k1.*(1-tt(ii)) +k2.*tt(ii);
    gamma = alpha1*(1-tt(ii)) + tt(ii)*alpha2;
    theta1 = cumsum(gamma);
    pp1(:,1) = p1(:,1);
    bc(:,1) = p1(:,1);
    
    for jj = 2:numel(p1)/2-4
        pp1(:,jj) = pp1(:,jj-1) + (sn2(jj) -sn2(jj-1))*[cos(theta1(jj)); sin(theta1(jj))];
        gg =jj +1;
        pp2(:,gg-2) = pp1(:,gg-1) + (sn1(gg) -sn1(gg-1))*[cos(theta1(gg)); sin(theta1(gg))];
        bc(:,jj) = pp1(:,jj) - pp2(:,jj-1)+bc(:,jj-1);
        %pp1(:,jj) = (2/(sn1(jj) -sn1(jj-1)))*(alpha1(jj-1)-alpha2(jj-1));
    end
    pp1 = [1 0 .5*ii; 0 1 0 ; 0 0 1]*[pp1;ones(1,numel(pp1(1,:)))];
    bc = [1 0 .5*ii; 0 1 0 ; 0 0 1]*[bc;ones(1,numel(bc(1,:)))];
    
    
    
    plot(pp1(1,:),pp1(2,:))
    %plot(bc(1,:),bc(2,:))
    clear pp1 bc
    %plot(kt,'-')
end

axis equal

% parametrizacion
ii = 2;

%%

clf
figure(1); hold on, axis equal
for  ii = 1:8
    aa = alpha(:,:,ii);
     hold on, plot(aa(1,:),aa(2,:),'.')
end


%% theta

%
% ii = 2;
% pp2(:,1) = p2(:,1)
% for ii = 2:numel(p1)/2-1
%     pp2(:,ii) = pp2(:,ii-1) + (sn2(ii) -sn2(ii-1))*[cos(theta2(ii-1)); sin(theta2(ii-1))];
% end
%
% %%
% figure(4);hold on
% tt = linspace(0,1,4)
% for ii = 1:numel(tt)
%     rt = pp1.*(1-tt(ii)) + pp2.*tt(ii);
%
%     plot(rt(1,:),rt(2,:),'.')
% end
%%












