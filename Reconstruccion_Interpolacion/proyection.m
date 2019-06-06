clear; close all; clc
%%
t = linspace(0,2*pi,100);
t2 = linspace(0,1,100);
c1x = 1.2*cos(t(1:end-1)); c1y = .9*sin(t(1:end-1));
c1 = [c1x;c1y]';
pointsC2 =[-0.5716    0.1771;-0.3528   -0.0458;0.1214   -0.1742;0.1795   -0.2175;0.3009   -0.5162;0.6160   -0.3615;0.2928   -0.0554;0.2527    0.0020;0.0996    0.5063;-0.2879    0.5342;-0.5716    0.1771];
c2 = interparc(t2,pointsC2(:,1),pointsC2(:,2));


%%
figure(2); hold on
plot(c1x,c1y,'.b')
plot(c2(:,1),c2(:,2),'.b')
plot(pointsC2(:,1),pointsC2(:,2),'+k')

% [xy,distance,t_a]= distance2curve(c2,c1,'linear');
% line([c1(:,1),xy(:,1)]',[c1(:,2),xy(:,2)]','color',[0 0 1])

%%  normal curve

c2 = c2(1:end-1,:)';
dx = gradient (c1x);
dy = gradient (c1y);

dx2 = gradient (c2(1,:));
dy2 = gradient (c2(2,:));
%quiver(c1x,c1y,-dx,dy)
dXY = cross([dx; dy;zeros(1,numel(dx))],repmat([0 0 1]',1,numel(dx)));
dXY2 = cross([dx2; dy2;zeros(1,numel(dx2))],repmat([0 0 1]',1,numel(dx2)));
dXY = dXY(1:2,:);
dXY2 = dXY2(1:2,:);
quiver(c1x,c1y,-dXY(1,:),-dXY(2,:))
quiver(c2(1,:),c2(2,:),dXY2(1,:),dXY2(2,:))

t= linspace(0,20,100);
xii = c1x(1) + t*(-dXY(1,1));
yii = c1y(1) + t*(-dXY(2,1));
aa = 7;
plot(pointsC2(aa,1), pointsC2(aa,2),'o')

for pp = 1:numel(pointsC2)/2
    dist = zeros(size(c1x));
    for ii = 1:numel(dXY)/2
        m = dXY(2,ii)/dXY(1,ii);
        d = c1y(ii) -m*c1x(ii);
        dist(ii) = abs((m*pointsC2(pp,1))-pointsC2(pp,2)+d)/sqrt(m.^2 +1);
    end
    [~,posD ] = min(dist)
    
%     xii = c1x(posD) - t*dXY(1,posD);
%     yii = c1y(posD) - t*dXY(2,posD);
%     plot(xii,yii,'-')
%     
end

for pp = 1:numel(pointsC2)/2
    dist = zeros(size(c2(1,:)));
    for ii = 1:numel(dXY2)/2
        m = dXY2(2,ii)/dXY2(1,ii);
        d = c2(1,ii) -m*c2(2,ii);
        dist(ii) = abs((m*pointsC2(pp,1))-pointsC2(pp,2)+d)/sqrt(m.^2 +1);
    end
    [~,posD ] = min(dist)
    
    xii = c2(1,posD) - t*dXY(1,posD);
    yii = c2(2,posD) - t*dXY(2,posD);
    plot(xii,yii,'-')
    
end
%%
% clf
% t = linspace(0,2*pi,100)';
% curvexy = [cos(t) - 1,3*sin(t) + cos(t) - 1.25];
%
% s = linspace(0,2*pi,100)';
% mapxy = 5*[cos(s),sin(s)];
%
% [xy,distance,t_a]= distance2curve(curvexy,mapxy);
% hold on
% plot(curvexy(:,1),curvexy(:,2),'k*')
% plot(mapxy(:,1),mapxy(:,2),'k.')
%
% for ii = 1:numel(curvexy)/2
%     posD = find(sum(ismember(xy,curvexy(ii,:)),2) == 2);
%     [~,pp]=min(distance(posD));
%     if ~isempty(posD)
%         line([mapxy(posD(pp(1)),1),xy(posD(pp(1)),1)]',...
%             [mapxy(posD(pp(1)),2),xy(posD(pp(1)),2)]','color',[0 0 1])
%     end
% end
% axis equal
% axis square



%%
%
% figure(1); clf;hold on, axis equal
% p1 = [1,1]
% p2 = [3,1]
% p3 = [1.22,2]
% plot([p1(1) p2(1) p3(1)],[p1(2) p2(2) p3(2)],'*k')
% plot([p1(1) p2(1)],[p1(2) p2(2) ])
%
% %R = (sqrt(sum((p1-p3).^2))+sqrt(sum((p1-p2).^2)))./(sqrt(sum((p2-p3).^2))+sqrt(sum((p2-p1).^2)));
%
% R = (sqrt(sum((p1-p3).^2)))./(sqrt(sum((p1-p2).^2)));
%
% xi = p3(1) +R*(p2(1)-p1(1));
% yi = p3(2) +R*(p2(2)-p1(2));
%
% plot(xi,yi,'*r')
% plot(yi,xi,'*r')
% line([yi p3(1)],[xi  p3(2)])