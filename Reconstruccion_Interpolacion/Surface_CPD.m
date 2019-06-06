clear; clc
%%

t = linspace(0,2*pi,100);
t2 = linspace(0,1,100);

Xc1  = [1*sin(t(1:end-1)); 1*cos(t(1:end-1))]';
pointsC2 =[-0.5716    0.1771;-0.3528   -0.0458;0.1214   -0.1742;0.1795   -0.2175;0.3009   -0.5162;0.6160   -0.3615;0.2928   -0.0554;0.2527    0.0020;0.0996    0.5063;-0.2879    0.5342;-0.5716    0.1771];
Xc2 = interparc(t2,pointsC2(:,1),pointsC2(:,2)); Xc2= Xc2(1:end-1,:);
Xc3 = [1*cos(t(1:end-1));1*sin(t(1:end-1))]';

Xc = {Xc1,Xc2,Xc3};
%% Parametros
addpath('mex/')
opt.method = 'nonrigid';
opt.normalize = 0;
opt.max_it = 100;
opt.tol = 1e-5;
opt.outliers = .01;
opt.sigma2 = 1;
opt.corresp = 1;
% reflejan la cantidad de regularización de la suavidad.
opt.beta = .1;
opt.lambda = 20;

%%
for  pp = 1:numel(Xc)-1
    X = Xc{pp};
    Y = Xc{pp+1};
    [M,D]=size(Y); [N, D2]=size(X);
    
    [Transform, C]=cpd_register(Y, X, opt);



    figure(pp+2);clf; hold on; axis equal
    T = Transform.Y;
    scatter(X(:,1),X(:,2),10)
    scatter(Y(:,1),Y(:,2),10)

    for ii = 1:numel(C)
        %plot([X(ii,1) T(ii,1)],[X(ii,2) T(ii,2)],'k')
        plot([X(ii,1) Y(C(ii),1)],[X(ii,2) Y(C(ii),2)],'r')
  
    end
    Cc(:,pp)= C;
end

Cc = [(1:1:size(Cc,1))',Cc];

%%
[~,bb]=sort(Cc(:,2));
Cc(:,3) = Cc(bb,3);
%%

tt = linspace(0,1,50);

figure(1);clf; axis equal; hold on
for pp = 1:size(Xc{1},1)
    disp(pp)
    px = zeros(1,numel(Xc));
    py = zeros(1,numel(Xc));
    for qq = 1:numel(Xc)
        px(qq) = Xc{qq}(Cc(pp,qq),1);
        py(qq) = Xc{qq}(Cc(pp,qq),2);
    end
    %pz = [tt(1)   tt(end)];
    pz = 1.5*[tt(1)  tt(round(numel(tt)/2))  tt(end)];
    splineP =  interparc(tt,px,py,pz);
    
    %scatter(px,py,'*b');
    scatter3(px,py,pz,20);
     scatter3(splineP(:,1),splineP(:,2),splineP(:,3),'.k')
    view([0 0 1])
    drawnow
    
    cont{pp,1} = splineP;
end

%%

clearvars -except cont
point_cloud  = vertcat(cont{:});
t = MyCrustOpen(point_cloud);
figure(4);
clf
hold on; title('Output Triangulation','fontsize',14); axis equal
% trisurf(t,point_cloud(:,1),point_cloud(:,2),point_cloud(:,3),'facecolor','c','edgecolor','b')
p = patch('Faces',t,'Vertices',point_cloud);
p.FaceColor = 'red';
p.EdgeColor = 'none';

daspect([1 1 1])
view(3);
axis vis3d
camlight