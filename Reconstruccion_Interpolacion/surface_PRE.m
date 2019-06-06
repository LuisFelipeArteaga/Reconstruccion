 clear; clc
%%

t = linspace(0,2*pi,200);

Xc1  = [1*sin(t(1:end-1)); 1.7*cos(t(1:end-1))]';
%Xc2 = .2*[3*cos(t)+cos(3*t); 3*sin(t)-sin(3*t)]'; Xc2 = Xc2(1:end-1,:);
Xc2 = [2.3*cos(t(1:end-1)); 1.3*sin(t(1:end-1))]';
Xc3 = [1.2*sin(t(1:end-1));1.9*cos(t(1:end-1))]';
%%
dx1 = gradient (Xc1(:,1)');
dy1 = gradient (Xc1(:,1)');

dx2 = gradient (Xc2(:,1)');
dy2 = gradient (Xc2(:,2)');

dx3 = gradient (Xc3(:,1)');
dy3 = gradient (Xc3(:,2)');

dXY1 = cross([dx1; dy1;zeros(1,numel(dx1))],repmat([0 0 1]',1,numel(dx1))); dXY1 = dXY1(1:2,:);
dXY2 = cross([dx2; dy2;zeros(1,numel(dx2))],repmat([0 0 1]',1,numel(dx2)));  dXY2 = dXY2(1:2,:);
dXY3 = cross([dx3; dy3;zeros(1,numel(dx3))],repmat([0 0 1]',1,numel(dx3)));  dXY3 = dXY3(1:2,:);

%XX{1} =  [Xc1,dXY1']; XX{2} =  [Xc2,dXY2']; XX{3} =  [Xc3,dXY3'];
XX{1} =  [Xc1]; XX{2} =  [Xc2]; XX{3} =  [Xc3];
%%
maxIter = 100;
for ii= 1:2
    X = XX{ii+1};
    Y = XX{ii};
    
    %%
    W0 = ones(numel(Xc2)/2,2);

   
    figure(1),
    
    gauss = @(X,Y,s)(exp(-pdist2(X,Y).^2/(2*s^2)));
    
    V = @(X,Y,s)mean(mean(gauss(X,Y,s)));
    
    s = median(pdist(Y));
    
    N = size(X,1);
    M = size(Y,1);
    
    Xi = X;
    
    b =50;
    
    cmap = parula(N);
    
    for i = 1:maxIter
        
        Vxy = V(Xi,Y,s);
        Vx = V(Xi,Xi,s);
        c = Vxy*M/Vx/N;
        
        Kxy = gauss(Xi,Y,s);
        Kx = gauss(Xi,Xi,s);
        
        X1 = c*(1-b)/b*Kx*Xi;
        X2 = Kxy*Y;
        X3 = bsxfun(@times,c*(1-b)/b*sum(Kx,2),Xi);
        Xi = bsxfun(@plus,X1,X2-X3);
        Xi = bsxfun(@times,Xi,1./sum(Kxy,2));
        
%         i
%         clf
%         hold on;
%         scatter(X(:,1),X(:,2),[],cmap,'filled');
%         scatter(Y(:,1),Y(:,2),'filled');
%         scatter(Xi(:,1),Xi(:,2),[],cmap,'filled');
%         
%         for j=1:N
%             plot([X(j,1) Xi(j,1)],[X(j,2) Xi(j,2)],'Color',cmap(j,:))
%         end
%         
%         hold off;
%         
%         pause(0.01);
    end
    Xii {ii}= Xi;
    for jj = 1:numel(Xi)/2
        [~,posM] = min(pdist2(Xi(jj,1:2),Y(:,1:2)));
        corresXY(1,jj) = posM;
        %line([X(jj,1) Y(posM,1)],[X(jj,2) Y(posM,2)])
    end
    Corres{ii} = corresXY;
    XX{ii+1} = XX{ii+1}(Corres{ii},:);
end

%%
figure(3)
numC = 100;
cont = cell(numC,1);
tt = linspace(0,1,numC);
hold on
for pp = 1:199
    pp
    splineP = interparc(tt,[XX{1}(pp,1) XX{2}(pp,1) XX{3}(pp,1) ],...
       [XX{1}(pp,2) XX{2}(pp,2) XX{3}(pp,2)], 4*[tt(1) tt(numC/2) tt(end)]);
    
    cont{pp,1} = splineP;
    plot3(splineP(:,1),splineP(:,2),splineP(:,3),'.k')
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