clc, clear,close all
 restoredefaultpath
%% Read files
 addpath('Functions','Files','ClosedCurves2D3D','tspline');
load Cont.mat
div = 50
numCon = 2;
pointsT = c4;

pointsT = cellfun(@(x)  [x(:,1),zeros(size(x(:,1))),x(:,2)],pointsT,'UniformOutput',0);
pointsT{2} = pointsT{2} + [0 .5 0];

pointsCi = cellfun(@(x)  mean(x),pointsT,'UniformOutput',0)'; % Centro contornos iniciales
d = zeros(numCon-1,1);
mC0 = mean(c0{1})
pointsCi{3}=  [mC0(1) 1 mC0(2)]
%
% load('Posiciones_Imag.mat') % Posiciones Imagenes
% load('Parametros.mat')% Parametros Registro
% fileNerve = 'CiaticoP02.xyz'; % Contornos segmentados nervios
% load('posBoundary3D.mat') % Posiciones imagenes
% pixels = load(fileNerve);
% 
% %Selecionar contorno y posicionar sobre imagen
% div =30;  % Division contornos totales
% numCon =4;  %Numero de contornos
% imag = 1:numCon;
% posImag = round(linspace(1,135,numCon));
% [pixelsROI,pixelsROIDownS] = cell_Cont(pixels,imag,posImag);
% vecPosCon=compute_Posicion_Pixels (Parametros,Posiciones_Imag,pixelsROIDownS);
% for qq = 1:numCon
%     vecPosCon{qq,2}(:,2) = qq*.5 + vecPosCon{qq,2}(:,2);
% end

%% Smooth trayectoria

vecPosCon = cell(2,2)
vecPosCon(:,2) = pointsT
% pointsCi = cellfun(@(x)  mean(x),vecPosCon(:,2),'UniformOutput',0); % Centro contornos iniciales
d = zeros(numCon-1,1);

% Distancia trayectoria
for kk = 1:size(pointsCi,1)-1
    d(kk,1) = sqrt(sum((pointsCi{kk} -pointsCi{kk+1}).^2));
end
dDelta = sum(d)/(div);

conS = round(d./dDelta); %Numero de interpolaciones entre contorno
mAuxC = [0 ;cumsum(conS(1:end))];
conS(2:end)= conS(2:end)+1;

posCon = hobbysplines3(pointsCi,'bezierpoints',conS); % interpolacion spline trayectoria
posConM = cell2mat(posCon);

%-

%% Interpolacion Contornos
[conT,RoT,qq] = Interpolation_Contours( vecPosCon,conS,numCon);

%% Interpolar Rotacion (Spline)
RoT{1}(2,2) = RoT{1}(2,2); 
for ii = 1:numCon-1
    q1(ii,:) = Matrix_to_quaternion(-RoT{ii});
end
q = [1 0 0 0;q1];
% Spline Rotation
[qOut]= spline_Quaternion(q,conS);


%%  Rotar and Transladar Contornos

[conF] = RT_contours(qOut,posCon,conT,vecPosCon) ;


%% Reescalar Contornos

[conF,sx]= Scale2(vecPosCon(:,2),conF(1),posCon(1)); 


%% Tunning

[conFT] = Tuning(conF,vecPosCon);
conFT1 = conFT;

for ii = 1:numel(conFT)-1
    conFT{ii,1}(end) = [];
end

conFT = vertcat( conFT{:});
%% Plot contour
% close all
figure(3);
clf
hold on
axis equal

mAuxC(1) = 1;
for qq = 1:numel(conFT1)
    aux = conFT1{qq};
    for pp = 1:numel(aux)
        if pp == numel(aux) ||pp == 1  
            plot3(aux{pp}(:,1),aux{pp}(:,2),aux{pp}(:,3),'.')
        else
             plot3(aux{pp}(:,1),aux{pp}(:,2),aux{pp}(:,3),'.b')
        end
    end
end
hold on 
plot3(posConM(:,1),posConM(:,2),posConM(:,3),'.b')
title('Trayectoria nuevos Contornos')
axis off

 hold on
for  ii = 1:numel(conF) +1
      plot3(vecPosCon{ii,2}(:,1), vecPosCon{ii,2}(:,2), vecPosCon{ii,2}(:,3),'+g')
 end
%% 
figure(5), clf
axis equal
hold on
for ii = 1:numel(conFT)
    
    plot3(conFT{ii}(:,1),conFT{ii}(:,2),conFT{ii}(:,3),'.b')
end


%%
%  close all
%  
%  Co1 = conF{1}{end};
%  Co2 = conF{2}{1};
%  [~,~,Ti,~]=icp(Co1(1:end-1,:),Co2(1:end-1,:)); % vector de ajuste
%  figure(1); hold on
%  Co3 = Co2 + [Ti(1) Ti(2) 0*Ti(3)];
%   plot3(Co1(:,1),Co1(:,2),Co1(:,3),'.b');
%   plot3(Co2(:,1),Co2(:,2),Co2(:,3),'.r');
%   plot3(Co3(:,1),Co3(:,2),Co3(:,3),'.');
 
%%
% Save point Cloud
conFT = cellfun(@(x) x(1:end-1,:),conFT,'UniformOutput',false);
point_cloud = cell2mat(conFT);
save PLY_Contours/c4.mat conFT
%save('Files/point_cloud.txt', 'point_cloud', '-ascii', '-double', '-tabs')

%% Plot Surface
% close all
t = MyCrustOpen(point_cloud);
figure(66);
%hold on; title('Output Triangulation','fontsize',14); axis equal
% trisurf(t,point_cloud(:,1),point_cloud(:,2),point_cloud(:,3),'facecolor','c','edgecolor','b')
p = patch('Faces',t,'Vertices',point_cloud);
p.FaceColor = 'red';
p.EdgeColor = 'none';

daspect([1 1 1])
view(3);
axis vis3d
camlight
% 
writePLY('PLY_Contours/c4.ply',point_cloud,t, 'ascii' )
% save PLY_Contours/c3.ply,point_cloud,t, 'ascii' )
% %


%% % Interpolacion escala
% % close all
% %
% %
% % figure(3); clf
% % title('Interpolacion escala')
% % hold on
% % plot(dt,intS(:,1),'--r',tt,sx(:,1),'*r')
% % plot(dt,intS(:,2),'--b',tt,sx(:,2),'*b')
% % plot(dt,intS(:,3),'--g',tt,sx(:,3),'*g')
% %
% % figure(4); clf
% % title('Interpolacion posicion')
% % t = linspace(1,100,numel(intS(:,1)));
% % hold on
% % plot(posCon(:,1),'.r')
% % plot(posCon(:,2),'.b')
% % plot(posCon(:,3),'.g')
% %
% % figure(5); clf
% % hold on
% % title('Interpolacion posicion')
% % plot3(pointsCim(:,1),pointsCim(:,2),pointsCim(:,3),'*g')
% % plot3(posCon(:,1),posCon(:,2),posCon(:,3),'.r')
% %
% % figure(6)
% % hold on
% % title('Interpolacion rotacion')
% % plot(qqq(:,1),'r')
% % plot(qqq(:,3),'b')
