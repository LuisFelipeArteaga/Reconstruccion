clc, clear,close all
div = 10;
%% Read files
addpath('Functions','icp','Files','ClosedCurves2D3D');
load('Posiciones_Imag.mat')
load('Parametros.mat')
fileNerve = 'CiaticoP02.xyz';
load('posBoundary3D.mat')
pixels = load(fileNerve);

%Seleccionar contorno y posicionar sobre imagen
imag = [7,7];
posImag = [1 135];
[pixelsROI,pixelsROIDownS] = cell_Cont(pixels,imag,posImag);
vecPosCon=compute_Posicion_Pixels (Parametros,Posiciones_Imag,pixelsROIDownS);
%% Contornos
mC1 = mean(vecPosCon{1,2});
mC2 = mean(vecPosCon{2,2});
C1 = (rotz(0)*rotz(0)*roty(0)*vecPosCon{1,2}')';
C2 = (rotz(0)*rotz(0)*roty(90)*vecPosCon{2,2}')';

C1 = Translation_matrix(C1',mC1);C1 = C1';
C2 = Translation_matrix(C2',mC2);C2 = C2';
 
%Plot contornos Iniciales
figure(1); clf;hold on ; axis equal
plot3(C1(:,1),C1(:,2),C1(:,3),'+b')
plot3(C2(:,1),C2(:,2),C2(:,3),'+b')
%%  Interpolacion
% Calcula matriz de rotacion PCA(normal plano)
C1mR = pca(C1);
C2mR = pca(C2);

%Pocisiona los contornos sobre el mismo plano
R1 = Rotation_matrix_2(C1mR(:,3),[ 0 1 0]');
R2 = Rotation_matrix_2(C2mR(:,3),[ 0 1 0]');
C1i = R1*C1';
C2i = R2*C2';

%Translada los contornos al origen
C1ii = Translation_matrix(C1i,[ 0 0 0]);
C2ii = Translation_matrix(C2i,[ 0 0 0]);

%Interpolacion de formas
[ftc,d,Geod,Ro] = GeodesicElasticClosed2(C1i,C2i,div-1,[],1);
view([ 0 1 0])

C3i = Ro'*C1i; 
C3ii = Translation_matrix(C3i,[ 0 0 0]);
figure(3); hold on; axis equal; view([ 0 1 0])
plot3(C1ii(1,:),C1ii(2,:),C1ii(3,:),'.');
plot3(C2ii(1,:),C2ii(2,:),C2ii(3,:),'.');
plot3(C3ii(1,:),C3ii(2,:),C3ii(3,:),'-');

%% Scalar contornos
cRS = Scale(ftc,C1i',C2i',div,Ro);
figure(12); axis equal, view([ 0 1 0])
for ii= 1:numel(cRS)
    c1 = cRS{ii}';
    hold on;plot3(c1(:,1),c1(:,2),c1(:,3),'.g');
end


%%
figure(1); hold on ; axis equal
cIn = cell(size(cRS));
normalP = cell(size(cRS));
t = linspace(0,1,div);

v1 = C1mR(:,3);
v2 = C2mR(:,3);
if abs(acos((v1'*v2)/(norm(v1)*norm(v2)))) > pi/4
    v2 = -v2;
end
for jj = 1:div
    cs = cRS{jj}; 
    vN = pca(cs');
    pos = (1-t(jj))*mC1 + t(jj)*(mC2);
    vecN = (1-t(jj))*v1+ t(jj)*v2;
    normalP{jj} = vecN;
%     cs = Translation_matrix(cs,mC1);
    Ro = Rotation_matrix_2(vN(:,3),vecN);  
    c1 = Ro*cs;
    c11 = Translation_matrix(c1,pos);
    cIn{jj} = c11;
    plot3(c11(1,:),c11(2,:),c11(3,:),'.r')     
    %quiver3(pos(1),pos(2),pos(3),vecN(1,1),vecN(2,1),vecN(3,1),'LineWidth',1,'color','g');
    
end

%% Registro Contour

figure(4);clf; hold on ; axis equal

% plot3(C2(:,1),C2(:,2),C2(:,3),'.b')
plot3(cIn{1}(1,:),cIn{1}(2,:),cIn{1}(3,:),'.r')
% plot3(C1(:,1),C1(:,2),C1(:,3),'.b')
plot3(cIn{end}(1,:),cIn{end}(2,:),cIn{end}(3,:),'.b')

[Ri1,TrV1,cOut1]=icp(C1,C2,200,100,3,1e-6);


Cif = pinv(Ri1)*C2'+ TrV1;
%Cif = pinv(Ri1)*cIn{end};
Cif = Translation_matrix(Cif,mC2);
plot3(cOut1(1,:),cOut1(2,:),cOut1(3,:),'.g')
plot3(Cif(1,:),Cif(2,:),Cif(3,:),'.m')
% plot3(cOut2(1,:),cOut2(2,:),cOut2(3,:),'.g')
% plot3(cOut3(1,:),cOut3(2,:),cOut3(3,:),'.m')


quiver3(mC1(1),mC1(2),mC1(3),C1mR(1,1),C1mR(2,1),C1mR(3,1),'LineWidth',2,'color','k');
quiver3(mC1(1),mC1(2),mC1(3),C1mR(1,2),C1mR(2,2),C1mR(3,2),'LineWidth',3,'color','k');
quiver3(mC1(1),mC1(2),mC1(3),C1mR(1,3),C1mR(2,3),C1mR(3,3),'LineWidth',4,'color','k');
quiver3(mC2(1),mC2(2),mC2(3),C2mR(1,1),C2mR(2,1),C2mR(3,1),'LineWidth',2,'color','k');
quiver3(mC2(1),mC2(2),mC2(3),C2mR(1,2),C2mR(2,2),C2mR(3,2),'LineWidth',3,'color','k');
quiver3(mC2(1),mC2(2),mC2(3),C2mR(1,3),C2mR(2,3),C2mR(3,3),'LineWidth',4,'color','k');
quiver3(mC1(1),mC1(2),mC1(3),C1mR(1,1),C1mR(2,1),C1mR(3,1),'LineWidth',2,'color','k');
quiver3(mC1(1),mC1(2),mC1(3),C1mR(1,2),C1mR(2,2),C1mR(3,2),'LineWidth',3,'color','k');
quiver3(mC1(1),mC1(2),mC1(3),C1mR(1,3),C1mR(2,3),C1mR(3,3),'LineWidth',4,'color','k');









