clc, clear,close all

%% Read files
addpath('Functions','Files','ClosedCurves2D3D');
load('Posiciones_Imag.mat') % Posiciones Imagenes
load('Parametros.mat')% Parametros Registro
fileNerve = 'CiaticoP02.xyz'; % Contornos segmentados nervios
load('posBoundary3D.mat') % Posiciones imagenes
pixels = load(fileNerve);

%Selecionar contorno y posicionar sobre imagen
div =40;  % Division contornos totales
numCon = 6 ;  %Numero de contornos
imag = 1:numCon;
posImag = round(linspace(1,135,numCon));
[pixelsROI,pixelsROIDownS] = cell_Cont(pixels,imag,posImag);
vecPosCon=compute_Posicion_Pixels (Parametros,Posiciones_Imag,pixelsROIDownS);
for qq = 1:numCon
    vecPosCon{qq,2}(:,2) = qq*.5 + vecPosCon{qq,2}(:,2);
end
%% Smooth trayectoria

pointsCi = cellfun(@(x)  mean(x),vecPosCon(:,2),'UniformOutput',0 ); % Centro contornos iniciales
d = zeros(numCon-1,1);

% Distancia trayectoria
for kk = 1:size(pointsCi,1)-1
    d(kk,1) = sqrt(sum((pointsCi{kk} -pointsCi{kk+1}).^2));
end
dDelta = sum(d)/(div);

conS = round(d./dDelta); %Numero entre interpolaciones entre contorno
mAuxC = [0 ;cumsum(conS(1:end))];
posCon = hobbysplines(pointsCi,'bezierpoints',conS); % interpolacion spline
posConM = cell2mat(posCon);

%% Contornos

conT = cell(numCon-1,1);
for pp = 1:1:numCon-1
   %Contornos
    C1 = vecPosCon{pp,2};
    C2 = vecPosCon{pp+1,2};
    %Matrix de rotacion contornos
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
    
    %%  Interpolacion
    
    %Interpolacion de formas
    [ftc,~,~,Ro] = GeodesicElasticClosed2(C1ii,C2ii,conS(pp)-1);
    %% Calcula nueva posicion de los contornos
    ftc= New_Contours(ftc,posCon{pp,1},Ro,C1,C2);
    
    conT{pp,1} = ftc;
end

conT= vertcat(conT{:});
conT = cellfun(@(x)  x(1:end-1,:),conT,'UniformOutput',0);

%%  Reescalar contornos interpolados

conF= cell(div,1);
intS = Scale2(vecPosCon(:,2),conT,mAuxC,conS); % vector [1x3] para escalar cada contorno

% Reescalar cada contorno
for qq = 1:div
    mD = diag(intS(qq,:));
    cAux = mD*conT{qq}';
    cAux = Translation_matrix(cAux,posConM(qq,:))'; % Posicionar cada contorno
    conF{qq,1} = cAux;
end
% Ajuste fino de trayactoria
for jj = 1:numCon-1
    Co1 = Translation_matrix(conF{mAuxC(jj+1),1}',[ 0 0 0])';
    Co2 = Translation_matrix(conF{mAuxC(jj+1)-1,1}',[ 0 0 0])';
    [~,T,~,~]=icp(Co1',Co2'); % vector de ajuste
    t = linspace(0,1,conS(jj)-1);
    for ii = mAuxC(jj)+1:mAuxC(jj+1)-1
        conF{ii} = conF{ii} +[T(1) 0 T(3)]*t(ii-mAuxC(jj)); % ajustar solo en X and Y
    end
end

%% Plot contornos Iniciales 
clf
figure(2)
hold on
axis equal
posCon1 = cell2mat(posCon );
pointsCim = cell2mat(pointsCi ) ;
for qq = 1:numCon
    plot3(vecPosCon{qq,2}(:,1),vecPosCon{qq,2}(:,2),vecPosCon{qq,2}(:,3),'.r');
   
end


hold on; title('Input Contornos','fontsize',14); axis equal

figure(2); hold on; axis equal
for pp =1:div
    cAux = conF{pp};
    plot3(cAux(:,1),cAux(:,2),cAux(:,3),'.b')
    conF{pp} = cAux;
end
%% Save point Cloud
point_cloud = cell2mat(conF);
save('Files/point_cloud.txt', 'point_cloud', '-ascii', '-double', '-tabs')

%% Plot Surface

t = MyCrustOpen(point_cloud);
figure(3)
hold on; title('Surface','fontsize',14); axis equal
%trisurf(t,point_cloud(:,1),point_cloud(:,2),point_cloud(:,3),'facecolor','r','edgecolor','b')
p = patch('Faces',t,'Vertices',point_cloud);
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1 1 1])
view(3); 
axis vis3d
camlight

