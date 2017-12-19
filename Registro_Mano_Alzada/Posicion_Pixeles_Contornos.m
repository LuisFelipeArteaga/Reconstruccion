close all; clear; clc;
addpath('Funciones'); % Carga la carpeta que contiene las funciones
addpath('Imagenes'); % Carga la carpeta que contiene las Imagenes
addpath('Archivos'); % Carga la carpeta que contiene los Archivos
%%  Carga los parametros basicos de las imagenes y sensores de la base de datos 
load Parametros.mat  
frames = Parametros(3,1); %Numero de frames
width = Parametros(3,2);  
height = Parametros(3,3); 
coor = Parametros(3,4); % Numero de coornadas [x y z]

%%  Crear matrices - posicion sensor-frame localizacion-isocentro(room)
% Se calcula las posiciones de cada pixel de cada frame
% Las posiciones se almacenan en un vector [ x, y, z, azimuth , elevation,  roll ] en angulos de Euler, Se realiza 
% la transformacio de las posicones expresadas en angulos de Euler a angulos Tait-Bryan. 
% Esto pude tomar demasiado tiempo, se puede calcular solo la frontera y
% ROI de cada frame

load Posiciones_Imag.mat % Posiciones Frames obtenidas de las base de datos
matrixPosFrame = cell(frames,2); %[frame,[x,y,z]]Cell que contiene las posiciones de cada pixel de cada imagen
pixels = cell(1,3); % Contiene los pixeles de cada frame

for i = 1:frames
    [xm,ym] = meshgrid(linspace(1,height,height),linspace(1,width,width));
    pixels{1,1} = i; pixels{1,2} = xm(:);pixels{1,3} = ym(:); 
    matrixPosFrame{i,1} = i; % Frame
    cellAux  = compute_Posicion_Pixels (Parametros,Posiciones_Imag,pixels); %Posiciones Frames
    vecAux = cellAux{1,2};
    matrixPos(:,:,1) = reshape(vecAux(:,1),width,height);
    matrixPos(:,:,2) = reshape(vecAux(:,2),width,height);
    matrixPos(:,:,3) = reshape(vecAux(:,3),width,height);
    matrixPosFrame{i,2} = matrixPos;
    fprintf('Frame %d de %d \n',i,frames);
end

save 'archivos\matrixPosFrame.mat' matrixPosFrame % Almacenar Posiciones


%% Posicion de Imagenes en el room;
 %  Dibuja las posiciones de cada imagen en el plano 3D de reconstruccion
figure(1)
hold on
for p = 1:frames
    frameAux = matrixPosFrame{p,2};
    vert = [frameAux(width,1,1) frameAux(width,1,2) frameAux(width,1,3) ;frameAux(width,height,1) frameAux(width,height,2) frameAux(width,height,3);...
              frameAux(1,height,1) frameAux(1,height,2) frameAux(1,height,3);  frameAux(1,1,1)  frameAux(1,1,2)  frameAux(1,1,3)]; % Vertices de las esquinas
    fac = [1 2 3 4];
    patch('Vertices',vert,'Faces',fac,'FaceColor','red');
end
title('Posicionamineto imagenes en el espacio 3D')
hold off
%% Graficar pixeles de un solo frame
%load matrixPos.mat
figure(2)
frame = 2; % Que desea graficar
frameAux = matrixPosFrame{frame,2};
hold on
vert = [frameAux(width,1,1) frameAux(width,1,2) frameAux(width,1,3) ;frameAux(width,height,1) frameAux(width,height,2) frameAux(width,height,3);...
              frameAux(1,height,1) frameAux(1,height,2) frameAux(1,height,3);  frameAux(1,1,1)  frameAux(1,1,2)  frameAux(1,1,3)];
fac = [1 2 3 4];
patch('Vertices',vert,'Faces',fac,'EdgeColor','red','FaceColor','none','LineWidth',2);
for i = 1:10:width
       for j = 1:10:height
        plot3(frameAux(i,j,1),frameAux(i,j,2),frameAux(i,j,3),'.b');
       end
end
title('Pixeles en el espacio 3D de una sola imagen ')
hold off

%% Calcular posiciones del contorno en el room.
% Calcula las posciones de los contorno, y los ubica en el plano 3D de
% reconstruccion
load contornos_Pos.mat % Vertices  de los contornos suministrado por la base de datos,  vertices(x,y) pixeles
figure(3)
[m,n] = size(contornos_Pos);
contorno_Pos3D = cell(m,n);
auxC = contornos_Pos{1,2};
% Busca las coordenas de los pixeles de los vertices en las posiciones  anteriormente almacenados
for j = 1:m
    frameAux = matrixPosFrame{j,2};
    Pix = contornos_Pos{j,2}; %
    v = zeros(numel(Pix)/2,3);
    for i = 1:numel(Pix)/2
        v(i,1) = frameAux(Pix(i,1),Pix(i,2),1);
        v(i,2) = frameAux(Pix(i,1),Pix(i,2),2);
        v(i,3) = frameAux(Pix(i,1),Pix(i,2),3);
    end
    auxD = contornos_Pos{j,1};
    contorno_Pos3D{j,1} = auxD;
    contorno_Pos3D{j,2} = v;
    % Dibujar recuadro
    if auxC ~= auxD 
        vert = [frameAux(width,1,1) frameAux(width,1,2) frameAux(width,1,3) ;frameAux(width,height,1) frameAux(width,height,2) frameAux(width,height,3);...
              frameAux(1,height,1) frameAux(1,height,2) frameAux(1,height,3);  frameAux(1,1,1)  frameAux(1,1,2)  frameAux(1,1,3)]; % Vertices de las esquinas
        fac = [1 2 3 4];
        auxC = auxD;    
    end
    hold on
    %Dibujar contornos con los vertices
    patch('Vertices',vert,'Faces',fac,'EdgeColor','blue','FaceColor','none','LineWidth',1);
    visx = [v(:,1) ;v(1,1)]; visy = [v(:,2); v(1,2)]; visz = [v(:,3); v(1,3)];
    plot3(visx,visy,visz,'g','LineWidth',1.5);
    title('Contornos  segmentados de la imagen de ultrasonido')
    %surf(v(:,1),v(:,2),v(:,3));
    hold off
end

