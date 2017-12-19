function [pixelsROI , pixelesBoundary]= Generar_ROI_Boundary (contornos_Pos,width,height,a ,t)
% Se genera los pixeles de la ROI con su respectiva frontera a partir de los vertices de
% los contornos de la base de datos, 
% load contornos_Pos.mat  Vertices  de los contornos suministrado por la base de datos,  vertices(x,y) pixeles
% width : ancho de la imagen
% heigth : altura de la imagen
% a : 1, si desea mostrar imagenes de las regiones y de los contornos
% t : tiempo por el que desea mostar las imagenes.
    vecCont =  cell2mat(contornos_Pos(:,1));
    imgCont = unique(vecCont);
    pixelsROI = cell(numel(imgCont),3);
    pixelesBoundary = cell(numel(imgCont),3);
    auxPos = 0;

for i = 1:numel(imgCont )
     BW = zeros(width,height);
     posCont = find(vecCont == imgCont(i));
     % Contorno del mismo frame
     for j = posCont(1):posCont(end)
        pixels = contornos_Pos{j,2};
        c = pixels(:,1);
        r = pixels(:,2);
        aux = roipoly(width,height,c,r); % Generar ROI
        auxPos = auxPos +1;
        BW = BW+aux; 
     end
     BW2 = boundarymask(BW);% Genera Boundary
     
     [xpix,ypix] = find(BW == 1);
     pixelsROI{i,1} = imgCont(i);
     pixelsROI{i,2} = xpix;
     pixelsROI{i,3} = ypix;
     
     [xpix2,ypix2] = find(BW2 == 1);
     pixelesBoundary{i,1} = imgCont(i);
     pixelesBoundary{i,2} = xpix2;
     pixelesBoundary{i,3} = ypix2;
     % Graficas las ROI y boundary
     if a == 1
         figure(1)
         subplot(1,3,1)
         imshow(BW);title('Region de Interes')
         subplot(1,3,2)
         imshow(BW2);title('Frontera')
         subplot(1,3,3)
         hold on
         imshow(imread([ num2str(imgCont(i)) '.png']))
         for j = posCont(1):posCont(end)
             pixels = contornos_Pos{j,2};
             c = pixels(:,1);
             r = pixels(:,2);
             plot(c,r,'m*')
         end
         title('Contorno')
         pause(t)
     end
end

end