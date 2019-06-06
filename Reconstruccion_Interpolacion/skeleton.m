%% Skeleton 
clear;clc
load /home/felipe/Desktop/Imagen_medic/PointCloud/slice.mat
load /home/felipe/Desktop/Imagen_medic/PointCloud/contsC.mat

figure(10)
for ii = 1:numel(slice)
    subplot(ceil(numel(slice)/2),2,ii)
    imagesc(slice{ii});hold on; yy = cell2mat(contsC{ii}); plot(yy(:,2),yy(:,1),'.r');
    title(num2str(ii))
end

