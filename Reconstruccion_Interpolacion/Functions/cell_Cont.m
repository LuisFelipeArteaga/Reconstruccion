function [pixelsROI,pixelsROIDownS] = cell_Cont(pixels,imag,posImag)
    %pixels = round(pixels);
    % Almacenar pixeles en cell
    [ai,ac,~] = unique(pixels(:,1));
    pixelsROI = cell(numel(ai),3);
    ac2 =[ac(2:end,:)-1;size(pixels,1)];
    for ii = 1:numel(ac)
        pixelsROI{ii,1} = ai(ii);
        pixelsROI{ii,2} = pixels(ac(ii):ac2(ii),2);
        pixelsROI{ii,3} = pixels(ac(ii):ac2(ii),3);
    end
    % Posicionar imagnes 
%     imag = [1,2];
%     posImag = [1 135];
    pixelsROIDownS = cell(numel(imag),3);
    for ii = 1:numel(imag)
        jj = imag(ii);
        pixelsROIDownS{ii,1} = posImag(ii);
        pixelsROIDownS{ii,2} = pixelsROI{jj,2}-mean(pixelsROI{jj,2});
        pixelsROIDownS{ii,3} = pixelsROI{jj,3}-mean(pixelsROI{jj,3});
    end
end