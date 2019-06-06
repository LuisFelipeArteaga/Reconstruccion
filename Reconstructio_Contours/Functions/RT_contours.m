function [conF] = RT_contours(qOut,posCon,conT,vecPosCon)

% figure(1); clf
% hold on, axis equal
% plot3(posConM(:,1),posConM(:,2),posConM(:,3),'*b')
conF = cell(size(posCon));
tmpca = cellfun(@(x) pca(x),vecPosCon(:,2),'UniformOutput',false);
tmpca = cellfun(@(x) x(:,2)',tmpca,'UniformOutput',false);
conS = cell2mat(cellfun(@(x) size(x,1),posCon,'UniformOutput',false));
vecN = hobbysplines3(tmpca,'bezierpoints',conS);
for jj = 1:numel(conT)
     intAux = conT{jj};

    for ii = 1:numel(intAux)
        conRT =  mult_M(qOut,conT,vecN,posCon,jj,ii) ;
        conF{jj,1}{ii,1}  =  conRT';
        %plot3(conRT(1,:),conRT(2,:),conRT(3,:),'.r')
    end
end
end
function conRT =  mult_M(qOut,conT,vecN,posCon,jj,ii) 
mR = quat2rotm(qOut{jj}(ii,:));
cAux = conT{jj}{ii};
conRT = mR*cAux';
vi = pca(conRT');
Ro = Rotation_matrix_N(vi(:,2),vecN{jj}(ii,:)');
Ro = eye(3);
conRT  = (Ro*conRT);
conRT = Translation_matrix(conRT,posCon{jj}(ii,:));
end