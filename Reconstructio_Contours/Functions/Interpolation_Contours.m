function [conT,RoT,qq] = Interpolation_Contours( vecPosCon,conS,numCon)
conT = cell(numCon-1,1);
RoT = cell(numCon-1,1);
Ro = eye(3,3);
q2 = [];
for pp = 1:numCon-1
    C1 = (Ro*vecPosCon{pp,2}')';
    %C1 = vecPosCon{pp,2};
    C2 = vecPosCon{pp+1,2};
    %%% Matrix de rotacion contornos
    C1mR = pca(C1);
    C2mR = pca(C2);
    
    %%% Pocisiona los contornos sobre el mismo plano
%     R1 = Rotation_matrix_2(C1mR(:,3),[ 0 1 0]');
%     R2 = Rotation_matrix_2(C2mR(:,3),[ 0 1 0]');
    R1 = eye(3);  
    R2 = eye(3);
    C1i = R1*C1';
    C2i = R2*C2';
    
    %%% Translada los contornos al origen
    C1ii = Translation_matrix(C1i,[ 0 0 0]);
    C2ii = Translation_matrix(C2i,[ 0 0 0]);
    
    %%% Interpolacion de formas
    [ftc,~,~,Ro,q2] = GeodesicElasticClosed2(C1ii,C2ii,conS(pp)-1,q2,pp);
    qq{pp} = q2;
   % [ftc,~,~,Ro] = GeodesicElasticClosed2(C1',C2',conS(pp)-1);
    
    RoT{pp,1} = Ro;
    conT{pp,1} = ftc;
end
end