function [conF,sx] =Scale2(Ci,Cf,posCon)

sx = cell(numel(Ci),1);
for ii = 1: numel(Cf)
    if ii == 1
        sx{ii,:} = deltaS(Ci{ii},Cf{ii}{ii});
    end
    sx{ii+1,:} = deltaS(Ci{ii+1},Cf{ii}{end});
end
 vAux = cell2mat(cellfun(@(x) size(x,1),posCon,'UniformOutput',false));
 cS = hobbysplines3(sx,'bezierpoints',vAux);

conF = R_C(Cf,posCon,cS);
end
%% Escala
function dxyz = deltaS(C1,C2)
n1 = C1-mean(C1);
n2 = C2-mean(C2);
dv1 = mean(abs(n1));
dv2 = mean(abs(n2));
if numel(dv2 < 1e-10) > 0
    pos = find(dv2 < 1e-10 == 1);
    dv2(pos) = 1;
end
dxyz = dv1./dv2;
end
%% Reescala contornos
function conF = R_C(cF,posCon,tmpS)
for pp = 1:numel(cF)
    for qq = 1:size(cF{pp,1},1)
        mD = diag(tmpS{pp,1}(qq,:));
        cAux = mD*(cF{pp,1}{qq,1})';
        cAux = Translation_matrix(cAux,posCon{pp,1}(qq,:))'; % Posicionar cada contorno
        conF{pp,1}{qq,1}= cAux;
    end
end
end