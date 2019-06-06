
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
