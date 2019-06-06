function cS =Scale2(Ci,Cf,cconS,conS)
        %conS(2:end)= conS(2:end)-1;%xxxxx
        cconS (1) = 1;
        sx = zeros(numel(cconS),3);
        for ii = 1: numel(cconS)
                sx(ii,:) = deltaS(Ci{ii},Cf{cconS(ii)}) ;
                %sx(ii,:) = deltaS(Ci{ii},Cf{ii}) ;%xxxxxxxx
        end
        t = linspace(0,1,numel(cconS));
        dt = linspace(0,1,sum(conS));
        cS = zeros(sum(conS),3);
        
         cS(:,1) = spline(t,sx(:,1),dt);
         cS(:,2) = spline(t,sx(:,2),dt);
         cS(:,3) = spline(t,sx(:,3),dt);
end
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