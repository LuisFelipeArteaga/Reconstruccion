function cRS = Scale(ftc,C1,C2,Ro)
    C2 = Ro*C2';C2 = C2';
    t = linspace(0,1,numel(ftc));
    n1 = C1-mean(C1);
    n2 = ftc{1}-mean(ftc{1});
    dv1 = mean(abs(n1));
    dv2 = mean(abs(n2));
    if numel(dv2 < 1e-10) > 0
        pos = find(dv2 < 1e-10 == 1);
        dv2(pos) = 1;
    end
    s1 = dv1./dv2;
    diag(s1);

    n4 = C2-mean(C2);
    n5 = ftc{end}-mean(ftc{end});
    dv3 = mean(abs(n4));
    dv4 = mean(abs(n5));
    if numel(dv4 < 1e-10) > 0
        pos = find(dv4 < 1e-10 == 1);
        dv4(pos) = 1;
    end
    s2 = dv3./dv4;

    cRS = cell(1,numel(ftc));

    for ii= 1:numel(ftc)
        c1 = ftc{ii}'-mean(ftc{ii})';
        st = (1-t(ii))*s1 + t(ii)*s2;
        mST = diag(st);
        c1 = mST*c1;
        cRS{ii} = c1;
    end
     
end