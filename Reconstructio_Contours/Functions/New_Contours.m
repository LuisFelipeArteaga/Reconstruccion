function [conF,qq,q2] = New_Contours(cRS,div,Ro,C1,C2)
    C1mR = pca(C1);
    C2mR = pca(C2);
    t = linspace(0,1,numel(cRS));
    q1 = [1 0 0 0];
    q2 = Matrix_to_quaternion(Ro');

    v1 = C1mR(:,3);
    v2 = C2mR(:,3);
    qq = [];
    conF = cell(numel(cRS),1);
    for jj = 1:numel(cRS)
        q = slerp( q1', q2', t(jj));
        qq = [qq;q'];
        pos = div(jj,:);
        Rq = Rotation_quaternion(q);
        Cnew = Rq*cRS{jj}';
        vecN = (1-t(jj))*v1+ t(jj)*v2;
        vi = pca(Cnew');
        Ro = Rotation_matrix_2(vi(:,3),vecN);
        Cnew = Ro*Cnew;
        Cnew = Translation_matrix(Cnew,pos)';
        conF{jj} = Cnew(1:end,:);
    end
end