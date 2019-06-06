function Cout = Translation_matrix(Cint,p)
    [m,n] = size(Cint);
    me = mean(Cint,2);
    mT = eye(m+1);
    mT(1:end-1,m+1) = p'-me;
    Cout = mT*[Cint; ones(1,n)]; 
    Cout = Cout(1:3,:);
end