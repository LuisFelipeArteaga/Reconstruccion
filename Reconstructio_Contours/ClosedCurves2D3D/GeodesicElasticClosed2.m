function [ftc,d,Geod,R,qOut] = GeodesicElasticClosed2(p1,p2,stp,qA,pp)
    if pp == 1
        q1 = curve_to_q(p1);
    else
        q1 = qA;
    end
    q2=curve_to_q(p2);

    [q2n,R] = Find_Rotation_and_Seed_unique(q1,q2,1);
    q2n = q2n/sqrt(InnerProd_Q(q2n,q2n));
    
    q2 = ProjectC(q2);
    p2=q_to_curve(q2);
    p2n=R*p2;
    
    qOut = q2n;
    d = acos(InnerProd_Q(q1,q2n));
    alpha = geodesic_sphere_Full(q1,q2n,stp);
    ftc = Path_Plot(alpha,p2n);
 
    [k] = size(alpha,3);

    for j=1:k
        Geod(:,:,j)=q_to_curve(alpha(:,:,j));
    end
end

