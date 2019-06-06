% This function calculates a geodesic on the sphere at the
% starting point x_init in the direction g
function X = geodesic_sphere_Full(q1,q2,stp)

theta = acos(InnerProd_Q(q1,q2));
if theta > 0.0001
    for t=1:stp+1
        tau = (t-1)/stp;
        X(:,:,t) = (sin((1-tau)*theta)*q1 + sin(tau*theta)*q2)/sin(theta);
        X(:,:,t) = ProjectC(X(:,:,t));  
    end
else
    for t=1:stp+1
        X(:,:,t) = q1;
    end
end
return;
