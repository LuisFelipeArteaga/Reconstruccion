function [ q3 ] = slerp( q1, q2, t )
%SLERP quaternion slerp
q1 = q1 ./ norm(q1);
q2 = q2 ./ norm(q2);

one = 1.0 - eps;
d = q1'*q2;
absD = abs(d);

if(absD >= one)
    scale0 = 1 - t;
    scale1 = t;
else 
    theta = acos(absD);
    sinTheta = sin(theta);
    
    scale0 = sin( ( 1.0 - t ) * theta) / sinTheta;
    scale1 = sin( ( t * theta) ) / sinTheta;
end
if(d < 0)
    scale1 = -scale1;
end

q3 = scale0 * q1 + scale1 * q2;
q3 = q3 ./ norm(q3);
end

