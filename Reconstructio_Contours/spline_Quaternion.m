function [qOut]= spline_Quaternion(q,int)
% q Nx4, N numero de quaterniones
% int N-1x1 numero de interpolaciones entre cuaterniones N-1

 figure(2)
[m,~] = size(q);
qOut = cell(m-1,1);
% Quaternion en las misma direccion
a = 0;
for  jj = 2:m
    
    for pp = 2:m
        d = dot(q(pp-1,:),q(pp,:));
        if (d < 0)
            q(pp,:) = -q(pp,:);
        end
    end
    qi = [];
    tt = linspace(0,1,int(jj-1));
    for qq =1:numel(tt)
        t = tt(qq);
        % Angulo entre quaternions
        if  t > 0
            C = dot(q(jj-1,:),q(jj,:));
            if ((1 - C) <= eps)  || a==0%quaternions paralelos Lerp
                qi(qq,:) = q(jj-1,:)*(1-t) + q(jj,:)*t;
            else
                % Calcular puntos intermedios a and b
                qa = Estimate_intermediate_Q(jj-1,q);
                qb = Estimate_intermediate_Q(jj,q);
                qpq = Slerp(q(jj-1,:),q(jj,:),t);
                qab = Slerp(qa,qb,t);
                qi(qq,:) = quatnormalize( Slerp(qpq,qab,2*t*(1-t)));
               
                 hold on
                 [yaw, pitch, roll] = quat2angle(qi(qq,:) );
                plot3(yaw, pitch, roll,'.')
                
            end
        end
    end
    
   qi(1,:) = q(jj-1,:);
   qOut{jj-1,1} = qi;
end
%
hold on
aa =q;
for ii = 1:numel(int)
    [yaw, pitch, roll] = quat2angle(aa(ii,:));
    plot3(yaw, pitch, roll,'*')
end
end

%% Estimar quaterniones between qn qn+1
function qa = Estimate_intermediate_Q(jj,q)
[m,~] = size(q);
if m== jj
    qa = q(jj,:);
    return;
elseif jj == 1
    qa = q(jj,:);
    return;
else
    qnI = quatinv(q(jj,:));
    qL1 = quatlog(quatmultiply(qnI,q(jj+1,:)));
    qL2 = quatlog(quatmultiply(qnI,q(jj-1,:)));
    qE = quatexp(-(qL1 +qL2)/4);
    qa = quatmultiply(q(jj,:),qE);
end
end
%% SLERP

function [ q3 ] = Slerp( q1, q2, t )
%SLERP quaternion slerp
q1 = q1 ./ norm(q1);
q2 = q2 ./ norm(q2);
one = 1.0 - eps;
d = dot(q1,q2);
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

%%

