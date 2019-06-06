classdef Quaternion
    
    methods (Static)
        %%  Suma entre quaternion q1 q2
        function q = Add(q1,q2)
            q = q1+q2;
        end
        
        %% Resta entre quaternion q1 q2
        function q = Subt(q1,q2)
            q = q1-q2;
        end
        
        %% Multipliacion entre quaternion q1 q2
        function q = Mult(q1,q2)
            q1 = transp(q1);
            q2 = transp(q2);
            
            s1 = q1(1);
            s2 = q2(1);
            v1 = q1(2:4);
            v2 = q2(2:4);
            
            s =s1*s2 - dot( v1,v2);
            v = s1*v2 + s2*v1 + cross( v1, v2 );
            v = reshape( v, 1, 3 );
            q = [s v];
        end
        %% Cojugada de un quaternion
        
        function q=Conj(q1)
            q = [q1(1) -q1(2) -q1(3) -q1(4)];
        end
        
        %% Norma de un quaternion
        function q = Norm( q1)
            q = sqrt(sum(q1.^2));
        end
        
        %% inverse de un quaternion
        
        function  q = Inv(q1)
            q = Quaternion.Conj(q1)./(Quaternion.Norm(q1)^2);
        end
        
        %% quaternion unitario
        function q = Unit(q1)
            q = q1./Quaternion.Norm(q1)+eps;
        end
        
        %% angulo quaternion
        
        function theta = Angle(q)
            theta = acos(q(1)/Quaternion.Norm(q));
        end
        
        %% log quaternion
        function q = Log(q2)
            q2 = transp(q2);
            theta = Quaternion.Angle(q2);
            q1 = Quaternion.Unit(q2);
            v = q1(2:4);
            q = [0 theta*v];
            q = q./norm(q1(2:4));
        end
        
        %% exp quaternion
        
        function q = Exp(q1)
            th = norm(q1(2:4));
            if th == 0
                q = [1 0 0 0];
            else
                q = exp(q1(1))*[cos(th) sin(th)*q1(2:4)/th];
            end
        end
        %% power quaternion
        function q = Power(q1,n)
            q = Quaternion.Exp(n*Quaternion.Log(q1));
        end
    end
end
function q = transp(q1)
[m,n] = size(q1);
if m>n
    q = q1';
else
    q = q1;
end
end