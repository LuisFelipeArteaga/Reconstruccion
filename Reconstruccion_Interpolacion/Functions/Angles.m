function [fi,theta,psi]  = Angles(R)

    if R(3,1) < +1  &&  R(3,1) > -1
        theta = asin(-R(3,1));
        psi = atan2(R(3,2)/cos(theta),R(3,3)/cos(theta));
        %psi = atan2(R(3,2),R(3,3));
        fi = atan2(R(2,1)/cos(theta),R(1,1)/cos(theta));
        %fi = atan2(R(2,1),R(1,1));
        disp(['newA = ' num2str(rad2deg([theta,psi,fi]))])

        if fi > pi/2 || fi < -pi/2
            theta = -theta;
            fi = pi-abs(fi);
        end
        if psi >= pi/2 || psi <= -pi/2
            psi = pi-psi;
        end
        if fi >= pi/4 && psi<=0
            psi = -psi;
        end
    else
        fi = 0;
        if R(3,1) == -1
            theta = pi/2;
            psi = atan2(R(1,2),R(1,3));
        else
            theta = -pi/2;
            psi = atan2(-R(1,2),-R(1,3));
        end
    end
    

%     if R(3,1) < 1
%         if R(3,1) > -1
%             theta =  asin(-R(3,1));
%             fi = atan2(R(2,1),R(1,1));
%             psi = atan2(R(3,2),R(3,3));
%         else
%             theta =  +pi/2;
%             fi = -atan2(-R(2,3),R(2,2));
%             psi = 0;
%         end
%     else
%         theta = -pi/2;
%        	fi = atan2(-R(2,3),R(2,2));
%         psi = 0;
% 
%     end


end