function R = Rotation_matrix_N(n1,n2)
    theta = acos((n1'*n2)/(norm(n1)*norm(n2)));
    if theta > pi/2
        theta = theta+pi;
    end
    k =  cross(n1,n2)/norm(cross(n1,n2));
    K = [ 0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
    
    R = eye(3) +(sin(theta)*K) +((1-cos(theta))*(K.^2));  
end
