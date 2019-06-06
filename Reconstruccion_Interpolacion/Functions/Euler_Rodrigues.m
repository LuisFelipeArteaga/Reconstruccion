function R = Euler_Rodrigues(n1,n2)
    theta = acos((n1'*n2)/(norm(n1)*norm(n2)));
    k =  cross(n1,n2)/norm(cross(n1,n2));
    
    a = cos(theta/2);
    b = k(1)*sin(theta/2);
    c = k(2)*sin(theta/2);
    d = k(3)*sin(theta/2);
    
    R = [(a^2+b^2-c^2-d^2) 2*(b*c-a*d)  2*(b*d+a*c);...
        2*(b*c+a*d) (a^2+c^2-b^2-d^2) 2*(c*d -a*b);...
        2*(b*d -a*c) 2*(c*d+a*b) (a^2+d^2-b^2-c^2)];
end