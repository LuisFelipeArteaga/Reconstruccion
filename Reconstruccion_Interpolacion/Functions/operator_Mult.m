function vector3D = operator_Mult( vector3Da,matrix16d )
%Mutilpevector by matriz 
    vector3D = zeros(3,1);
    if vector3Da (3) == 0
        vector3D(1) = vector3Da(1)*matrix16d(1,1)+vector3Da(2)*matrix16d(2,1)+matrix16d(4,1);
        vector3D(2) = vector3Da(1)*matrix16d(1,2)+vector3Da(2)*matrix16d(2,2)+matrix16d(4,2);
        vector3D(3) = vector3Da(1)*matrix16d(1,3)+vector3Da(2)*matrix16d(2,3)+matrix16d(4,3);
    else
        vector3D(1) = vector3Da(1)*matrix16d(1,1)+vector3Da(2)*matrix16d(2,1)+vector3Da(3)*matrix16d(3,1)+matrix16d(4,1);
        vector3D(2) = vector3Da(1)*matrix16d(1,2)+vector3Da(2)*matrix16d(2,2)+vector3Da(3)*matrix16d(3,2)+matrix16d(4,2);
        vector3D(3) = vector3Da(1)*matrix16d(1,3)+vector3Da(2)*matrix16d(2,3)+vector3Da(3)*matrix16d(3,3)+matrix16d(4,3);
    end
end

