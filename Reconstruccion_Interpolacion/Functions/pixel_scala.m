function  vector3D = pixel_scala( x,y,z,xs,xy,xz )
%Escala el vector para cada componente
    x_out = x*xs;
    y_out = y*xy;
    z_out = z*xz;
    vector3D = [x_out y_out z_out];
end

