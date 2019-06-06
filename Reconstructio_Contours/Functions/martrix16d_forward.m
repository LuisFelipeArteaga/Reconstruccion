
function  m = martrix16d_forward (x , y, z , a ,e ,r)
% Crea una translacion y roacion de Euler angles 
% L a notacion de los angulos Tait-Bryan in ZYX, es:
% alpha = azimuth = yaw = a
% beta = elevation = pitch = e
%  gamma = roll = r
    m = zeros(4,4);
    ra = (a*pi)/180;
    rb = (e*pi)/180;
    rg = (r*pi)/180;
    m(1,1) = cos(ra)*cos(rb); 
    m(1,2) = sin(ra)*cos(rb); 
    m(1,3)=-sin(rb); 
    m(1,4) = 0;
    m(2,1) = cos(ra)*sin(rb)*sin(rg)-sin(ra)*cos(rg);
    m(2,2) = sin(ra)*sin(rb)*sin(rg)+cos(ra)*cos(rg);
    m(2,3)= cos(rb)*sin(rg);
    m(2,4) = 0; 
    m(3,1) = cos(ra)*sin(rb)*cos(rg)+sin(ra)*sin(rg);
    m(3,2) = sin(ra)*sin(rb)*cos(rg)-cos(ra)*sin(rg);
    m(3,3) = cos(rb)*cos(rg); 
    m(3,4) = 0;
    m(4,1) = x; 
    m(4,2)=y; 
    m(4,3)=z;
    m(4,4)=1;
end
