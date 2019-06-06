function  vecPosCon  = compute_Posicion_Pixels (Parametros,Posiciones_Imag,pixelsROI)
    %Parametros: Imagenes y sensor 
    %Posiciones_Imag : [x, y, z, azimuth , elevation,  roll]
    % pixelsROI = pixeles que desea calcular las posiciones
    %Calcula la posicion de cada pixel de cada frame, en el espacio
    %3D(X,Y,Z) de reconstruccion.
    vecCal = Parametros(1,:); % Posicion del sensor 
    vecIso = Parametros(2,:); % Posicion inicio del room(Espacio de reconstruccion)
    cal = martrix16d_forward(vecCal(1),vecCal(2),vecCal(3),vecCal(4),vecCal(5),vecCal(6)); % Convertir angulos
    iso = martrix16d_forward(vecIso(1),vecIso(2),vecIso(3),vecIso(4),vecIso(5),vecIso(6)); % Convertir angulos
 
    vecPosCon = cell(size(pixelsROI,1),2);
  

    for p = 1:size(pixelsROI,1)
        %Posicionese Frame
        vecPos = Posiciones_Imag(pixelsROI{p,1},:); % Frame
        frame = martrix16d_forward(vecPos(1),vecPos(2),vecPos(3),vecPos(4),vecPos(5),vecPos(6)); %Frame
        posAuxX = pixelsROI{p,2}; % Posicion X
        posAuxY = pixelsROI{p,3}; % Posicion Y
        vectorCont = zeros(3);
        for i = 1:numel(posAuxX)
                %vectores 3D [x,y,z]
                x = posAuxX(i,1); y = posAuxY(i,1); z = 0;
                pos = operator_Mult(pixel_scala(x,y,z,Parametros(3,5),Parametros(3,6),0),cal);
                world = operator_Mult(pos,frame);
                room =  operator_Mult(world,iso);
                vectorCont(i,1) = room(1); %x
                vectorCont(i,2) = room(2); %y
                vectorCont(i,3) = room(3); %z
        end
        vecPosCon{p,1} = pixelsROI{p,1}; % Almacena el frame
        vecPosCon{p,2} = vectorCont; % Almacena el vector posicion [x y z]
    end

end