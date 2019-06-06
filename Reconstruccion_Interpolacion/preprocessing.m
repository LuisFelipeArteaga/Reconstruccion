function [e,a,pendenza,n]=preprocessing(coord)
    %pendenza mi dice come deve essere messo il primo edge, di quanto si scosta dall'asse x
    
    n=size(coord,1);
    % a: angoli con segno
    a=zeros(n,1);
    e=zeros(n,1);
    
    v=coord(2,:)-coord(1,:);
    v_norm=norm(v);
    asse=[1;0];
    phi=acos((v*asse)/v_norm);
    matrix=[coord(2,:), 1; coord(1,:), 1; coord(2,1) coord(1,2) 1];
    if det(matrix) > 0
        segno=1;
    else
        segno=-1;
    end
    pendenza=segno*phi;
    
    for i=1:1:n
        if i==n
            e(i)=sqrt((coord(i,1)-coord(1,1))^2 + (coord(i,2)-coord(1,2))^2);
        else
            e(i)=sqrt((coord(i,1)-coord(i+1,1))^2 + (coord(i,2)-coord(i+1,2))^2);
        end
    end
    
    %angoli con segno;
    c=coord(1,:)-coord(n,:);
    b=coord(2,:)-coord(1,:);
    val_ass=acos((c*b')/(norm(c)*norm(b)));
    matrix=[coord(n,:), 1; coord(1,:), 1; coord(2,:), 1];
    if det(matrix) > 0
        segno=1;
    else 
        segno=-1;
    end
    a(1)=segno*val_ass;
    
    for i=2:1:n-1
        c=coord(i,:)-coord(i-1,:);
        b=coord(i+1,:)-coord(i,:);
        val_ass=acos((c*b')/(norm(c)*norm(b)));
        matrix=[coord(i-1,:), 1; coord(i,:), 1; coord(i+1,:), 1];
        if det(matrix) > 0
            segno=1;
        else 
            segno=-1;
        end
        a(i)=segno*val_ass;
    end
    
    c=coord(n,:)-coord(n-1,:);
    b=coord(1,:)-coord(n,:);
    val_ass=acos((c*b')/(norm(c)*norm(b)));
    matrix=[coord(n-1,:), 1; coord(n,:), 1; coord(1,:), 1];
    if det(matrix) > 0
        segno=1;
    else 
        segno=-1;
    end
    a(n)=segno*val_ass;
