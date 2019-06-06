%camel;
tu;
%easy_polygon;

[l, a, teta1, n1]=preprocessing(coord1);
[m, b, teta2, n2]=preprocessing(coord2);
% a, l:  angles (curvatures) and edges of the first polygon
% b, m:  angles (curvatures) and edges of the  second polygon

len1=sum(l)
len2=sum(m)
param1=zeros(n1+1,1);
param2=zeros(n2+1,1);

% param1 and param2 are the parameter points in [0,len1] and [0,len2]
% respectively; par1 and par2 on the interval [0,1].
for i=2:1:(n1+1);
    param1(i)= param1(i-1)+l(i-1);
end
par1=param1./len1;

for i=2:1:(n2+1);
    param2(i)= param2(i-1)+m(i-1);
end
par2=param2./len2;

% the vector p contains the warping of parameters par1 and par2;
[angles1, angles2, p] =warping(par1, par2, n1, n2, a, b);

n=size(p,1)-1

step=input('insert step: ');
for t=0:step:1
    teta_start= (1-t)*teta1 + t*teta2;
    k=(1-t)*angles1+t*angles2;
    L=(1-t)*len1 + t*len2; % Length of the interpolated curve

    % teta is the vector of angles that the i-th edges creates with the x-axis
    teta=zeros(n,1); 
    teta(1)= teta_start;
    for i=2:1:n
        teta(i)=teta(i-1)+k(i);
    end

    % Open curve obtained from interpolation
    aa=zeros(n+1,1);
    bb=zeros(n+1,1);
    for i=2:1:n+1
        aa(i)=aa(i-1)+L*(p(i)-p(i-1))*cos(teta(i-1));
        bb(i)=bb(i-1)+L*(p(i)-p(i-1))*sin(teta(i-1));
    end
    plot(aa+4*round(t*10),bb,'b');
%     hold on;
%     pause;
%     plot(aa,bb,'c*');
%     hold on;
    pause;

    % closure constraints:
    Aeq=zeros(3,2*n-1);
    
    for i=1:1:(n-1)
        Aeq(1,2*i-1)= sin(teta(i));
        Aeq(1,2*i)= sin(teta(i));
        Aeq(2,2*i-1)= cos(teta(i));
        Aeq(2,2*i)=cos(teta(i));        
    end
    Aeq(1,2*n-1)=sin(teta(n));
    Aeq(2,2*n-1)=cos(teta(n));
    Aeq(3,:)=ones(1, 2*n-1);
    
    beq=[0; 0; 1];
    [x,fval,exitflag,output] = linprog(minimiz(k,n), [], [], Aeq, beq, zeros(2*n-1,1))
    
    % The new parameters of the closed curve are
    q=zeros(n+1,1);
    q(1)=0;
    for i=1:1:(n-1)
        q(i+1)=q(i)+x(2*i)+x(2*i-1);
    end
    q(n+1)=1;
    
    % The reconstructed curve:
    xx=zeros(n+1,1);
    yy=zeros(n+1,1);
    for i=2:1:n+1
        xx(i)=xx(i-1)+L*(q(i)-q(i-1))*cos(teta(i-1));
        yy(i)=yy(i-1)+L*(q(i)-q(i-1))*sin(teta(i-1));
    end
    %plot(xx,yy,'r');
    hold on;
    pause;

%----------------------------------
    % I imposed that the edges q(i)-q(i-1) must be positive
    A=zeros(n, 2*n-1);
    for i=1:1:n-1
        A(i,2*i)=-1;
        A(i,2*i-1)=-1;
    end
    A(n,2*n-1)=-1;
    
    eps=0.01;
    b=-eps.*ones(n,1);
    
    [x,fval,exitflag,output] = linprog(minimiz(k,n), A, b, Aeq, beq, zeros(2*n-1,1));

    % The new parameters:
    q=zeros(n+1,1);
    q(1)=0;
    for i=1:1:(n-1)
        q(i+1)=q(i)+x(2*i)+x(2*i-1);
    end
    q(n+1)=1;

    % The reconstructed curve:
    xx=zeros(n+1,1);
    yy=zeros(n+1,1);
    for i=2:1:n+1
        xx(i)=xx(i-1)+L*(q(i)-q(i-1))*cos(teta(i-1));
        yy(i)=yy(i-1)+L*(q(i)-q(i-1))*sin(teta(i-1));
    end
    %plot(xx,yy,'g');
    pause;
    
    %We observe that increasing the value of eps the result becomes always better
end