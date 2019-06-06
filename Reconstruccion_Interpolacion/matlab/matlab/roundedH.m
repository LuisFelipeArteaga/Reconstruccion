function v=roundedH(T)

x=zeros(size(T));
y=zeros(size(T));
L=4*pi+22;

index=1;
for t=T
    if t>=0 && t<=((pi/2)/L)
        s=getS(t);
        x(index)=cos(-s-pi/2);
        y(index)=sin(-s-pi/2)+1;
    elseif t>=((pi/2)/L) && t<=(((pi/2+3)/L));
        s=getS(t);
        x(index)=-1;
        y(index)=s-pi/2+1;
    elseif t>=(((pi/2+3)/L)) && t<=((pi+3)/L)
        s=getS(t);
        x(index)=0.5*cos(2*s-pi-6)-3/2;
        y(index)=0.5*sin(2*s-pi-6)+4;
    elseif t>=((pi+3)/L) && t<=((pi+10)/L)
        s=getS(t);
        x(index)=-2;
        y(index)=pi+7-s;
    elseif t>=((pi+10)/L) && t<=((3/2*pi+10)/L)
        s=getS(t);
        x(index)=0.5*cos(2*s-20-pi)-3/2;
        y(index)=0.5*sin(2*s-20-pi)-3;
    elseif t>=((3/2*pi+10)/L) && t<=((3/2*pi+11)/L)
        s=getS(t);
        x(index)=-1;
        y(index)=s-3*pi/2-13;
    elseif t>=((3*pi/2+11)/L) && t<=((5/2*pi+11)/L)
        s=getS(t);
        x(index)=cos(pi/2+11-s);
        y(index)=sin(pi/2+11-s)-2;
    elseif t>=((5/2*pi+11)/L) && t<=((5/2*pi+12)/L)
        s=getS(t);
        x(index)=1;
        y(index)=5*pi/2-s+9;
    elseif t>=((5/2*pi+12)/L) && t<=((3*pi+12)/L)
        s=getS(t);
        x(index)=0.5*cos(2*s-4*pi-24)+3/2;
        y(index)=0.5*sin(2*s-4*pi-24)-3;
    elseif t>=((3*pi+12)/L) && t<=((3*pi+19)/L)
        s=getS(t);
        x(index)=2;
        y(index)=-3*pi+s-15;
    elseif t>=((3*pi+19)/L) && t<=((7/2*pi+19)/L)
        s=getS(t);
        x(index)=0.5*cos(2*s-6*pi-38)+3/2;
        y(index)=0.5*sin(2*s-6*pi-38)+4;
    elseif t>=((7/2*pi+19)/L) && t<=((7/2*pi+22)/L)
        s=getS(t);
        x(index)=1;
        y(index)=-s+7*pi/2+23;
    elseif t>=((7/2*pi+22)/L) && t<=((4*pi+22)/L)
        s=getS(t);
        x(index)=cos(7*pi/2+22-s);
        y(index)=sin(7*pi/2+22-s)+1;
    end

    index=index+1;
end

v=[x; y];

return

function s=getS(t)
s=t*(4*pi+22);
return