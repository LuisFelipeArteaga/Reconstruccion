function v=roundedCroce(T)

x=zeros(size(T));
y=zeros(size(T));
L=4*pi+22;

index=1;
for t=T
    if t>=0 && t<=((pi/4)/L)
        s=getS(t);
        x(index)=1/2+ (1/2)*cos(pi/2+2*s);
        y(index)=3+(1/2)*sin(pi/2+2*s);
    elseif t>=((pi/4)/L) && t<=((pi/4+2)/L);
        s=getS(t);
        x(index)=0;
        y(index)=pi/4+3-s;
    elseif t>=((pi/4+2)/L) && t<=((3*pi/4+2)/L)
        s=getS(t);
        x(index)=-1+cos(pi/4+2-s);
        y(index)=1+sin(pi/4+2-s);
    elseif t>=((3*pi/4+2)/L) && t<=((3*pi/4+4)/L)
        s=getS(t);
        x(index)=3*pi/4+1-s;
        y(index)=0;
    elseif t>=((3*pi/4+4)/L) && t<=((5*pi/4+4)/L)
        s=getS(t);
        x(index)=-3+(1/2)*cos(2*s-8-pi);
        y(index)=-1/2+(1/2)*sin(2*s-8-pi);
    elseif t>=((5*pi/4+4)/L) && t<=((5*pi/4+6)/L)
        s=getS(t);
        x(index)=s-7-5*pi/4;
        y(index)=-1;
    elseif t>=((5*pi/4+6)/L) && t<=((7*pi/4+6)/L)
        s=getS(t);
        x(index)=-1+cos(-s+6-pi/4);
        y(index)=-2+sin(-s+6-pi/4);
    elseif t>=((7*pi/4+6)/L) && t<=((7*pi/4+11)/L)
        s=getS(t);
        x(index)=0;
        y(index)=7*pi/4+4-s;
    elseif t>=((7*pi/4+11)/L) && t<=((9*pi/4+11)/L)
        s=getS(t);
        x(index)=1/2+(1/2)*cos(2*s-5*pi/2-22);
        y(index)=-7+(1/2)*sin(2*s-5*pi/2-22);
    elseif t>=((9*pi/4+11)/L) && t<=((9*pi/4+16)/L)
        s=getS(t);
        x(index)=1;
        y(index)=s-9*pi/4-18;
    elseif t>=((9*pi/4+16)/L) && t<=((11*pi/4+16)/L)
        s=getS(t);
        x(index)=2+cos(5*pi/4+16-s);
        y(index)=-2+sin(5*pi/4+16-s);
    elseif t>=((11*pi/4+16)/L) && t<=((11*pi/4+18)/L)
        s=getS(t);
        x(index)=-14+s-11*pi/4;
        y(index)=-1;
    elseif t>=((11*pi/4+18)/L) && t<=((13*pi/4+18)/L)
        s=getS(t);
        x(index)=4+(1/2)*cos(-6*pi+2*s-36);
        y(index)=-1/2+(1/2)*sin(-6*pi+2*s-36);
    elseif t>=((13*pi/4+18)/L) && t<=((13*pi/4+20)/L)
        s=getS(t);
        x(index)=22+13*pi/4-s;
        y(index)=0;
    elseif t>=((13*pi/4+20)/L) && t<=((15*pi/4+20)/L)
        s=getS(t);
        x(index)=2+cos(11*pi/4+20-s);
        y(index)=1+sin(11*pi/4+20-s);
    elseif t>=((15*pi/4+20)/L) && t<=((15*pi/4+22)/L)
        s=getS(t);
        x(index)=1;
        y(index)=s-15*pi/4-19;
    elseif t>=((15*pi/4+22)/L) && t<=((4*pi+22)/L)
        s=getS(t);
        x(index)=1/2+(1/2)*cos(2*s-44-15*pi/2);
        y(index)=3+(1/2)*sin(2*s-44-15*pi/2);
    end

    index=index+1;
end

v=[x; y];

return

function s=getS(t)
s=t*(4*pi+22);
return
