function [aa,ab,ac]  =Angulo_T(a,b,c)
    da = b - a; da = sqrt(sum(da.^2));
    db = b - c; db = sqrt(sum(db.^2));
    dc = c - a; dc = sqrt(sum(dc.^2));

    aa = acos((dc^2+db^2-da^2)/(2*dc*db)); 
    ab = acos((da^2+dc^2-db^2)/(2*da*dc));
    ac = acos((da^2+db^2-dc^2)/(2*da*db));
   
end