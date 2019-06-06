function [angoli1, angoli2, p] =warping(par1, par2, n1, n2, a, b)
eps=0.000001;
angoli1 =a(1);
angoli2 =b(1);
p=0; % tutti i parametri ordinati e in [0,1]
i=2;
for j=2:1:n2
    found_interval=0;
    while found_interval == 0
        if (par2(j) > par1(i-1)+eps) && (par2(j) < par1(i)-eps)  
            p= [p; par2(j)];
            angoli1=[angoli1; 0];
            angoli2=[angoli2; b(j)];
            found_interval =1;
        elseif (par2(j) > par1(i)-eps) && (par2(j) < par1(i)+eps)
            p=[p; par2(j)];
            angoli1=[angoli1; a(i)];
            angoli2=[angoli2; b(j)];
            i=i+1;
            found_interval = 1;
        else
            p=[p; par1(i)];
            angoli1=[angoli1; a(i)];
            angoli2=[angoli2; 0];
            i=i+1;
        end
    end
end

% Potrei non aver finito di inserire tutti gli elementi della prima curva,
% e se non sono stati inseriti è perchè erano più grandi del più grande
% valore della curva b.
if i <= n1
    for s=i:1:n1
        p=[p; par1(s)];
        angoli1=[angoli1; a(s)];
        angoli2=[angoli2; 0];
    end
end
p=[p; 1];

