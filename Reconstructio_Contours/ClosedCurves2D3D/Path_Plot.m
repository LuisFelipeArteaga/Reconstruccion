function ftc = Path_Plot(alpha,p2n)

[n,T,k] = size(alpha);
stp = k-1;
if(n == 3)
    dt = 0.3;
else
    dt = 0.25;
end
ftc = {};
nn=1;
for tau = 1:stp+1
    if (tau==stp+1)
        %p = p2n;
        p = q_to_curve(alpha(:,:,tau));
    else
        p = q_to_curve(alpha(:,:,tau));
    end
%     ft = p;
%     ft(1,:) = p(1,:) + dt*tau;
    ft(1,:) = p(1,:) ;
    %ft(1,:) = ft(1,:) - mean(ft(1,:)) +dt*tau;
    ft(1,:) = ft(1,:) - mean(ft(1,:))+dt*tau;
    for i = 2:n
        ft(i,:) = p(i,:) - mean(p(i,:));
    end
    ftc{nn}=ft';
    nn = nn+1;
end

return;