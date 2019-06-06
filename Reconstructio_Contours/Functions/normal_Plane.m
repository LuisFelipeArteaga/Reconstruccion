function n = normal_Plane(data)
    warning('off','all')
    can_solve = zeros(1,3);
    
    for ii = 1:3
        X = data;
        X(:,ii) = 1;
        X_m = X'*X;
        if det(X_m) == 0
			can_solve(ii) = 0;
            continue
        end
        can_solve(ii) = 1;
        coeff = (X_m)^-1 * X' * data(:,ii);
        c_neg = -coeff;
        c_neg(ii) = 1;
        coeff(ii) = 1;
        n(:,ii) = c_neg / norm(coeff);
    end
	center = mean(data);
	off_center = [data(:,1)-center(1) data(:,2)-center(2) data(:,3)-center(3)];
    for jj = 1:3
        if can_solve(jj) == 0
            residual_sum(jj) = NaN;
            continue
        end
        residuals = off_center * n(:,jj);
        residual_sum(jj) = sum(residuals .* residuals);	
    end
	best_fit = find(residual_sum == min(residual_sum));
	n = n(:,best_fit(1));
    warning
end