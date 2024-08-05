function zn = g_update_me_admm(rho_l, x, N, mu)
    etta = x;
    zn = zeros(1,N);
    for j = 1:N
        if etta(j) >= mu + 1 / (rho_l)
            zn(j) = etta(j) - 1 / (rho_l);
        elseif etta(j) <= -mu - 1 / ( rho_l)
            zn(j) = etta(j) + 1 / ( rho_l);
        else
            zn(j) = etta(j) / (1 + 1 / (rho_l * mu));
        end
    end
end