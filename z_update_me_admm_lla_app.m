function zn = z_update_me_admm_lla_app(rho_l, x, y, lambda_dual, tau, w, N, mu)
    etta = (y - (x * w')' - (lambda_dual / rho_l)) + (1 - 2 * tau) / (2 * rho_l);
    zn = zeros(N, 1);
    for j = 1:N
        if etta(j) >= mu + 1 / (2 * rho_l)
            zn(j) = etta(j) - 1 / (2 * rho_l);
        elseif etta(j) <= -mu - 1 / (2 * rho_l)
            zn(j) = etta(j) + 1 / (2 * rho_l);
        else
            zn(j) = etta(j) / (1 + 1 / (2 * rho_l * mu));
        end
    end
end