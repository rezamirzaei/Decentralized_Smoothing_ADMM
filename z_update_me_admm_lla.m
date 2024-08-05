function zn = z_update_me_admm_lla(rho_l, x, y, lambda_dual, tau, w, N)
    etta = (y - (x * w) - (lambda_dual / rho_l));
    zn = sign(etta + (1 - 2 * tau) / (2 * rho_l)) .* max(abs(etta + (1 - (2 * tau)) / (2 * rho_l)) - 1 / (2 * rho_l), zeros(N, 1));
end