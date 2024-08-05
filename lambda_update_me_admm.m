function lambda = lambda_update_me_admm(lambda_o, x, y, w, z, rho_l)
    lambda = lambda_o + rho_l * (z - y + x * w);
end