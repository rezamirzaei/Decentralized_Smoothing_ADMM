function wn = w_update_me_admm_lla_prox(rho_l, qf, z, x, y, lambda_dual, wp, Mx, sig, w, N)
    u = (((-lambda_dual - rho_l * (z - y))' * x + w' * Mx) / sig);
    w_temp = w;
    for j = 1:N
        qap = u(j);
        cqf = wp(j) * qf;
        if cqf == 0
            w_temp(j) = qap;
        elseif abs(qap) >= cqf
            w_temp(j) = ((-sign(qap) * cqf) + qap);
        else
            w_temp(j) = 0;
        end
    end
    wn = (w_temp);
end