function wn = w_update_me_admm_lla(rho_l, qf, z, x, y, lambda_dual, wp, A, w, co_bl, N)
    u = (-lambda_dual / rho_l - (z-y))'*x ;
    w_temp = w;
    for i = 1:co_bl
        w_temp1 = w_temp;
        for j = 1:N
            ap = u(j) -  A(j,:)*w_temp + w_temp(j) * A(j, j);
            cc = A(j, j);
            qap = ap / cc;
            cqf = wp(j) * qf / cc;
            if cqf == 0
                w_temp(j) = qap;
            elseif abs(qap) >= cqf
                w_temp(j) = (-sign(qap)*cqf) + qap;
            else
                w_temp(j) = 0;
            end
        end
        if ((w_temp - w_temp1)' * (w_temp - w_temp1)) < 10^-6
            break;
        end
    end
    wn = (w_temp);
end