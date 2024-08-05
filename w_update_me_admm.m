function wn = w_update_me_admm(rho_l, qf, z, x, y, lambda_dual, lambda, gamma, A, w, co_bl, N, penalty,g,xi,rho_xi,n,node,num)
    [nodeNumber,~,~]=size(xi);
  
    g_node=sum(rho_xi*(reshape(g(node,:,:),nodeNumber,N)));
    %g_node=w;
    xi_node=-sum(reshape(xi(node,:,:),nodeNumber,N));
    u = ((-lambda_dual) - rho_l*(z-y))*x +g_node+xi_node;
    %u = ((-lambda_dual) - rho_l*(z-y))*x;
    w_temp = w;
    % w_temp(1)= (u(1) -  A(1,:)*w_temp + w_temp(1) * A(1, 1))/A(1,1);
    for i = 1:co_bl 
        if penalty==0
            for j = 1:N
                ap = u(j) +  (-A(j,:)*w_temp' + w_temp(j) * A(j, j))*rho_l;
                cc = (A(j, j)*rho_l)+(num*rho_xi);
                qap = ap / cc;
                cqf = (n) / cc;
                aqap = abs(qap);
                if j~=N
                if (aqap <= ((cqf+1)*lambda)) && (aqap >= cqf*lambda)
                    w_temp(j) = sign(qap) * (aqap - cqf * lambda);
                  % w_temp(j) = sign(qap)*(abs(qap)-(cqf*lambda));
                else
                    if (aqap >= ((cqf+1)*lambda)) && (aqap <= (gamma*lambda))
                        w_temp(j) = ((gamma-1)*qap - (sign(qap)*cqf*gamma*lambda)) / (-cqf + gamma - 1);
                    else
                        if aqap >= gamma*lambda
                            w_temp(j) = qap;
                        else
                            w_temp(j) = 0;
                        end
                    end
                end
                else
                    w_temp(j) = qap;
                end    
            
            end
        else
            for j = 1:N
                 ap = u(j) +  (-A(j,:)*w_temp' + w_temp(j) * A(j, j))*rho_l;
                cc = (A(j, j)*rho_l)+(num*rho_xi);
                qap = ap / cc;
                cqf = (n) / cc;
                aqap = abs(qap);
                if j~=N
                if aqap <= (cqf * lambda)
                    w_temp(j) = 0;
                else
                    if aqap >= (gamma*lambda)
                        w_temp(j) = qap;
                    else
                        w_temp(j) = (qap - sign(qap) * cqf * lambda) / (1 - (cqf / gamma));
                    end
                end
                else
                   w_temp(j) = qap; 
                end    
            end
        end
    end
    wn = (w_temp);
end

