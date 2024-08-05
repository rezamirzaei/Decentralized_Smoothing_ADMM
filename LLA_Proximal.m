function [beta,dist,acc_rec] = LLA_Proximal(X,y,tau,a,lambda,x_true,active,option)
rho_index = option.rho_index;
beta_init = option.beta_init;
max_iters_outer = option.max_iters_outer;
max_iters_inner = option.max_iters_inner;
max_it = option.max_it;
dist=zeros(max_it,1);
acc_rec=zeros(max_it,1);
rho_l=option.rho;
% W = undirected_graph_generator(Num_Nodes);
%x = cell(Num_Nodes,1);
[n, p] = size(X); % number of samples and predictors
 A=X'*X;


beta = beta_init; % initialize beta
%old_beta = beta;  % old_beta is used to check for convergence
z = zeros(n, 1); % auxiliary variable for ADMM
u = zeros(n, 1); % dual variables for ADMM
w = zeros(p, 1);
om=0;

sig = rho_l*(max(eigs(A))+n);
 Mx = sig*eye(p)-rho_l*A;
 qf = n/sig;
 counter=5;
 it=0;
while true 
    it=it+1;
    if om>=max_it
        break;
    end
    old_beta = beta; 
    penalty = zeros(p,1);
    
     for j=1:p
         if rho_index == 0 % SCAD
           
         if abs(beta(j))>= a*lambda
            penalty(j)=0;
         elseif lambda<=abs(beta(j)) && abs(beta(j))<= a*lambda
            penalty(j) = -(beta(j)-sign(beta(j))*a*lambda)/((a-1));
         elseif  0<abs(beta(j)) && abs(beta(j))<= lambda
            penalty(j) = sign(beta(j)) * lambda;
         elseif beta(j)== 0 
            penalty(j) =  lambda;
         else 
            penalty(j)=0;
         end
         
        elseif rho_index == 1 % MCP
            
            if abs(beta(j))> a*lambda 
                penalty(j)=0;
            elseif abs(beta(j))<=a*lambda && abs(beta(j))>0 
                penalty(j) = sign(beta(j))*lambda -(beta(j)/a);
            else 
                penalty(j) = lambda;
            end  
            
   
        else % LASSO
            penalty(j) = lambda;
         end
    end  
     penalty=abs(penalty);
   z = zeros(n, 1); % auxiliary variable for ADMM
u = zeros(n, 1); % dual variables for ADMM
w = zeros(p, 1);
 counter=5;
    for iter_inner = 1:max_iters_inner
        z_o=z;
        om=om+1;
        dist(om)=norm(beta-x_true)^2;
          if it>1
         if om>1
             if  dist(om)>norm(old_beta-x_true)^2
                 dist(om)=norm(old_beta-x_true)^2;
             end   
         end    
         end
    acc_rec(om)= sum(abs(sign(beta))==active);
    if it>1
    if om>1
             if  acc_rec(om)<sum(abs(sign(old_beta))==active)
                 acc_rec(om)=sum(abs(sign(old_beta))==active);
             end   
    end 
    end
        old_beta_in = beta; 
        %w update step
        w = w_update_me_admm_lla_prox(rho_l, qf, z, X, y, u, penalty, Mx, sig, w, p);
        %z update step
        z =z_update_me_admm_lla(rho_l, X, y, u, tau, w, n);
        %lambda update step
        %u = lambda.update.me.admm(u,x.admm,y,w,z,rho_l,coo);
        u = u + rho_l * (z - y + X * w);
        
        
         
    
   beta=w;
    if norm(z + X * w - y) <= (10^-3 * sqrt(n) + 10^-3 * max(norm(X * w), max(norm(z), norm(y))))
    if rho_l * norm(X' * (z - z_o)) < 10^-3 * sqrt(p) + 10^-3 * norm(X' * u)
        counter = counter - 1;
        if counter <= 0
            break;
        end
    end
    end
     if norm(beta - old_beta_in) < 1e-8
            break;
    end
    if om>=max_it
        break;
    end
    end
    if norm(beta - old_beta, 2) < 1e-6
            break;
    end
       
    
end    
if om<max_it
    dist(om:max_it)=dist(om);
    acc_rec(om:max_it)=acc_rec(om);
end
dist(max_it)=norm(old_beta-x_true)^2;
end