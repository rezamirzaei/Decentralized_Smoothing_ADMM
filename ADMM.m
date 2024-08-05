function [beta,dist,acc_rec] = ADMM(X,y,tau,a,lambda,x_true,active,option)
rho_index = option.rho_index;
beta_init = option.beta_init;
max_iters_outer = option.max_iters_outer;
max_iters_inner = option.max_iters_inner;
rho_l = option.rho;
max_it = option.max_it;
dist=zeros(max_it,1);
acc_rec=zeros(max_it,1);
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
co_bl=1;

for om=1:max_it
    dist(om)=norm(beta-x_true)^2;
    acc_rec(om)= sum(abs(sign(beta))==active);
   % rho_l=(sqrt(om))/3;
   % mu=sqrt(3)/rho_l;
    qf=n/rho_l;
        
        %w update step
        w = w_update_me_admm(rho_l, qf, z, X, y, u, lambda, a, A, w, co_bl, p, rho_index);
        %z update step
        z =z_update_me_admm_lla(rho_l, X, y, u, tau, w, n);
        %lambda update step
        %u = lambda.update.me.admm(u,x.admm,y,w,z,rho_l,coo);
        u = u + rho_l * (z - y + X * w);
        
        
         
    
    beta=w;

  
   
       
    
end    

end