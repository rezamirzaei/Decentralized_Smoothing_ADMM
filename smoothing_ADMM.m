function [beta,dist,dist_mean,acc_rec] = smoothing_ADMM(X,y,tau,a,lambda,x_true,active,maxe,option)
 
rho_index = option.rho_index;
beta_init = option.beta_init;
max_iters_outer = option.max_iters_outer;
max_iters_inner = option.max_iters_inner;
rho_l = option.rho;
betta = option.beta;
c=option.c;
d=option.d;
type=option.type;
max_it = option.max_it;
dist=zeros(max_it,1);
dist_mean=zeros(max_it,1);
acc_rec=zeros(max_it,1);
% W = undirected_graph_generator(Num_Nodes);
%x = cell(Num_Nodes,1);
[Num_Nodes,~]=size(y);
net = option.net;
bta=option.w;%10000
[n, p] = size(X{1}); % number of samples and predictors
 A = cell(Num_Nodes,1);
% A=X'*X;
 

    




beta = beta_init; % initialize beta
%old_beta = beta;  % old_beta is used to check for convergence
z = zeros(Num_Nodes, n); % auxiliary variable for ADMM
u = zeros(Num_Nodes, n); % dual variables for ADMM
w = zeros(Num_Nodes, p);
g = zeros(Num_Nodes,Num_Nodes, p);
xi = zeros(Num_Nodes,Num_Nodes, p);
xii = zeros(Num_Nodes,Num_Nodes, p);
om=0;
co_bl=1;
counter=100;
rho_xi=rho_l;
for node=1:Num_Nodes
    A{node}=X{node}'*X{node};
end    
for om=1:max_it
        
    
    
    rho_l=(om^(type))*c;%cc=1/3
    rho_xi=d*(om^(type));%20
    mu=betta/(om^(type));%beta=sqrt(3)
    
    qf=n/rho_l;
        mean_w=sum(w)/Num_Nodes;
        for node=1:Num_Nodes 
        dist_mean(om)=dist_mean(om)+norm(w(node,:)-mean_w)^2/Num_Nodes;    
        dist(om)=dist(om)+norm(w(node,:)-x_true')^2/Num_Nodes;
        acc_rec(om)= acc_rec(om)+(sum(sign(abs(w(node,:)))==active')/p)/Num_Nodes;
        x_local=X{node};
        A_local=A{node};
        %w update step
        num_neigh=sum(net(node,:));
        w(node,:) = w_update_me_admm(rho_l, qf, z(node,:), x_local, y(node,:), u(node,:), lambda, a, A_local, w(node,:), co_bl, p, rho_index,g,xi,rho_xi,n,node,num_neigh);
        %z update step
        z(node,:) =z_update_me_admm_lla_app(rho_l, x_local, y(node,:), u(node,:), tau, w(node,:), n,mu);
        %lambda update step
        %u = lambda.update.me.admm(u,x.admm,y,w,z,rho_l,coo);
        u(node,:) = u(node,:) + rho_l * (z(node,:) - y(node,:) + (x_local * w(node,:)')');
        end
        
     for node = 1:Num_Nodes
        for neig=1:Num_Nodes
            if net(node,neig)==1
                if(neig>node)
                  Ab = (-reshape(w(node,:),1,1,p) + reshape(w(neig,:),1,1,p) - xi(node,neig,:)/rho_xi + xi(neig,node,:)/rho_xi);
                  eta =   rho_xi/(2 * bta);
                  muu=sqrt(2)/eta;
                  E =  reshape(g_update_me_admm(eta, reshape(Ab,1,p), p,muu),1,1,p);
                  %E =  reshape(zeros(1,p),1,1,p);
                  core = reshape(w(node,:),1,1,p) + reshape(w(neig,:),1,1,p) + xi(node,neig,:)/rho_xi + xi(neig,node,:)/rho_xi;
                  
                  g(node,neig,:) = (core - E) / 2;
                  g(neig,node,:) = (core + E) / 2;   
                  
                  xi(node,neig,:) = xi(node,neig,:) + rho_xi*(reshape(w(node,:),1,1,p)-g(node,neig,:));
                  xi(neig,node,:) = xi(neig,node,:) + rho_xi*(reshape(w(neig,:),1,1,p)-g(neig,node,:));
                
                end                  
            end    
        end    
    end

        
         
    
%     beta=w;
%     if norm(z + X * w - y) <= (10^-3 * sqrt(n) + 10^-3 * max(norm(X * w), max(norm(z), norm(y))))
%     if rho_l * norm(X' * (z - z_o)) < 10^-3 * sqrt(p) + 10^-3 * norm(X' * u)
%         counter = counter - 1;
%         if counter <= 0
%             break;
%         end
%     end
%     end

  
   
       
    
end   
beta=w;

end