function [beta_new,omega_new,dist] = gradient_descent_knn(X,tau,x_true,option)

%beta_init = option.beta_init;
kn = option.k;
max_it = option.max_it;
dist=zeros(max_it,1);
% W = undirected_graph_generator(Num_Nodes);
%x = cell(Num_Nodes,1);
[n, p] = size(X); % number of samples and predictors
%inz= ones(p-1,1);
%inz=[0;inz];
beta = zeros(n, p); % initialize beta
dQR_a = 100;
dQR_b = 0.7;
dQR_c = 100;
dQR_d = 0.8;
omega=zeros(n, p);
beta_new=beta;
omega_new=omega;


L=n;
while true
r1 = unifrnd(0,2.5,L,1);
r2 = unifrnd(0,2.5,L,1);
r3 = [r1,r2];
G_rand_graph = zeros(L,L);
flag = true;
for i=1:L
    for j=1:L
        
        if i~=j && norm(r3(i,:)-r3(j,:))<0.8
            G_rand_graph(i,j)=1;
        end
        
    end
end
 % G_rand_graph = round(rand(L));
% G_rand_graph = triu(G_rand_graph) + triu(G_rand_graph,1)';
% G_rand_graph = G_rand_graph - diag(diag(G_rand_graph));
  adj = G_rand_graph;
% % 
  
deg_a = sum(adj,2);
C = zeros(L,L);
for i=1:L
    for j=1:L
        if(i==j)
            C(i,j)=0;
        else
            C(i,j) = adj(i,j)/(max(deg_a(i),deg_a(j)));
        end
    end
end


max(deg_a);
min(deg_a);
metra_wt = C + diag(1-sum(C,2));
for v=1:L
if 0==sum(adj(v,:))||1==sum(adj(v,:))||2==sum(adj(v,:))
    flag=false;
end
end
if flag==true
    break;
end    
end
 
   
for om = 1:max_it
    for j=1:n

       
    
       loss = X(j) - beta_new(j);
      
      gq =(sign(abs(loss)).*(-0.5*sign(loss) + 0.5 - tau));
      
      
       %Delta_i(:,i) = (1/(M_i*L))*s_bar(:,i) + SCAP(dQR_w_h(:,i),lam_f,a_f);
  %     Delta_i(:,i) = (1/(M_i*L))*s_bar(:,i) + (1/(L*sqrt(M_i)))*inz.*MCP(dQR_w_h_MCP(:,i),lam_f_mcp,a_f_mcp);
  
       g = (1/(n))*gq;
        
       
       learning_rate= (dQR_a/(om^dQR_b));
       %learning_rate= 0.01;
        % Update rule: gradient descent
       beta(j) = beta_new(j) -learning_rate * g;
       
       learning_rate1=(dQR_c/(om^dQR_d));
       if X(j)<=(beta(j)+0.001 )
           g= (omega_new(j)-X(j))/kn;
           omega(j)= omega_new(j) -  g*0.01;
       else
           omega(j)= omega_new(j);
       end    
        % Project onto positive orthant (if necessary)
    end
    
    for j=1:n
       %qQR_phi = dQR_w(i,:)- (dQR_a/(dQR_b^k))*Delta_i(:,i);
       beta_new(j) =  beta'*metra_wt(:,j);
       omega_new(j)= omega'*metra_wt(:,j);
      % dQR_w_L1(:,i) = dQR_phi(:,i) + (dQR_a/(L* k^dQR_c )) *(dQR_w_h_L1* adj(i,:)'-degi(i)* dQR_w_h_L1(:,i));
    end
   
    dist(om) = sum(( omega_new - ones(n,1).*x_true).^2)/n;
    
end

end