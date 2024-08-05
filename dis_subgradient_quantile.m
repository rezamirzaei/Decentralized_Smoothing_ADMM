function [beta,dist,dist_mean,acc_rec] = dis_subgradient_quantile(X,y,tau,a,lambda,x_true,active,option)
rho_index = option.rho_index;
beta_init = option.beta_init;
max_iters_outer = option.max_iters_outer;
max_iters_inner = option.max_iters_inner;
rho_l = option.rho;
betta = option.beta;
cc=option.c;
type=option.type;
max_it = option.max_it;
net=option.net;
dist=zeros(max_it,1);
dist_mean=zeros(max_it,1);
acc_rec=zeros(max_it,1);
% W = undirected_graph_generator(Num_Nodes);
%x = cell(Num_Nodes,1);
[Num_Nodes,~]=size(y);
dQR_a = 3.5;
dQR_b = 0.51;
[n, p] = size(X{1}); % number of samples and predictors
% A=X'*X;
 

     
deg_a = sum(net,2);
C = zeros(Num_Nodes,Num_Nodes);
for i=1:Num_Nodes
    for j=1:Num_Nodes
        if(i==j)
            C(i,j)=0;
        else
            C(i,j) = net(i,j)/(max(deg_a(i),deg_a(j)));
        end
    end
end


max(deg_a);
min(deg_a);

metra_wt = C + diag(1-sum(C,2));


beta = beta_init; % initialize beta
%old_beta = beta;  % old_beta is used to check for convergence
z = zeros(Num_Nodes, n); % auxiliary variable for ADMM
u = zeros(Num_Nodes, n); % dual variables for ADMM
w = zeros(Num_Nodes, p);
w_temp = zeros(Num_Nodes, p);
om=0;
co_bl=1;
counter=100;
rho_xi=rho_l;
for om=1:max_it
        mean_w=sum(w)/Num_Nodes;
    for node=1:Num_Nodes
       dist_mean(om)=dist_mean(om)+norm(w(node,:)-mean_w)^2/Num_Nodes;    
       dist(om)=dist(om)+norm(w(node,:)-x_true')^2/Num_Nodes;
       dist(om)=dist(om)+norm(w(node,:)-x_true')^2/Num_Nodes;
       acc_rec(om)= acc_rec(om)+(sum(sign(abs(w(node,:)))==active')/p)/Num_Nodes;
       x_local=X{node};
       temp = y(node,:) - w(node,:)*x_local';
       temp1 = (sign(abs(temp)).*(-0.5*sign(temp) + 0.5 - tau))*x_local;
       if rho_index==0
       temp2= (1/(n*Num_Nodes))*temp1 + (1/(Num_Nodes))*SCAP(w(node,:),lambda,a);
       elseif rho_index==1
       temp2= (1/(n*Num_Nodes))*temp1 + (1/(Num_Nodes))*MCP(w(node,:),lambda,a);
       else
       temp2= (1/(n*Num_Nodes))*temp1 + (1/(Num_Nodes))*L1(w(node,:),lambda);    
       end    
           
       w_temp(node,:) = w(node,:)- (dQR_a/(om^dQR_b))*temp2;
   end    
   for node=1:Num_Nodes
       w(node,:) = metra_wt(node,:)*w_temp;    
   end
   
end  
beta=w;

end