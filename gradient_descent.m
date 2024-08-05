function [beta,dist,acc_rec] = gradient_descent(X,y,tau,a,lambda,x_true,active,option)
rho_index = option.rho_index;
beta_init = option.beta_init;
max_iters_outer = option.max_iters_outer;
max_iters_inner = option.max_iters_inner;
max_it = option.max_it;
dist=zeros(max_it,1);
acc_rec=zeros(max_it,1);
% W = undirected_graph_generator(Num_Nodes);
%x = cell(Num_Nodes,1);
[n, p] = size(X); % number of samples and predictors
%inz= ones(p-1,1);
%inz=[0;inz];
beta = beta_init; % initialize beta
dQR_a = 3.5;
dQR_b = 0.51;
for om = 1:max_it
    dist(om)=norm(beta-x_true)^2;
       acc_rec(om)= sum(abs(sign(beta))==active);
        % Compute loss and gradient
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
    
       loss = y - X*beta;
       gq = X'*(sign(abs(loss)).*(-0.5*sign(loss) + 0.5 - tau));
       %Delta_i(:,i) = (1/(M_i*L))*s_bar(:,i) + SCAP(dQR_w_h(:,i),lam_f,a_f);
  %     Delta_i(:,i) = (1/(M_i*L))*s_bar(:,i) + (1/(L*sqrt(M_i)))*inz.*MCP(dQR_w_h_MCP(:,i),lam_f_mcp,a_f_mcp);
  
       g = (1/(n))*gq + penalty;
        
       
       learning_rate= (dQR_a/(om^dQR_b));
        % Update rule: gradient descent
       beta = beta - learning_rate * g;
       
        % Project onto positive orthant (if necessary)
     
end

end