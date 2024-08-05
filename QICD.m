function [beta,dist,acc_rec] = QICD(X,y,tau,a,lambda,x_true,active,option)
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

beta = beta_init; % initialize beta
%old_beta = beta;  % old_beta is used to check for convergence

om=0;
counter_o=1;
it=0;
while true 
    it=it+1;
    if om>=max_it
        break;
    end
    old_beta = beta; 
    penalty = zeros(p,1);
    counter=5;
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
            elseif beta(j)== 0 
                penalty(j) =  lambda;
            else
                penalty(j) = 0;
            end  
            
   
        else % LASSO
            penalty(j) = lambda;
         end
    end  
    penalty=abs(penalty);
    beta=beta_init;
    for iter_inner = 1:max_iters_inner
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
    weights = zeros(n+1, p);
    u = zeros(n+1, p);
    
    for j = 1:p
        % Calculate residuals and weights
        for i = 1:n
            x_is = X(i, :);
            x_is(j) = 0; % exclude j-th predictor
            u(i, j) = (y(i) - x_is * beta) / X(i, j);
            weights(i, j) = abs(X(i, j) * (tau - double(u(i, j) * X(i, j) < 0))) / n;
        end
        
        % Calculate penalty and weights
        u(n+1, j) = 0;
       
        weights(n+1, j) = penalty(j);
        
        % Update beta_j
        weights_temp = weights(:, j);
        u_temp = u(:, j);
        weightPairs = [u_temp(:), weights_temp(:)];
        sortedWeightPairs = sortrows(weightPairs, 1);
        cumulativeWeights = cumsum(sortedWeightPairs(:, 2));
        medianIndex = find(cumulativeWeights >= (sum(weights_temp) / 2), 1, 'first');
        beta(j) = sortedWeightPairs(medianIndex, 1);
    end
   
    if norm(beta - old_beta_in) < 1e-6
        counter=counter-1;
        if counter <=0
            break;
        end    
    end
    
    if om>=max_it
        break;
    end
    end
   
     if norm(beta - old_beta) < 1e-6
         counter_o=counter_o-1;
         if counter_o<=0
            break;
         end   
    end  
    
end    
if om<max_it
    dist(om:max_it)=dist(om);
    acc_rec(om:max_it)=acc_rec(om);
end
end







