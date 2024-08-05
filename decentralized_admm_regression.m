function [X,rhol,x_average,dist] = decentralized_admm_regression(Num_Nodes,A,b,x_true,option)
[n,m] = size(A{1});
% W = undirected_graph_generator(Num_Nodes);
%x = cell(Num_Nodes,1);
X = zeros(1,Num_Nodes,n);
x_true= reshape(kron(ones(Num_Nodes,1),x_true'),1,Num_Nodes,n);

eta=option.eta;
rhol=option.mu;
gamm=option.gamma;
cof=option.cof;
zeta=option.zeta;
beta = zeros(Num_Nodes, Num_Nodes,n); 
lambda = zeros(Num_Nodes, Num_Nodes,n);
z = zeros(Num_Nodes, Num_Nodes,n); 
seed=option.seed;

for node = 1:Num_Nodes
    X(1,node,:) = option.x0;
end
v = X;

mu = option.mu;
if option.timevary == false 
    if strcmp(option.top,'star')
     W = sign(undirected_star_generator(Num_Nodes))-eye(Num_Nodes);
    lam = svds(W,2);
    disp(['sigma: ', num2str(min(lam(2)))]);
    elseif strcmp(option.top,'ring')
        W = sign(undirected_ring_generator(Num_Nodes))-eye(Num_Nodes);
        lam = svds(W,2);
        disp(['sigma: ', num2str(min(lam(2)))]);
    else
        W = sign(undirected_graph_generator(Num_Nodes, option.p, 10))-eye(Num_Nodes);
        lam = svds(W,2);
        disp(['sigma: ', num2str(min(lam(2)))]);
    end    
end



% maxa=0;
% flag=false;
% if(flag)
%     for i = 1:Num_Nodes
%         if(maxa<A{i}'*A{i})
%             maxa = A{i}'*A{i};
%         end    
%     end 
% 
%     if or( option.sen,(min(sum(W))*rhol>maxa))
%         1;
%     else
%         0;
%         rhol=maxa/min(sum(W))+1/min(sum(W));
%     end
% end

for i = 1:option.max_it_d
    Xt =X;
%     if mod(i,1000)==0
%     i 
%     rhol
%     end
   
    if option.timevary
        seed = option.seed(i);
        W = sign(undirected_graph_generator(Num_Nodes, option.p, seed))-eye(Num_Nodes);
        lam = svds(W,2);
        sigma(i) = lam(2);
    end
    

     for node = 1:Num_Nodes
      zt =  z(node,:,:);
      lambdat = lambda(node,:,:);
      nl = sum(W(node,:));
      ai = reshape(A{node}(:,i),n,1);
      bi = b(node,i);
      w = (rhol*sum(zt) -sum(lambdat)+eta*(Xt(1,node,:)))/(eta+rhol*nl);
      w= reshape(w,n,1)+ (2*bi*ai/(eta+rhol*nl));
      gam = 1/(rhol*nl + eta);
      X(1,node,:) = reshape(proxr(ai,bi,gam,w),1,1,n);
     end

    
    for node = 1:Num_Nodes
        for neig=1:Num_Nodes
            if W(node,neig)==1 
%                  z(node,neig,:)=(rhol*(X(1,node,:)+X(1,neig,:)) + beta(node,neig,:)/gamm + lambda(node,neig,:)+ lambda(neig,node,:))/(2*rhol+ 1/gamm);
                    z(node,neig,:)=(rhol*(X(1,node,:)+X(1,neig,:)) + lambda(node,neig,:)+ lambda(neig,node,:))/(2*rhol);
            end    
        end    
%         z(node,:,:)=(rhol*(reshape(kron(W(:,node),reshape(X(1,node,:),1,n)),1,Num_Nodes,n)+W(node,:).*X(1,:,:)) + beta(node,:,:)/gamm + reshape(reshape(lambda(:,node,:),Num_Nodes,n)+ reshape(lambda(:,node,:),Num_Nodes,n),1,Num_Nodes,n))/(2*rhol+ 1/gamm);
        
    end

%        beta = beta - zeta*(beta - z);
    for node = 1:Num_Nodes
%         for neig=1:Num_Nodes
%       if W(node,neig)==1 
%               lambda(node,neig,:)=lambda(node,neig,:)+rhol*(X(1,node,:)-z(node,neig,:));
        lambda(node,:,:)=lambda(node,:,:)+rhol*(reshape(kron(W(:,node),reshape(X(1,node,:),1,n)),1,Num_Nodes,n)-z(node,:,:));
 
%       end  
%         end
    end
%     lambda(node,:,:)=lambda(node,:,:)+rhol*(reshape(kron(W(:,node),reshape(X(1,node,:),1,n)),1,Num_Nodes,n)-z(node,:,:));
    


%            rhol=rhol*cof;
if cof>1
            eta=eta+0.5;
           rhol=rhol+0.5;
end   
         
   
%  
%  gamm=gamm*1.001;
    
    %    x_average = zeros(n,1);
    %f(i) = norm(AX_b,1);
    x_average = X;
    dist(i) =   min( sum(sum((x_true-X).^2))/Num_Nodes, sum(sum((x_true+X).^2))/Num_Nodes);
    %flag(i) =   ( norm(x_true-x_average)< norm(x_true+x_average));
end

%semilogy(f(i) - min(f));
disp(['Distance: ', num2str(dist(i))]);

if option.timevary && option.cnt == 1
    disp(['min_sigma: ', num2str(min(sigma))]);
    plot(sigma,   'mo-', 'linewidth',1 );
    set(gca, 'FontSize', 12);
    xlabel('time'); ylabel('\sigma_2(k)');
end