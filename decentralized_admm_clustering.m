function [X,rhol,x_average,dist] = decentralized_admm_clustering(Num_Nodes,A,b,x_true,option)
[n,m] = size(A{1});
% W = undirected_graph_generator(Num_Nodes);
%x = cell(Num_Nodes,1);
X = zeros(1,Num_Nodes,n);
x_true= reshape(cell2mat(x_true),n,Num_Nodes);
x_true=x_true';
x_true=reshape(x_true,1,Num_Nodes,n);
eta=option.eta;
rhol=option.mu;
gamm=option.gamma;
cof=option.cof;
zeta=option.zeta;
beta = zeros(Num_Nodes, Num_Nodes,n); 
lambda = zeros(Num_Nodes, Num_Nodes,n);
z = zeros(Num_Nodes, Num_Nodes,n); 
seed=option.seed;
flagu=zeros(1,n);
bta=option.eta;
norm_type=option.norm_type;
cl=option.cl;
phi=0.1;
% for node = 1:Num_Nodes
%     X(1,node,:) = option.x0;
% end
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
        W = sign(undirected_graph_generator_cluster(Num_Nodes,cl ,option.p, 10))-eye(Num_Nodes);
        lam = svds(W,2);
        disp(['sigma: ', num2str(min(lam(2)))]);
    end    
end



 minn=inf;
% flag=false;
% if(flag)
     for i = 1:Num_Nodes
         if(minn>min(eigs(A{i}'*A{i},n)))
             minn = min(eigs(A{i}'*A{i},n));
         end    
     end 
disp(minn)
% 
%     if or( option.sen,(min(sum(W))*rhol>maxa))
%         1;
%     else
%         0;
%         rhol=maxa/min(sum(W))+1/min(sum(W));
%     end
% end
Wa=W;
for i = 1:option.max_it_d
    flagu=zeros(1,Num_Nodes);
    Xt =X;
    if(norm_type == 0)
    rhol= sqrt(i)/3;
    end
%     if mod(i,1000)==0
%     i 
%     rhol
%     end
   
%     if option.timevary
%         seed = option.seed(i);
%         W = sign(undirected_graph_generator(Num_Nodes, option.p, seed))-eye(Num_Nodes);
%         lam = svds(W,2);
%         sigma(i) = lam(2);
%     end
      for node = 1:Num_Nodes
        for neig=1:Num_Nodes
            if W(node,neig)==1
                if(neig>node)
                  Ab = (X(1,neig,:) - X(1,node,:) - lambda(node,neig,:)/rhol + lambda(neig,node,:)/rhol);
                  eta =   rhol/(2 * bta);
            if(norm_type == 0)
                  E= Ab;
            elseif(norm_type == 1) % L1 norm
                 % muu = 0.00001;
                 muu = sqrt(3)/eta;
                  E=reshape(z_update_me_admm_lla_app(Ab,eta, muu),1,1,n);
                 % E=shrinkage(Ab,eta);
           %E = prox_SCAD(A,eta,3.7,0.001);
            elseif(norm_type == 2) % L2 norm
                  E = prox_SCAD(Ab,1/eta,2.1,zeta);
            elseif(norm_type == 3) % Laplacian Norm
                  E =prox_MCP(Ab,1/eta,1.4,zeta);
            elseif(norm_type == 6) % New: Nuclear Norm
                  [U, S, V] = svd(Ab);
                  E = U * shrinkage(S, eta) * V';
            end
       
            core = X(1,node,:) + X(1,neig,:) + lambda(node,neig,:)/rhol + lambda(neig,node,:)/rhol;
            z(node,neig,:) = (core - E) / 2;
            z(neig,node,:) = (core + E) / 2;   
                
                end  
%                  z(node,neig,:)=(rhol*(X(1,node,:)+X(1,neig,:)) + beta(node,neig,:)/gamm + lambda(node,neig,:)+ lambda(neig,node,:))/(2*rhol+ 1/gamm);
                 %   z(node,neig,:)=(rhol*(X(1,node,:)+X(1,neig,:)) + lambda(node,neig,:)+ lambda(neig,node,:))/(2*rhol);
                   % z(neig,node,:)=(rhol*(X(1,node,:)+X(1,neig,:)) + lambda(node,neig,:)+ lambda(neig,node,:))/(2*rhol);
            end    
        end    

%         z(node,:,:)=(rhol*(reshape(kron(W(:,node),reshape(X(1,node,:),1,n)),1,Num_Nodes,n)+W(node,:).*X(1,:,:)) + beta(node,:,:)/gamm + reshape(reshape(lambda(:,node,:),Num_Nodes,n)+ reshape(lambda(:,node,:),Num_Nodes,n),1,Num_Nodes,n))/(2*rhol+ 1/gamm);
        
    end

     for node = 1:Num_Nodes
      %if sum(Wa(node,:))>=phi*sum(W(node,:))
      
      zt =  z(node,:,:);
      lambdat = lambda(node,:,:);
      nl = sum(W(node,:));
      ai = reshape(A{node}(:,:),n,m);
      bi = b(node,:);
      w = (rhol*sum(zt) -sum(lambdat))/(rhol*nl);
      w= reshape(w,n,1)+ (2*ai*bi'/(rhol*nl));
      gam = 1/(rhol*nl);
      X(1,node,:) = reshape(proxr(ai,bi,gam,w),1,1,n);
    %  end
     end
     

  
     

    
  

%        beta = beta - zeta*(beta - z);
    for node = 1:Num_Nodes
         for neig=1:Num_Nodes
       if W(node,neig)==1 
               lambda(node,neig,:)=lambda(node,neig,:)+rhol*(X(1,node,:)-z(node,neig,:));
      
 %       lambda(node,:,:)=lambda(node,:,:)+rhol*(reshape(kron(W(:,node),reshape(X(1,node,:),1,n)),1,Num_Nodes,n)-z(node,:,:));

       end  
         end
    end
%     lambda(node,:,:)=lambda(node,:,:)+rhol*(reshape(kron(W(:,node),reshape(X(1,node,:),1,n)),1,Num_Nodes,n)-z(node,:,:));
    


%            rhol=rhol*cof;
% if cof>1
%             eta=eta+0.5;
%            rhol=rhol+0.5;
% end   
         
   
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