function W = undirected_graph_generator_cluster(n,cl, p, seed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code generates an undirected graph.
%
% I: Number of agents over the graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(seed)
LL = cell(cl,1);
L = zeros(n,n);
G= zeros(n,n);
while true
for i=1:cl
while true
    
    G_1 = rand(floor(n/cl),floor(n/cl)) < p;
    G_1 = triu(G_1,1);
    Adj_1 = G_1 + G_1';
%     p_data_gen = 1 - sqrt(1 - pErdosRenyi);
%     Adj = rand(I);
%     idx1 = (Adj >= p_data_gen);
%     idx2 = (Adj < p_data_gen);
%     Adj(idx1) = 0;
%     Adj(idx2) = 1;
    
%     NotI = ~eye(I);
%     Adj = Adj.*NotI;        % set diagonal entries to 0's
%     Adj = or(Adj,Adj');     % symmetrize, undirected
    degree_1=diag(sum(Adj_1));  % degree matrix
    L_1 = degree_1 - Adj_1;       % standard Laplacian matrix
    lambda = sort(eig(L_1));
 
    if lambda(2) > 0.05
        %         fprintf(['The Erdos-Renyi graph is generated. Algebraic Connectivity: ',...
        %             num2str(lambda(2)),'\n']);
        break;
    end
end
G((i-1)*floor(n/cl)+1:i*floor(n/cl),(i-1)*floor(n/cl)+1:i*floor(n/cl))=G_1;
end
G=sign(G+(rand(n,n) < p/cl));
G = triu(G,1);
Adj = G + G';
degree=diag(sum(Adj));  % degree matrix
    L = degree - Adj;       % standard Laplacian matrix
    lambda = sort(eig(L));
 
    if lambda(2) > 0.05
        break;
    end

end



%%%% generate MH weight matrix %%%%
A = zeros(n);
for i = 1:n
    i_link = find(Adj(i,:)>0);
    for j=1:n
        if i~=j && sum(find(j == i_link))>0
            A(i,j) = 1 / (max(degree(i,i), degree(j,j)) + 1);
        end
    end
end

W = eye(n)-diag(sum(A))+ A; % row stochastic matrix for x