function W = undirected_graph_generator(n, p, seed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code generates an undirected graph.
%
% I: Number of agents over the graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(seed)
 
while true
    
    G = rand(n,n) < p;
    G = triu(G,1);
    Adj = G + G';
%     p_data_gen = 1 - sqrt(1 - pErdosRenyi);
%     Adj = rand(I);
%     idx1 = (Adj >= p_data_gen);
%     idx2 = (Adj < p_data_gen);
%     Adj(idx1) = 0;
%     Adj(idx2) = 1;
    
%     NotI = ~eye(I);
%     Adj = Adj.*NotI;        % set diagonal entries to 0's
%     Adj = or(Adj,Adj');     % symmetrize, undirected
    degree=diag(sum(Adj));  % degree matrix
    L = degree - Adj;       % standard Laplacian matrix
    lambda = sort(eig(L));
 
    if lambda(2) > 0.05
        %         fprintf(['The Erdos-Renyi graph is generated. Algebraic Connectivity: ',...
        %             num2str(lambda(2)),'\n']);
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