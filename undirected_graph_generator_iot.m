function W = undirected_graph_generator_iot(n)
while true
r1 = unifrnd(0,2.5,n,1);
r2 = unifrnd(0,2.5,n,1);
r3 = [r1,r2];
G_rand_graph = zeros(n,n);
flag = true;
for i=1:n
    for j=1:n
        
        if i~=j && norm(r3(i,:)-r3(j,:))<0.8
            G_rand_graph(i,j)=1;
        end
        
    end
end

adj = G_rand_graph;
flag=true;
deg_a = sum(adj,2);
max(deg_a);
min(deg_a);
 
   
%for v=1:L
% if 0==sum(adj(v,:))||1==sum(adj(v,:))||2==sum(adj(v,:))
%     flag=false;
% end
if max(deg_a)~=10||min(deg_a)~=2
        flag=false;
end
[con,~]=CheckConnected(adj);
if con==0
     flag=false;
end     
if flag==true
    break;
end    
end
W=adj;
end
