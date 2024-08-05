function [a,b,x_true] = generate_data_cluster_federated_learning(N,cl,cofn,m,n,s,snr)
a = cell(N,1);
x_true = cell(N,1);
x_true1 = cell(N,1);
b = zeros(N,m);

 sup=randperm(n,s);
 xx= mvnrnd(zeros(n,1),eye(n),cl);
 for i=1:cl
     for j=1:floor(N/cl)
        x_true1{(floor(N/cl)*(i-1))+j} = xx(i,:)'+randn(n,1)*cofn;
        x_true{(floor(N/cl)*(i-1))+j}=xx(i,:)';
     end    
 end    
 
 data= mvnrnd(zeros(n,1),eye(n),m*N);
 
 
for i = 1:N
    for j = 1:m

        a{i}(:,j) = data((i-1)*m+j,1:n)'; 
      b(i,j)= (a{i}(:,j)'*x_true1{i})+snr*randn(1);
      
    end
end
end