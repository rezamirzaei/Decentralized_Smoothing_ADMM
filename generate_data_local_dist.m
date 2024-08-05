function [a,b,x_true,active,maxe] = generate_data_local_dist(m,n,N,s,cof,tau)

a = cell(N,1);
b = zeros(N,m);
%m                 = 1;   % number of local sample
%n                 = 3;    % problem dimension

%b = zeros(m,1);
Sigma=zeros(n+1,n+1);
for i=1:n+1
  for j=1:n+1
    Sigma(i,j) = 0.5^(abs(i-j));
  end
end

Sigma(:,n+1)=0;
Sigma(n+1,:)=0;
Sigma(n+1,n+1)=1;
 % sup=[6,12,15,20];
  active=zeros(n+1,1);
  if tau~=0.5
  active(n+1)=1;
  end
 
 
 sup=randperm(n,s);
 active(sup)=1;
 x_true= zeros(n+1,1);
 x_true(sup)=1;
% x_true = x_true+randn(n,1)*10^(-3);
 maxe=-inf;
for i=1:N
 data = mvnrnd(zeros(n+1,1),Sigma,m);
 
 
% if s~=0
 %sup=randperm(n-1,s);
 %sup=sup+1;
 %else
 %    sup=[];
 %end    
  %x_true = x_true;
 eps=data(:,n+1);
 data(:,n+1)=0;
 b(i,:) = (data*x_true)+cof*(eps);
 data(:,n+1)=1;
 a{i}=data;
 if max((sum(abs(a{i}'))))>maxe
    maxe=max((sum(abs(a{i}'))));
end  
end
x_true(n+1)=x_true(n+1)+cof*norminv(tau);
end