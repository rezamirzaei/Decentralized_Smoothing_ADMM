function [a,b,x_true,active] = generate_data_local_fixed(m,n,s,cof,tau)
%m                 = 1;   % number of local sample
%n                 = 3;    % problem dimension

b = zeros(m,1);
Sigma=zeros(n+1,n+1);
for i=1:n+1
  for j=1:n+1
    Sigma(i,j) = 0.5^(abs(i-j));
  end
end

Sigma(:,n+1)=0;
Sigma(n+1,:)=0;
Sigma(n+1,n+1)=1;
 data = mvnrnd(zeros(n+1,1),Sigma,m);
 a=data(:,1:n);
% if s~=0
 %sup=randperm(n-1,s);
 %sup=sup+1;
 %else
 %    sup=[];
 %end    
 sup=[6,12,15,20];
 active=zeros(n,1);
 if tau~=0.5
 active(1)=1;
 end
 active(sup)=1;
 
 %temp=randn(s,1);
 x_true= zeros(n,1);
 x_true(sup)=1;
 x_true = x_true+randn(n,1)*10^(-3);
  %x_true = x_true;
 x_true(1)=0;
 eps=data(:,n+1);
 a(:,1)=normcdf(a(:,1));
  b = (a*x_true)+cof*(eps.*a(:,1));
  x_true(1)=x_true(1)+cof*norminv(tau);
end