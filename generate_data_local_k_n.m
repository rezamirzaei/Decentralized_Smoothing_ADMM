function [a,kt,x_true] = generate_data_local_k_n(m,n,k_n)
%m                 = 1;   % number of local sample
%n                 = 3;    % problem dimension


 %a=abs(randn(m,n))*10+1;
 sup=randperm(m,k_n);
 a=zeros(m,n);
 a(sup)=1;
 
 sort_a=sort(a);
 x_true=sum(sort_a(1:k_n))/k_n;
 kt=sort_a(k_n);

end