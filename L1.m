function p = L1(u,lamda)
size_u = size(u);
s = size_u(2);
p = zeros(size_u);
for i=1:s-1
   if abs(u(i))>= 0
       p(i)=sign(u(i))*lamda;
   else
       p(i)=lamda;
   end
end  
end