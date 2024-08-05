function dw = SCAP(w,lamda, a)
size_w = size(w);
s = size_w(2);
dw = zeros(size_w);
for i=1:s-1
   if abs(w(i))>= a*lamda
       dw(i)=0;
   elseif lamda<=abs(w(i)) && abs(w(i))<= a*lamda
       dw(i) = -(w(i)-sign(w(i))*a*lamda)/((a-1));
   elseif  0<abs(w(i)) && abs(w(i))<= lamda
       dw(i) = sign(w(i)) * lamda;
   elseif w(i)== 0 
       dw(i) =  lamda;
   else 
       dw(i)=0;
    end  
end