function dw = MCP(w, lambda, a)
size_w = size(w);
s = size_w(2);
dw = zeros(size_w);
for i=1:s-1
   if abs(w(i)) > a*lambda 
       dw(i) = 0;
   elseif abs(w(i)) <= a*lambda && abs(w(i)) > 0 
       dw(i) = sign(w(i)) * lambda - (w(i) / a);
   else 
       dw(i) =  lambda;
   end
end
end