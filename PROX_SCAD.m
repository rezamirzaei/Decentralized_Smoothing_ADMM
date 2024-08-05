function [dw] = PROX_SCAD(x,gamma,lambda)
p=size(x,2)
    dw = zeros(1, p);
    for i = 1:p
        if abs(w(i)) >= a*lamda
            dw(i) = 0;
        elseif lamda <= abs(w(i)) && abs(w(i)) <= a*lamda
            dw(i) = abs(-(w(i) - sign(w(i)) * a * lamda) / (a - 1));
        elseif 0 < abs(w(i)) && abs(w(i)) <= lamda
            dw(i) = abs(sign(w(i)) * lamda);
        elseif w(i) == 0
            dw(i) = abs(lamda);
        else
            dw(i) = 0;
        end
    end
end

