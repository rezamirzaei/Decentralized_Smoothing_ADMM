function [idx,Group] = CheckConnected(E)
L=length(E);
T=find(E(:,1)~=0);
T=[T;1]';
M=sum(E(:,T),2);
T1=union(find(M~=0),T);
while length(T1)>length(T)
    T=T1;
    M=sum(E(:,T),2);
    T1=union(find(M~=0),T);
end
if length(T)==L
    idx=1;
    Group=T;
else
    idx=0;
    Group=T;
end
end