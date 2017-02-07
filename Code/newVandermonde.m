function [Vmn,b]=newVandermonde(m,n)

b=zeros(m,1);
t=zeros(m,1);
Vmn=zeros(m,n);
for i=1:m
    t(i)=i*.02;
    for j=1:n
        Vmn(i,j)=t(i)^(j-1);
    end
    b(i)=cos(4*t(i));
end

