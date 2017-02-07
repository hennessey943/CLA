%Vandermonde Matrix creator
function [Vm]=Vandermonde(m)
for i=1:m
    x(i)=(i-1)/(m-1);
    for j=1:m
        Vm(i,j)=x(i)^(j-1);
    end
end
