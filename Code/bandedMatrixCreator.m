function A=bandedMatrixCreator(m,p)

A=zeros(m,m);
for i=1:m
    for j=1:m
        if abs(i-j)<=p
           A(i,j)=rand;
        end
    end
end
