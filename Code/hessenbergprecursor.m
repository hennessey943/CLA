%hessenberg matrix precursor
m=5;
A=zeros(m,m);
for i=1:m
    for j=1:m
        if i==j
            A(i,j)=9;
        else
            A(i,j)=1/(i+j);
        end
    end
end
