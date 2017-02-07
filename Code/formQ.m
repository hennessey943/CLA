function Q=formQ(W)

[m,n]=size(W);

for i=1:n
    x=zeros(m,1);
    x(i)=1;
    for k=n:-1:1
        x(k:m)=x(k:m)-2*W(k:m,k)*(W(k:m,k)'*x(k:m));
        
    end
    Q(1:m,i)=x;
end
    