function Q=formQh(W)

[m,~]=size(W);
Q=zeros(m,m);
for i=1:m
    x=zeros(m,1);
    x(i)=1;
    for k=m-2:-1:1
        x(k:m)=x(k:m)-2*W(k:m,k)*(W(k:m,k)'*x(k:m));
        
    end
    Q(1:m,i)=x;
end