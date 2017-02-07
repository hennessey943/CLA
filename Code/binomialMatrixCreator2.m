function P=binomialMatrixCreator2(M,r,p)

P=zeros(M+r+1,M+r+1);
P(1,M+r+1)=1;
for i=2:M+r+1
    for j=1:M+r+1
        if i<=M+1
            P(i,j)=binopdf(M+r+1-j,i-1,p/(i-1));
        else
            P(i,j)=binopdf(M+r+1-j,M,p/M);
        end
        
    end
end
t=1;
x=1;
[V,D,W]=eig(P);
while abs(x)>=.0001 
    x=1-D(t,t);
    t=t+1;
end
pi=W(:,t-1)/sum(W(:,t-1));
expcost=0;
for i=0:M-r-1
    expcost=(M-i)*pi(i+1)+expcost;
end
for i=M-r:M-1
    expcost=2*(M-i)*pi(i+1)+expcost;
end
for i=M+1:M+r
    expcost=(i-M)*pi(i+1)+expcost;
end
expcost


