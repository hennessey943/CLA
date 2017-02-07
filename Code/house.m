function [W,R]=house(A)
%householder algorithm
[m,n]=size(A);

for k=1:n
    I=k:m;
    x=A(k:m,k);
    e=zeros(m-k+1,1);
    e(1)=1;
    if x(1)==0
        V(I,k)=norm(x,2)*e+x;
    else
        V(I,k)=sign(x(1))*norm(x,2)*e+x;
    end
    V(I,k)=V(I,k)/norm(V(I,k),2);
    A(k:m,k:n)=A(k:m,k:n)-2*V(I,k)*(V(I,k)'*A(k:m,k:n));
end
W=V;
R=A(1:n,1:n);