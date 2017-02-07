function [Qh,Rh]=householder(A)

[m,n]=size(A);

for k=1:n
    I=1:(m-k);
    x=A(I,k);
    B=eye(m-k+1);
    e=B(I,1);
    V(I,k)=sign(x(1))*norm(x,2)*e+x;
    V(I,k)=V(I,k)/norm(V(I,k),2);
    A(k:m,k:n)=A(k:m,k:n)-2*V(I,k)*(V(I,k)'*A(k:m,k:n));
end
