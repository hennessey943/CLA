function [W,H]=hessenberg(A)

[~,m]=size(A);

for k=1:m-2
    I=(k+1):m;
    %V=zeros(I,1);
    x=A(k+1:m,k);
    e=zeros(m-k,1);
    e(1)=1;
    if x(1)==0
        V(I,k)=norm(x,2)*e+x;
    else
        V(I,k)=sign(x(1))*norm(x,2)*e+x;
    end
    V(I,k)=V(I,k)/norm(V(I,k),2);
    A(k+1:m,k:m)=A(k+1:m,k:m)-2*V(I,k)*(V(I,k)'*A(k+1:m,k:m));
    A(1:m,k+1:m)=A(1:m,k+1:m)-2*(A(1:m,k+1:m)*V(I,k))*V(I,k)';
end
W=V;
H=A(1:m,1:m);