% Gaussian Elimination with Partial Pivoting
function [L,U,P]=lufactor(A)
[~,m]=size(A);
U=A;
L=eye(m,m);
P=eye(m,m);
for k=1:(m-1)
    max=abs(U(k,k));
    maxindex=k;
    for i=k+1:m
        if abs(U(i,k))>max
            max=abs(U(i,k));
            maxindex=i;
        end
    end
    b=U(k,k:m);
    U(k,k:m)=U(maxindex,k:m);
    U(maxindex,k:m)=b;
    c=L(k,1:(k-1));
    L(k,1:(k-1))=L(maxindex,1:(k-1));
    L(maxindex,1:(k-1))=c;
    d=P(k,1:m);
    P(k,1:m)=P(maxindex,1:m);
    P(maxindex,1:m)=d;
    for j=k+1:m
        L(j,k)=U(j,k)/U(k,k);
        U(j,k:m)=U(j,k:m)-L(j,k)*U(k,k:m);
    end
end
