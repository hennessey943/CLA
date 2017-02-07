

function [Qm,Rm]=mgs(A)
%Modified Gram-Schmidt
%Compute the reduced QR factorization

[m,n]=size(A);
Rm=zeros(n,n);
Qm=zeros(m,n);
I=1:m;
v=A;
for i=1:n
    Rm(i,i)=norm(v(I,i),2);
    Qm(I,i)=v(I,i)/Rm(i,i);
    for j=(i+1):n
        Rm(i,j)=dot(Qm(I,i),v(I,j));
        v(I,j)=v(I,j)-Rm(i,j)*Qm(I,i);
    end
end

