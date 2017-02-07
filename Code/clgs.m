
function [Qc,Rc]=clgs(A)
%classical G-S
% Compute the reduced QR factorization

[m,n]=size(A); %dimensions of A
Rc=zeros(n,n);
Qc=zeros(m,n);
I=1:m; %q0 index range 1,2,...,m
for j=1:n
    vj=A(I,j); %column j
    for i=1:j-1;
        Rc(i,j)=dot(Qc(I,i),A(I,j));
        vj=vj-Rc(i,j)*Qc(I,i);
    end
    Rc(j,j)=norm(vj,2); %2-norm
    Qc(I,j)=vj/Rc(j,j);
end