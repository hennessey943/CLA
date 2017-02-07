function [xa,xb,xc,xd,xe,xf,xg]=leastSquaresSolver(A,b)
[m,n]=size(A);
%Normal Equations

B=A'*A;
R=chol(B);
w=(R')^(-1)*A'*b;
xa=R^(-1)*w;


%QR factorization and Classical Gram-Schmidt

[Qc,Rc]=clgs(A);
xb=Rc^(-1)*Qc'*b;

%QR factorization and Modified Gram-Schmidt

[Qm,Rm]=mgs(A);
xc=Rm^(-1)*Qm'*b;

%QR factorization and Householder triangularization

[W,Rh]=house(A);
Qh=formQ(W);
xd=Rh^(-1)*Qh'*b;

%QR factorization and Matlab qr

[Q,R]=qr(A,0);
xe=R^(-1)*Q'*b;

%Matlab backslash

xf=A\b;

%Matlab SVD

[U,S,V]=svd(A,0);
y=S^(-1)*U'*b;
xg=V*y;

fprintf('    Normal                CLGS                    MGS                   HOUSE\n');
for i=1:n
    fprintf(' %22.15e %22.15e %22.15e %22.15e\n', xa(i),xb(i),xc(i),xd(i))
end;

fprintf('    Matlab QR             Matlab backslash        SVD\n');
for i=1:n
    fprintf(' %22.15e %22.15e %22.15e\n',xe(i),xf(i),xg(i))
end;