function [lambda,x,r ] = powerEigenvalueIteration( A,v0 )
%Applies the power iteration method to approximate the largest eigenvalue
%of the matrix A
delta=1;
tol=10^(-5);
v=v0;
[~,D]=eig(A);
k=1;
deltaOld=1;
detlambda=max(max(D));
while delta>=tol
    w=A*v;
    v=w/norm(w);
    lambda=v'*A*v;
    delta=abs(detlambda-lambda);
    r(k)=delta/deltaOld;
    k=k+1;
    deltaOld=delta;
end
x=v;
