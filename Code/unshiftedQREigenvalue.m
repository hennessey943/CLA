function [schurA,eigA,maxerror]=unshiftedQREigenvalue(A)
%This algorithm takes an arbitrary matrix A and iterates the pure QR
%algorithm to produce a Schur form for A
tol=10^(-5);
k=1;
[~,m]=size(A);
delta=1;
deltaOld=1;
while delta>=tol
    [Q,R]=qr(A);
    A=R*Q;
    deltaM=A;
    
    for i=1:m
        deltaM(i,i)=0;
    end
    delta=max(max(abs(deltaM))); 
            
    if mod(k,10)==0
        fprintf('QR:k=%d : delta=%8.2e, ratio=%5.3f\n',k,delta,delta/deltaOld);
    end
    deltaOld=delta;
    k=k+1;
end
schurA=A;
eigA=diag(A);
E=eig(A);
x=eigA-E;
maxerror=max(abs(x));