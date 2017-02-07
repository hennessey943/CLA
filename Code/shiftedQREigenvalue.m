function [eigA,maxerror]=shiftedQREigenvalue(A)
%This algorithm takes an arbitrary matrix A and iterates the practical QR
%algorithm with shifts to produce a Schur form for A

%Initialize tolerance and variables
tol=10^(-5);
k=1;
deltaOld=1;
[~,m]=size(A);
E=eig(A);
%delta=zeros(14,1);
%loop the practical QR algorithm until the tolerance is satisfied
while m>1
    %Choose the Wilkinson Shift for this step
    del=(A(m-1,m-1)-A(m,m))/2;
    if del==0
        mu=A(m,m)+A(m,m-1)^2/(abs(del)+sqrt(del^2+A(m,m-1)^2));
    else
        mu=A(m,m)-sign(del)*A(m,m-1)^2/(abs(del)+sqrt(del^2+A(m,m-1)^2));
    end
    
    %QR factorize A-mu I
    [Q,R]=qr(A-mu*eye(m));
    %Recombine factors in reverse order
    A=R*Q+mu*eye(m);
    
    delta(k)=abs(A(m,m-1));
            
    fprintf('QR:k=%d : shift=%5.3f, dimension=%d, delta=%8.2e, ratio=%5.3f\n',k,mu,m,delta(k),delta(k)/deltaOld);
    
    
    %Test the convergence of the eigenvalue in the row m, save the
    %converged eigenvalue and deflate the matrix
    if delta(k)<tol
        eigA(m)=A(m,m);
        if m==2
            eigA(m-1)=A(m-1,m-1);
        end
        
        A=A(1:m-1,1:m-1);
    end
    [~,m]=size(A);
    
    deltaOld=delta(k);
    k=k+1;
end
eigA=sort(eigA);
x=eigA'-E;
maxerror=max(abs(x));
t=1:k-1;
plot(t,delta)