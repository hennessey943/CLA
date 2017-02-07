function [x,nit]=gmres0(A,b,x0,maxit,tol)
%This function implements the GMRES algorithm to determine an approximate
%solution to a linear system Ax=b. 
%maxit defines the maximum number of iterations allowed
%tol is the residual error tolerance
%x0 is the initial guess of the solution
%The solution is returned as x along with the number of iteration nit
%required
x=x0;
[~,m]=size(A);
I=1:m;
Q(I,1)=b/norm(b);
residOld=1;
for n=1:maxit
    %Calculate the residual from step n-1
    resid=norm(b-A*x);
    if resid<tol
        break;
    end
    
    %Arnoldi iteration
    v=A*Q(I,n);
    for j=1:n
        H(j,n)=Q(I,j)'*v;
        v=v-H(j,n)*Q(I,j);
    end
    H(n+1,n)=norm(v);
    if norm(v)==0
        break;
    end
    Q(I,n+1)=v/H(n+1,n);
    %minimize
    [~,R]=qr(H,0);
    [Omega,~]=qr(H);
    [q,~]=size(Omega');
    e=eye(q,1);
    g=norm(b)*Omega'*e;
    [p,~]=size(R);
    y=R\g(1:p);
    
    x=Q(I,1:n)*y;
    
    fprintf('GMRES: n=%3d, ||r_n||_2 = %8.2e, ratio=%8.2e\n',n,resid,resid/residOld);
    residOld=resid;
end

nit=n;

%Check error
xTrue=A\b;

fprintf('FINAL: Solution from GMRES: (nit=%d)\n',nit);
fprintf(' 2-norm error=%8.2e\n',norm(x-xTrue,2));
fprintf('max-norm error=%8.2e\n',norm(x-xTrue,inf));

        
        
    
    