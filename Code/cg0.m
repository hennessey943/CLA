function [x,nit]=cg0(A,b,x0,maxit,tol)
%This function applies the conjugate gradient iteration to determine an 
%approximate solution to the system Ax=b
%maxit defines the maximum number of iterations allowed
%tol is the residual error tolerance
%x0 is the initial guess of the solution
%The solution is returned as x along with the number of iteration nit
%required


[~,m]=size(A);
I=1:m;
x=x0;
r(I,1)=b-A*x0;
p(I,1)=r(I,1);
residOld=norm(b-A*x0,2);
%note we start at 2 so as to allow a more convenient notation
for n=2:maxit+1
    alpha(n-1)=(r(I,n-1)'*r(I,n-1))/(p(I,n-1)'*A*p(I,n-1)); %step length
    x=x+alpha(n-1)*p(I,n-1); %approximate solution
    r(I,n)=r(I,n-1)-alpha(n-1)*A*p(I,n-1); %residual
    resid=norm(r(I,n),2);
    fprintf('CG: n=%3d, || r_n ||_2 = %8.2e, ratio=%8.2e\n',n-1,resid,resid/residOld);
    if resid<tol
        break;
    end
    residOld=resid;
    beta(n-1)=(r(I,n)'*r(I,n))/(r(I,n-1)'*r(I,n-1)); %improvement
    p(I,n)=r(I,n)+beta(n-1)*p(I,n-1); %search direction
end
nit=n-1;

xTrue=A\b;
fprintf('FINAL: Solution from CG: (nit=%d)\n',nit);
fprintf(' 2-norm error=%8.2e\n',norm(x-xTrue,2));
fprintf('max-norm error=%8.2e\n',norm(x-xTrue,inf));
