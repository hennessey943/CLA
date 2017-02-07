function x=lusolve(b,L,U,P)
[~,m]=size(L);
 z=P*b;
 
 
 y=zeros(m,1);
 x=zeros(m,1);
 
 for i=1:m
     y(i)=(z(i)-dot(y,L(i,1:m)))/L(i,i);
 end
 
 for k=m:-1:1
     x(k)=(y(k)-dot(x,U(k,1:m)))/U(k,k);
 end
 