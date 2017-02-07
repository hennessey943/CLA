m=50;
R=triu(randn(m)); %upper triagnular matrix with random entries
[Q,X]=qr(randn(m)); %generates some Q
y=norm(Q'*Q-eye(m),2)
A=Q*R; %form QR

x=norm(A-Q*R,2)
[Q2,R2]=qr(A);
a=norm(Q2'*Q2-eye(m),2)
b=norm(Q2-Q,2)
c=norm(R2-R,2)
d=norm(A,2)
e=norm(A-Q2*R2,2)/norm(A,2)
Q3=Q+eps*randn(m);
R3=R+eps*randn(m);
f=norm(A-Q3*R3,2)