% Exp 1: Discrete Legendre
m=257;
mh=(m-1)/2;
x=(-mh:mh)'/mh;
A=[x.^0,x.^1,x.^2,x.^3,x.^4,x.^5,x.^6,x.^7,x.^8,x.^9,x.^10,x.^11,x.^12];
[Q,R]=qr(A);
scale=Q(m,:);
Q=Q*diag(1./scale);
plot(Q,'LineWidth',2)
axis([0,m,-1,1]);
grid on;
title('Lecture 9 discrete Legendre Poly');