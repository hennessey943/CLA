%Exper. 2 CLGS vs MGS
m=80;
[U,X]=qr(rand(m)); %form U
[V,X]=qr(rand(m)); %form V*
S=diag(2.^(-1:-1:-m)); %diag(1/2,1/4,1/8,...,1/2^j,...1/2^m)
A=U*S*V;
[Qc,Rc]=clgs(A); %classical GS
[Qm,Rm]=mgs(A); %Modified GS
x=(1:m)';
yc=diag(Rc); %diagonalizes rii
ym=diag(Rm);
semilogy(x,yc,'ro',x,ym,'bx');
grid on;
title('Experiment 2. CLGS vs. MGS');