function [Qp,Rp]=qr_plus(Q,R)
%make all diagonal elements of R and corresponding columns of Q positive

[m,m]=size(R);
I=1:m;
for i=1:m
        if R(i,i)<0
            Rp(i,I)=-1*R(i,I);
            Qp(I,i)=-1*Q(I,i);
        else
            Rp(i,I)=R(i,I);
            Qp(I,i)=Q(I,i);
        end
end