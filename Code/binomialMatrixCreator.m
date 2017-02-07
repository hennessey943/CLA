% P=zeros(6,6);
% P(1,6)=binopdf(0,0,.5);
% P(2,6)=binopdf(0,1,.5);
% P(2,5)=binopdf(1,1,.5);
% P(3,6)=binopdf(0,2,.5/2);
% P(3,5)=binopdf(1,2,.5/2);
% P(3,4)=binopdf(2,2,.5/2);
% P(4,6)=binopdf(0,3,.5/3);
% P(4,5)=binopdf(1,3,.5/3);
% P(4,4)=binopdf(2,3,.5/3);
% P(4,3)=binopdf(3,3,.5/3);
% P(5,6)=binopdf(0,4,.5/4);
% P(5,5)=binopdf(1,4,.5/4);
% P(5,4)=binopdf(2,4,.5/4);
% P(5,3)=binopdf(3,4,.5/4);
% P(5,2)=binopdf(4,4,.5/4);
% P(6,6)=binopdf(0,5,.5/5);
% P(6,5)=binopdf(1,5,.5/5);
% P(6,4)=binopdf(2,5,.5/5);
% P(6,3)=binopdf(3,5,.5/5);
% P(6,2)=binopdf(4,5,.5/5);
% P(6,1)=binopdf(5,5,.5/5);

function P=binomialMatrixCreator(M,p)

P=zeros(M+1,M+1);
P(1,M+1)=1;
for i=2:M+1
    for j=1:M+1
        P(i,j)=binopdf(M+1-j,i-1,p/(i-1));
    end
end

[V,D,W]=eig(P);
pi=W(:,1)/sum(W(:,1));
for t=1:M+1
    expected=(t-1)*pi(t)+expected;
end
expected=M-expected
