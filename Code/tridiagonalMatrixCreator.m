%create an m by m tridiagonal matrix 
m=100;
A=zeros(m,m);
A(1,1)=4; %avoids the issue of accessing A(1,0)
A(1,2)=1;
for i=2:m
    A(i,i)=4; %Diagonal entries
    A(i,i+1)=1; %Superdiagonal entries
    A(i,i-1)=1; %subdiagonal entries
end
A=A(1:m,1:m); %Gets rid of accidental m+1 column

%form b
b=zeros(m,1);
for i=1:m
    b(i)=1+i/m;
end
