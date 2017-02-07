%Experiment 3
A=[.7000,.70711;.70001,.70711];
[Q,R]=qr(A); %built-in
norm(Q'*Q-eye(2)) %should give 2.3382e-16

%MGS result should be 2.3014e-11 same as CLGS

