stime=0.0025;
tN=1;
deltat=stime/tN;
N=10;
I=speye(2*N^3,2*N^3);
A1=I-deltat*A/2;
A2=I+deltat*A/2;
rhs=A2*c;
x=A1\rhs;