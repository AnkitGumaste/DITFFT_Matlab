close all;clear all;clc;

%input sequence
x=[0,1,2,3,4,5,6,7]
N=length(x);

%Twiddle factor
w0=1.0000 + 0.0000i;         
w1=0.7071 - 0.7071i;
w2=0.0000 - 1.0000i;
w3=-0.7071 - 0.7071i;

%Computing A
A0=x(1)+w0*x(5);
A1=x(1)-w0*x(5);

%Computing B
B0=x(3)+w0*x(7);
B1=x(3)-w0*x(7);

%Computing C
C0=x(2)+w0*x(6);
C1=x(2)-w0*x(6);

%Computing D
D0=x(4)+w0*x(8);
D1=x(4)-w0*x(8);

%Computing G
G0=A0+w0*B0;
G1=A1+w2*B1;
G2=A0-w0*B0;
G3=A1-w2*B1;

%Computing H
H0=C0+w0*D0;
H1=C1+w2*D1;
H2=C0-w0*D0;
H3=C1-w2*D1;

%Computing X
X0=G0+w0*H0;
X1=G1+w1*H1;
X2=G2+w2*H2;
X3=G3+w3*H3;
X4=G0-w0*H0;
X5=G1-w1*H1;
X6=G2-w2*H2;
X7=G3-w3*H3;

X=[X0 X1 X2 X3 X4 X5 X6 X7]

X_Verify=fft(x)