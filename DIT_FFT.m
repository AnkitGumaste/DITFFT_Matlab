close all;clear all;clc;
%input sequence
x=[0,1,2,3,4,5,6,7]
N=length(x);


%dividing into even sequence and odd sequence
%Stage 1
xe=0;
x0=0;
for n=0:N-1
    if mod(n,2)==0
        xe(n/2+1)=x(n+1);
    else
        xo((n+1)/2+1)=x(n+1);
    end
end
g=xe
h=xo(2:N/2+1)
%stage 2
for n=0:N/2-1
    if mod(n,2)==0
        ge(n/2+1)=g(n+1);
    else
        go((n+1)/2+1)=g(n+1);
    end
end
a=ge
b=go(2:N/4+1)

for n=0:N/2-1
    if mod(n,2)==0
        he(n/2+1)=h(n+1);
    else
        ho((n+1)/2+1)=h(n+1);
    end
end
c=he
d=ho(2:N/4+1)







