clc
clear
close all
Max_iter=15;
delta_f=4;
L=90;
threshold_1=0.1 % ratio between the maximum and minimum amplitude
n=-128:127;
N=256;
Sig=1*exp(1i*n.^3/N^2+1i*3*pi*n.^2/(4*N))+1*exp(1i*n.^4/(70*N^2)+1i*pi*n/(8))+1*exp(1i*n.^4/(81*N^2)-1i*pi*n.^2/N);

IF(1,:)=(3*n.^2/N^2+6*pi*n/(4*N))/(2*pi);
IF(2,:)=(4*n.^3/(70*N^2)+pi/8)/(2*pi);
IF(3,:)=(4*n.^3/(81*N^2)-2*pi*n/(N))/(2*pi);

tic
[findexmult] = FAST_IF_c(Sig,121, Max_iter, delta_f,L,0.1,0.0);

toc
figure;plot((findexmult.'),'b','linewidth',3)
hold on;plot(IF.','r:','linewidth',3)
axis([1  256  -0.5  0.5])
title('Estimated IF (blue) vs original IF (red)');