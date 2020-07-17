clc
clear all;
close all
SampFreq = 256/2;
addpath('D:\tfsa_5-5\windows\win64_bin');
t = 0:1/SampFreq:1-1/SampFreq;


Sig1 = 1*exp(1i*(1*pi*(30*t.^3))+1i*(2*pi*(0*t))); %300t或者150t
Sig2 = 1*exp(1i*(-1*pi*(30*t.^3))+1i*(1*pi*(90*t))); %300t或者150t

Sig3 = exp(1i*(1*pi*(20*t +30*t.^3)));
Sig =1*Sig1 +1*Sig3+1*Sig2;
SigO =Sig;
IF_O(:,1)=90*t.^2/2;
IF_O(:,2)=-90*t.^2/2+90/2;
IF_O(:,3)=90*t.^2/2+10;


%Sig=Sig.*([1:128 128:-1:1]);
num=3;

% HADTFD BASED
[fidexmult,A] = FAST_IF(Sig,length(Sig)/(2)-1, num, 2,100,0,0);
plot(t,IF_O,'-',t,SampFreq*fidexmult,':','linewidth',3);
axis([0 1 0 64]);
xlabel('Time (s)');
ylabel('Instantaneous Frequency (Hz)')