clc
clear
close all
SampFreq = 128*1;
addpath('D:\tfsa_5-5\windows\win64_bin');
t = 0:1/SampFreq:1-1/SampFreq;

A=[zeros(1,32) ones(1,64) zeros(1,32)];

Sig=exp(1*1i*(2*pi*(SampFreq*t/4 +2*sin(4*pi*t))))+exp(1*1i*(2*pi*(SampFreq*t/4 -2*sin(4*pi*t))))+1*A.*exp(1i*2*pi*SampFreq*t/4);
%Sig=awgn(Sig,10,'measured');
%IF_O(:,1)=100*t/1;
%IF_O(:,2)=100*t/1+10;
IF_O(:,1)=90*t.^2/2;
%IF_O(:,4)=-90*t.^2/2+110;
%IF_O(:,3)=90*t.^2/2+15;
IF_O(:,2)=-90*t.^2/2+90/2;

%Sig =1*Sig1 +1*Sig3+1*Sig2;
%Sig=Sig.*([1:128 128:-1:1]);
num=3;


%plot(findex1.')





%plot(findex1.')

%findex1= non_tfd_IF_new_display(Sig,25, num, 3,50,0.1,0.4);

%[I,O]=HTFD_neww(Sig,3,5,30);


[findex1]= FAST_IF(Sig,25+0*26, num, 3,30,0.2,0.4);
figure;plot(findex1.')


%figure;imagesc(I)
%c=findridges_new_viterbi_adtfd(I,O);
%plot(c)
%I=HTFD_neww(Sig,3,4,30);
%figure;imagesc(I)