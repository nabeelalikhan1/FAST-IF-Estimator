clc
clear
close all
mul=1;
SampFreq = 128*mul;
addpath('D:\tfsa_5-5\windows\win64_bin');
t = 0:1/SampFreq:1-1/SampFreq;
Sig1 = 1*exp(mul*1i*(1*pi*(30*t.^3))+mul*1i*(2*pi*(0*t))); %300t或者150t
Sig2 = 1*exp(mul*1i*(-1*pi*(30*t.^3))+mul*1i*(1*pi*(90*t))); %300t或者150t

Sig3 = exp(mul*1i*(1*pi*(30*t +30*t.^3)));
Sig4 =1*exp(mul*1i*(1*pi*(120*t -30*t.^3)));
A=[zeros(1,32) ones(1,64) zeros(1,32)];

Sig =1*Sig1.*1 +0*Sig4 +0*Sig3+1*Sig2.*1;
%Sig=awgn(Sig,10,'measured');
%IF_O(:,1)=100*t/1;
%IF_O(:,2)=100*t/1+10;
IF_O(:,1)=mul*90*t.^2/2;
%IF_O(:,4)=-mul*90*t.^2/2+mul*120/2;
%IF_O(:,3)=mul*90*t.^2/2+mul*30/2;
IF_O(:,2)=-mul*90*t.^2/2+mul*90/2;

%Sig=Sig.*([1:128 128:-1:1]);
num=2;
tic
findex1= Proposed_IF_estimation_fast(Sig, num,5);
findex1=RPRG(findex1,0);
toc

tic
%[fidexmult] = non_tfd_IF(Sig,61, num, 3,20);
%[fidexmult,A] = non_tfd_IF_new_display(Sig,length(Sig)/2-1, num, 2,100/10,0.2,0);
%[fidexmult] = non_tfd_IF_new(Sig,length(Sig)/2-1, num, 2,100,0.2,0);


[fidexmult,A] = non_tfd_IF_new_display(Sig,length(Sig)/(2)-1, num, 2,100,0,0);
%[fidexmult,A] = non_tfd_IF_new_fast(Sig,length(Sig)/(2)-1, num, 2,32*4,0,0,1);
                             [findex_viterbi] = ADTFD_IF_estimation_viterbi_modified(Sig, num);

%[fidexmult] = non_tfd_IF_new(Sig,length(Sig)/2-1, num, 2,100,0,0);
toc
figure; 
plot(t,IF_O,':',t,findex1/2,'-','linewidth',3);
axis([0 1 0 64]);
xlabel('Time (s)');
ylabel('Instantaneous Frequency (Hz)')
figure;
plot(t,IF_O,'-',t,SampFreq*fidexmult,':','linewidth',3);
axis([0 1 0 64]);
xlabel('Time (s)');
ylabel('Instantaneous Frequency (Hz)')

figure;
plot(t,IF_O,':',t,findex_viterbi/2,'--','linewidth',3);
axis([0 1 0 64]);
xlabel('Time (s)');
ylabel('Instantaneous Frequency (Hz)')

figure;idx = find(A>=1);
[X, Y, Z] = ind2sub(size(A), idx);
pointsize = 30;
%scatter3(X(:), Y(:), Z(:), pointsize, A(idx));
scatter3(Z(:), X(:), Y(:), pointsize, A(idx));

xlabel('Frequency (Hz)');
ylabel('Time(s)');
zlabel('Chirp rate ');



% HADTFD BASED
