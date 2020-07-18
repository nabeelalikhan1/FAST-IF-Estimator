clc
clear
close all
 load seizure;
 SampFreq = 16;
addpath('D:\tfsa_5-5\windows\win64_bin');

t=0:1/SampFreq:8-1/SampFreq;
f=(0:8/128:8-8/128);
n=1:128;
num=3;
%Sig=awgn(Sig,10,'measured');

Sig=hilbert(Sig);
%[fidexmult,A] = non_tfd_IF_new_display_real(Sig,length(Sig)/(2)-1, 3, 2,30,0,0);
[fidexmult,A] = FAST_IF_EEG(Sig,length(Sig)/(2)-1, num+10, 3,64,0.01,0.1);
figure;
plot(t(1:4:end),SampFreq*fidexmult(:,1:4:end).','-','linewidth',2);
axis([0 8 0 8]);
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
title('IF of seizure signal using FAST-IF method');


I=HTFD_neww(Sig,2,20,64);
figure;imagesc(t,f,I);

set(gcf,'Position',[20 100 640 500]);	    
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
title('ADTFD of the seizure signal');


[f_ridge_tracking] = Proposed_IF_estimation(Sig, 3, 5)/(2*length(Sig));

figure;

plot(t(1:4:end),SampFreq*f_ridge_tracking(:,1:4:end).','-','linewidth',2);
axis([0 8 0 8]);
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
title('IF estimated of seizure using ADTFD(ridge tracking)');


[f_viterbi] = ADTFD_IF_estimation_viterbi_modified_eeg(Sig, 3)/(2*length(Sig));
figure;
plot(t(1:4:end),SampFreq*f_viterbi(:,1:4:end).','-','linewidth',2);
axis([0 8 0 8]);
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
title('IF estimated of seizure +synthetic signal using Viterbi');

%%%%%%%%%%%%%  WITH ADDED NLF CHIRP

IF_O=fidexmult;
IF_O(4,:)=0.01+3*0.000007*n.^2;
IF_O(5,:)=0.48-3*0.000009*n.^2;

n=1:128;
 Sig=Sig+1*0.04*cos(2*pi*(0.01*n+0.000007*n.^3))+0.1*cos(2*pi*(0.48*n-0.000009*n.^3)); 
%Sig=awgn(Sig,15,'measured');
num=5;

Sig=hilbert(Sig);

[fidexmult,A] =FAST_IF_EEG(Sig,length(Sig)/(2)-1, num+10, 3,64,0.01,0.1);

%[fidexmult,A] = FAST_IF(Sig,length(Sig)/(2)-1, num+10, 3,32,0.05,0.1);
figure;
plot(t(1:4:end),SampFreq*fidexmult(:,1:4:end).','--',t(1:4:end),SampFreq*IF_O(:,1:4:end),'o','linewidth',2);
axis([0 8 0 8]);
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
title('FAST_IF of seizure + synthetic signal');


I=HTFD_neww(Sig,2,20,64);
figure;imagesc(t,f,I);
set(gcf,'Position',[20 100 640 500]);	    
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
title('ADTFD of seizure + synthetic signal');



IF_O=f_ridge_tracking;
IF_O(4,:)=0.01+3*0.000007*n.^2;
IF_O(5,:)=0.48-3*0.000009*n.^2;
[f_ridge_tracking] = Proposed_IF_estimation(Sig, 5, 5)/(2*length(Sig));
figure;

plot(t(1:4:end),SampFreq*f_ridge_tracking(:,1:4:end).','--',t(1:4:end),SampFreq*IF_O(:,1:4:end),'o','linewidth',2);
axis([0 8 0 8]);
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
title('IF estimated of seizure +synthetic signal using ADTFD (ridge tracking)');

IF_O=f_viterbi;
IF_O(4,:)=0.01+3*0.000007*n.^2;
IF_O(5,:)=0.48-3*0.000009*n.^2;

[f_viterbi] = ADTFD_IF_estimation_viterbi_modified_eeg(Sig, 5)/(2*length(Sig));
figure;
plot(t(1:4:end),SampFreq*f_viterbi(:,1:4:end).','--',t(1:4:end),SampFreq*IF_O(:,1:4:end),'o','linewidth',2);
axis([0 8 0 8]);
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
title('IF estimated of seizure +synthetic signal using Viterbi');


