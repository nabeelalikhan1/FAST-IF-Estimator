clc
clear
close all
 load seizure;

SampFreq = 16;

n=1:128;
addpath('D:\tfsa_5-5\windows\win64_bin');
%Sig=awgn(Sig,10,'measured');
num=4;

Sig=hilbert(Sig);
tic
[fidexmult,A] = FAST_IF_EEG(Sig,length(Sig)/(2)-1, num+10, 3,64,0.05,0.1);

toc
t=0:1/SampFreq:8-1/SampFreq;
f=(0:8/128:8-8/128);
TF=zeros(length(Sig),length(Sig));
[N,M]=size(fidexmult);
for i=1:N
    for j=1:length(Sig)
        TF(round(2*length(Sig)*fidexmult(i,j))+1,j)=1;
    end
end
figure;imagesc(t,f,TF);
set(gcf,'Position',[20 100 640 500]);	    
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
tic
I=HTFD_neww(Sig,2,20,84);
toc
figure;imagesc(t,f,I);
set(gcf,'Position',[20 100 640 500]);	    
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
figure;
plot(t(1:4:end),SampFreq*fidexmult(:,1:4:end).','-','linewidth',2);
axis([0 8 0 8]);
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');

%%%%%%%%%%%%%  WITH ADDED NLF CHIRP

IF_O=fidexmult;
IF_O(4,:)=0.01+3*0.000007*n.^2;
IF_O(5,:)=0.48-3*0.000009*n.^2;

 Sig=Sig+1*0.04*cos(2*pi*(0.01*n+0.000007*n.^3))+1*0.1*cos(2*pi*(0.48*n-0.000009*n.^3)); 
%Sig=awgn(Sig,15,'measured');
num=4;

Sig=hilbert(Sig);
tic
[fidexmult,A] = FAST_IF_EEG(Sig,length(Sig)/(2)-1, num+10, 3,64,0.01,0.1);
toc
t=0:1/SampFreq:8-1/SampFreq;
TF=zeros(length(Sig),length(Sig));
[N,M]=size(fidexmult);
for i=1:N
    for j=1:length(Sig)
        TF(round(2*length(Sig)*fidexmult(i,j))+1,j)=1;
    end
end
figure;imagesc(t,f,TF);
set(gcf,'Position',[20 100 640 500]);	    
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
I=HTFD_neww(Sig,2,20,64);
figure;imagesc(t,f,I);
set(gcf,'Position',[20 100 640 500]);	    
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
figure;
plot(t(1:4:end),SampFreq*fidexmult(:,1:4:end).','--',t(1:4:end),SampFreq*IF_O(:,1:4:end),'o','linewidth',2);
axis([0 8 0 8]);
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');



