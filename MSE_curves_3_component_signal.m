clc
clear all;
close all
SampFreq = 256/2;
addpath('D:\tfsa_5-5\windows\win64_bin');
t = 0:1/SampFreq:1-1/SampFreq;
%Sig1 = 1*exp(1i*(2*pi*(50*t.^1))+1i*(2*pi*(0*t))); %300t或者150t
%Sig2 = 1*exp(1i*(2*pi*(50*t.^2))+1i*(2*pi*(10*t))); %300t或者150t

%Sig3 = exp(1i*(1*pi*(10*t +35*t.^3)));
%Sig4 =1*exp(1i*(1*pi*(110*t -35*t.^3)));
%SigO =1*Sig1 +1*Sig4 +Sig3;
%IF_O(:,1)=100*t/1;
%IF_O(:,2)=100*t/1+10;
%IF_O(:,1)=105*t.^2/2+5;
%IF_O(:,2)=-105*t.^2/2+55;
%IF_O(:,3)=50;


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
NS=200;
IF_O=2*IF_O/length(IF_O);
% HADTFD BASED

%for snr=-10:2:10
iiii=0;
for snr=0:2:10
    iiii=iiii+1;
    
    for k1=1:NS
        
        Sig=awgn(SigO,snr,'measured');
        
        for kkkkk=0:4
            
            % ORIGINAL
            delta=5;
            alpha = 5;
            if kkkkk==0   %ADTFD+VITERBI
                             [findex] = ADTFD_IF_estimation_viterbi_modified(Sig, num);
               
            elseif kkkkk==1 %SPEC+RPGP
          
                %ADTFD+ridge tracking
                findex= Proposed_IF_estimation(Sig, num,5);
            elseif kkkkk==2 %the new algorithm
                [findex] = FAST_IF(Sig,length(Sig)/2-1, num, 2,100,0,0)*2*SampFreq;

            elseif kkkkk==3 %the new algorithm
                [fidexmult] = MB_IF_estimation(Sig, num, delta);
                fidexmult(fidexmult>128)=128;
                fidexmult(fidexmult<1)=1;
                [findex,~] = RPRG(fidexmult,5);            
            
            elseif kkkkk==4 %the new algorithm
               [fidexmult] = Proposed_IF_estimation_spec_RPGP(Sig, num,5);
               fidexmult(fidexmult>128)=128;
                fidexmult(fidexmult<1)=1;
                [findex,interset] = RPRG(fidexmult,5);           
              
            end
            
            msee=0.1*ones(1,num);
            IF=zeros(1,length(Sig));
            dis=0;
            clear c;
            
            for ii=1:num
                
                t=1:SampFreq;
                IF=findex(ii,:)/length(Sig);
                t=t(5:end-5);
                for i=1:num
                    c(i)=sum(abs(IF(t)'-IF_O(t,i)).^2);
                end
                [a1 b1]=min(c);
                if msee(b1)>=a1(1)/length(Sig)
                    msee(b1)=a1(1)/length(Sig);
                end
                if dis==1
                    figure;
                    plot(t,IF(t),'-',t,IF_O(t,b1),'d');
                end
            end
            if kkkkk==0
                mse_adtfd_viterbi_modified(k1)=mean(msee);
            elseif kkkkk==1
                mse_adtfd_ridgetracking(k1)=mean(msee);
            elseif kkkkk==2
                mse_non_tfd(k1)=mean(msee);
            
            elseif kkkkk==3
                mse_Spec_rpgp(k1)=mean(msee);
            elseif kkkkk==4
                mse_mb_rpgp(k1)=mean(msee);
            end
            
            
        end
        
        
    end
    mse_viterbi(iiii)=mean(mse_adtfd_viterbi_modified)
    mse_spec(iiii)=mean(mse_Spec_rpgp)
    mse_adtfd(iiii)=mean(mse_adtfd_ridgetracking)
        mse_mb(iiii)=mean(mse_mb_rpgp)

    mse_proposed(iiii)=mean(mse_non_tfd)
end


snr=0:2:10;
plot(snr, 10*(log10(mse_viterbi)),'-rh','linewidth',4);


hold on;

plot(snr, 10*(log10(mse_adtfd)),'-bh','linewidth',4);

hold on;
plot(snr, 10*(log10(mse_proposed)),'-.k+','linewidth',4);

hold on;
plot(snr, 10*(log10(mse_mb)),'-.y+','linewidth',4);

hold on;
plot(snr, 10*(log10(mse_spec)),'-.g+','linewidth',4);




xlabel('Signal to Noise Ratio');
ylabel('Mean Square Error (dB)');
legend('Viterbi based on ADTFD','Ridge detection and tracking using ADTFD','The Proposed Algorithm', 'Ridge path regrouping using MBD', 'Ridge path regrouping using Spectrogram');
% axis([min(snr) max(snr)  -50  0])
