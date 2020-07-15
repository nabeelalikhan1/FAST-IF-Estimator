function [fidexmult] = ADTFD_IF_estimation_viterbi_modified_eeg(Sig, num)

% code to reproduce results of the paper, "A modified Viterbi
% algorithm-based if estimation algorithm for adaptive directional time–frequency distributions"

Spec12=zeros(length(Sig),length(Sig));
for i = 1:num
    
    %[Spec,orienttfd]=HTFD_neww(Sig,2,18,54);
    [Spec,orienttfd]=HTFD_neww(Sig,2,20,64);
Spec=Spec/max(Spec(:));
    Spec(Spec<0.05)=0;
    %figure; imagesc(Spec)
    c = findridges_new_viterbi_adtfd(Spec,orienttfd);%(Spec,orienttfd,delta);
    %c=findridges_new_viterbi(Spec);
    %c = findridges_neww(Spec,orienttfd,delta);
    
    
    IF=(c)/(2*length(Sig));
    
    Phase=2*pi*filter(1,[1 -1],IF);
    s_dechirp=exp(-1i*Phase);
    
    
    %im_label2=bwmorph(im_label2,'dilate',3);
    
    % For each sensor do the following steps
    
    L=2;
    %TF filtering for each sensor
    s1 = Sig.*(s_dechirp);
    s2=fftshift(fft(s1));
    PPP=length(s2)/2;
    s3=zeros(1,length(Sig));
    s3(PPP-L:PPP+L)=s2(PPP-L:PPP+L);
    s2(PPP-L:PPP+L)=0;
    extr_Sig=ifft(ifftshift(s3)).*conj(s_dechirp);
    s2=ifft(ifftshift(s2)).*conj(s_dechirp);
    
    %Sig(iii)=Sig(iii)-extr_Sig(iii);
    Sig=s2;%-extr_Sig(iii);
    [Spec11,orienttfd]=HTFD_neww(extr_Sig,2,8,48);
    Spec12=Spec12+Spec11;
    % extr_Sig1=extr_Sig1+extr_Sig;
    %hold on; quiver(1:256,1:256,cos(orienttfd*pi/180),sin(orienttfd*pi/180));
    fidexmult(i,:) = c;
    
end
figure;tfsapl( Sig, Spec12,'SampleFreq',16, 'GrayScale','on' );
end