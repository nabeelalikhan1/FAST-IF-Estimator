function [fidexmult] = MB_IF_estimation(Sig, num, delta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  output   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fidexmult: obtained ridge index for multi-component signals

for i = 1:num
%Spec=tfr_stft_high(Sig);


    Spec = quadtfd(Sig, length(Sig)/4-1, 1, 'mb',0.2,length(Sig));

c = findridges(Spec,delta);
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
                s3=zeros(1,length(Sig));
                s3(length(Sig)/2-L:length(Sig)/2+L)=s2(length(Sig)/2-L:length(Sig)/2+L);
                s2(length(Sig)/2-L:length(Sig)/2+L)=0;
                extr_Sig=ifft(ifftshift(s3)).*conj(s_dechirp);
                s2=ifft(ifftshift(s2)).*conj(s_dechirp);
                
                %Sig(iii)=Sig(iii)-extr_Sig(iii);
                Sig=s2;%-extr_Sig(iii);

           % extr_Sig1=extr_Sig1+extr_Sig;
%hold on; quiver(1:256,1:256,cos(orienttfd*pi/180),sin(orienttfd*pi/180));
fidexmult(i,:) = c;

end

end