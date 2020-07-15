function [fidexmult] = Proposed_IF_estimation_spec_RPGP(Sig, num,delta)

for i = 1:num
Spec=quadtfd(Sig,length(Sig)-1,1,'specx',31,'hamm');
%Spec=quadtfd(Sig,length(Sig)/4-1,1,'mb',0.05,128);

c = findridges(Spec,delta);%(Spec,orienttfd,delta);
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

           % extr_Sig1=extr_Sig1+extr_Sig;
%hold on; quiver(1:256,1:256,cos(orienttfd*pi/180),sin(orienttfd*pi/180));
fidexmult(i,:) = c;

end

end