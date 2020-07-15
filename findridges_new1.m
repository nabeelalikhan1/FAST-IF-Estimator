function index = findridges_new1(Spec,Specangle,delta)
% IF estimation algorithm published in "Instantaneous frequency estimation
% of intersecting and close multi-component signals with varying amplitudes"

Spec = abs(Spec);
T=15;
II(:,:,1)=zeros(size(Spec));
II(:,:,2)=zeros(size(Spec));
II(:,:,3)=zeros(size(Spec));

[M,N] = size(Spec);
index = zeros(1,N);
index1=index;
Spec=Spec/max(Spec(:));
[fmax,tmax] = find(Spec == max(Spec(:)));
fmax = fmax(1);
tmax = tmax(1);
index(tmax) = fmax;
Specangle(Specangle>90)=-180+Specangle(Specangle>90);
theta=Specangle(fmax,tmax);
f0 = fmax;
jj=-1;
for j = (min(tmax+1,N)):N
    f0=round(f0);
    
    low = max(1,f0-delta);
    up = min(M,f0+delta);
 
[~,f0] = max(Spec(low:up,j));
    f0 = f0 + low - 1;
    
    if abs(Specangle(round(f0),j)-theta)>T
        if jj==-1
            jj=j-1;
        end
        
        f0=index(jj)-(j-jj)*tan(pi*theta/180);
        if f0>length(Spec)
            f0=length(Spec);
        elseif f0<1
            f0=1;
        end
        II(round(f0),j,1)=1;

    else
        jj=-1;
        
            theta=0.8*theta+0.2*Specangle(round(f0),j);
                       theta=Specangle(round(f0),j);

    if f0<1
        f0=1;
    end
    if f0>length(Spec)
        f0=length(Spec);
    end
   
             II(f0,j,2)=1;
        
        
    end
    
    
    index(j) = f0;
    index1(j)=theta;
    %f0
end
theta=Specangle(fmax,tmax);
jj=-1;

f1 = fmax;
for j = (max(1,tmax-1)):-1:1
    f1=round(f1);
    low = max(1,f1-delta);
    up = min(M,f1+delta);
    [~,f1] = max(Spec(low:up,j));
    f1 = f1 + low - 1;
    if abs(Specangle(round(f1),j)-theta)>T
        if jj==-1
            jj=j+1;
        end

        

        f1=index(jj)-(j-jj)*tan(pi*theta/180);
    if f1>length(Spec)
        f1=length(Spec);
    end
    if f1<1
        f1=1;
    end

        II(round(f1),j,1)=1;

    else
        jj=-1;
            theta=0.8*theta+0.2*Specangle(round(f1),j);
                        theta=Specangle(round(f1),j);

    if f1>length(Spec)
        f1=length(Spec);
    end
    if f1<1
        f1=1;
    end
   
                        II(f1,j,2)=1;
    
    end
    
        
    index(j) = f1;
     II(round(f1),j)=1;
   
end
index=round(index);
%figure;imagesc(II)
end

