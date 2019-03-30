clear 
clc
close all

%Read in the original image
img=imread('cameraman.tif');
subplot(4,3,1),imshow(img);
title('original')
[m,n]=size(img);

%transform the image from spaital domain to frequency domain
G=fftshift(fft2(img));

%get the motion blur frequency domain filter
H=zeros(m,n);
for u=1:m
    for v=1:n
        if (pi*(0.15*(u-m/2)+0.1*(v-n/2)))==0
            H(u,v)=1;
        else
            H(u,v)=(1*sin(pi*(0.15*(u-m/2)+0.1*(v-n/2)))*exp(-1i*pi*(0.15*(u-m/2)+0.1*(v-n/2))))/(pi*(0.15*(u-m/2)+0.1*(v-n/2)));
        end
    end
end

%apply the filter the image 
F=G.*H;

%transform the image from frequency domain back to spatial domain
f=real(ifft2(ifftshift(F)));
f=(mat2gray(f));
subplot(4,3,2),imshow(f);
title('Motion Blur')

%apply 3 kinds of noice to the blurred image
noise5=imnoise(f,'gaussian',0,0.5);
subplot(4,3,4),imshow(noise5);
title('corrupted 1')
noise002=imnoise(f,'gaussian',0,0.002);
subplot(4,3,5),imshow(noise002);
title('corrupted 2')
noise00065=imnoise(f,'gaussian',0,0.00065);
subplot(4,3,6),imshow(noise00065);
title('corrupted 3')

%inverse filter
Hinv=zeros(m,n);
for u=1:m
    for v=1:n
        if (abs(H(u,v)))>0.003
            Hinv(u,v)=1/H(u,v);
        else
            Hinv(u,v)=0;
        end
    end
end
Gapprox=fftshift(fft2(noise5));
Fapprox=zeros(u,v);
for u=1:m
    for v=1:n
        Fapprox(u,v)=Gapprox(u,v)*Hinv(u,v);
    end
end
fapprox=real(ifft2(ifftshift(Fapprox)));
%fapprox = fapprox./(max(fapprox(:)));
subplot(4,3,7),imshow(fapprox);
title('inverse filter-1')

Gapprox2=fftshift(fft2(noise002));
Fapprox2=Gapprox2.*Hinv;
fapprox2=real(ifft2(ifftshift(Fapprox2)));
subplot(4,3,8),imshow(fapprox2);
title('inverse filter-2')

Gapprox3=fftshift(fft2(noise00065));
Fapprox3=Gapprox3.*Hinv;
fapprox3=real(ifft2(ifftshift(Fapprox3)));
subplot(4,3,9),imshow(fapprox3);
title('inverse filter-3')


%Wiener 1 - using abs(H).^2
K = 0.01;
H_Wiener = ((abs(H).^2)./((abs(H).^2)+K)).*(1./H);
F_Wiener = H_Wiener .*  Gapprox;  
f_Wiener = real(ifft2(F_Wiener)); 
%f_Wiener = f_Wiener ./(max(f_Wiener(:)));
subplot(4,3,10),imshow(f_Wiener);
title('Winener filter(1)-1')

F_Wiener2 = H_Wiener .*  Gapprox2;  
f_Wiener2 = real(ifft2(F_Wiener2)); 
%f_Wiener = f_Wiener ./(max(f_Wiener(:)));
subplot(4,3,11),imshow(f_Wiener2);
title('Winener filter(1)-2')

F_Wiener3 = H_Wiener .*  Gapprox3;  
f_Wiener3 = real(ifft2(F_Wiener3)); 
%f_Wiener = f_Wiener ./(max(f_Wiener(:)));
subplot(4,3,12),imshow(f_Wiener3);
title('Winener filter(1)-3')

%Wiener 2 - using Hconj*H - basically Wiener 1 and 2 are the same 
Hconj=conj(H);
Habs2=Hconj.*H;
F_Wiener21=(Hinv.*Habs2./(Habs2+K)).*Gapprox;
f_Wiener21 = real(ifft2(F_Wiener21)); 
%subplot(5,3,13),imshow(f_Wiener21);
%title('Winener filter(2)-1')

F_Wiener22=(Hinv.*Habs2./(Habs2+K)).*Gapprox2;
f_Wiener22 = real(ifft2(F_Wiener22)); 
%subplot(5,3,14),imshow(f_Wiener22);
%title('Winener filter(2)-2')

F_Wiener23=(Hinv.*Habs2./(Habs2+K)).*Gapprox3;
f_Wiener23 = real(ifft2(F_Wiener23)); 
%subplot(5,3,15),imshow(f_Wiener23);
%title('Winener filter(2)-3')

