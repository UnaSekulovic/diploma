function psnr = PSNR(Xfull,Xrecover,maxP)
%
% Written by Canyi Lu (canyilu@gmail.com)
% 
Xrecover = max(0,Xrecover);
Xrecover = min(maxP,Xrecover);
[n1,n2,n3] = size(Xrecover);
%norm(double(Xfull(:))-double(Xrecover(:)))
MSE = norm(Xfull(:)-Xrecover(:))^2/(n1*n2*n3);
%maxP^2/MSE
psnr = 10*log10(maxP^2/MSE);
