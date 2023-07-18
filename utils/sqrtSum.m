function res = sqrtSum(x,isfft,isflipabs)
%sqrt of sum of abs
%input
%x: 3D k-space matrix, nx,ny,ncoil
%isfft: FFT option, 1 fft
%isflipabs: abs option

if nargin<3
    isflipabs=1;
end

if isflipabs==1
    if nargin<=1
        res = sqrt(sum(abs(ifft2c(x)).^2,3));%rssq(ifft2c(x),3);%
    else
        if isfft==0
            res = sqrt(sum(abs(x).^2,3));%rssq(x,3);%
        else
            res = sqrt(sum(abs(ifft2c(x)).^2,3));%rssq(ifft2c(x),3);%
        end
    end
else
    if nargin<=1
        res = sqrt(sum((ifft2c(x)).^2,3));%rssq(ifft2c(x),3);%
    else
        if isfft==0
            res = x;%rssq(x,3);%
        else
            res = sqrt(sum((ifft2c(x)).^2,3));%rssq(ifft2c(x),3);%

        end
    end
end

