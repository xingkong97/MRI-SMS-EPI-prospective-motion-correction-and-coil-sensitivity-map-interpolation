function res = SMS_CAIPIshift(Kdata,CAIPIshifts, shiftMode)
% Kdata: kx ky coil slice
% shiftMode:0 for even, 1 for odd
% CAIPIshifts is a 1D vector of length Nslices (=size(K,4))
% RFphases (optional) is the same dimension as CAIPIshifts (not sure if needed)
%
% Shift image space corresponding to K according to CAIPI shifts
%
% Stan Rapacchi

Ksz = size(Kdata);


RFphases = 0*CAIPIshifts;


Nslices = size(Kdata,4);%Ksz(4);

if Ksz(2)==1
    [xx,~]=meshgrid(1-Ksz(1)/2:Ksz(1)/2,1);  
    xx=xx';
else
    [xx,~]=meshgrid(1-Ksz(2)/2:Ksz(2)/2,1:Ksz(1));
end

 xx = repmat(xx,[1 1 Ksz(3)]);

%kspace=kspace*exp(i*ky*Angle)
res=zeros(size(Kdata));
for s=1:Nslices
    res(:,:,:,s)=Kdata(:,:,:,s).*exp(1i*CAIPIshifts(s)*xx) * exp(1i*RFphases(s));
end

if nargin>=3
    if shiftMode==0 %for even PE lines
        res(:,1:2:end,:,:)=Kdata(:,1:2:end,:,:);   
    else %for odd lines
        res(:,2:2:end,:,:)=Kdata(:,2:2:end,:,:); 
    end
end

end
