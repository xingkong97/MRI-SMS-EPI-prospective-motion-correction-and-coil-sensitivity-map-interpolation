function recon_slice=DoSplitSliceGrappa(kdata_sliceimg, KSpaceDATA, sliceOrderSMS, prot)
%this function is to recon SMS aliased images by split slice-GRAPPA method
%author: Bo Li, 4-28-2021
%input
%kdata_sliceimg: single-slice reference k-sapce data
%KSpaceDATA: k-space of SMS imaging
%sliceOrderSMS: SMS acquisition order
%prot: imaging protocol

%output
%recon_slice: recon images

    nx=size(KSpaceDATA,1);ny=size(KSpaceDATA,2);
    recon_slice=zeros(nx,ny,prot.chn,prot.OriNslice, 'double');
    kspaceSingleSlice=zeros(size(kdata_sliceimg,1),size(kdata_sliceimg,2),size(kdata_sliceimg,3),prot.lMultiBandFactor);
    
    for iSMS=1:prot.OriNslice / prot.lMultiBandFactor 
        kspaceSingleSlice(:,:,:,1:end)=kdata_sliceimg(:,:,:,sliceOrderSMS((iSMS-1)*prot.lMultiBandFactor+1:iSMS*prot.lMultiBandFactor));
       if prot.lMultiBandFactor==2
            if prot.Nslice==1
                CAIPIshifts=[0,pi/2];
            else
                CAIPIshifts=[0,pi/2];    
            end
        elseif prot.lMultiBandFactor==3
            CAIPIshifts= [0,pi*2/3,pi*4/3];
        elseif prot.lMultiBandFactor==4
            CAIPIshifts= [0,pi*2/3,pi*4/3,pi*2];
        elseif prot.lMultiBandFactor==5
            CAIPIshifts= [0,pi*2/3,pi*4/3,pi*2,pi*2/3]; 
        elseif prot.lMultiBandFactor==6
            CAIPIshifts= [0,pi*2/3,pi*4/3,pi*2,pi*2/3,pi*4/3]; 
        end

        %k-space shift along ky direction
        kspaceSingleSlice=SMS_CAIPIshift(kspaceSingleSlice,CAIPIshifts);
        kspace_sms=KSpaceDATA(:,:,:,iSMS);
        calib=kspaceSingleSlice;
        
        %img_sms=sqrtSum(kspace_sms);
        clear kData_ori
        nt=1; %nt is time frame (nreps)
        kspace_sms=repmat(kspace_sms, [1,1,1,nt]);

        kernel_size=[7,7];
        lamda=0.1;
        
        res_k=SplitSliceGrappa(kspace_sms, calib, kernel_size, lamda); 
        res_k=squeeze(res_k(:,:,:,1,:));%nx,ny,ncoil,nreps,nslices
        
        res_k=SMS_CAIPIshift(res_k,-CAIPIshifts);
        
        imgIndex=sliceOrderSMS((iSMS-1)*prot.lMultiBandFactor+1:iSMS*prot.lMultiBandFactor);
        
        recon_slice(:,:,:,imgIndex)=ifft2c(res_k);
    end
end