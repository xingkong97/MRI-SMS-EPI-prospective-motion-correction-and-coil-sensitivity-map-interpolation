function recon_slice=DoSENSE_SMS(coilSen, KSpaceDATA, sliceOrderSMS, prot,iSMS_index)
%this function is to reconstruct SMS aliased images by SENSE method
%author: Bo Li 4-28-2021
%input
%coilSen: coilSen maps for all single slices. [nx,ny,nslc_ori,nchn]
%KSpaceDATA: SMS kspace. [nx,ny,nchn,nslc]
%sliceOrderSMS: the slice acquisition order
%iSMS_index: the index of SMS acquisition. 0 means for all acquisitions
    nx=size(KSpaceDATA,1);ny=size(KSpaceDATA,2);
    recon_slice=zeros(nx,ny,prot.OriNslice, 'double');
    coilSen_SingleSlice=zeros(size(coilSen,1),size(coilSen,2),size(coilSen,3),prot.lMultiBandFactor);
    
    if nargin<=5
        iSMS_index=0;
    end
    
    for iSMS=1:prot.OriNslice / prot.lMultiBandFactor 
        if iSMS_index~=0
            if iSMS~=iSMS_index
                continue;
            end
        end
        if size(coilSen,5)>1
            coilSen_SingleSlice(:,:,:,1:end)=coilSen(:,:,:,sliceOrderSMS((iSMS-1)*prot.lMultiBandFactor+1:iSMS*prot.lMultiBandFactor),iSMS);
        else
            coilSen_SingleSlice(:,:,:,1:end)=coilSen(:,:,:,sliceOrderSMS((iSMS-1)*prot.lMultiBandFactor+1:iSMS*prot.lMultiBandFactor));
        end
        
        if prot.lMultiBandFactor==2
            coilSen_SingleSlice(:,:,:,2)=(circshift(coilSen_SingleSlice(:,:,:,2),[0,-ny/4,0]));
        elseif prot.lMultiBandFactor==3
            %CAIPIshifts= [0,pi*2/3,pi*4/3];
            coilSen_SingleSlice(:,:,:,2)=(circshift(coilSen_SingleSlice(:,:,:,2),[0,round(-ny/3),0]));
            coilSen_SingleSlice(:,:,:,3)=(circshift(coilSen_SingleSlice(:,:,:,3),[0,round(ny/3),0]));
        elseif prot.lMultiBandFactor==4
            coilSen_SingleSlice(:,:,:,2)=(circshift(coilSen_SingleSlice(:,:,:,2),[0,round(-ny/3),0]));
            coilSen_SingleSlice(:,:,:,3)=(circshift(coilSen_SingleSlice(:,:,:,3),[0,round(ny/3),0]));
            coilSen_SingleSlice(:,:,:,4)=(circshift(coilSen_SingleSlice(:,:,:,4),[0,0,0]));
        elseif prot.lMultiBandFactor==5
            coilSen_SingleSlice(:,:,:,2)=(circshift(coilSen_SingleSlice(:,:,:,2),[0,ceil(-ny/3),0]));
            coilSen_SingleSlice(:,:,:,3)=(circshift(coilSen_SingleSlice(:,:,:,3),[0,ceil(ny/3),0]));
            coilSen_SingleSlice(:,:,:,4)=(circshift(coilSen_SingleSlice(:,:,:,4),[0,0,0]));
            coilSen_SingleSlice(:,:,:,5)=(circshift(coilSen_SingleSlice(:,:,:,5),[0,ceil(-ny/3),0]));
        end
        kspace_sms=KSpaceDATA(:,:,:,iSMS);
        
        imgIndex=sliceOrderSMS((iSMS-1)*prot.lMultiBandFactor+1:iSMS*prot.lMultiBandFactor);
        
        %SMS recon
        recon_slice_1=SENSE_SMS((ifft2c(kspace_sms)),coilSen_SingleSlice,prot.lMultiBandFactor);
        
        if prot.lMultiBandFactor==2
            recon_slice_1(:,:,2)=(circshift(recon_slice_1(:,:,2),[0,ny/4,0]));
        elseif prot.lMultiBandFactor==3
            recon_slice_1(:,:,2)=(circshift(recon_slice_1(:,:,2),[0,round(ny/3),0]));
            recon_slice_1(:,:,3)=(circshift(recon_slice_1(:,:,3),[0,round(-ny/3),0]));
        elseif prot.lMultiBandFactor==4
            recon_slice_1(:,:,2)=(circshift(recon_slice_1(:,:,2),[0,round(ny/3),0]));
            recon_slice_1(:,:,3)=(circshift(recon_slice_1(:,:,3),[0,round(-ny/3),0]));
            recon_slice_1(:,:,4)=(circshift(recon_slice_1(:,:,4),[0,0,0]));
        end
        recon_slice(:,:,imgIndex)=recon_slice_1(:,:,:);
    end
end