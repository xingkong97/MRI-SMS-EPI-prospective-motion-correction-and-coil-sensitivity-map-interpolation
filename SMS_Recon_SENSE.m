function recon_slice=SMS_Recon_SENSE(coilSen, prot, out)
%SMS SENSE recon
%INPUT: coilSen, coil sensitivity map
%       prot, image parameters
%       out, single-slice reference k-space data and SMS collapsed data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT: recon_slice, reconstruction slice images
%
%%%%%% Written by: Bo Li, University of Maryland, Baltimore
%%%%%% for manuscript "SMS-EPI prospective motion correction 
%%%%%% by real-time phase compensation and coil sensitivity map interpolation"
%%%%%% Created on Sep. 22, 2022

    if prot.lMultiBandFactor==1
        KSpaceDATA=out.kdata;
        if prot.grappa~=0 %if IPAT is on
            recon_slice=zeros(size(out.kdata,1), size(out.kdata,2), size(out.kdata,4));
            sliceOrderSMS=out.sliceOrderSMS;
            for islc=1:prot.OriNslice
                res=DoGrappa(out.kdata(:,:,:,islc),out.refscan(:,:,:,sliceOrderSMS(islc)));
                recon_slice(:,:,sliceOrderSMS(islc))=sqrtSum(res);
            end
            recon_slice=permute(recon_slice,[1,2,4,3]);
            showimagefft(recon_slice,0);
        else
            recon_slice=zeros(size(KSpaceDATA,1), size(KSpaceDATA,2), size(KSpaceDATA,4));
            sliceOrderSMS=out.sliceOrderSMS;

            for islc=1:prot.Nslice
                recon_slice(:,:,sliceOrderSMS(islc))=sqrtSum(KSpaceDATA(:,:,:,islc));
            end
            recon_slice=permute(recon_slice,[1,2,4,3]);
            showimagefft(recon_slice,0); %x,y,chn,slice   
        end
    else
        KSpaceDATA=out.kdata;
        
        recon_slice=DoSENSE_SMS(coilSen, KSpaceDATA, out.sliceOrderSMS, prot);
        recon_slice=permute(recon_slice,[1,2,4,3]);
        
    end
end

function recon_slice=DoSENSE_SMS(coilSen, KSpaceDATA, sliceOrderSMS, prot)%
%SMS-SENSE reconstruction
% coilSen: coil sensitivity for all single slices. [nx,ny,nch,nslc]
% KSpaceDATA: SMS kspace. [nx,ny,nchn,nslc]
    nx=size(KSpaceDATA,1);ny=size(KSpaceDATA,2);
    recon_slice=zeros(nx,ny,prot.OriNslice, 'double');
    coilSen_SingleSlice=zeros(size(coilSen,1),size(coilSen,2),size(coilSen,3),prot.lMultiBandFactor);
   
    
    for iSMS=1:prot.OriNslice / prot.lMultiBandFactor 
        if size(coilSen,5)>1
            coilSen_SingleSlice(:,:,:,1:end)=coilSen(:,:,:,sliceOrderSMS((iSMS-1)*prot.lMultiBandFactor+1:iSMS*prot.lMultiBandFactor),iSMS);
        else
            coilSen_SingleSlice(:,:,:,1:end)=coilSen(:,:,:,sliceOrderSMS((iSMS-1)*prot.lMultiBandFactor+1:iSMS*prot.lMultiBandFactor));
        end
        if prot.lMultiBandFactor==2
            coilSen_SingleSlice(:,:,:,2)=(circshift(coilSen_SingleSlice(:,:,:,2),[0,-ny/2,0]));
        elseif prot.lMultiBandFactor==3
            %CAIPIshifts= [0,pi*2/3,pi*4/3];
              coilSen_SingleSlice(:,:,:,2)=(circshift(coilSen_SingleSlice(:,:,:,2),[0,-ny/2,0]));
              coilSen_SingleSlice(:,:,:,3)=(circshift(coilSen_SingleSlice(:,:,:,3),[0,0,0]));
        elseif prot.lMultiBandFactor==4
            %CAIPIshifts= [0,pi*2/3,pi*4/3];
            coilSen_SingleSlice(:,:,:,2)=(circshift(coilSen_SingleSlice(:,:,:,2),[0,round(-ny/3),0]));
            coilSen_SingleSlice(:,:,:,3)=(circshift(coilSen_SingleSlice(:,:,:,3),[0,round(ny/3),0]));
            coilSen_SingleSlice(:,:,:,4)=(circshift(coilSen_SingleSlice(:,:,:,4),[0,0,0]));
        elseif prot.lMultiBandFactor==5
            %CAIPIshifts= [0,pi*2/3,pi*4/3];
            coilSen_SingleSlice(:,:,:,2)=(circshift(coilSen_SingleSlice(:,:,:,2),[0,ceil(-ny/3),0]));
            coilSen_SingleSlice(:,:,:,3)=(circshift(coilSen_SingleSlice(:,:,:,3),[0,ceil(ny/3),0]));
            coilSen_SingleSlice(:,:,:,4)=(circshift(coilSen_SingleSlice(:,:,:,4),[0,0,0]));
            coilSen_SingleSlice(:,:,:,5)=(circshift(coilSen_SingleSlice(:,:,:,5),[0,ceil(-ny/3),0]));
        elseif prot.lMultiBandFactor==6
            %CAIPIshifts= [0,pi*2/3,pi*4/3];
            coilSen_SingleSlice(:,:,:,2)=(circshift(coilSen_SingleSlice(:,:,:,2),[0,ceil(-ny/3),0]));
            coilSen_SingleSlice(:,:,:,3)=(circshift(coilSen_SingleSlice(:,:,:,3),[0,ceil(ny/3),0]));
            coilSen_SingleSlice(:,:,:,4)=(circshift(coilSen_SingleSlice(:,:,:,4),[0,0,0]));
            coilSen_SingleSlice(:,:,:,5)=(circshift(coilSen_SingleSlice(:,:,:,5),[0,ceil(-ny/3),0]));
            coilSen_SingleSlice(:,:,:,6)=(circshift(coilSen_SingleSlice(:,:,:,6),[0,ceil(ny/3),0]));
        end
        
        kspace_sms=KSpaceDATA(:,:,:,iSMS);
        imgIndex=sliceOrderSMS((iSMS-1)*prot.lMultiBandFactor+1:iSMS*prot.lMultiBandFactor);
        recon_slice_1=SENSE_SMS((ifft2c(kspace_sms)),coilSen_SingleSlice,prot.lMultiBandFactor);
        if prot.lMultiBandFactor==2
            recon_slice_1(:,:,2)=circshift(recon_slice_1(:,:,2),[0,ny/2,0]);%ny/4 (circshift(recon_slice_1(:,:,2),[0,ny/4,0]));
        elseif prot.lMultiBandFactor==3
              recon_slice_1(:,:,2)=(circshift(recon_slice_1(:,:,2),[0,round(ny/2),0]));
              recon_slice_1(:,:,3)=(circshift(recon_slice_1(:,:,3),[0,0,0]));
        elseif prot.lMultiBandFactor==4
            recon_slice_1(:,:,2)=(circshift(recon_slice_1(:,:,2),[0,round(ny/3),0]));
            recon_slice_1(:,:,3)=(circshift(recon_slice_1(:,:,3),[0,round(-ny/3),0]));
            recon_slice_1(:,:,4)=(circshift(recon_slice_1(:,:,4),[0,0,0]));
        elseif prot.lMultiBandFactor==5
            recon_slice_1(:,:,2)=(circshift(recon_slice_1(:,:,2),[0,round(ny/3),0]));
            recon_slice_1(:,:,3)=(circshift(recon_slice_1(:,:,3),[0,round(-ny/3),0]));
            recon_slice_1(:,:,4)=(circshift(recon_slice_1(:,:,4),[0,0,0]));
            recon_slice_1(:,:,5)=(circshift(recon_slice_1(:,:,5),[0,round(ny/3),0]));
        elseif prot.lMultiBandFactor==6
            recon_slice_1(:,:,2)=(circshift(recon_slice_1(:,:,2),[0,round(ny/3),0]));
            recon_slice_1(:,:,3)=(circshift(recon_slice_1(:,:,3),[0,round(-ny/3),0]));
            recon_slice_1(:,:,4)=(circshift(recon_slice_1(:,:,4),[0,0,0]));
            recon_slice_1(:,:,5)=(circshift(recon_slice_1(:,:,5),[0,round(ny/3),0]));
            recon_slice_1(:,:,6)=(circshift(recon_slice_1(:,:,6),[0,round(-ny/3),0]));
        end
        recon_slice(:,:,imgIndex)=recon_slice_1(:,:,:);
    end

end

function res=rm_ghostWZ(kdata,nslice,prot,phasecor,lAccelFactPE)
%kdata: nx,ncoil,ny,nslc phasecor:nx,ncoil,nslc,3
        res=permute(kdata,[3,1,2,4]);
        if nargin<=4
            res=rm_ghost(res,nslice,prot,phasecor);
        else
            res=rm_ghost(res,nslice,prot,phasecor,lAccelFactPE);
        end
        res=permute(res,[2,1,3,4]);
end

function res=DoGrappa(kdata, refdata)%kdata:nx,ny,ncoil,refdata:nx,ny,ncoil
    kdata=permute(kdata,[2,1,3]);
    R=[1,2];%the undersampling rate of readout(nx) and phase(ny) directions
    kernel=[3,2];%one is odd and one is even
    kcabliData=permute(refdata,[2,1,3]);
    res=GrappaRecon(kdata, kcabliData, kernel, R);
    res=permute(res,[2,1,3]);
end
