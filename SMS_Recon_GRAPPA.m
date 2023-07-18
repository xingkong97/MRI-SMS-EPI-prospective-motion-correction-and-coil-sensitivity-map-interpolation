function recon_slice=SMS_Recon_GRAPPA(coilSen,prot,out)
%Split Slice-GRAPPA and Slice-GRAPPA recon
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

    coilSen=permute(coilSen,[1,2,4,3]);
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
            showimagefft(recon_slice,0);    
        end
    else
        isSplit=1;%1 Split Slice-GRAPPA; 0 Slice-GRAPPA
        lamda=1;
        %single-slice reference data acquired in pre-scan
        kdata_sliceimg=out.kdata_sliceimg;
        %SMS collapsed data
        KSpaceDATA=out.kdata;

        %%update single slice by updated coil sensitivity
        img_csm=zeros(prot.Nread,prot.Nphase,prot.chn,prot.OriNslice);
            img_slc_chn=zeros(prot.Nread,prot.Nphase,prot.chn);
            for iSMS=1:prot.Nslice
                for islc=(iSMS-1)*prot.lMultiBandFactor+1:iSMS*prot.lMultiBandFactor
                    iIndex_slc=out.sliceOrderSMS(islc);
                    img_slc=sqrtSum((kdata_sliceimg(:,:,:,iIndex_slc)),1);
                    for ichn=1:prot.chn
                        img_slc_chn(:,:,ichn)=img_slc.*coilSen(:,:,iIndex_slc,ichn);
                        kdata_sliceimg(:,:,ichn,iIndex_slc)=fft2c(img_slc_chn(:,:,ichn));
                    end
                    img_csm(:,:,:,iIndex_slc)=img_slc_chn;
                end
            end
        %do Split Slice-GRAPPA or Slice-GRAPPA recon
        recon_slice=DoSplitSliceGrappa(kdata_sliceimg, KSpaceDATA, out.sliceOrderSMS, prot,lamda,isSplit);
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

