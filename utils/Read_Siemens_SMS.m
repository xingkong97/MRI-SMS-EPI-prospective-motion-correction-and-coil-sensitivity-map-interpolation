function [prot,out]=Read_Siemens_SMS(dname, RepIndex)
%read Siemens SMS raw data
%input: dname is the raw data
%       RepIndex is the Repetition index
%       
%output: prot: imaging parameters
%        out : phase corrected k-space data
%
%Written by: Bo Li, University of Maryland, Baltimore
%Created on Dec. 20, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    twix_obj=mapVBVD(dname);
    
    %IsMotion=1;

    % Is the data single-RAID or multi-RAID?
    % struct - single-RAID
    % cell - multi-RAID, with info in the last cell element
    if isstruct(twix_obj)
        disp('Loading single-RAID file...' )
    elseif iscell(twix_obj)
        disp('Loading multi-RAID file...')
        twix_obj = twix_obj{end};
    end

    %% read in some useful header info
    prot.SiemensVersion       = twix_obj.image.softwareVersion; % Siemens software version (VA,VB,VC,VD,VE?)
    prot.sequenceString       = twix_obj.hdr.Config.SequenceString; % Short sequence name
    prot.mrprot               = twix_obj.hdr.MeasYaps;
    prot.sqzSize              = twix_obj.image.sqzSize; % dimensions (data points, averages, number of coils, dynamics (ON and OFF))
    prot.sqzDims              = twix_obj.image.sqzDims; % variable names for dimensions
    prot.dPhaseFOV            = twix_obj.hdr.Config.PhaseFoV;
    prot.dReadoutFOV          = twix_obj.hdr.Config.ReadFoV;
    
    prot.Nread                = twix_obj.hdr.Config.NImageCols;
    prot.Nphase               = twix_obj.hdr.Config.NImageLins;
    prot.NRep                  = twix_obj.image.NRep; %how many measurements
    prot.lAccelFactPE        = twix_obj.hdr.Dicom.lAccelFactPE;
    prot.BaseResolution     =twix_obj.hdr.Config.BaseResolution;
    
    %slice 1 pos: prot.slicePos{1,1}.sPosition, slice 2 pos: prot.slicePos{1,2}.sPosition
    prot.slicePos           = twix_obj.hdr.Phoenix.sSliceArray.asSlice;%slice position

    prot.sliceThickness     = prot.slicePos{1,1}.dThickness;
    Tra_slc1=prot.slicePos{1,1}.sPosition.dTra;
    Tra_slc2=prot.slicePos{1,2}.sPosition.dTra;
    prot.sliceDist          = uint16(100*(abs(Tra_slc1-Tra_slc2)-prot.sliceThickness)/prot.sliceThickness);
    prot.sliceCetr2Cetr     = abs(Tra_slc1-Tra_slc2);%neiboring slices center to center distance
    
    prot.OriginPositionZ    = twix_obj.hdr.Config.SBCSOriginPositionZ;
    
    prot.TR=twix_obj.hdr.Config.TR;
    prot.ScanTime=twix_obj.hdr.Config.TotalScanTimeFromUIinms; %second
    
    %twix_obj.hdr.Config.NSlc=twix_obj.hdr.Config.NSlc+1;
    %the number of SMS slices plus single slice 
    prot.Nslice               = twix_obj.hdr.Config.NSlc;
    prot.NAve                 = twix_obj.hdr.Config.NAve;
    lMultiBandFactor=twix_obj.hdr.MeasYaps.sSliceAcceleration.lMultiBandFactor;
    %the number of slices
    prot.lMultiBandFactor=lMultiBandFactor;
    if lMultiBandFactor>1 
        prot.OriNslice=twix_obj.hdr.Config.OriNSlc;
    else
        prot.OriNslice=prot.Nslice;
    end
    
    prot.SliceFOV   = prot.OriNslice*prot.sliceThickness + prot.sliceDist/100*prot.sliceThickness*(prot.OriNslice-1);
    
    slicePos_Matrix=cell(prot.OriNslice,1);
    %slice center point coordinate: [img_pos_x,img_pos_y,img_pos_z]
    %for transversal slice
    for islc=1:prot.OriNslice
        if isfield(prot.slicePos{1,islc},'sPosition')==1
            if isfield(prot.slicePos{1,islc}.sPosition,'dSag')==1
                img_pos_x = prot.slicePos{1,islc}.sPosition.dSag;
            else
                img_pos_x = 0;
            end
            if isfield(prot.slicePos{1,islc}.sPosition,'dCor')==1
                img_pos_y = prot.slicePos{1,islc}.sPosition.dCor;
            else
                img_pos_y = 0;
            end

            img_pos_z = prot.slicePos{1,islc}.sPosition.dTra;
        else
            img_pos_x = 0;
            img_pos_y = 0;
            img_pos_z = 0;
        end
        lx = 0.5*linspace(-1,1,prot.Nread+1).*prot.dReadoutFOV+1;lx=lx(1:end-1);
        ly = 0.5*linspace(-1,1,prot.Nphase+1).*prot.dPhaseFOV +1;ly=ly(1:end-1);
        %u(1,1)=1;u(2,1)=0;u(3,1)=0;
        %v(1,1)=0;v(2,1)=1;v(3,1)=0;
        %[X,Y] = ndgrid(lx,ly);
        %slice all points coordinates: [img_pos_X_Cor,img_pos_Y_Cor,img_pos_Z_Cor]
        img_pos_X_Cor = img_pos_x + lx;%*u(1) + Y*v(1);
        img_pos_Y_Cor = img_pos_y + ly;%*u(2) + Y*v(2);
        img_pos_Z_Cor = img_pos_z;
        tmp_PosMatrix_info=struct('img_pos_X_Cor', img_pos_X_Cor,...
                                  'img_pos_Y_Cor', img_pos_Y_Cor,...
                                  'img_pos_Z_Cor', img_pos_Z_Cor);
        
        slicePos_Matrix{islc,1}=tmp_PosMatrix_info;
    end
    prot.slicePos_Matrix=slicePos_Matrix;
    
    prot.rampup=twix_obj.hdr.Meas.alRegridRampupTime(1);
    prot.rampdown=twix_obj.hdr.Meas.alRegridRampdownTime(1);
    prot.ramptop=twix_obj.hdr.Meas.alRegridFlattopTime(1);
    prot.sampdelay=twix_obj.hdr.Meas.alRegridDelaySamplesTime(1);
    prot.iNoOfFourierColumns=twix_obj.hdr.Meas.iNoOfFourierColumns;
    prot.lAccelFactPE=twix_obj.hdr.Meas.lAccelFactPE;
    prot.iNoOfFourierLines=twix_obj.hdr.Meas.iNoOfFourierLines;
    if prot.lAccelFactPE>1
        prot.CenterLineNo=twix_obj.refscan.dataSize(3);
        prot.CenterLineIndex=twix_obj.refscan.Lin(1):twix_obj.refscan.Lin(prot.CenterLineNo);
    end
    prot.lPhaseEncodingLines=twix_obj.hdr.Meas.lPhaseEncodingLines;
    
    isrm_ghost_NZ=0;
    if isfield(twix_obj.hdr.Config,'ADCDuration')
        prot.adcduration=twix_obj.hdr.Config.ADCDuration;
        isrm_ghost_NZ=0;
    end
    %isrm_ghost_NZ=1;
    
    prot.chronSliceIndices=twix_obj.hdr.Config.chronSliceIndices(1:prot.OriNslice)+1;
    
    prot.iNoOfFourierLines    = twix_obj.hdr.Config.NoOfFourierLines;
    prot.chn                  = prot.sqzSize(2) ; %how many channels
    prot.seq                  = twix_obj.hdr.Config.SequenceFileName; % Full sequence name
    % if IPAT is on, read in grappy info
    if isfield(twix_obj.hdr.Config,'GrappaPSize') %if the grappa is turned on
        prot.grappa               = twix_obj.hdr.Config.GrappaPSize/prot.lMultiBandFactor; 
        prot.sliceACC             = twix_obj.hdr.Config.GrappaSegments; % slice acceleration number?
        prot.lFirstFourierLine    = twix_obj.hdr.Config.NFirstLin;% N.L. ?
        prot.lFirstRefLine        = twix_obj.hdr.Config.NFirstRefLin; % 35
    else
        prot.grappa=0;
    end
    
    twix_obj.image.flagDoAverage = true;
    twix_obj.image.flagRemoveOS = true;
    twix_obj.image.flagIgnoreSeg = true;

    %% if IPAT is off
    imagekData = twix_obj.image();
    noisekData = permute((twix_obj.noise()),[2,1]);
    %Noise covariance matrix
    NSamples=size(noisekData,2);
    NoisePsi = (1/(NSamples-1))*(noisekData*noisekData');
    out.NoisePsi=NoisePsi;
    %figure;imagesc(abs(NoisePsi)); axis equal; axis off; colormap(jet);
    
    if lMultiBandFactor>1 
        twix_obj.sliceimgSMS.flagDoAverage = twix_obj.image.flagDoAverage;
        twix_obj.sliceimgSMS.flagRemoveOS = twix_obj.image.flagRemoveOS;
        twix_obj.sliceimgSMS.flagIgnoreSeg = twix_obj.image.flagIgnoreSeg;

        twix_obj.phasecorSMS.flagRemoveOS = twix_obj.image.flagRemoveOS;
    end

    %twix_obj.phasecor.flagDoAverage = twix_obj.image.flagDoAverage;
    twix_obj.phasecor.flagRemoveOS = twix_obj.image.flagRemoveOS;

    twix_obj.refscan.flagRemoveOS = twix_obj.image.flagRemoveOS;
    twix_obj.refscan.flagIgnoreSeg = twix_obj.image.flagIgnoreSeg;

    twix_obj.refscanPC.flagRemoveOS = twix_obj.image.flagRemoveOS;
    
    sliceAcqOrder=twix_obj.hdr.Config.chronSliceIndices(1:prot.OriNslice)+1;
    
    %imaging order for SMS
    sliceOrderSMS_tmp=calSliceOrderSMS(prot.OriNslice,lMultiBandFactor);
    if mod(prot.OriNslice,2*lMultiBandFactor)==0 %if even, then 2,4,6,...,N,1,3,5,...,N-1. if odd, then 1,3,5,...,N,2,4,6,...,N-1.
        if prot.OriNslice==2
            sliceOrderSMS = sliceAcqOrder;
        else
            tt=prot.OriNslice/2/prot.lMultiBandFactor;
            for k=2:-1:1
                slice_even=k:2:prot.OriNslice;
                slice_used=slice_even;
                
                for i=2:prot.lMultiBandFactor
                    slice_tmp=circshift(slice_even,-(i-1)*tt);
                    slice_used=[slice_used;slice_tmp];
                    
                end
                slice_used=slice_used(1:prot.lMultiBandFactor,1:tt);
                if k==2
                    sliceOrderSMS=slice_used(:)';
                else
                    sliceOrderSMS=[sliceOrderSMS,slice_used(:)'];
                end
            end
        end
    else
        %for lMultiBandFactor=2
        if mod(prot.OriNslice,lMultiBandFactor)==0 
            %if mod(idivide(prot.OriNslice,int32(2*lMultiBandFactor), 'fix'),2)==0
            if mod(lMultiBandFactor,2)==0
                loopk=2:-1:1;
            else
                loopk=1:1:2;
            end
            
            slice_half=prot.OriNslice/lMultiBandFactor;
                for k=loopk %1:1:2
                    slice_odd=k:2:prot.OriNslice/lMultiBandFactor;                    
                    slice_odd_1=zeros(lMultiBandFactor-1,size(slice_odd,2),class(slice_odd));
                    slice_odd_1(1,:)=slice_odd+slice_half;
                    slice_used=cat(1,slice_odd,slice_odd_1(1,:));
                    for mk=2:lMultiBandFactor-1
                        slice_odd_1(mk,:)=slice_odd_1(mk-1,:)+slice_half;
                        slice_used=cat(1,slice_used,slice_odd_1(mk,:));
                    end
                    
                    if k==loopk(1)
                        sliceOrderSMS=slice_used(:)';
                    else
                        sliceOrderSMS=[sliceOrderSMS,slice_used(:)'];
                    end
                end
        end
    end

    out.sliceOrderSMS=sliceOrderSMS;%[1, 6, 11,3,8,13, 5,10,15, 2, 7, 12, 4,9,14];%sliceOrderSMS;%[2,7,4,9,1,6,3,8,5,10];%
    
    %kdata_sliceimg=twix_obj.sliceimgSMS();
    if lMultiBandFactor>1  
        kdata_sliceimg_tmp=twix_obj.sliceimgSMS();%single slice k-space data in SMS scan
    end
    if prot.grappa~=0
        if lMultiBandFactor>1 
            kdata_sliceimg_tmp=padarray(kdata_sliceimg_tmp,[0 0 prot.lAccelFactPE-1],'pre');
        end
        refscanPC_tmp=squeeze(twix_obj.refscanPC());
        refscanPC=zeros(size(refscanPC_tmp,1), size(refscanPC_tmp,2), size(refscanPC_tmp,3), 3, class(refscanPC_tmp));
        refscanPC(:,:,sliceAcqOrder(1:prot.OriNslice),1)=refscanPC_tmp(:,:,1:prot.OriNslice,1,2);
        refscanPC(:,:,sliceAcqOrder(1:prot.OriNslice),2)=refscanPC_tmp(:,:,1:prot.OriNslice,1,1);
        refscanPC(:,:,sliceAcqOrder(1:prot.OriNslice),3)=refscanPC_tmp(:,:,1:prot.OriNslice,2,2);
%         clear refscanPC_tmp
%         refscan_tmp=squeeze(twix_obj.refscan());
%         refscan(:,:,:,sliceAcqOrder(1:prot.OriNslice))=refscan_tmp(:,:,:,1:prot.OriNslice);

%         refscanPC(:,:,:,1)=refscanPC_tmp(:,:,:,1,2);
%         refscanPC(:,:,:,2)=refscanPC_tmp(:,:,:,1,1);
%         refscanPC(:,:,:,3)=refscanPC_tmp(:,:,:,2,2);
        clear refscanPC_tmp
        refscan_tmp=squeeze(twix_obj.refscan());
        refscan=zeros(size(refscan_tmp));
        refscan(:,:,:,sliceAcqOrder(1:prot.OriNslice))=refscan_tmp(:,:,:,1:prot.OriNslice);
        refscan=rm_ghost_NZorWZ(isrm_ghost_NZ,refscan,prot,prot.OriNslice,refscanPC,prot.lAccelFactPE);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         load('kdata.mat');
%         clear refscan
%         refscan=zeros(size(refscan_tmp));
%         kdata=permute(kdata,[1,3,2,4]);
%         refscan(:,:,:,sliceAcqOrder(1:prot.OriNslice))=kdata(:,:,prot.CenterLineIndex,:);
%         refscan=permute(refscan,[1,3,2,4]);
%         clear kdata
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         if isrm_ghost_NZ==1
%             NZ(refscan,prot.OriNslice,prot.chn);    
%         else
%             refscan=rm_ghostWZ(refscan,prot.OriNslice,prot,refscanPC,1);
%         end
        clear refscan_tmp
    end
    if lMultiBandFactor>1
        kdata_sliceimg=zeros(size(kdata_sliceimg_tmp));
        kdata_sliceimg(:,:,:,:,sliceAcqOrder(1:prot.OriNslice))=kdata_sliceimg_tmp(:,:,:,:,1:prot.OriNslice);
    end
    %clear kdata_sliceimg_tmp
    if lMultiBandFactor>1 
        phasecorSMS_tmp=squeeze(twix_obj.phasecorSMS());
        %phasecorSMS is phasecor for single slice sliceimgSMS
        phasecorSMS=zeros(size(phasecorSMS_tmp,1), size(phasecorSMS_tmp,2), size(phasecorSMS_tmp,3), 3, class(phasecorSMS_tmp));
        phasecorSMS(:,:,sliceAcqOrder(1:prot.OriNslice),1)=phasecorSMS_tmp(:,:,1:prot.OriNslice,1,2);
        phasecorSMS(:,:,sliceAcqOrder(1:prot.OriNslice),2)=phasecorSMS_tmp(:,:,1:prot.OriNslice,1,1);
        phasecorSMS(:,:,sliceAcqOrder(1:prot.OriNslice),3)=phasecorSMS_tmp(:,:,1:prot.OriNslice,2,2);
    end
    clear phasecorSMS_tmp

    phasecor_tmp=squeeze(twix_obj.phasecor());
    phasecor=zeros(size(phasecor_tmp,1), size(phasecor_tmp,2), prot.Nslice, 3, class(phasecor_tmp));
    if prot.NRep==1
        phasecor(:,:,:,1)=phasecor_tmp(:,:,end,:,1,2);%the 3rd index has to be set "end". It means which phase was acquired in phase correction. 
        phasecor(:,:,:,2)=phasecor_tmp(:,:,end,:,1,1);
        phasecor(:,:,:,3)=phasecor_tmp(:,:,end,:,2,2);
    else
        %phasecor=zeros(size(phasecor_tmp,1), size(phasecor_tmp,2), size(phasecor_tmp,4), size(phasecor_tmp,6), 3, class(phasecor_tmp));
        if prot.Nslice==1
            phasecor(:,:,1,1)=phasecor_tmp(:,:,end,1,RepIndex,2); 
            phasecor(:,:,1,2)=phasecor_tmp(:,:,end,1,RepIndex,1);
            phasecor(:,:,1,3)=phasecor_tmp(:,:,end,2,RepIndex,2);
        else
            phasecor(:,:,:,1)=phasecor_tmp(:,:,end,:,1,RepIndex,2); 
            phasecor(:,:,:,2)=phasecor_tmp(:,:,end,:,1,RepIndex,1);
            phasecor(:,:,:,3)=phasecor_tmp(:,:,end,:,2,RepIndex,2);
        end
    end
    
    clear phasecor_tmp

    %% read motion timestamps stored in each k-space line
    realres=size(imagekData, 3)/prot.lAccelFactPE;
    %the No.5 iceparam save the high 16 bit of timestamps
    ulHigh16=twix_obj.image.iceParam(5, :);
    %the No.6 iceparam save the low 16 bit of timestamps
    ulLow16=twix_obj.image.iceParam(6, :);
    %get all timestamps for whole scan
    out.mcTimeStamp_all=bitor(bitshift(ulHigh16,16),ulLow16)';
    
    %get timestamps only for one repetition
    ulHigh16=twix_obj.image.iceParam(5, realres*prot.Nslice*(RepIndex-1)+1:realres*prot.Nslice*RepIndex);
    ulLow16=twix_obj.image.iceParam(6, realres*prot.Nslice*(RepIndex-1)+1:realres*prot.Nslice*RepIndex);
    out.mcTimeStamp=bitor(bitshift(ulHigh16,16),ulLow16)';

    if prot.grappa==0 
        if prot.NRep==1
            kdatatmp=squeeze(imagekData(:,:,:,1,:));
            kdata=rm_ghost_NZorWZ(isrm_ghost_NZ,kdatatmp,prot,prot.Nslice,phasecor);
        else
            imagekData=squeeze(imagekData);
            if prot.Nslice==1
                imagekData=permute(imagekData,[1,2,3,5,4]);    
            end
            kdatatmp=squeeze(imagekData(:,:,:,:,RepIndex));
            
           kdata=rm_ghost_NZorWZ(isrm_ghost_NZ,kdatatmp,prot,prot.Nslice,phasecor);
           
           %add additional phase to phase corrected k-space data
%            if prot.lMultiBandFactor==2
%                FOV_factor=2.0;
%                FOV_sms=double(FOV_factor*(out.sliceOrderSMS(1,2)-out.sliceOrderSMS(1,1))*prot.sliceThickness*(1+prot.sliceDist/100));
%                delta=360.0/FOV_sms; %degree per mm 
%                FI2=[0,delta,2.0*delta,-delta]/180;%pi*1/2; [0,-27,-54,27]/180
%                FI2=repmat(FI2, prot.Nread, prot.Nphase/4 * prot.Nslice);
%                tmpmcMatrix=repmat(out.mcMatrix(:, 4)', prot.Nread, 1);
%                FI2=FI2.*(tmpmcMatrix);%for motion injection, FI2=FI2.*(-tmpmcMatrix);
%                FI2=exp(-1i*FI2*pi);%FI2*pi;
%                FI2=reshape(FI2, prot.Nread, prot.Nphase, prot.Nslice);
%                FI2=repmat(FI2, 1, 1, 1, prot.chn);
%                FI2=permute(FI2, [1,2,4,3]);
%                kdata=kdata.*FI2;
%            end
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              
        end
        %for single slice kspace data
        if lMultiBandFactor>1 
            kdata_sliceimg=squeeze(kdata_sliceimg);
            kdata_sliceimg=rm_ghost_NZorWZ(isrm_ghost_NZ,kdata_sliceimg,prot,prot.OriNslice,phasecorSMS);
        end
        
    else %if IPAT is on
        %twix_obj.refscan() returns all reference lines (including the lines that
        %are also read in by twix_obj.image(). fill in the reference
        %scan at the appropriate position (duplicate imaging lines are overwritten
        
        isSliceGrappa=0;
        
        if lMultiBandFactor>1 
            kdata_sliceimg=squeeze(kdata_sliceimg);
            kdata_sliceimg=rm_ghost_NZorWZ(isrm_ghost_NZ,kdata_sliceimg,prot,prot.OriNslice,phasecorSMS);
            
            if size(kdata_sliceimg,2)<size(kdata_sliceimg,1)
                kdata_sliceimg_tmp=zeros(size(kdata_sliceimg,1), size(kdata_sliceimg,1), size(kdata_sliceimg,3), size(kdata_sliceimg,4));
                kdata_sliceimg_tmp(:,1:size(kdata_sliceimg,2),:,:)=kdata_sliceimg;
                kdata_sliceimg=kdata_sliceimg_tmp;
            end
        end
        
        if lMultiBandFactor==1
            imagekData=padarray(imagekData,[0 0 prot.lAccelFactPE],'post');
        end
        
        if prot.NRep==1
            kdatatmp=squeeze(imagekData(:,:,:,1,:));
            kdata=rm_ghost_NZorWZ(isrm_ghost_NZ,kdatatmp,prot,prot.Nslice,phasecor);
        else
            kdata=squeeze(imagekData);
            kdata=squeeze(kdata(:,:,:,:,RepIndex));
            phasecor_tmp=squeeze(phasecor(:,:,:,:));
            kdata=rm_ghost_NZorWZ(isrm_ghost_NZ,kdata,prot,prot.Nslice,phasecor_tmp);
            
            if size(kdata,2)<size(kdata,1)
                kdata_tmp=zeros(size(kdata,1), size(kdata,1), size(kdata,3), size(kdata,4));
                kdata_tmp(:,1:size(kdata,2),:,:)=kdata;
                kdata=kdata_tmp;
            end 
        end

        %the image of grappa recon for under-sampling dataset prior to SMS scan
        %is not applicable. The result is still aliasing, meaning refscan is
        %not for the dataset. refscan is just for grappa scan after SMS. 
        %if isSliceGrappa set 1, then do grappa for all slices acquisition
        %prior to SMS reconstruction.
        if isSliceGrappa
            %sliceimg=DoGrappaWZ(kdata_sliceimg, refscan, prot);
            sliceimg=zeros(size(kdata,1),size(kdata,2),prot.OriNslice);%nx,ny,nslice
            for islc=1:prot.OriNslice
                kdata_sliceimg_tmp_1=squeeze(kdata_sliceimg(:,:,:,islc));%nx,ny,ncoil,nslice
                refdata=squeeze(refscan(:,:,:,islc));
                kdata_sliceimg_tmp_1=DoGrappa(kdata_sliceimg_tmp_1,refdata);
                sliceimg(:,:,islc)=sqrtSum(kdata_sliceimg_tmp_1);
            end
        end
    end
    
 %%%%%%%%%%%%%%%%%%%%%%return Parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    out.kdata_ori=squeeze(imagekData(:,:,:,:,RepIndex));
    out.kdata=kdata; %nx,ny,ncoil,nslc
    if lMultiBandFactor>1
        %kdata_sliceimg is the single-slice reference data acquired in
        %pre-scan
        out.kdata_sliceimg=kdata_sliceimg;%nx,ny,ncoil,orinslc
        %SMS phase correction three lines
        out.phasecorSMS=phasecorSMS;
    end
    %phase correction three lines
    out.phasecor=phasecor;
    if prot.grappa~=0
        out.refscan=refscan;
        out.refscanPC=refscanPC;
    end
end

function res=rm_ghost_NZorWZ(rm_ghost_mode,kdata,prot,nslc,phasecor,lAccelFactPE)
    if nargin<=3
        res=rm_ghostNZ(kdata,nslc,prot.chn);
        return;
    end

    if rm_ghost_mode==1
        res=rm_ghostNZ(kdata,nslc,prot.chn);    
    else
        if nargin<=5
            res=rm_ghostWZ(kdata,nslc, prot,phasecor);
        else
            res=rm_ghostWZ(kdata,nslc, prot,phasecor,lAccelFactPE);
        end
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

function res=rm_ghostNZ(kdata,nslice,chn)
%kdata: nx,ncoil,ny,nslc
        res=permute(kdata,[1,3,4,2]);
        res=rm_ghost_NZ(res,nslice,chn);
        res=permute(res,[1,2,4,3]);
end

function res=DoGrappa(kdata, refdata)%kdata:nx,ny,ncoil,refdata:nx,ny,ncoil
    kdata=permute(kdata,[2,1,3]);
    R=[1,2];%the undersampling rate of readout(nx) and phase(ny) directions
    kernel=[3,2];%one is odd and one is even
    kcabliData=permute(refdata,[2,1,3]);
    res=GrappaRecon(kdata, kcabliData, kernel, R);%;%GRAPPA(kdata,kcabliData,[5,5],0.05)
    res=permute(res,[2,1,3]);
end

function res=DoGrappaWZ(kdata, refdata,prot)%kdata:nx,ny,ncoil,nslice,refdata:nx,ny,ncoil,nslice
    kdata=permute(kdata,[3,4,2,1]);%ncoil, nslice, ny,nx
    R=2;
    kcabliData=permute(refdata,[3,4,2,1]);
    res=GrappaWZ(kdata, kcabliData, prot, R);%ncoil,ny,nx
    %res=permute(res,[3,2,1]);%nx,ny,ncoil
end

function showimagefft(kspace, fftparam) %kspace:nx,ny,ncoil,nslice
    
    if nargin<=1
        isfft=1;
    else
        if fftparam==0
            isfft=0;
        else
            isfft=1;
        end
    end
    NSlice=size(kspace,4);
    if NSlice==1
        imagesc(flip(sqrtSum(kspace(:,:,:),isfft)'));colormap gray;
        return;
    end
    
    if mod(NSlice,4)==0
        nrow=NSlice/4;
    elseif mod(NSlice,3)==0
        nrow=NSlice/3;
    elseif mod(NSlice,5)==0
        nrow=NSlice/5;
    elseif mod(NSlice,2)==0
        nrow=NSlice/2;
    end
    ncol=NSlice/nrow;%ncoil nslice
    ha = tight_subplot(nrow,ncol,   [   .07         .005    ],  [ 0.05  .015], [ 0.002   .002]);
    %                  row, col, [gap_height, gap_width],  [bottom  top], [left    right]
    %figure;
    for irow=1:nrow
        for icol=1:ncol
            iIndex=(irow-1)*ncol+icol;
            axes(ha(iIndex));
            imagesc(flip(sqrtSum(kspace(:,:,:,iIndex),isfft)'));colormap gray;%
            %imagesc((sqrtSum(kspace(:,:,:,iIndex),isfft)));colormap gray;%
            %imagesc(abs((kspace(:,:,:,iIndex))));colormap gray;%
        end
    end
end
