close all

ksize = [6,6]; % ESPIRiT kernel-window-size
eigThresh_k = 0.02; %defaukt 0.02 threshold of eigenvectors in k-space
eigThresh_im = 0.95; %default 0.95
%CSM_from_1stMeas=1;

prot.ksize=ksize;
prot.eigThresh_k=eigThresh_k;
prot.eigThresh_im=eigThresh_im;
kdata_sliceimg_sens=out.kdata_sliceimg;

IsUpdateCoilSen=0;
[tukey_window,tukey_window_red]=filter_2D(prot.Nphase);%the tukey_window size won't change recon results
IsMask=1;
IsRef=0;

mask_thresh=0.05;%0.03;coilSen_1, mask_slice
coilSen_thresh=0.08;%0.03;coilSen_matrix, mask_coilSen

IsSmooth=1;
IsSaveIntp=0;
IsNormalization=1;

x_off=0; y_off=0; z_off=0; %mm

if IsUpdateCoilSen==1
     
    isSENSE=1;
    %interpolation Approch 1: real&imaginary, 2: magnitude&phase
    InterpApp=1;
    
    %for phantom, 1:prot.OriNslice-3 will show artifacts.
    prot.InterpSlice=1:prot.OriNslice-3;
    
    tmp_MotionMatrix=zeros(6,1);
    
    if IsRef==1
        %Mat_filename=strcat(raw_data_filename(1:(strfind(raw_data_filename, '.dat')-1)), '_ref.mat');
        %coilSen_ori=CalCoilSens(permute(kdata_sliceimg_sens,[1,2,4,3]),prot.ksize,prot.eigThresh_k,prot.eigThresh_im);%nx,ny,nslc,ncoil
        Mat_filename=strcat(raw_data_filename(1:(strfind(raw_data_filename, '.dat')-1)), '_ref.mat');
        %save(Mat_filename,'coilSen_ori');
        load(Mat_filename)
       %[coilSen_matrix,coilSen_OriCor,coilSen_1,slcPos_matZ,mask_coilSen,mask_slice]=getCoilSen(mask_thresh,coilSen_thresh,tukey_window,x_off,y_off,z_off,coilSen_ori,kdata_sliceimg_sens,10,prot,IsSmooth);%prot.lMultiBandFactor+10
       %for Fig. 2
       %[coilSen_matrix,coilSen_OriCor,coilSen_1,coilSen_ori,slcPos_matZ,mask_coilSen,mask_slice]=getCoilSen_7_8(mask_thresh,coilSen_thresh,tukey_window,x_off,y_off,z_off,kdata_sliceimg_sens,prot,IsSmooth);
    else
        Mat_filename=strcat(raw_data_filename(1:(strfind(raw_data_filename, '.dat')-1)), '_1st.mat');
        SMS_Factor=prot.lMultiBandFactor;
        load(Mat_filename)
        [coilSen_matrix,coilSen_OriCor,coilSen_1,slcPos_matZ,mask_coilSen,mask_slice]=getCoilSen(mask_thresh,coilSen_thresh,tukey_window,x_off,y_off,z_off,coilSen_ori,fft2c(recon_slice_1st),SMS_Factor,prot,IsSmooth);
    end
    
    %[coilSen_matrix,coilSen_OriCor,coilSen_1,coilSen_ori,slcPos_matZ,mask_coilSen,mask_slice]=getCoilSen_7_8(mask_thresh,coilSen_thresh,tukey_window,x_off,y_off,z_off,fft2c(recon_slice_1st),prot,IsSmooth);
    
%% 
    
%     [tmp_tukey_window,tmp_tukey_window_red]=filter_2D(prot.Nphase,32);
%     tmp_kdata_sliceimg=kdata_sliceimg_sens;%.*tmp_tukey_window;
    


    
    
%     if IsNormalization==1
%         coilSen_ori=smoothCoilSen(coilSen_ori,prot);%
%             %coilSen_Update_recon=smoothCoilSen(coilSen_Update_recon,mask_thresh,tukey_window,kdata_sliceimg_sens,prot,mask_coilSen);
%     end
%     out.sliceimg_coilSen_ori=coilSen_ori;%nx,ny,nslc,nch
%     out.kdata_sliceimg_coilSen_ori=permute(fft2c(coilSen_ori),[1,2,4,3]);
    %SMS_coilSen_ori=CalCoilSens(permute(out.kdata,[1,2,4,3]),prot.ksize,prot.eigThresh_k,prot.eigThresh_im);
%     if IsNormalization==1
%         SMS_coilSen_ori=smoothCoilSen(SMS_coilSen_ori,prot);%
%             %coilSen_Update_recon=smoothCoilSen(coilSen_Update_recon,mask_thresh,tukey_window,kdata_sliceimg_sens,prot,mask_coilSen);
%     end
%     out.SMS_coilSen_ori=SMS_coilSen_ori;
%     out.kdata_SMS_coilSen_ori=permute(fft2c(SMS_coilSen_ori),[1,2,4,3]);
    
    
    %kdata_sliceimg_sens usually show unexpected artifact of coil sensitivity map
    
%     SMS_Factor=prot.lMultiBandFactor;
%     [coilSen_matrix,coilSen_OriCor,coilSen_1,slcPos_matZ,mask_coilSen,mask_slice]=getCoilSen(mask_thresh,coilSen_thresh,tukey_window,x_off,y_off,z_off,coilSen_ori,fft2c(recon_slice_1st),SMS_Factor,prot,IsSmooth);
    %[coilSen_matrix,coilSen_OriCor,coilSen_1,slcPos_matZ,mask_coilSen,mask_slice]=getCoilSen(mask_thresh,coilSen_thresh,tukey_window,x_off,y_off,z_off,coilSen_ori,fft2c(coilSen_ori),SMS_Factor,prot,IsSmooth);
%     
    coilSen_OriCor_R=coilSen_OriCor;
    
    %
%     coilSen_pred=permute(prediction,[1,2,4,3]);
%     coilSen_pred=smoothCoilSen(coilSen_pred.*mask_slice,prot);
%     coilSen_1(:,:,2:end,:)=coilSen_pred(:,:,1:end-1,:);

    coilSen_matrix_U=coilSen_ori;%coilSen_ori;%coilSen_matrix;%coilSen_1;%;coilSen_ori 
    
    %coilSen_matrix_U(27:38,38:57,16,:)=coilSen_matrix_U(27:38,38:57,16,:)*(-1);

    %results are compromised with Poly interpolation expandation
%     for islc=1:prot.OriNslice
%         for jch=1:prot.chn
%             coilSen_matrix_U(:,:,islc,jch)=CoilSenFill2DPoly(coilSen_matrix_U(:,:,islc,jch),0.02); 
%             %coilSen_matrix_U(:,:,islc,jch)=(1-mask_slice(:,:,islc,jch)).*coilSen_matrix_U(:,:,islc,jch)+mask_slice(:,:,islc,jch).*coilSen_matrix(:,:,islc,jch);
%         end
%     end
    
%     coilSen_matrix_U=smoothCoilSen(coilSen_matrix_U,prot);%.*mask_slice
    %Phase Correction
    %[coilSen_matrix_U,coilSen_sign]=PhaCor(coilSen_matrix_U,prot);
    
    SMS_MotionCell_1st=cell(prot.Nslice,1);
    SMS_MotionCell=cell(prot.Nslice,1);
    %%
    %coilSen_Update=coilSen_ori;
    
    coilSen_Update=zeros(prot.Nread,prot.Nphase,prot.OriNslice,prot.chn);
    
    for iSMS=1:prot.Nslice
        %%1st measurement translation and rotation
        x_trans_arr=zeros(prot.Nphase,1); y_trans_arr=zeros(prot.Nphase,1); z_trans_arr=zeros(prot.Nphase,1);
        x_rot_arr=zeros(prot.Nphase,1); y_rot_arr=zeros(prot.Nphase,1); z_rot_arr=zeros(prot.Nphase,1);
        for ikline=1:prot.Nphase
            iIndex=find(array_timestamp(:,1)==out.mcTimeStamp_all((iSMS-1)*prot.Nphase+ikline));
            tmp_RotMatrix=MotionData{iIndex(1,1),1}.RotMatrix;
            tmp_RotMatrix_info=struct('timestamp', array_timestamp(iIndex(1,1),1), 'RotMatrix', tmp_RotMatrix,...
                                      'transX', MotionData{iIndex(1,1),1}.transX,...
                                      'transY', MotionData{iIndex(1,1),1}.transY,...
                                      'transZ', MotionData{iIndex(1,1),1}.transZ,...
                                      'Euler_deg', MotionData{iIndex(1,1),1}.Euler_deg);
            SMS_MotionCell_1st{iSMS, 1}=tmp_RotMatrix_info;

            x_trans_arr(ikline,1)=SMS_MotionCell_1st{iSMS,1}.transX; y_trans_arr(ikline,1)=SMS_MotionCell_1st{iSMS,1}.transY; z_trans_arr(ikline,1)=SMS_MotionCell_1st{iSMS,1}.transZ;
            x_rot_arr(ikline,1)=SMS_MotionCell_1st{iSMS,1}.Euler_deg(1,1);y_rot_arr(ikline,1)=SMS_MotionCell_1st{iSMS,1}.Euler_deg(1,2);z_rot_arr(ikline,1)=SMS_MotionCell_1st{iSMS,1}.Euler_deg(1,3);
        end
        x_trans_1st=median(x_trans_arr(:)); y_trans_1st=median(y_trans_arr(:)); z_trans_1st=median(z_trans_arr(:));
        x_rot_1st=median(x_rot_arr(:)); y_rot_1st=median(y_rot_arr(:)); z_rot_1st=median(z_rot_arr(:));
        %%
        x_trans_arr=zeros(prot.Nphase,1); y_trans_arr=zeros(prot.Nphase,1); z_trans_arr=zeros(prot.Nphase,1);
        x_rot_arr=zeros(prot.Nphase,1); y_rot_arr=zeros(prot.Nphase,1); z_rot_arr=zeros(prot.Nphase,1);
        for ikline=1:prot.Nphase
            iIndex=find(array_timestamp(:,1)==out.mcTimeStamp((iSMS-1)*prot.Nphase+ikline));
            tmp_RotMatrix=MotionData{iIndex(1,1),1}.RotMatrix;
            tmp_RotMatrix_info=struct('timestamp', array_timestamp(iIndex(1,1),1), 'RotMatrix', tmp_RotMatrix,...
                                      'transX', MotionData{iIndex(1,1),1}.transX,...
                                      'transY', MotionData{iIndex(1,1),1}.transY,...
                                      'transZ', MotionData{iIndex(1,1),1}.transZ,...
                                      'Euler_deg', MotionData{iIndex(1,1),1}.Euler_deg);
            SMS_MotionCell{iSMS, 1}=tmp_RotMatrix_info;

            x_trans_arr(ikline,1)=SMS_MotionCell{iSMS,1}.transX; y_trans_arr(ikline,1)=SMS_MotionCell{iSMS,1}.transY; z_trans_arr(ikline,1)=SMS_MotionCell{iSMS,1}.transZ;
            x_rot_arr(ikline,1)=SMS_MotionCell{iSMS,1}.Euler_deg(1,1);y_rot_arr(ikline,1)=SMS_MotionCell{iSMS,1}.Euler_deg(1,2);z_rot_arr(ikline,1)=SMS_MotionCell{iSMS,1}.Euler_deg(1,3);
        end
        x_trans=median(x_trans_arr(:))-x_trans_1st; y_trans=median(y_trans_arr(:))-y_trans_1st; z_trans=median(z_trans_arr(:))-z_trans_1st;
        x_rot=median(x_rot_arr(:))-x_rot_1st; y_rot=median(y_rot_arr(:))-y_rot_1st; 
        if iSMS==3
            z_rot=median(z_rot_arr(:))-z_rot_1st;
        else
            z_rot=median(z_rot_arr(:))-z_rot_1st;
        end
%         x_trans=median(x_trans_arr(:)); y_trans=median(y_trans_arr(:)); z_trans=median(z_trans_arr(:));
%         x_rot=median(x_rot_arr(:)); y_rot=median(y_rot_arr(:)); 
%         if iSMS==3
%             z_rot=median(z_rot_arr(:));
%         else
%             z_rot=median(z_rot_arr(:));
%         end
        %%
        x_rot=x_rot; y_rot=y_rot; z_rot=-z_rot;
        x_trans=x_trans; y_trans=y_trans; z_trans=-z_trans;
        
        x_rot_matrix=[1 0 0;0 cosd(x_rot) sind(x_rot);0 -sind(x_rot) cosd(x_rot)];
        y_rot_matrix=[cosd(y_rot) 0 -sind(y_rot);0 1 0;sind(y_rot) 0 cosd(y_rot)];
        %z_rot_matrix=[cosd(-z_rot) sind(-z_rot) 0;-sind(-z_rot) cosd(-z_rot) 0;0 0 1];
        z_rot_matrix=[cosd(z_rot) sind(z_rot) 0;-sind(z_rot) cosd(z_rot) 0;0 0 1];
        %xyz_RotMatrix=zeros(3,3,6);
        f_RotMatrix=@(x,y,z)(x*y*z);
        tmp_xyz_RotMatrix=f_RotMatrix(z_rot_matrix,y_rot_matrix,x_rot_matrix);
        
        %x,y,-z translation
        tmp_MotionMatrix(1,1)=x_trans;
        tmp_MotionMatrix(2,1)=y_trans;
        %tmp_MotionMatrix(3,1)=-z_trans;
        tmp_MotionMatrix(3,1)=z_trans;
        %x,y,-z rotation
        tmp_MotionMatrix(4,1)=x_rot;
        tmp_MotionMatrix(5,1)=y_rot;
        %tmp_MotionMatrix(6,1)=-z_rot;
        tmp_MotionMatrix(6,1)=z_rot;
 %%
        tmp_MotionMatrix=round(tmp_MotionMatrix,1);
        
        [coilSen_matrix_Update,coilSen_OriCor_R,coilSen_Update_append]=getUpdatedCoilSen(tukey_window,tmp_MotionMatrix,coilSen_matrix_U,coilSen_OriCor_R,slcPos_matZ,x_off,y_off,z_off,prot,IsSaveIntp,IsSmooth,InterpApp,tmp_xyz_RotMatrix);
        
        coilSen_Update_recon=coilSen_matrix_Update;%.*coilSen_sign;
        
        
        iCount=0;
        for islc=(iSMS-1)*prot.lMultiBandFactor+1:iSMS*prot.lMultiBandFactor
            iCount=iCount+1;
            iIndex_slc=out.sliceOrderSMS(islc);
            
%             if (iIndex_slc==6)% || (iIndex_slc==5)
%                 %coilSen_Update(:,:,iIndex_slc,:)=imrotate(coilSen_Update_recon(:,:,iIndex_slc,:),-9, 'bicubic','crop').*mask_coilSen(:,:,iIndex_slc,:);
%                coilSen_Update(:,:,iIndex_slc,:)=(coilSen_Update_recon(:,:,iIndex_slc,:)).*mask_slice(:,:,iIndex_slc,:);
%                 %coilSen_Update(40:55,
%                 %32:45,iIndex_slc,:)=(coilSen_Update_recon(40:55,
%                 %32:45,iIndex_slc,:)-0.05); se=5
%                 %coilSen_Update(23:47, 16:36,iIndex_slc,:)=(coilSen_Update_recon(23:47, 16:36,iIndex_slc,:)-0.1);%se=12
% %                 coilSen_Update(41:58, 31:45,iIndex_slc,:)=(coilSen_Update_recon(41:58, 31:45,iIndex_slc,:)-0.1);%se=3
%                 %coilSen_Update(23:40,21:33,iIndex_slc,:)=(coilSen_Update_recon(23:40,21:33,iIndex_slc,:)+0.01+0.01i);%se=14
%                  %coilSen_Update(28:40,20:50,iIndex_slc,:)=(coilSen_Update_recon(28:40,20:50,iIndex_slc,:)-0.1);%se=3
%                  %coilSen_Update(34:57,27:36,iIndex_slc,:)=(coilSen_Update_recon(34:57,27:36,iIndex_slc,:)*0.05);%se=3
%                  coilSen_Update(30:end,17:end,iIndex_slc,:)=(coilSen_Update_recon(30:end,17:end,iIndex_slc,:)-1);%se=5
%           else

%             if iIndex_slc==13 || iIndex_slc==14
%                 coilSen_Update(:,:,iIndex_slc,:)=(coilSen_matrix(:,:,iIndex_slc,:)).*mask_slice(:,:,iIndex_slc,:); 
%             else
              if IsMask==1
                coilSen_Update(:,:,iIndex_slc,:)=(coilSen_Update_recon(:,:,iIndex_slc,:)).*mask_slice(:,:,iIndex_slc,:);%mask_slice(:,:,iIndex_slc,:); 
              else
                coilSen_Update(:,:,iIndex_slc,:)=(coilSen_Update_recon(:,:,iIndex_slc,:));%mask_slice(:,:,iIndex_slc,:);
              end
%            end
        end
    end
    %Phase Correction
    %[coilSen_Update,coilSen_sign]=PhaCor(coilSen_Update,prot);
    
    if IsNormalization==1
        coilSen_Update=smoothCoilSen(coilSen_Update,prot,0,tukey_window);
            %coilSen_Update_recon=smoothCoilSen(coilSen_Update_recon,mask_thresh,tukey_window,kdata_sliceimg_sens,prot,mask_coilSen);
    end
else
%     coilSen_ori=CalCoilSens(permute(kdata_sliceimg_sens,[1,2,4,3]),prot.ksize,prot.eigThresh_k,prot.eigThresh_im);%nx,ny,nslc,ncoil
%     coilSen_Update=coilSen_ori;
    %Mat_filename=strcat(raw_data_filename(1:(strfind(raw_data_filename, '.dat')-1)), '.mat');
    %load(Mat_filename)
    %coilSen_Update=coilSen_ori;
    %coilSen_Update(:,:,2:end,:)=prediction_10mm(:,:,1:end-1,:);
    %coilSen_Update(:,:,4:end,:)=prediction_30mm(:,:,1:end-3,:);
    %coilSen_Update(:,:,1:end,:)=coilSen_Update_InterP(:,:,1:end,:);
    
    if IsRef==1
        Mat_filename=strcat(raw_data_filename(1:(strfind(raw_data_filename, '.dat')-1)), '_ref.mat');
    %coilSen_ori=CalCoilSens(permute(kdata_sliceimg_sens,[1,2,4,3]),prot.ksize,prot.eigThresh_k,prot.eigThresh_im);%nx,ny,nslc,ncoil
   % Mat_filename=strcat(raw_data_filename(1:(strfind(raw_data_filename, '.dat')-1)), '_ref.mat');
    %save(Mat_filename,'coilSen_ori');
       load(Mat_filename)
       %[coilSen_matrix,coilSen_OriCor,coilSen_1,slcPos_matZ,mask_coilSen,mask_slice]=getCoilSen(mask_thresh,coilSen_thresh,tukey_window,x_off,y_off,z_off,coilSen_ori,kdata_sliceimg_sens,10,prot,IsSmooth);%prot.lMultiBandFactor+10
       %for Fig. 2
       %coilSen_matrix,coilSen_OriCor,coilSen_1,coilSen_ori,slcPos_matZ,mask_coilSen,mask_slice]=getCoilSen_7_8(mask_thresh,coilSen_thresh,tukey_window,x_off,y_off,z_off,kdata_sliceimg_sens,prot,IsSmooth);
    else
        Mat_filename=strcat(raw_data_filename(1:(strfind(raw_data_filename, '.dat')-1)), '_1st.mat');
        SMS_Factor=prot.lMultiBandFactor;
        load(Mat_filename)
        [coilSen_matrix,coilSen_OriCor,coilSen_1,slcPos_matZ,mask_coilSen,mask_slice]=getCoilSen(mask_thresh,coilSen_thresh,tukey_window,x_off,y_off,z_off,coilSen_ori,fft2c(recon_slice_1st),SMS_Factor,prot,IsSmooth);
    end
    coilSen_Update=coilSen_ori;
    if IsMask==1
        %[coilSen_matrix,coilSen_OriCor,coilSen_1,slcPos_matZ,mask_coilSen,mask_slice]=getCoilSen(mask_thresh,coilSen_thresh,tukey_window,x_off,y_off,z_off,coilSen_ori,fft2c(recon_slice_1st),SMS_Factor,prot,IsSmooth);
        coilSen_Update=coilSen_Update.*mask_slice;
%         for islc=1:prot.OriNslice
%             coilSen_Update(:,:,iIndex_slc,:)=(coilSen_Update_recon(:,:,iIndex_slc,:)).*mask_slice(:,:,iIndex_slc,:);
%         end
    end
    if IsNormalization==1
        coilSen_Update=smoothCoilSen(coilSen_Update,prot,0,tukey_window);%smoothCoilSen(coilSen_Update,prot,out,1,tukey_window);
    end
end

if prot.lMultiBandFactor==1
    KSpaceDATA=out.kdata;
    ncoil=size(KSpaceDATA,3);
    if prot.grappa~=0 %if IPAT is on
        recon_slice=zeros(size(out.kdata,1), size(out.kdata,2), size(out.kdata,4));
        sliceOrderSMS=out.sliceOrderSMS;
        for islc=1:prot.OriNslice
            res=DoGrappa(out.kdata(:,:,:,islc),out.refscan(:,:,:,sliceOrderSMS(islc)));
            recon_slice(:,:,sliceOrderSMS(islc))=sqrtSum(res);
        end
        recon_slice=permute(recon_slice,[1,2,4,3]);
        showimagefft(recon_slice,0);
        %         figure;
        %         subplot(1,2,1); imshow(sqrtSum(permute(out.kdata,[2,1,3])),[]); title('aliasing image');
        %         subplot(1,2,2); imshow(sqrtSum(res),[]); title('Grappa recon image')
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
    %coilSen_Update recon_slice_coilSens
    %coilSen=permute(coilSen_1,[1,2,4,3]);%permute(coilSen_Update(:,:,:,:,1),[1,2,4,3]);
%     for ichn=1:prot.chn
%         tmp_coilSen=squeeze(coilSen(:,:,ichn,:));  
%         tmp_coilSen=CoilSenInterPzPoly(tmp_coilSen);
%         coilSen(:,:,ichn,:)=tmp_coilSen;
%     end
%     mask_coilSen=zeros(size(coilSen,1), size(coilSen,2), size(coilSen,4));
%     for islc=1:prot.OriNslice
% %         mask=sqrtSum(out.kdata_sliceimg(:,:,:,islc),1);
% %         mask(mask>0.3*max(mask(:))) = 1;
% %         level = graythresh(mask);
% %         mask = imbinarize(mask,level);
% %         se90=strel('line',2,90);
% %         se0=strel('line',2,0);
% %         mask=imdilate(mask,[se90 se0]);
% %         se = strel('disk',4);
% %         mask=imclose(mask,se);
% %         mask=imopen(mask,se);
% %         mask = imfill(mask, 'holes');
% %         mask_coilSen(:,:,islc)=mask;
%         mask=mask_coilSen(:,:,islc);
%         coilSen(:,:,:,islc)=coilSen(:,:,:,islc).*mask;
%     end

    %recon_slice_all=zeros(prot.Nread,prot.Nphase,1,prot.OriNslice);
    coilSen=permute(coilSen_Update,[1,2,4,3]);%nx,ny,nch,nslc,coilSen_1 
    recon_slice=DoSENSE_SMS(coilSen, KSpaceDATA, out.sliceOrderSMS, prot);
    recon_slice=permute(recon_slice,[1,2,4,3]);
    %showimagefft(recon_slice,0);
    showimagefft(recon_slice(:,:,:,1:end),0,1,gray,1,0,0,'haha');
    
%     sliceimg_Sense=zeros(size(coilSen));
%     for islc=1:prot.OriNslice
%         tmp_slc=recon_slice(:,:,islc);
%         for ichn=1:prot.chn
%             sliceimg_Sense(:,:,ichn,islc)=tmp_slc.*coilSen(:,:,ichn,islc);%coilSen_ori(:,:,islc,ichn);
%         end
%     end
%     sliceimg_Sense=fft2c(sliceimg_Sense);
end

function [recon_slice, kdata_sliceimg_up]=DoSENSE_SMS(coilSen, KSpaceDATA, sliceOrderSMS, prot,iSMS_index)%
% coilSen: coil sensitivity for all single slices. [nx,ny,nch,nslc]
% KSpaceDATA: SMS kspace. [nx,ny,nchn,nslc]
  %if the NCO phase was input by programmer, then set isInputPhase to 1. 
    nx=size(KSpaceDATA,1);ny=size(KSpaceDATA,2);ncoil=size(KSpaceDATA,3);
    isInputPhase=0;
    recon_slice=zeros(nx,ny,prot.OriNslice, 'double');
    coilSen_SingleSlice=zeros(size(coilSen,1),size(coilSen,2),size(coilSen,3),prot.lMultiBandFactor);
    
    if nargin<=4
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
        
        if prot.lAccelFactPE>1 
            if prot.lMultiBandFactor==2
                CAIPIshifts=[0,pi*5/4];  
            elseif prot.lMultiBandFactor==3
                CAIPIshifts=[0,pi*4/3, pi*2/3];
            elseif prot.lMultiBandFactor==4
                CAIPIshifts=[0,pi*4/3, pi*2/3,2*pi];
            elseif prot.lMultiBandFactor==5
                CAIPIshifts= [0,pi*4/3,pi*2/3,pi,pi*4/3];
            elseif prot.lMultiBandFactor==6
                CAIPIshifts= [0,pi*2/3,pi*4/3,pi*2,pi*2/3,pi*4/3]; 
            end
        else
            if isInputPhase==0
                if prot.lMultiBandFactor==2
                    if prot.Nslice==1
                        CAIPIshifts=[0,pi];%[pi/2,0];%[-3*pi/2,-2*pi];%[3*pi/2,pi];%% [-3*pi/4,-5*pi/4];%%[13*pi/6, 5*pi/3];%
                    else
                        CAIPIshifts=[0,pi];%[0,pi/2];    
                    end
                elseif prot.lMultiBandFactor==3
                    CAIPIshifts= [0,pi,2*pi];%[0,pi*2/3,pi*4/3];
                elseif prot.lMultiBandFactor==4
                    CAIPIshifts= [0,pi*2/3,pi*4/3,pi*2];
                elseif prot.lMultiBandFactor==5
                    CAIPIshifts= [0,pi*2/3,pi*4/3,pi*2,pi*2/3]; 
                elseif prot.lMultiBandFactor==6
                    CAIPIshifts= [0,pi*2/3,pi*4/3,pi*2,pi*2/3,pi*4/3]; 
                end
            else
                if prot.lMultiBandFactor==2
                    if prot.Nslice==1
                        CAIPIshifts=[0,pi/2];%[][0,3*pi/2][-3*pi/2,-2*pi];%[3*pi/2,pi];%% [-3*pi/4,-5*pi/4];%%[13*pi/6, 5*pi/3];%
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
            end
        end
        coilSen_SingleSlice_1=coilSen_SingleSlice;
        k_calib=(fft2c(coilSen_SingleSlice_1));
        k_calib=SMS_CAIPIshift(k_calib,CAIPIshifts);%CAIPIshifts
        coilSen_SingleSlice_1=(ifft2c(k_calib));%permute(ifft2c(k_calib),[2,1,3,4]);
        
        if prot.lMultiBandFactor==2
            coilSen_SingleSlice(:,:,:,2)=(circshift(coilSen_SingleSlice(:,:,:,2),[0,-ny/2+1,0]));%-ny/4+1 (circshift(coilSen_SingleSlice(:,:,:,2),[0,-ny/4+1,0]));
        elseif prot.lMultiBandFactor==3
            %CAIPIshifts= [0,pi*2/3,pi*4/3];
%             coilSen_SingleSlice(:,:,:,2)=(circshift(coilSen_SingleSlice(:,:,:,2),[0,round(-ny/3),0]));
%             coilSen_SingleSlice(:,:,:,3)=(circshift(coilSen_SingleSlice(:,:,:,3),[0,round(ny/3),0]));
              %CAIPIshifts= [0,pi,2*pi];
              coilSen_SingleSlice(:,:,:,2)=(circshift(coilSen_SingleSlice(:,:,:,2),[0,-ny/2+1,0]));
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
        %coilSen_SingleSlice=permute(coilSen_SingleSlice,[2,1,3,4]);
        kspace_sms=KSpaceDATA(:,:,:,iSMS);%permute(KSpaceDATA(:,:,:,iSMS),[2,1,3,4]);
        %kspace_sms=permute(kspace_sms, [2,1,3,4]);
        
        imgIndex=sliceOrderSMS((iSMS-1)*prot.lMultiBandFactor+1:iSMS*prot.lMultiBandFactor);
%         SNS=ESPIRiT(coilSen_SingleSlice);
%         nIterCG=12;
%         [reskSENSE, resSENSE] = cgESPIRiT(kspace_sms, SNS, nIterCG,0.01 ,0);
        aliasCoilSen=permute(CalCoilSens(permute(kspace_sms,[1,2,4,3]),[6 6],0.06,0.95),[1,2,4,3]);
        
        recon_slice_1=SENSE_SMS((ifft2c(kspace_sms)),coilSen_SingleSlice,prot.lMultiBandFactor);
        %recon_slice_1=SENSE_SMS(aliasCoilSen,coilSen_SingleSlice,prot.lMultiBandFactor);
        if prot.lMultiBandFactor==2
            recon_slice_1(:,:,2)=circshift(recon_slice_1(:,:,2),[0,ny/2,0]);%ny/4 (circshift(recon_slice_1(:,:,2),[0,ny/4,0]));
        elseif prot.lMultiBandFactor==3
%             recon_slice_1(:,:,2)=(circshift(recon_slice_1(:,:,2),[0,round(ny/3),0]));
%             recon_slice_1(:,:,3)=(circshift(recon_slice_1(:,:,3),[0,round(-ny/3),0]));
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
    %imshow(squeeze(res0(1:nx*2,:,1)),[]);
    %showimagefft(recon_slice_1);
end

% function showimagefft(kspace, fftparam) %kspace:nx,ny,ncoil,nslice
%     
%     if nargin<=1
%         isfft=1;
%     else
%         if fftparam==0
%             isfft=0;
%         else
%             isfft=1;
%         end
%     end
%     NSlice=size(kspace,4);
%     if NSlice==1
%         imagesc(flip(sqrtSum(kspace(:,:,:),isfft)'));colormap gray;
%         set(gca,'xtick',[],'xticklabel',[]);
%         set(gca,'ytick',[],'yticklabel',[]);
%         return;
%     end
%     
%     if mod(NSlice,4)==0
%         nrow=NSlice/4;
%     elseif mod(NSlice,3)==0
%         nrow=NSlice/3;
%     elseif mod(NSlice,5)==0
%         nrow=NSlice/5;
%     elseif mod(NSlice,2)==0
%         nrow=NSlice/2;
%     end
%     ncol=NSlice/nrow;%ncoil nslice
%     ha = tight_subplot(nrow,ncol,   [   .07         .005    ],  [ 0.05  .015], [ 0.002   .002]);
%     %                  row, col, [gap_height, gap_width],  [bottom  top], [left    right]
%     %figure;
%     for irow=1:nrow
%         for icol=1:ncol
%             iIndex=(irow-1)*ncol+icol;
%             axes(ha(iIndex));
%             imagesc(flip(sqrtSum(kspace(:,:,:,iIndex),isfft)'));colormap gray;
%             %for single slice images
%             %imagesc(fliplr(sqrtSum(kspace(:,:,:,iIndex),isfft)'));colormap gray;
%             set(gca,'xtick',[],'xticklabel',[]);
%             set(gca,'ytick',[],'yticklabel',[]);
%             %imagesc((sqrtSum(kspace(:,:,:,iIndex),isfft)));colormap gray;%
%             %imagesc(abs((kspace(:,:,:,iIndex))));colormap gray;%
%         end
%     end
% end

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

