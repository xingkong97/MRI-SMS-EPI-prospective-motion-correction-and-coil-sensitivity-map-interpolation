%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Main function 
%%%%%%
%%%%%% An example for the reconstruction of SMS-EPI with motion correction
%%%%%%
%%%%%% Written by: Bo Li in Thomas Ernst's Lab 
%%%%%%             University of Maryland, Baltimore
%%%%%%             
%%%%%% Citation: 
%%%%%% Created on 4-28-2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%for matlab version 2016a or above
addpath(strcat(pwd,'/utils'));
addpath(strcat(pwd,'/ESPIRIT/utils'));
addpath(strcat(pwd,'/ESPIRIT/ESPIRIT'));
%clear;
close all;
load 'data_1.mat'
load 'data_2.mat'
out=out_1;
out.kdata_sliceimg=out_2.kdata_sliceimg;
out.phasecorSMS=out_2.phasecorSMS;
out.phasecor=out_2.phasecor;

%the k-space data is  for No. 15 measurement with SMS acceleration factor 2
RepIndex=15;

%interpolation Approch 1: real&imaginary, 2: magnitude&phase
InterpApp=1;

%mask_thresh: the threshold for masking single slice for reconstruction
%coilSen_thresh: the mask threshold for coil sensitivity interpolation.
mask_thresh=0.22;
coilSen_thresh=0.03;
ksize = [6,6]; %ESPIRiT kernel-window-size
eigThresh_k = 0.02; % threshold of eigenvectors in k-space
eigThresh_im = 0.95;

%prot is SMS imaging protocol
prot.ksize=ksize;
prot.eigThresh_k=eigThresh_k;
prot.eigThresh_im=eigThresh_im;

%IsSmooth: 1, smooth for coil sensitivity maps, 2, no smooth
IsSmooth=1;
%Normalization for all coil sensitivity maps
IsNormalization=1;
%the matrix for six motion parameters. 
MotionMatrix=zeros(6,1);

x_off=0; y_off=0; z_off=0; 

%tukey_window is the low-pass filter for the smooth purpose
[tukey_window,tukey_window_red]=filter_2D(prot.Nphase);%the tukey_window size won't change recon results

%out is the class which includes SMS k-space data, acquisition order...
%out.kdata_sliceimg is single-slice reference k-space data
%out.sliceOrderSMS is the slice acquisition order
kdata_sliceimg_sens=out.kdata_sliceimg;
kdata_sliceimg_R=kdata_sliceimg_sens;

%get coil sensitivity from single slice reference data
[coilSen_matrix,coilSen_OriCor,coilSen_1,coilSen_ori,slcPos_matZ,mask_coilSen,mask_slice]=getCoilSen(mask_thresh,coilSen_thresh,tukey_window,x_off,y_off,z_off,kdata_sliceimg_sens,prot,IsSmooth);
coilSen_OriCor_R=coilSen_OriCor;

%updated coil sensitivity
coilSen_Update=zeros(prot.Nread,prot.Nphase,prot.OriNslice,prot.chn);

%motion struct for SMS imaging
SMS_MotionCell=cell(prot.Nslice,1);

%generate RotMatrix for each acquisition
for iSMS=1:prot.Nslice
    %search the index of timestamp for the acquisition
    iIndex=find(array_timestamp(:,1)==out.mcTimeStamp((iSMS-1)*prot.Nphase+1));
    tmp_RotMatrix=MotionData{iIndex(1,1),1}.RotMatrix;
    tmp_RotMatrix_info=struct('timestamp', array_timestamp(iIndex(1,1),1), 'RotMatrix', tmp_RotMatrix,...
                              'transX', MotionData{iIndex(1,1),1}.transX,...
                              'transY', MotionData{iIndex(1,1),1}.transY,...
                              'transZ', MotionData{iIndex(1,1),1}.transZ,...
                              'Euler_deg', MotionData{iIndex(1,1),1}.Euler_deg);
    SMS_MotionCell{iSMS, 1}=tmp_RotMatrix_info;

    %x,y,z translations
    x_trans=SMS_MotionCell{iSMS,1}.transX; 
    y_trans=SMS_MotionCell{iSMS,1}.transY; 
    z_trans=SMS_MotionCell{iSMS,1}.transZ;
    
    tmp_RotMatrix=SMS_MotionCell{iSMS,1}.RotMatrix;

    %x,y,z rotations
    x_rot=SMS_MotionCell{1,1}.Euler_deg(1,1);
    y_rot=SMS_MotionCell{1,1}.Euler_deg(1,2);
    z_rot=SMS_MotionCell{1,1}.Euler_deg(1,3);
    
    x_rot_matrix=[1 0 0;0 cosd(x_rot) sind(x_rot);0 -sind(x_rot) cosd(x_rot)];
    y_rot_matrix=[cosd(y_rot) 0 -sind(y_rot);0 1 0;sind(y_rot) 0 cosd(y_rot)];
    z_rot_matrix=[cosd(-z_rot) sind(-z_rot) 0;-sind(-z_rot) cosd(-z_rot) 0;0 0 1];
    
    f_RotMatrix=@(x,y,z)(x*y*z);
    tmp_xyz_RotMatrix=f_RotMatrix(z_rot_matrix,y_rot_matrix,x_rot_matrix);

    %x,y,-z rotation
    MotionMatrix(1,1)=x_trans;
    MotionMatrix(2,1)=y_trans;
    MotionMatrix(3,1)=-z_trans;
    %x,y,-z translation
    MotionMatrix(4,1)=x_rot;
    MotionMatrix(5,1)=y_rot;
    MotionMatrix(6,1)=-z_rot;
    
    %get updated coilSen maps
    coilSen_matrix_Update=getUpdatedCoilSen(tukey_window,MotionMatrix,coilSen_ori,coilSen_OriCor_R,slcPos_matZ,z_off,prot,IsSmooth,tmp_xyz_RotMatrix);
    
    if IsNormalization==1
        coilSen_matrix_Update=smoothCoilSen(coilSen_matrix_Update.*mask_slice,prot);
    end

    for islc=(iSMS-1)*prot.lMultiBandFactor+1:iSMS*prot.lMultiBandFactor
        iIndex_slc=out.sliceOrderSMS(islc);
        coilSen_Update(:,:,iIndex_slc,:)=coilSen_matrix_Update(:,:,iIndex_slc,:);
    end

end

KSpaceDATA=out.kdata;
kdata_sliceimg=out.kdata_sliceimg;
%SENSE recon
%recon by original coilSen maps
coilSen=permute(coilSen_1,[1,2,4,3]);
recon_slice_Ori=DoSENSE_SMS(coilSen, KSpaceDATA, out.sliceOrderSMS, prot);
recon_slice_Ori=permute(recon_slice_Ori,[1,2,4,3]);
img_S_1=flip(abs(recon_slice_Ori(:,:,1,4))');
img_S_2=flip(abs(recon_slice_Ori(:,:,1,6))');
img_S_3=flip(abs(recon_slice_Ori(:,:,1,7))');

%recon by updated coilSen maps
coilSen=permute(coilSen_Update,[1,2,4,3]);
recon_slice_U=DoSENSE_SMS(coilSen, KSpaceDATA, out.sliceOrderSMS, prot);
recon_slice_U=permute(recon_slice_U,[1,2,4,3]);
img_U_1=flip(abs(recon_slice_U(:,:,1,4))');
img_U_2=flip(abs(recon_slice_U(:,:,1,6))');
img_U_3=flip(abs(recon_slice_U(:,:,1,7))');

%single-slice reference images
img_R_1=flip(sqrtSum(kdata_sliceimg(:,:,:,4),1,1)');
img_R_2=flip(sqrtSum(kdata_sliceimg(:,:,:,6),1,1)');
img_R_3=flip(sqrtSum(kdata_sliceimg(:,:,:,7),1,1)');

%show all recon images
%showimagefft(recon_slice_U,0);

subplot(2,1,1);
imagesc(cat(2,img_S_1,img_S_2,img_S_3));colormap jet;
title('SENSE images recon with Original coilSen maps');
subplot(2,1,2);
imagesc(cat(2,img_U_1,img_U_2,img_U_3));colormap jet;
title('SENSE images recon with Updated coilSen maps');
% subplot(3,1,3);
% imagesc(cat(2,img_R_1,img_R_2,img_R_3));colormap jet;
% title('Single-slice reference images');

%Split slice-GRAPPA recon
%recon by original coilSen maps
recon_slice_L=DoSplitSliceGrappa(kdata_sliceimg, KSpaceDATA, out.sliceOrderSMS, prot);
img_L_1=flip(sqrtSum(recon_slice_L(:,:,:,4),0,1)');
img_L_2=flip(sqrtSum(recon_slice_L(:,:,:,6),0,1)');
img_L_3=flip(sqrtSum(recon_slice_L(:,:,:,7),0,1)');

%recon by updated coilSen maps
kdata_sliceimg=out.kdata_sliceimg;
img_csm=zeros(prot.Nread,prot.Nphase,prot.chn,prot.OriNslice);
img_slc_chn=zeros(prot.Nread,prot.Nphase,prot.chn);
for iSMS=1:prot.Nslice
    for islc=(iSMS-1)*prot.lMultiBandFactor+1:iSMS*prot.lMultiBandFactor
        iIndex_slc=out.sliceOrderSMS(islc);
        img_slc=sqrtSum((kdata_sliceimg(:,:,:,iIndex_slc)),1);
        for ichn=1:prot.chn
            img_slc_chn(:,:,ichn)=img_slc.*coilSen_Update(:,:,iIndex_slc,ichn);
            kdata_sliceimg(:,:,ichn,iIndex_slc)=fft2c(img_slc_chn(:,:,ichn));
        end
        img_csm(:,:,:,iIndex_slc)=img_slc_chn;
    end
end

recon_slice_T=DoSplitSliceGrappa(kdata_sliceimg, KSpaceDATA, out.sliceOrderSMS, prot);
%recon_slice_T=permute(recon_slice_T,[1,2,4,3]);
img_T_1=flip(sqrtSum(recon_slice_T(:,:,:,4),0,1)');
img_T_2=flip(sqrtSum(recon_slice_T(:,:,:,6),0,1)');
img_T_3=flip(sqrtSum(recon_slice_T(:,:,:,7),0,1)');

figure;
subplot(2,1,1);
imagesc(cat(2,img_L_1,img_L_2,img_L_3));colormap jet;
title('Split slice-GRAPPA images recon with Original coilSen maps');
subplot(2,1,2);
imagesc(cat(2,img_T_1,img_T_2,img_T_3));colormap jet;
title('Split slice-GRAPPA images recon with Updated coilSen maps');
% subplot(3,1,3);
% imagesc(cat(2,img_R_1,img_R_2,img_R_3));colormap jet;
% title('Single-slice reference images')

%plot motion parameters. x axis means the abstract time
%measMotionInfo is the motion parameters for all measurements
figure;plotTransAndRot(measMotionInfo, RepIndex);



