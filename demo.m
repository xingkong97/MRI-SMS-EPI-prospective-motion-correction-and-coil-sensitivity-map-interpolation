%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Main function 
%%%%%%
%%%%%% An example for SMS-EPI motion correction
%%%%%%
%%%%%% Written by: Bo Li, University of Maryland, Baltimore
%%%%%% for paper "Simultaneous multislice EPI prospective motion correction 
%%%%%% by real‐time receiver phase correction and coil sensitivity map
%%%%%% interpolation. Magn Reson Med. 2023;1-17. doi:10.1002/mrm.29789"
%%%%%% Created on Sep. 22, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%load SMS-EPI k-space rawdata and motion data
clear all;
close all;
addpath([pwd,'/mapVBVD']);
addpath([pwd,'/utils']);
addpath([pwd,'/ESPIRIT']);
addpath([pwd,'/SMSRecon']);

load('data/SMS4.mat');
%plot current repetition motion data
%RepIndex is the index of Rep
%this data is for RepIndex repetition
figure(100);
plotTransAndRot(measMotionInfo, RepIndex,1);

%%
%calculate coil sensitivity map (CSM)
ksize = [6,6]; % ESPIRiT kernel-window-size
eigThresh_k = 0.02; %defaukt 0.02 threshold of eigenvectors in k-space
eigThresh_im = 0.95; %default 0.95
%kdata_sliceimg_sens is the single slice raw data acquired in prescan
kdata_sliceimg_sens=out.kdata_sliceimg;
%calculate CSM by ESPIRiT algorithm
coilSen_ori=CalCoilSens(permute(kdata_sliceimg_sens,[1,2,4,3]),ksize,eigThresh_k,eigThresh_im);

%%
%interpolating (extrapolating) CSM with motion parameters to get updated CSM
IsSmooth=1; %smooth the coil sensitivity
IsUpdateCoilSen=1;
coilSen_Update = UpdateCSM(coilSen_ori,prot,out,IsSmooth,array_timestamp,MotionData,IsUpdateCoilSen);
%smooth coilSen_ori
IsSmooth=1;
IsUpdateCoilSen=0;
coilSen_ori=UpdateCSM(coilSen_ori,prot,out,IsSmooth,array_timestamp,MotionData,IsUpdateCoilSen);

%%
%SMS-SENSE 

%SMS-SENSE Recon with original CSM (oCSM)
recon_slice_oCSM=SMS_Recon_SENSE(coilSen_ori,prot,out);
figure(10);
showimagefft(recon_slice_oCSM,0);%kspace, fftparam, flipparam, cmap_or_jet, flipabs, isSameScale,IsSaveImageToDisk, folder_name
%showimagefft(recon_slice_oCSM,0,1,gray,1,0,0,'');
title('SMS-SENSE Recon with oCSMs', 'FontSize',20);

%%
%SMS-SENSE Recon with updated CSM (uCSM)
recon_slice_uCSM=SMS_Recon_SENSE(coilSen_Update,prot,out);
figure(20);
showimagefft(recon_slice_uCSM,0);
%showimagefft(recon_slice_uCSM,0,1,gray,1,0,0,'');
title('SMS-SENSE Recon with uCSMs', 'FontSize',20);

%%
%Split Slice-GRAPPA Recon 
%SP-SG Recon with original CSM (oCSM)
recon_slice_oCSM_SSG=SMS_Recon_GRAPPA(coilSen_ori,prot,out);
figure(30);
showimagefft(recon_slice_oCSM_SSG,0);
%showimagefft(recon_slice_oCSM,0,1,gray,1,0,0,'');
title('Split Slice-GRAPPA Recon with oCSMs', 'FontSize',20);
%%
%SP-SG Recon with updated CSM (uCSM)
recon_slice_uCSM_SSG=SMS_Recon_GRAPPA(coilSen_Update,prot,out);
figure(40);
showimagefft(recon_slice_uCSM_SSG,0);
%showimagefft(recon_slice_oCSM,0,1,gray,1,0,0,'');
title('Split Slice-GRAPPA Recon with uCSMs', 'FontSize',20);

%%
%single-slice images acquired in pre-scan, used for SMS reconstruction
figure(50);
showimagefft(kdata_sliceimg_sens,1);
%showimagefft(recon_slice_oCSM,0,1,gray,1,0,0,'');
title('Single-slice reference images', 'FontSize',20);



