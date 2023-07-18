function [coilSen_matrix,coilSen_OriCor,coilSen_1,slcPos_matZ,mask_coilSen,mask_slice]=getCoilSen(mask_thresh,coilSen_thresh,tukey_window,x_off,y_off,z_off,coilSen_ori,kdata_sliceimg_sens,SMS_factor,prot,IsSmooth)
%get smoothed and masked coil sensitivity map from the original CSM (coilSen_ori) calculated by ESPIRIT
%INPUT: mask_thresh, slice-mask threshold
%       coilSen_thresh, a specific mask threshold 
%       tukey_window, smooth filter
%       x_off=0,y_off=0,z_off=0
%       coilSen_ori, the original CSM calculated by ESPIRIT
%       kdata_sliceimg_sens, the single-slice reference k-space
%       SMS_factor, SMS acceleration factor
%       prot, image parameters
%       IsSmooth, CSM smooth option
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT: coilSen_matrix, the CSM with a specific mask threshold
%        coilSen_OriCor, x,y,z Coordinates on CSM
%        coilSen_1, the CSM with slice-mask threshold
%        slcPos_matZ, slice Z coordinates
%        mask_coilSen, the specific mask matrix
%        mask_slice, the slice mask matrix
%%%%%% Written by: Bo Li, University of Maryland, Baltimore
%%%%%% for manuscript "SMS-EPI prospective motion correction 
%%%%%% by real-time phase compensation and coil sensitivity map interpolation"
%%%%%% Created on Sep. 22, 2022
    
    single_sliceimg=ifft2c(kdata_sliceimg_sens); %nx,ny,ncoil,nslc
    single_sliceimg=permute(single_sliceimg,[1,2,4,3]); %nx,ny,nslc,ncoil
    mask_coilSen=zeros(prot.Nread,prot.Nphase,prot.OriNslice);
    mask_slice=mask_coilSen;
    nx=prot.Nphase; ny=prot.Nread;
     
    for islc=1:prot.OriNslice
            if SMS_factor==2
                mask=sqrtSum(squeeze(coilSen_ori(:,:,islc,:)),0);
            else 
                mask=sqrtSum(squeeze(single_sliceimg(:,:,islc,:)),0);
            end
            mask(mask>coilSen_thresh*max(mask(:))) = 1;
            level = graythresh(mask);
            mask = imbinarize(mask,level);
            if SMS_factor>2 
                se90=strel('line',4,90);
                se0=strel('line',4,0);
                mask=imdilate(mask,[se90 se0]);
                se = strel('disk',4);
                mask=imclose(mask,se);
                mask=imopen(mask,se);
            end
            mask = imfill(mask, 'holes');
            mask_coilSen(:,:,islc)=mask;
        
            if SMS_factor==2
                mask=sqrtSum(squeeze(coilSen_ori(:,:,islc,:)),0);
            else
                mask=sqrtSum(squeeze(single_sliceimg(:,:,islc,:)),0);
            end
            mask(mask>mask_thresh*max(mask(:))) = 1;%
            level = graythresh(mask);
            mask = imbinarize(mask,level);
            if SMS_factor>2 
                se90=strel('line',4,90);
                se0=strel('line',4,0);
                mask=imdilate(mask,[se90 se0]);
                se = strel('disk',4);
                mask=imclose(mask,se);
                mask=imopen(mask,se);
            end
            mask = imfill(mask, 'holes');
            mask_slice(:,:,islc)=mask;
    end
     
     
    mtotal_Zsize=prot.OriNslice;%
    
    x_off_mat=x_off;%ceil(x_off/x_reso); 
    y_off_mat=y_off;%ceil(y_off/y_reso);
    z_off_mat=z_off;%ceil(z_off/z_reso);

    coilSen_matrix=zeros(size(coilSen_ori,1)+x_off_mat*2,size(coilSen_ori,2)+y_off_mat*2,mtotal_Zsize+2*z_off_mat, size(coilSen_ori,4));
    slcPos_matZ=zeros(1,prot.OriNslice);

    for islc=1:prot.OriNslice
        tmp_z=islc;
        slcPos_matZ(1,islc)=tmp_z;
        %fill all element in coilSen(:,:,islc,:)
        if IsSmooth==1
            tmp=ifft2c(fft2c(coilSen_ori(:,:,islc,:)).*tukey_window).*mask_coilSen(:,:,islc);
        else
            tmp=coilSen_ori(:,:,islc,:).*mask_coilSen(:,:,islc);
        end
        coilSen_matrix(x_off_mat+1:x_off_mat+nx,y_off_mat+1:y_off_mat+ny,tmp_z,:)=tmp;
    end

    coilSen_1=coilSen_matrix(x_off_mat+1:x_off_mat+nx,y_off_mat+1:y_off_mat+ny,slcPos_matZ,:);

    for islc=1:prot.OriNslice
            mask=mask_slice(:,:,islc);
            tmp=coilSen_ori(:,:,islc,:);
            if IsSmooth==1
                tmp=ifft2c(fft2c(tmp).*tukey_window).*mask;%
            else
                tmp=tmp.*mask;
            end
            coilSen_1(:,:,islc,:) = tmp;
    end    
    
    %coilSen x,y,z locations at expanded coil sensitivity
    coilSen_OriCor=cell(prot.OriNslice,1);
    tmp_z_cor=zeros(1,prot.OriNslice);
    for islc=1:prot.OriNslice
        tmp_z_cor(1,islc)=prot.slicePos_Matrix{islc,1}.img_pos_Z_Cor;
    end
    for islc=1 : prot.OriNslice
        tmp_y_cor=prot.slicePos_Matrix{islc,1}.img_pos_Y_Cor;
        tmp_x_cor=prot.slicePos_Matrix{islc,1}.img_pos_X_Cor;
        tmp_x_cor_all=tmp_x_cor;
        tmp_y_cor_all=tmp_y_cor;
        tmp_z_cor_all=tmp_z_cor;
        tmp_z_cor_allSlc=tmp_z_cor;

        tmp_coilSen_OriCor_info=struct('x_cor', tmp_x_cor,'y_cor', tmp_y_cor, ...
                                  'z_cor', tmp_z_cor,...
                                  'x_cor_all', tmp_x_cor_all,...
                                  'y_cor_all', tmp_y_cor_all,...
                                  'z_cor_all', tmp_z_cor_all,...
                                  'z_cor_allSlc', tmp_z_cor_allSlc,'x_cor_U',tmp_x_cor,'y_cor_U',tmp_y_cor,'z_cor_allSlc_append',tmp_z_cor_allSlc);
        coilSen_OriCor{islc, 1}=tmp_coilSen_OriCor_info;                    
    end
end