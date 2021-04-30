function [coilSen_matrix,coilSen_OriCor,coilSen_1,coilSen_ori,slcPos_matZ,mask_coilSen,mask_slice]=getCoilSen(mask_thresh,coilSen_thresh,tukey_window,x_off,y_off,z_off,kdata_sliceimg_sens,prot,IsSmooth)
%get coil sensitivity from single slice reference data and prot (SMS imaging protocol, such as te, tr, matrix size, etc.)
%author: Bo Li, 4-28-2021
%%input: 
%mask_thresh: the mask threshold for single slice used for reconstruction
%coilSen_thresh: the mask threshold for coil sensitivity interpolation.
%The interpolation mask is larger than reconstruction mask to get better
%results on image edges. 
%tukey_window is the low-pass filter
%x_off,y_off,z_off are 0
%kdata_sliceimg_sens is single-slice reference k-space data
%prot is SMS imaging protocol
%IsSmooth is smooth option for coil sensitivity maps

%%output: coilSen_ori is the original coil sensitivity maps
%coilSen_matrix: (smoothed) masked coilSen maps with coilSen_thresh threshold
%coilSen_OriCor: a struct consisting of the motion information of
%acquisition
%coilSen_1: masked coilSen with mask_thresh threshold 
%coilSen_ori: masked coilSen with coilSen_thresh threshold
%slcPos_matZ: slice position along z (slice) direction
%mask_slice: the binary masks for single slices with mask_thresh threshold 
%mask_coilSen: the binary masks for coilSen with coilSen_thresh threshold 
%%
    %low-rank coil sensitivity
    coilSen_ori=CalCoilSens(permute(kdata_sliceimg_sens,[1,2,4,3]),prot.ksize,prot.eigThresh_k,prot.eigThresh_im);
    mask_coilSen=zeros(prot.Nread,prot.Nphase,prot.OriNslice);
    mask_slice=mask_coilSen;
    nx=prot.Nphase; ny=prot.Nread;
     
    for islc=1:prot.OriNslice
        tmp=coilSen_ori(:,:,islc,:);        
        mask=sqrtSum(kdata_sliceimg_sens(:,:,:,islc),1);
        mask(mask>coilSen_thresh*max(mask(:))) = 1;
        level = graythresh(mask);
        mask = imbinarize(mask,level);
        se90=strel('line',4,90);
        se0=strel('line',4,0);
        mask=imdilate(mask,[se90 se0]);
        se = strel('disk',6);
        mask=imclose(mask,se);
        mask=imopen(mask,se);
        mask = imfill(mask, 'holes');
        mask_coilSen(:,:,islc)=mask;
        coilSen_ori(:,:,islc,:) = tmp;
    end
     
    mtotal_Zsize=prot.OriNslice;
    %x and y directions allow for maximum 32/0.5=64mm offset
    %z direction allows 50mm offset, 50mm/0.5=100
    x_off_mat=x_off;%
    y_off_mat=y_off;%
    z_off_mat=z_off;%

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
        
        tmp=squeeze(coilSen_matrix(:,:,islc,:));
        mask=sqrtSum(kdata_sliceimg_sens(:,:,:,islc),1);
        mask(mask>mask_thresh*max(mask(:))) = 1;%
        level = graythresh(mask);
        mask = imbinarize(mask,level);
        se90=strel('line',4,90);
        se0=strel('line',4,0);
        mask=imdilate(mask,[se90 se0]);
        se = strel('disk',6);
        mask=imclose(mask,se);
        mask=imopen(mask,se);
        mask = imfill(mask, 'holes');
        mask_slice(:,:,islc)=mask;

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