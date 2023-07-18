function coilSen_Update=getUpdatedCoilSen(tukey_window,MotionMatrix, coilSen_matrix,coilSen_OriCor,prot, IsSmooth, InterpApp,tmp_xyz_RotMatrix)
%get CSM interpolation
%INPUT: tukey_window, smooth filter
%       MotionMatrix, translations and rotations
%       coilSen_thresh, a specific mask threshold 
%       coilSen_matrix, CSM used for the update
%       coilSen_OriCor, x,y,z Coordinates on CSM
%       x_off=0,y_off=0,z_off=0
%       prot, image parameters
%       IsSmooth, CSM smooth option
%       InterpApp, interpolation type 1: real&imaginary, 2: magnitude&phase
%       kdata_sliceimg_sens, the single-slice reference k-space
%       SMS_factor, SMS acceleration factor
%       tmp_xyz_RotMatrix, motion RotMatrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT: coilSen_Update, Updated CSM
%
%%%%%% Written by: Bo Li, University of Maryland, Baltimore
%%%%%% for paper "Simultaneous multislice EPI prospective motion correction 
%%%%%% by realtime receiver phase correction and coil sensitivity map interpolation. 
%%%%%% Magn Reson Med. 2023;1-17. doi:10.1002/mrm.29789"
%%%%%% Created on Sep. 22, 2022

    coilSen_Update=zeros(prot.Nread,prot.Nphase,prot.OriNslice,prot.chn);
    
    T_R_matrix=zeros(4,4,'double');
    T_R_matrix(1,1)=1;
    T_R_matrix(2,2)=1;
    T_R_matrix(3,3)=1;

    x_trans=MotionMatrix(1,1); 
    y_trans=MotionMatrix(2,1); 
    z_trans=MotionMatrix(3,1);
    
    x_rot=MotionMatrix(4,1);
    y_rot=MotionMatrix(5,1);
    z_rot=MotionMatrix(6,1);%degree
    
    x_rot_matrix=[1 0 0;0 cosd(x_rot) sind(x_rot);0 -sind(x_rot) cosd(x_rot)];
    y_rot_matrix=[cosd(y_rot) 0 -sind(y_rot);0 1 0;sind(y_rot) 0 cosd(y_rot)];
    z_rot_matrix=[cosd(z_rot) sind(z_rot) 0;-sind(z_rot) cosd(z_rot) 0;0 0 1];
    
    if size(tmp_xyz_RotMatrix,1)==0
        T_R_matrix(1:3,1:3)=x_rot_matrix*y_rot_matrix*z_rot_matrix;
    else
        T_R_matrix(1:3,1:3)=tmp_xyz_RotMatrix;
    end
    T_R_matrix(4,1:3)=[x_trans, y_trans, z_trans];
    T_R_matrix(4,4)=1;
    tform=affine3d(T_R_matrix);
    z_start_ind=prot.InterpSlice(1);
    z_end_ind=prot.InterpSlice(end);
    z_cor_start=coilSen_OriCor{1,1}.z_cor_allSlc(z_start_ind);
    z_cor_end=coilSen_OriCor{1,1}.z_cor_allSlc(z_end_ind);
    [R_cor, C_cor, Z_cor]=ndgrid(coilSen_OriCor{1,1}.x_cor, coilSen_OriCor{1,1}.y_cor, coilSen_OriCor{1,1}.z_cor_allSlc(z_start_ind:z_end_ind));%z_cor_start:z_cor_end
    [R_mot, C_mot, Z_mot]=transformPointsForward(tform,R_cor, C_cor, Z_cor);
    R_mot=round(R_mot,3);C_mot=round(C_mot,3);Z_mot=round(Z_mot,3);
    [RR, CC, ZZ]=ndgrid(coilSen_OriCor{1,1}.x_cor, coilSen_OriCor{1,1}.y_cor, coilSen_OriCor{1,1}.z_cor_allSlc_append(z_start_ind:z_end_ind));
    for ich=1:prot.chn        
        if InterpApp==1 %real&imaginary
            tmp_coilSen =interpn(RR,CC,ZZ,coilSen_matrix(:,:,z_start_ind:z_end_ind,ich),R_mot,C_mot,Z_mot,'makima');
        elseif InterpApp==2 %magnitude&phase
            mag_coilSen_matrix=abs(coilSen_matrix(:,:,:,ich));
            uw_ph_coilSen_matrix=angle(coilSen_matrix(:,:,:,ich));  
            tmp_mag_coilSen =interpn(RR,CC,ZZ,mag_coilSen_matrix,R_mot,C_mot,Z_mot,'makima');
            tmp_uw_ph_coilSen =interpn(RR,CC,ZZ,uw_ph_coilSen_matrix,R_mot,C_mot,Z_mot,'makima');
            tmp_coilSen = tmp_mag_coilSen .* exp(tmp_uw_ph_coilSen .* 1i);
        elseif InterpApp==3
            mag_coilSen_matrix=flip(abs(coilSen_matrix(:,:,:,ich)),3);
            uw_ph_coilSen_matrix=angle(coilSen_matrix(:,:,:,ich));
            uw_ph_coilSen_matrix=unwrap(uw_ph_coilSen_matrix,[],1);
            uw_ph_coilSen_matrix=unwrap(uw_ph_coilSen_matrix,[],2);
            uw_ph_coilSen_matrix=unwrap(uw_ph_coilSen_matrix,[],3);
            uw_ph_coilSen_matrix=flip(uw_ph_coilSen_matrix,3);
            tmp_mag_coilSen =interpn(R_mot,C_mot,Z_mot,mag_coilSen_matrix,RR,CC,ZZ,'makima');
            tmp_uw_ph_coilSen =interpn(R_mot,C_mot,Z_mot,uw_ph_coilSen_matrix,RR,CC,ZZ, 'nearest');
            
            tmp_uw_ph_coilSen_1 =interpn(R_mot,C_mot,Z_mot,uw_ph_coilSen_matrix,RR,CC,ZZ,'makima');
            iIndex=isnan(tmp_uw_ph_coilSen(:));
            iIndex=find(iIndex(:)==1);
            tmp_uw_ph_coilSen(iIndex)=tmp_uw_ph_coilSen_1(iIndex);

            tmp_mag_coilSen=flip(tmp_mag_coilSen,3);
            tmp_uw_ph_coilSen=flip(tmp_uw_ph_coilSen,3);
            tmp_coilSen = tmp_mag_coilSen .* exp(tmp_uw_ph_coilSen .* 1i);
            
        end
        coilSen_Update(:,:,z_start_ind:z_end_ind,ich)=tmp_coilSen(:,:,:);
        if z_start_ind>1
            coilSen_Update(:,:,1:z_start_ind-1,ich)=coilSen_matrix(:,:,1:z_start_ind-1,ich);    
        end
        if z_end_ind<prot.OriNslice
            coilSen_Update(:,:,z_end_ind+1:prot.OriNslice,ich)=coilSen_matrix(:,:,z_end_ind+1:prot.OriNslice,ich);    
        end
        if IsSmooth==1
            coilSen_Update(:,:,:,ich)=ifft2c(fft2c(coilSen_Update(:,:,:,ich)).*tukey_window);
        else
            coilSen_Update(:,:,:,ich)=coilSen_Update(:,:,:,ich);
        end
    end
    
end