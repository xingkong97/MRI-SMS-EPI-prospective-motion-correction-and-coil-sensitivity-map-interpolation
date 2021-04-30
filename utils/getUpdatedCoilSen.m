function coilSen_Update=getUpdatedCoilSen(tukey_window,MotionMatrix, coilSen_matrix,coilSen_OriCor,slcPos_matZ,z_off,prot, IsSmooth, tmp_xyz_RotMatrix)
%this function is to generate a new coil sensitivity maps by interpolation
%method
%author: Bo Li, 4-28-2021
%input
%tukey_window: smooth filter
%MotionMatrix: the matrix for six motion parameters
%coilSen_matrix: the coilSen matrix used for interpolation 
%coilSen_OriCor: a struct consisting of the motion information of acquisition
%slcPos_matZ: slice position along z (slice) direction
%z_off: 0
%prot: SMS imaging protocol
%IsSmooth: Smooth option
%tmp_xyz_RotMatrix

%output
%coilSen_Update: interpolated coil sensitivity maps

    coilSen_Update=zeros(prot.Nread,prot.Nphase,prot.OriNslice,prot.chn);
    z_off_mat=z_off;
    
    T_R_matrix=zeros(4,4,'double');
    T_R_matrix(1,1)=1;
    T_R_matrix(2,2)=1;
    T_R_matrix(3,3)=1;

    x_trans=MotionMatrix(1,1); 
    y_trans=MotionMatrix(2,1); 
    z_trans=MotionMatrix(3,1);
    
    T_R_matrix(1:3,1:3)=tmp_xyz_RotMatrix;
    T_R_matrix(4,1:3)=[x_trans, y_trans, z_trans];
    T_R_matrix(4,4)=1;
    tform=affine3d(T_R_matrix);
    
    [R_cor, C_cor, Z_cor]=ndgrid(coilSen_OriCor{1,1}.x_cor, coilSen_OriCor{1,1}.y_cor, coilSen_OriCor{1,1}.z_cor_allSlc); 

    [R_mot, C_mot, Z_mot]=transformPointsForward(tform,R_cor, C_cor, Z_cor);

    [RR, CC, ZZ]=ndgrid(coilSen_OriCor{1,1}.x_cor, coilSen_OriCor{1,1}.y_cor, coilSen_OriCor{1,1}.z_cor_allSlc_append);
    
    for ich=1:prot.chn        
        %real&imaginary interpolation
        tmp_coilSen =interpn(RR,CC,ZZ,coilSen_matrix(:,:,:,ich),R_mot,C_mot,Z_mot,'makima');
        coilSen_Update(:,:,:,ich)=tmp_coilSen(:,:,slcPos_matZ-z_off_mat);
        if IsSmooth==1
            coilSen_Update(:,:,:,ich)=ifft2c(fft2c(coilSen_Update(:,:,:,ich)).*tukey_window);%.*mask_coilSen
        else
            coilSen_Update(:,:,:,ich)=coilSen_Update(:,:,:,ich);%.*mask_coilSen
        end
    end
end