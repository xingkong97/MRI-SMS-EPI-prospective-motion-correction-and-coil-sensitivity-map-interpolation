function coilSen_Update = UpdateCSM(coilSen,prot,out,IsSmooth,array_timestamp,MotionData,IsUpdateCoilSen)
%Update coil sensitivity map from the coilSen_ori
%INPUT: mask_thresh, slice-mask threshold
%       coilSen, the CSM used for interpolation
%       prot, image parameters
%       out, image k-space data
%       IsSmooth, 1 smooth CSM, 0 not
%       array_timestamp, timestamps for the whole scan 
%       MotionData, motion parameters for the whole scan
%       IsUpdateCoilSen, 1 update CSM, 0 not update, but may smooth CSM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT: coilSen_Update, the new CSM by interpolating the old one (coilSen)
%
%%%%%% Written by: Bo Li, University of Maryland, Baltimore
%%%%%% for manuscript "SMS-EPI prospective motion correction 
%%%%%% by real-time phase compensation and coil sensitivity map interpolation"
%%%%%% Created on Sep. 22, 2022

    mask_thresh=0.03;%for slice mask 
    coilSen_thresh=0.03;%for cloi mask
    [tukey_window,tukey_window_red]=filter_2D(prot.Nphase);
    [coilSen_matrix,coilSen_OriCor,coilSen_1,slcPos_matZ,mask_coilSen,mask_slice]=getCoilSen(mask_thresh,coilSen_thresh,tukey_window,0,0,0,coilSen,out.kdata_sliceimg,prot.lMultiBandFactor,prot,IsSmooth);
    %updated CSM multiply mask option
    IsMask=1;
    %CSM normalization
    IsNormalization=1;

    if IsUpdateCoilSen==1
        isSENSE=1;
        %interpolation Approch 1: real&imaginary, 2: magnitude&phase
        InterpApp=1;
        %for human images, prot.InterpSlice=1:prot.OriNslice-3
        %for phantom, prot.InterpSlice=1:prot.OriNslice
        prot.InterpSlice=1:prot.OriNslice-3;

        tmp_MotionMatrix=zeros(6,1);

        coilSen_OriCor_R=coilSen_OriCor;
        %which CSM is used for updating CSM
        %coilSen_ori, coilSen_matrix or coilSen_1
        coilSen_matrix_U=coilSen;

        SMS_MotionCell_1st=cell(prot.Nslice,1);
        SMS_MotionCell=cell(prot.Nslice,1);

        coilSen_Update=zeros(prot.Nread,prot.Nphase,prot.OriNslice,prot.chn);

        %update CSM according to motion parameters
        for iSMS=1:prot.Nslice
            %1st repetition translation and rotation
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
            x_rot=median(x_rot_arr(:))-x_rot_1st; y_rot=median(y_rot_arr(:))-y_rot_1st; z_rot=median(z_rot_arr(:))-z_rot_1st;

            %Coil movement translation x,y,-z, rotation x,y,-z
            x_rot=x_rot; y_rot=y_rot; z_rot=-z_rot;
            x_trans=x_trans; y_trans=y_trans; z_trans=-z_trans;

            x_rot_matrix=[1 0 0;0 cosd(x_rot) sind(x_rot);0 -sind(x_rot) cosd(x_rot)];
            y_rot_matrix=[cosd(y_rot) 0 -sind(y_rot);0 1 0;sind(y_rot) 0 cosd(y_rot)];
            z_rot_matrix=[cosd(z_rot) sind(z_rot) 0;-sind(z_rot) cosd(z_rot) 0;0 0 1];
            f_RotMatrix=@(x,y,z)(x*y*z);
            tmp_xyz_RotMatrix=f_RotMatrix(z_rot_matrix,y_rot_matrix,x_rot_matrix);
            %x,y,-z translation
            tmp_MotionMatrix(1,1)=x_trans;
            tmp_MotionMatrix(2,1)=y_trans;
            tmp_MotionMatrix(3,1)=z_trans;
            %x,y,-z rotation
            tmp_MotionMatrix(4,1)=x_rot;
            tmp_MotionMatrix(5,1)=y_rot;
            tmp_MotionMatrix(6,1)=z_rot;

            tmp_MotionMatrix=round(tmp_MotionMatrix,1);
            coilSen_matrix_Update=getUpdatedCoilSen(tukey_window,tmp_MotionMatrix,coilSen_matrix_U,coilSen_OriCor_R,prot,IsSmooth,InterpApp,tmp_xyz_RotMatrix);
            coilSen_Update_recon=coilSen_matrix_Update;
            iCount=0;
            for islc=(iSMS-1)*prot.lMultiBandFactor+1:iSMS*prot.lMultiBandFactor
                iCount=iCount+1;
                iIndex_slc=out.sliceOrderSMS(islc);
                  if IsMask==1
                    coilSen_Update(:,:,iIndex_slc,:)=(coilSen_Update_recon(:,:,iIndex_slc,:)).*mask_slice(:,:,iIndex_slc,:);
                  else
                    coilSen_Update(:,:,iIndex_slc,:)=(coilSen_Update_recon(:,:,iIndex_slc,:));
                  end
            end
        end

        if IsNormalization==1
            coilSen_Update=smoothCoilSen(coilSen_Update,prot,0,tukey_window);
        end
    else
        coilSen_Update=coilSen;
        if IsMask==1
            coilSen_Update=coilSen_Update.*mask_slice;
        end
        if IsNormalization==1
            coilSen_Update=smoothCoilSen(coilSen_Update,prot,0,tukey_window);
        end
    end
    coilSen_Update=permute(coilSen_Update,[1,2,4,3]);
end

