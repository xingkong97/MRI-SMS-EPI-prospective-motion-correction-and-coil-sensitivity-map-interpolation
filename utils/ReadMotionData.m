function [res,array_timestamp]=ReadMotionData(MotionLogfile)
%read motion data
%input: MotionLogfile is the log file save in host
%output: res, res{:,1}.transX;
%             res{:,1}.transY;
%             res{:,1}.transZ;
%             res{:,1}.Euler_deg(1,1); x rotation
%             res{:,1}.Euler_deg(1,2); y rotation
%             res{:,1}.Euler_deg(1,3); z rotation
%        array_timestamp, all timesstamp for imaging, the first one is
%        baseline
%Written by: Bo Li, University of Maryland, Baltimore
%Created on Dec. 20, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    MLogfiletext = textread(MotionLogfile,'%s');
    for iLine=1:size(MLogfiletext,1)
        pos_rotMat=strfind(MLogfiletext{iLine,1},'rotMat');
        if (size(pos_rotMat,2)==9)
            break;
        end
    end
    iLine=iLine+1;
    TotalMotion=size(MLogfiletext,1)-iLine+1;
    %start reading motion data
    res=cell(TotalMotion,1);
    rotMatrix=zeros(3,3);
    array_timestamp=zeros(TotalMotion,1,'double');
    for iMotion=iLine:size(MLogfiletext,1)
        
        motionCell=textscan(MLogfiletext{iMotion,1},'%s','delimiter',',');
        MotionLine=motionCell{1,1};
        timestamp_s=MotionLine{1,1};
        
        if strcmp(timestamp_s,'26795979')
            aaa=1;
        end
        
        transX_mm=str2double(MotionLine{2,1});
        transY_mm=str2double(MotionLine{3,1});
        transZ_mm=str2double(MotionLine{4,1});
        array_timestamp(iMotion-iLine+1,1)=str2double(timestamp_s);
        
          if size(MotionLine,1)==20 %VE11C motion log file
              istart=12;
          elseif size(MotionLine,1)==22 %VE11E motion log file
              istart=14;
          end
          rotMatrix(1,1)=str2double(MotionLine{istart,1});
          rotMatrix(1,2)=str2double(MotionLine{istart+1,1});
          rotMatrix(1,3)=str2double(MotionLine{istart+2,1});
          rotMatrix(2,1)=str2double(MotionLine{istart+3,1});
          rotMatrix(2,2)=str2double(MotionLine{istart+4,1});
          rotMatrix(2,3)=str2double(MotionLine{istart+5,1});
          rotMatrix(3,1)=str2double(MotionLine{istart+6,1});
          rotMatrix(3,2)=str2double(MotionLine{istart+7,1});
          rotMatrix(3,3)=str2double(MotionLine{istart+8,1});

        Euler_deg=RotateMatrixToEuler(rotMatrix);
        
        motion_info=struct('timestamp', timestamp_s, 'transX', transX_mm, 'transY', transY_mm, 'transZ', transZ_mm, ...
                           'Euler_deg', Euler_deg, 'RotMatrix', rotMatrix);
        res{iMotion-iLine+1, 1}=motion_info;
    end
end
