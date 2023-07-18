function plotTransAndRot(MotionMatrix, measiIndex, IsPlotInOne)
%plot six motion parameters: x,y,z-translation and x,y,z-rotation
%IsPlotInOne: plot six motion parameters to one figure or not
    isOnePlot=0;
    if nargin==1 %for all measurements
        timestamp_star=MotionMatrix(1,1);
        isOnePlot=1;
    else %for one measurement
        measMotionInfo=MotionMatrix; 
        resol=measMotionInfo.BaseResolution;
        Nslc=measMotionInfo.Nslice;
        mcTimeStamp=measMotionInfo.mcTimeStamp_all((measiIndex-1)*Nslc*resol+1:measiIndex*Nslc*resol, 1);
        mcTimeStamp_uniq=unique(mcTimeStamp);
        MotionMatrix=zeros(size(mcTimeStamp_uniq,1), 7, 'double');
        for i=1:size(mcTimeStamp_uniq,1)
            tmpTimeStape=mcTimeStamp_uniq(i,1);
            iIndex=find(measMotionInfo.MotionMatrix_all(:,1)==tmpTimeStape);
            MotionMatrix(i, :)=measMotionInfo.MotionMatrix_all(iIndex(1,1),1:7);
        end
        timestamp_star=measMotionInfo.mcTimeStamp_all(1, 1);
    end
    if nargin==3 
        isOnePlot=IsPlotInOne;
    end
    
    
    XTime=2.5*(MotionMatrix(:,1)-timestamp_star)/1000;
    MotionMatrix(:,2:7)=round(MotionMatrix(:,2:7),1);
    MotionMatrix(:,1)=round(XTime(:),1);
    if isOnePlot==0
        subplot(3,2,1);plot(XTime, MotionMatrix(:,2));title('X translation (mm)');
        xlabel('time (s)');
        subplot(3,2,3);plot(XTime, MotionMatrix(:,3));title('Y translation (mm)');
        xlabel('time (s)');
        subplot(3,2,5);plot(XTime, MotionMatrix(:,4));title('Z translation (mm)');
        xlabel('time (s)');
        subplot(3,2,2);plot(XTime, MotionMatrix(:,5));title('X rotation (°)');
        xlabel('time (s)');
        subplot(3,2,4);plot(XTime, MotionMatrix(:,6));title('Y rotation (°)');
        xlabel('time (s)');
        subplot(3,2,6);plot(XTime, MotionMatrix(:,7));title('Z rotation (°)');
        xlabel('time (s)');
    else
        plot(XTime, MotionMatrix(:,2),'LineWidth',2);hold on;
        plot(XTime, MotionMatrix(:,3),'LineWidth',2);hold on
        plot(XTime, MotionMatrix(:,4),'LineWidth',2);hold on
        plot(XTime, MotionMatrix(:,5),'LineWidth',2);hold on
        plot(XTime, MotionMatrix(:,6),'LineWidth',2);hold on
        plot(XTime, MotionMatrix(:,7),'LineWidth',2);hold on
        hgx=xlabel('time (s)');
        set(hgx,'FontSize',20);
        hgy=ylabel('mm or degree');
        set(hgy,'FontSize',20);
        hg1 = legend('x tran (mm)', 'y tran (mm)','z tran (mm)','x rot (°)','y rot (°)','z rot (°)');
        set(hg1,'FontSize',18);
        set(gca,'FontName','arial','FontSize',20);
        title('Motion Parameters', 'FontSize',20);
    end
    
end
