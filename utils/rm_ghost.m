function kdataout=rm_ghost(kdatain,NSlices, prot,phasecor,lAccelFactPE)

kdataout=zeros(size(kdatain));
if nargin<=4
    rf=prot.lAccelFactPE;
else
    rf=lAccelFactPE;
end


if rf==1
%     acquiredlines=size(kdataout,2);
%     kxmax=size(kdataout,1);
%     dim_y=size(kdataout,2);
%     kymax=dim_y;
    acquiredlines=size(kdataout,1);
    kxmax=size(kdataout,2);
    dim_y=size(kdataout,1);
    kymax=dim_y;
else
    centerline=prot.CenterLineNo;
    kxmax=size(kdataout,2);
    %acquiredlines=prot.Nphase;%prot.iNoOfFourierLines;
    acquiredlines=size(kdataout,1);%centerline;
    kymax=centerline*2;
end
centy=kymax/2;
%acquiredlines=size(kdataout,1);%prot.iNoOfFourierLines;
IsInterp=1;

[Ksn_ramp,Ksn_dramp,Ksn_top]=regridpar(kxmax,prot.rampup,prot.ramptop,prot.rampdown,prot.sampdelay,prot.adcduration);
ROftLen=Ksn_ramp+Ksn_dramp+Ksn_top;
kleftEdge=fix( (kxmax-ROftLen)/2);

ktmp_ph=zeros(acquiredlines,kxmax);
for sl=1:size(kdatain,4)
    for ch=1:prot.chn
        ktmp=kdatain(:,:,ch,sl);
%         ktmp_even=fliplr(ktmp(2:2:end,:));
%         ktmp_1=ktmp;
%         ktmp_1(2:2:end,:)=ktmp_even;
        reftmp=permute(squeeze(phasecor(:,ch,sl,:)),[2,1]);
        %ktmp=ktmp.';
        [kout]=kspcorrectV6(ktmp,reftmp,prot.rampup,prot.ramptop,prot.rampdown,prot.sampdelay,prot.adcduration,0,1,IsInterp);
        ktmp_ph(1:acquiredlines,kleftEdge+(1:ROftLen))=kout;
        kdataout(:,:,ch,sl)=ktmp_ph;
        ktmp_ph(:)=0;
        
    end
end