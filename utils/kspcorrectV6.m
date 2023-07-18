function [out,sn_top]=kspcorrectV6(tmim,ref,rampup,ramptop,rampdown,sampdelay,adcduration,fstdir,fstphdir,IsInterpolation)
% Function for Kspace phase correction and RO interpolation (regridding).
% This function is especially for epi raw data, we assumed there are 3
% lines of reference scan. Inputs are row vectors.
% Made by:
% Ze Wang at center for neuroimaging, Upenn
% zewang@mail.med.upenn.edu 
% All rights reserved, 4-23-2004
% fstphdir:  0 forward 1 reverse
% fstdir:  0 forward 1 reverse  
% tmim is the kspace data to be corrected
% ref is the 3 lines 
% rampup (hdr.Config.RampupTime), top, down(hdr.Config.RampdownTime) are the points sampled during the three duration (readout has a linear rampup, then constant magnitude, then a rampdown time); 
% fstdir means the direction of the first ADC acquisition (right to left or the opposite); 
% fstphdir means the direction of the first line in the phase correction ref scans (3 lines); 
% adcduration (hdr.Config.ADCDuration) is the sampling time (used for calculating the kspace position for the rampup and ramdown points). 
% sampdelay (hdr.Config.DelaySamplesTime) might be related to the time difference between the onset of gradient and ACD (you can read the header if there is or check the code to make sure). 
% IsInterpolation: 1-interpolation for downramp part

clear interp_im
% if nargin<10
type=1;
% end
if (size(ref,1)~=3) 
    return; sprintf('This function only works for 3 reference lines\n'); 
end
[r,c]=size(tmim);

IsInterpolation=0;
%-----------------------------------
% Nonlinear phase correction
fref=fft(ref,[],2);
switch type
    case 1
         revline=(ref(1,:)+ref(3,:))/2.0;
         revline=fft(revline);
         norline=fref(2,:);
         if fstphdir   % the 1st and 3rd are reverse scan
            mtffilter=norline./revline;
        else
            mtffilter=revline./norline;
        end
         mtffilter=mtffilter./abs(mtffilter);
     case 2
        revline=fref(1,:);
        norline=fref(2,:);
        if fstphdir
            mtffilter=norline./revline;
        else
            mtffilter=revline./norline;
        end
        mtffilter=mtffilter./abs(mtffilter);
    otherwise
        if fstphdir
            mtffilter=exp(-j* ( angle(fref(3,:).*conj(fref(2,:)))+angle(fref(1,:).*conj(fref(2,:))))/2);
        else
            mtffilter=exp(-j* ( angle(fref(2,:).*conj(fref(3,:)))+angle(fref(2,:).*conj(fref(1,:))))/2);
        end
end
tmim=fft(tmim,[],2);  % perform fft across column, i.e., row fft

evenindex=2:2:r;
ehr=fix(r/2);    % hr means one half of r.
ohr=ehr;
if mod(r,2)
    oddindex(1,1:ehr)=evenindex-1;
    oddindex(1,ehr+1)=r;
    ohr=ehr+1;
else
    oddindex=evenindex-1;
end
%tmim(evenindex,:)=tmim(evenindex,:).*repmat(ph,row/2,1);
%tmim=tmim.*repmat(freq,row,1);
%tmim(evenindex,:)=tmim(evenindex,:).*repmat(revfilter,row/2,1);
%tmim(oddindex,:)=tmim(oddindex,:).*repmat(norfilter,row/2,1);

%LB note 2-5-2020
if fstdir   %1 reverse
    tmim(oddindex,:)=tmim(oddindex,:).*repmat(mtffilter,ohr,1);
else
    tmim(evenindex,:)=tmim(evenindex,:).*repmat(mtffilter,ehr,1);
end

% revline=revline.*mtffilter;
% norfilter=norline./mnorline;
% norfilter=conj(norfilter)./abs(norfilter);
% revfilter=revline./mnorline;
% revfilter=conj(revfilter)./abs(revfilter);
% if fstdir
%     tmim(oddindex,:)=tmim(oddindex,:).*repmat(revfilter,ohr,1);
%     tmim(evenindex,:)=tmim(evenindex,:).*repmat(norfilter,ehr,1);
% else
%     tmim(evenindex,:)=tmim(evenindex,:).*repmat(revfilter,ehr,1);
%     tmim(oddindex,:)=tmim(oddindex,:).*repmat(norfilter,ohr,1);
% 
% end


tmim=ifft(tmim,[],2);
% --------------------------------------------------------------
% readout interpolation
% The readout gradient for this sequence is assumed to be like:
%          /--------------------\                                   /-----
%        /                       \                                 /
%      /                          \                               /
%-----                             \                             /
%                                   \                           /
%                                    \_________________________/     
readoutscans=c;
rampscantime=rampup-sampdelay;
repetition=adcduration/(readoutscans-1);
delT=repetition;
G=1.0;
delK=delT*G;
tup=rampup;
tdp=rampdown;
tus=sampdelay;
ttp=ramptop;
% Uniform time domain sampling points we get.
%rampscan=round(rampscantime*readoutscans/(rampscantime+ramptop+rampscantime))+1; 
sn_ramp=fix(rampscantime/delT)+1;  % sampling steps during the upgrade ramp.
sn_top=fix( (rampscantime-(sn_ramp-1)*delT+ramptop)/delT);
rsn_top=sn_top;
sn_dramp=readoutscans-sn_ramp-rsn_top;
tts=sampdelay+sn_ramp*delT;  % time of the first sampling point during the plateau period
Stts=G*(tts-tup); % area from trp to tts
Strue=delK-Stts; % area needed in the rampup period for one Kspace sample.
Sru=1/2.0*G*(tup-tus^2/tup);
rsn_ramp=fix((Sru-Strue)/delK)+1;  % sample numbers after regridding

tde=tus+adcduration;
stde=tup+ttp+tdp-tde;
tds=stde;   % start time of rampdown area after flipping
tte=tus+(sn_ramp+sn_top-1)*delT;
Sdtts=G*(tup+ttp-tte);
Sdtrue=delK-Sdtts;
Srd=1/2.0*(tdp-stde^2/tdp);
rsn_dramp=fix((Srd-Sdtrue)/delK)+1;
Kramp_idx=1:rsn_ramp;
Kdramp_idx=1:rsn_dramp;
Srtn=(rsn_ramp-Kramp_idx)*delK+Strue;
trun=sqrt(tup^2-2*Srtn*tup/G);  % true Kspace sampling points in the rampup area
Ssrtn=(rsn_dramp-Kdramp_idx)*delK+Sdtrue;
strdn=sqrt(tdp^2-2*Ssrtn*tdp/G);  % true Kspace sampling points in the rampdown area

interp_im=zeros(r,rsn_ramp+sn_top+rsn_dramp);
interp_im(:,rsn_ramp+1:rsn_ramp+sn_top)=tmim(:,sn_ramp+1:sn_ramp+sn_top);

kinterp=(trun-sampdelay)/delT;
m0=round(kinterp)+1;         % The first index is 1.
% Create interpolation coefficients lookup tables
%m0=round(kinterp);
a=zeros(4,rsn_ramp);
indicator=(kinterp>m0-1);
% I use 4 neighbor points to get the interpolated data.
for ii=1:rsn_ramp
    if m0(ii)==1
        a(:,ii)=m0(ii)+[0 1 2 3]';
        %aa(:,ii)=readoutscan-topscan-rampscan-m0(ii)+[0 1 2 3]';
    elseif m0(ii)==2
        a(:,ii)=m0(ii)+[-1 0 1 2]';
    elseif m0(ii)==sn_ramp
        a(:,ii)=m0(ii)+[-2 -1 0 1]';
    else
        if indicator(ii)>0     %tm-m0+1<0.5
            a(:,ii)=m0(ii)+[-1 0 1 2]';
        else
            a(:,ii)=m0(ii)+[-2 -1 0 1]';
        end
    end
end
b=sinc((repmat(kinterp,[4 1])-(a-1)));
%interp_im=tmim; %interp_im=zeros(r,sn_ramp+sn_top+sn_dramp+1:readoutscans
%interp_im(:,1:sn_ramp-Ksn_ramp)=0;

for ii=1:r
      for jj=1:rsn_ramp
        interp_im(ii,jj)=tmim(ii,a(:,jj))*b(:,jj);
          %interp_im(ii,sn_ramp-Ksn_ramp+jj)=tmim(ii,a(:,jj))*b(:,jj);  %need to be modified
        %interp_dnsamp(ii,jj)=dnsamp(ii,a(:,jj))*b(:,jj);
      end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% interpolation for downramp part
if IsInterpolation==1
    clear a dnsamp b interp_dnsamp
    kinterp=(strdn-tds)/delT;
    m0=round(kinterp)+1;         % The first index is 1.
    a=zeros(4,rsn_dramp);
    % Create interpolation coefficients lookup tables
    indicator=(kinterp>m0-1);
    % I use 4 neighbor points to get the interpolated data.
    for ii=1:rsn_dramp
        if m0(ii)==1
            a(:,ii)=m0(ii)+[0 1 2 3]';
            %aa(:,ii)=readoutscan-topscan-rampscan-m0(ii)+[0 1 2 3]';
        elseif m0(ii)==2
            a(:,ii)=m0(ii)+[-1 0 1 2]';
        elseif m0(ii)==sn_dramp
            a(:,ii)=m0(ii)+[-2 -1 0 1]';
        elseif m0(ii)>sn_dramp
            a(:,ii)=sn_dramp+[-2 -1 0 0]';
        else
            if indicator(ii)>0     %tm-m0+1<0.5
                a(:,ii)=m0(ii)+[-1 0 1 2]';
            else
                a(:,ii)=m0(ii)+[-2 -1 0 1]';
            end
        end
    end
    b=sinc((repmat(kinterp,[4 1])-a+1));
    dnsamp=tmim(:,readoutscans :-1:readoutscans-sn_dramp-4);
    for ii=1:r
      for jj=1:rsn_dramp
            interp_dnsamp(ii,jj)=dnsamp(ii,a(:,jj))*b(:,jj);
      end
    end

%---------------------------------------------
%interp_im(:,sn_ramp+sn_top+1:readoutscans)=0;
%interp_im(:,sn_ramp+sn_top+1:sn_ramp+sn_top+Ksn_dramp)=interp_dnsamp(:,1:Ksn_dramp);
interp_im(:,rsn_ramp+sn_top+1:rsn_ramp+sn_top+rsn_dramp)=interp_dnsamp(:,rsn_dramp:-1:1);
end
out=interp_im;
return;