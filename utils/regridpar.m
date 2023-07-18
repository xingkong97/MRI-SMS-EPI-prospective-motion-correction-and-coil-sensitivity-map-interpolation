% The MRI image recon toolkit is a collection of matlab files for reconstructing MRI images acquired using either Cartesian or arbitrary trajectory.
% The toolkit has only been provided to internal Upenn users thus far. Please contact me if any part of it will be delivered to the other parties.

function [rsn_ramp,rsn_dramp,rsn_top]=regridpar(RoLen,rampup,ramptop,rampdown,sampdelay,adcduration)
% Function for calculating the regridding parameters for epi
% reconstruction.
G=1.0;
readoutscans=RoLen;
rampscantime=rampup-sampdelay;
repetition=adcduration/(readoutscans-1);
delT=repetition;
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
tte=tus+(sn_ramp+sn_top-1)*delT;
Sdtts=G*(tup+ttp-tte);
Sdtrue=delK-Sdtts;
Srd=1/2.0*(tdp-stde^2/tdp);
rsn_dramp=fix((Srd-Sdtrue)/delK)+1;

% t_drampend=(readoutscans-1)*delT-rampup-ramptop+sampdelay;
% t_drampstart=rampup+ramptop-(sn_ramp+Ksn_top-1)*delT-sampdelay;
% t_rampend=sampdelay+(sn_ramp-1)*delT;
% 
% % K-space sampling duration corresponding to the ramp duration.
% Ksd_ramp=(t_rampend^2-sampdelay*sampdelay)/(rampup*2); 
% Ksn_ramp=floor(Ksd_ramp/delK)+1;  %The true Kspace samples needed for this duration.
% 
% % K-space sampling duration corresponding to the down ramp duration.
% Ksd_dramp=(-t_drampend^2+t_drampstart^2)/(rampdown*2)+t_drampend-t_drampstart;  
% Ksn_dramp=floor(Ksd_dramp/delK)+1;
return;
