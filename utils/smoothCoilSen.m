function coilSen_Update=smoothCoilSen(coilSen,prot,IsSmooth,tukey_window)%IsSmooth        
%smooth the coil sensitivity maps
%Bo Li 9/22/2022
if nargin<3
    IsSmooth=0;
end
coilSen_Update=coilSen;    
%normalized sensitivity map
    for islc=1:prot.OriNslice
        tmp_coilSen_abs=sqrtSum(squeeze(coilSen_Update(:,:,islc,:)),0);
        tmp_coilSen=squeeze(coilSen_Update(:,:,islc,:));
        tmp_coilSen=tmp_coilSen./tmp_coilSen_abs;
        tmp_coilSen(isnan(tmp_coilSen)) = 0;
        if IsSmooth==1
            coilSen_Update(:,:,islc,:)=ifft2c(fft2c(tmp_coilSen).*tukey_window);
        else
            coilSen_Update(:,:,islc,:)=tmp_coilSen;
        end
    end
end