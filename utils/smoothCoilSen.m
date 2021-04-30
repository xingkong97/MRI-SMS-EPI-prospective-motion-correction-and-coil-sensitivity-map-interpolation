function coilSen_Update=smoothCoilSen(coilSen,prot)        
%normalized coil sensitivity maps
%author: Bo Li, 4-28-2021
    coilSen_Update=coilSen;    
    for islc=1:prot.OriNslice
        tmp_coilSen_abs=sqrtSum(squeeze(coilSen_Update(:,:,islc,:)),0);
        tmp_coilSen=squeeze(coilSen_Update(:,:,islc,:));
        tmp_coilSen=tmp_coilSen./tmp_coilSen_abs;
        tmp_coilSen(isnan(tmp_coilSen)) = 0;
        coilSen_Update(:,:,islc,:)=tmp_coilSen;
    end
end