function showimage(kspace) %kspace:nx,ny,ncoil,nslice
    
    NSlice=size(kspace,4);
    if NSlice==1
        imagesc(flip(sqrtSum(kspace(:,:,:),isfft)'));colormap gray;
        set(gca,'xtick',[],'xticklabel',[]);
        set(gca,'ytick',[],'yticklabel',[]);
        return;
    end
    
    if mod(NSlice,4)==0
        nrow=NSlice/4;
    elseif mod(NSlice,3)==0
        nrow=NSlice/3;
    elseif mod(NSlice,5)==0
        nrow=NSlice/5;
    elseif mod(NSlice,2)==0
        nrow=NSlice/2;
    end
    ncol=NSlice/nrow;%ncoil nslice
    ha = tight_subplot(nrow,ncol,   [   .07         .005    ],  [ 0.05  .015], [ 0.002   .002]);
    %                  row, col, [gap_height, gap_width],  [bottom  top], [left    right]
    %figure;
    for irow=1:nrow
        for icol=1:ncol
            iIndex=(irow-1)*ncol+icol;
            axes(ha(iIndex));
            %imagesc(flip(sqrtSum(kspace(:,:,:,iIndex),isfft)'));colormap gray;
            imagesc(flip(fliplr(squeeze(kspace(:,:,:,iIndex)))));colormap gray;
            %for single slice images
            %imagesc(fliplr(sqrtSum(kspace(:,:,:,iIndex),isfft)'));colormap gray;
            set(gca,'xtick',[],'xticklabel',[]);
            set(gca,'ytick',[],'yticklabel',[]);
            %imagesc((sqrtSum(kspace(:,:,:,iIndex),isfft)));colormap gray;%
            %imagesc(abs((kspace(:,:,:,iIndex))));colormap gray;%
        end
    end
end