function showimagefft(kspace, fftparam, flipparam, cmap, flipabs) 
%this function is to show 4d image matrix
%author: Bo Li 4-28-2021
%input
%kspace: nx,ny,ncoil,nslice
%fftparam: fft option
%flipparam: flip option
%cmap: gray or jet
%flipabs: flip option in abs
    
    if nargin<=1
        isfft=1;
        isflip=0;
        cmap=gray;
        isflipabs=1;
    elseif nargin==2
            if fftparam==0
                isfft=0;
            else
                isfft=1;
            end
            isflip=0;
            cmap=gray;
            isflipabs=1;
    elseif nargin==3
            isflip=flipparam;
            isfft=fftparam;
            cmap=gray;
            isflipabs=1;
    elseif nargin==4
            isflip=flipparam;
            isfft=fftparam;
            isflipabs=1;
    elseif nargin==5
            isflip=flipparam;
            isfft=fftparam;
            isflipabs=flipabs;
    end
    NSlice=size(kspace,4);
    if NSlice==1
        if isflip==1
            imshow(flip(sqrtSum(kspace(:,:,:),isfft,isflipabs)'),[],'colormap', cmap);
        else
            imshow(sqrtSum(kspace(:,:,:),isfft,isflipabs),[], 'colormap', cmap);
        end
        set(gca,'xtick',[],'xticklabel',[]);
        set(gca,'ytick',[],'yticklabel',[]);
        return;
    end
    
    if mod(NSlice,7)==0
        nrow=ceil(NSlice/7);
    elseif mod(NSlice,4)==0
        nrow=NSlice/4;
    elseif mod(NSlice,3)==0
        nrow=NSlice/3;
    elseif mod(NSlice,5)==0
        nrow=NSlice/5;
    elseif mod(NSlice,2)==0
        nrow=NSlice/2;
    end
    ncol=NSlice/nrow;%ncoil nslice
    ha = tight_subplot(nrow,ncol,   [   .01         .005    ],  [ 0.05  .015], [ 0.002   .002]);
    %                  row, col, [gap_height, gap_width],  [bottom  top], [left    right]
    %figure;
    for irow=1:nrow
        for icol=1:ncol
            iIndex=(irow-1)*ncol+icol;
            axes(ha(iIndex));
            if isflip==1
                if cmap==gray
                    imagesc(flip(sqrtSum(kspace(:,:,:,iIndex),isfft, isflipabs)'));colormap gray;
                else
                    imagesc(flip(sqrtSum(kspace(:,:,:,iIndex),isfft,isflipabs)'));colormap jet;
                end
                %imshow(flip(sqrtSum(kspace(:,:,:,iIndex),isfft)'),[],'colormap', cmap);
            else
                if cmap==gray
                    imagesc(flip(sqrtSum(kspace(:,:,:,iIndex),isfft,isflipabs)'));colormap gray;
                else
                    imagesc(flip(sqrtSum(kspace(:,:,:,iIndex),isfft,isflipabs)'));colormap jet;
                end
                %imshow(sqrtSum(kspace(:,:,:,iIndex),isfft),[],'colormap', cmap);
            end
            %imagesc(squeeze(kspace(:,:,:,iIndex)));colormap gray;
            %for single slice images
            %imagesc(fliplr(sqrtSum(kspace(:,:,:,iIndex),isfft)'));colormap gray;
            set(gca,'xtick',[],'xticklabel',[]);
            set(gca,'ytick',[],'yticklabel',[]);
            %imagesc((sqrtSum(kspace(:,:,:,iIndex),isfft)));colormap gray;%
            %imagesc(abs((kspace(:,:,:,iIndex))));colormap gray;%
        end
    end
end