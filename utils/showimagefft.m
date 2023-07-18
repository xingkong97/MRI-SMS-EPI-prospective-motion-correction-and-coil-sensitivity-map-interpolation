function showimagefft(kspace, fftparam, flipparam, cmap_or_jet, flipabs, isSameScale,IsSaveImageToDisk, folder_name) 
%this function is to show 4d image matrix
%input
%kspace: nx,ny,ncoil,nslice
%fftparam: fft option, default 1
%flipparam: flip option, default 1
%cmap: gray or jet, default gray
%flipabs: flip option in abs, default 1
%isSameScale: sameScale for all images, default 1
%IsSaveImageToDisk: save images to disk, default 0
%folder_name: the folder save images
%written by: Bo Li 4-28-2021

    isSaveImage=0;
    if nargin<=1
        isfft=1;
        isflip=1;
        cmap=0;
        isflipabs=1;
        isSmScale=1;
    elseif nargin==2
            if fftparam==0
                isfft=0;
            else
                isfft=1;
            end
            isflip=1;
            cmap=0;
            isflipabs=1;
            isSmScale=1;
    elseif nargin==3
            isflip=flipparam;
            isfft=fftparam;
            cmap=0;
            isflipabs=1;
            isSmScale=1;
    elseif nargin==4
            isflip=flipparam;
            isfft=fftparam;
            isflipabs=1;
            isSmScale=1;
            if cmap_or_jet==gray
                cmap=0;
            else
                cmap=1;
            end
    elseif nargin==5
            isflip=flipparam;
            isfft=fftparam;
            isflipabs=flipabs;
            isSmScale=1;
            if cmap_or_jet==gray
                cmap=0;
            else
                cmap=1;
            end
    elseif nargin==6
            isflip=flipparam;
            isfft=fftparam;
            isflipabs=flipabs;
            isSmScale=isSameScale;
            if cmap_or_jet==gray
                cmap=0;
            else
                cmap=1;
            end
    elseif nargin==8
            isflip=flipparam;
            isfft=fftparam;
            isSmScale=isSameScale;
            isflipabs=flipabs;
            isSaveImage=IsSaveImageToDisk;
            if cmap_or_jet==gray
                cmap=0;
            else
                cmap=1;
            end
    end
    
    NSlice=size(kspace,4);
    nx=size(kspace,1);
    ny=size(kspace,2);
    
    if isSmScale==1
        sqrtSum_min_max=zeros(size(kspace,1),size(kspace,2),NSlice);
        for islc=1:NSlice
            sqrtSum_min_max(:,:,islc)=sqrtSum(kspace(:,:,:,islc),isfft,isflipabs);
        end
        minV = min(sqrtSum_min_max(:));
        maxV = max(sqrtSum_min_max(:));
        rangeV = maxV - minV;
%         if cmap==0
%             thismap = colormap(gray(256));
%         else
%             thismap = colormap(jet);
%         end
        thismap = colormap(gray(256));
        maxcol = size(thismap, 1) - 1;
        sqrtSum_min_max=uint8(floor((sqrtSum_min_max(:) - minV) ./ rangeV .* maxcol));
        sqrtSum_min_max=reshape(sqrtSum_min_max,[size(kspace,1),size(kspace,2),NSlice]);
    end
    
    
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
    ncol=NSlice/nrow;
    
    if isSmScale==1
        showMatrix=zeros(nrow*nx,ncol*ny);
        for irow=1:nrow
            for icol=1:ncol
                iIndex=(irow-1)*ncol+icol;
                if isflipabs==1
                    showMatrix((irow-1)*nx+1:irow*nx,(icol-1)*ny+1:icol*ny)=flip((sqrtSum_min_max(:,:,iIndex))');
                else
                    showMatrix((irow-1)*nx+1:irow*nx,(icol-1)*ny+1:icol*ny)=sqrtSum_min_max(:,:,iIndex);
                end
                if iIndex==-1
                    img_cur=sqrtSum(kspace(:,:,:,iIndex),isfft,isflipabs)';
                    tmp_img_cur=(img_cur(:)-min(img_cur(:)))/(max(img_cur(:))-min(img_cur(:)));
                    tmp_img_cur=reshape(tmp_img_cur,[size(img_cur,1),size(img_cur,2)]);
                    figure;imshow(flipud(tmp_img_cur),[]);
                end
                if isSaveImage==1
                    filepath=pwd; 
                    cd(folder_name)
                    img_cur=sqrtSum(kspace(:,:,:,iIndex),isfft,isflipabs)';
                    tmp_img_cur=(img_cur(:)-min(img_cur(:)))/(max(img_cur(:))-min(img_cur(:)));
                    tmp_img_cur=reshape(tmp_img_cur,[size(img_cur,1),size(img_cur,2)]);
                    imwrite(flipud(tmp_img_cur),strcat(num2str(iIndex),'.jpg'))
                    cd(filepath)    
                end
            end
        end
        if (cmap==1)
            imshow(showMatrix,[]);colormap jet;colorbar;
        else
            imshow(showMatrix,[]);colormap gray;colorbar;
        end
%         for irow=1:nrow
%             for icol=1:ncol
%                 iIndex=(irow-1)*ncol+icol;
%                 text((icol-1)*ny+10,(irow-1)*nx+10,num2str(iIndex),'FontSize',12,'Color','white')
%             end
%         end
    else
    
    if isSaveImage==0
        ha = tight_subplot(nrow,ncol,   [   .01         .005    ],  [ 0.05  .015], [ 0.002   .002]);
    end
    for irow=1:nrow
        for icol=1:ncol
            iIndex=(irow-1)*ncol+icol;
            if isSaveImage==0
                axes(ha(iIndex));
            end
            if isflip==1
                img_cur=flip(sqrtSum(kspace(:,:,:,iIndex),isfft, isflipabs)');
                if isSaveImage==1
                    filepath=pwd;           %????????
                    cd(folder_name)
                    tmp_img_cur=(img_cur(:)-min(img_cur(:)))/(max(img_cur(:))-min(img_cur(:)));
                    tmp_img_cur=reshape(tmp_img_cur,[size(img_cur,1),size(img_cur,2)]);
                    imwrite(tmp_img_cur,strcat(num2str(iIndex),'.jpg'))
                    cd(filepath)            %???????
                else
                    if cmap==0
                        if isSmScale==0
                            imagesc(img_cur);colormap gray;
                        else
                            imagesc(flip((sqrtSum_min_max(:,:,iIndex))'));colormap(ha(iIndex),gray(256));colorbar;
                        end
                    else
                        imagesc(img_cur);colormap jet;
                    end
                end
            else
                img_cur=sqrtSum(kspace(:,:,:,iIndex),isfft,isflipabs)';
                if isSaveImage==1
                    filepath=pwd; 
                    cd(folder_name)
                    imwrite(img_cur,strcat(num2str(iIndex),'.jpg'))
                    cd(filepath)            
                else
                    if cmap==0
                        imagesc(img_cur);colormap gray;
                    else
                        imagesc(img_cur);colormap jet;
                    end
                end
            end
            if isSaveImage==0
                set(gca,'xtick',[],'xticklabel',[]);
                set(gca,'ytick',[],'yticklabel',[]);
            end
        end
    end
    end
end