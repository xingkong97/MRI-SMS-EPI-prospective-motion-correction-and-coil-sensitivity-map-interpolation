function CoilSens=CalCoilSens(kDATA,ksize,eigThresh_k,eigThresh_im)
%this function is to calculte coilSen by ESPIRIT method
%kDATA, nx,ny,nslc,ncoil

if nargin<=1
    ksize = [6,6]; % ESPIRiT kernel-window-size
    eigThresh_k = 0.02; % threshold of eigenvectors in k-space
    eigThresh_im = 0.95;
end

    CoilSens=zeros(size(kDATA));
    for islc=1:size(kDATA,3)
        kDATA_slc=squeeze(kDATA(:,:,islc,:));
        kDATA_slc = kDATA_slc/max(max(max(abs(ifft2c(kDATA_slc))))) + eps;

        [sx,sy,Nc] = size(kDATA_slc);

        % create a sampling mask to simulate x2 undersampling with autocalibration
        % lines
        %mask = mask_randm_x4;
        mask=ones(sx,sy);
        mask = repmat(mask,[1,1,Nc]);
        ncalib = getCalibSize(mask);



        kDATAc = kDATA_slc.*mask;
        calib = crop(kDATAc,[ncalib,Nc]);

        % compute Calibration matrix, perform 1st SVD and convert singular vectors
        % into k-space kernels

        [k,S] = dat2Kernel(calib,ksize);

        idx = max(find(S >= S(1)*eigThresh_k));

        %%
        % Display the singular vectors and values of the calibration matrix


        %%
        % crop kernels and compute eigen-value decomposition in image space to get
        % maps


        [M,W] = kernelEig(k(:,:,:,1:idx),[sx,sy]);

        CoilSens(:,:,islc,:)=M(:,:,:,end).*repmat(W(:,:,end)>eigThresh_im,[1,1,Nc]);
    end    
end


