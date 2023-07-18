function CoilSens=CalCoilSens(kDATA,ksize,eigThresh_k,eigThresh_im)
%kDATA, nx,ny,nslc,ncoil
%ksize, ESPIRiT kernel-window-size
%eigThresh_k, threshold of eigenvectors in k-space
%eigThresh_im, threshold of eigenvectors in image domain
%code extract from ESPIRiT code, Uecker et. al, MRM 2013 DOI 10.1002/mrm.24751
%Bo Li, University of Maryland, Baltimore
%12.20.2020

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
        % crop kernels and compute eigen-value decomposition in image space to get
        % maps
        [M,W] = kernelEig(k(:,:,:,1:idx),[sx,sy]);

        CoilSens(:,:,islc,:)=M(:,:,:,end).*repmat(W(:,:,end)>eigThresh_im,[1,1,Nc]);
    end    
end


