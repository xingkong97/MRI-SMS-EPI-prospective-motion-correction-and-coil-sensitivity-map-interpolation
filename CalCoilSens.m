function CoilSens=CalCoilSens(kDATA,ksize,eigThresh_k,eigThresh_im)
%kDATA, nx,ny,nslc,ncoil
%ksize, ESPIRiT kernel-window-size
%eigThresh_k, threshold of eigenvectors in k-space
%eigThresh_im, threshold of eigenvectors in image domain
%code extract from the demo of ESPIRIT — An Eigenvalue Approach to 
%Autocalibrating Parallel MRI: Where SENSE meets GRAPPA 
%Bo Li 9/22/2022

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

        %% Compute Eigen-Value Maps
        % Maps are computed in two steps. 


        % compute Calibration matrix, perform 1st SVD and convert singular vectors
        % into k-space kernels

        [k,S] = dat2Kernel(calib,ksize);

        idx = max(find(S >= S(1)*eigThresh_k));

        [M,W] = kernelEig(k(:,:,:,1:idx),[sx,sy]);

        CoilSens(:,:,islc,:)=M(:,:,:,end).*repmat(W(:,:,end)>eigThresh_im,[1,1,Nc]);%squeeze(abs(M(:,:,:,end)));
    end    
end


