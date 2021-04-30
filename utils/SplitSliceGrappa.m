function res=SplitSliceGrappa(kspace, calib, kernel_size, lamda, IsSplit)
%(Split-)slice grappa reconstruction for SMS scan. 
%the MATLAB code was translated from python version split slice-GRAPPA (pygrappa)
%author: Bo Li 4-28-2021

%input
%kspace: aliasing kspace data
%calib: calibration data
%knernel_size: [5 5] default
%lamda: 0.01 default
%IsSplit: 0 slice grappa; 1 split-slice grappa

%output
%res: recon images

    if nargin<5
        IsSplit=1;
    end

    [nx,ny,ncoil,nt]=size(kspace);
    kx=kernel_size(1);ky=kernel_size(2);kx2=floor(kx/2);ky2=floor(ky/2);
    [cx,cy,ncoil,cslice]=size(calib);
    
    kspace=padarray(kspace,[kx2 ky2]);
    calib=padarray(calib,[kx2 ky2]);
    
    W=zeros(cslice,kx*ky*ncoil, ncoil,'double');
    
    for sl=1:cslice
        %T coil sensitivity map, and reshape(nx*ny, ncoil)
        T=calib(kx2+1:end-kx2, ky2+1:end-ky2,:,sl);
        T=reshape(T, nx*ny, ncoil);
        
        if IsSplit==0
            %generate 5x5 coil sensitivity map for each point and reshape(kx*ky*nc)
            S=zeros(nx,ny,kx,ky,ncoil,'double'); 
            for irow=kx2+1:nx+kx2
                for icol=ky2+1:ny+ky2
                        S(irow-kx2,icol-ky2,:,:,1:end)=kspace(irow-kx2:irow+kx2,icol-ky2:icol+ky2,1:end);
                 end
            end
            S=reshape(S,[nx*ny,kx*ky*ncoil]);
            ShS=S'*S;
            ShT=S'*T;
            lamda0=lamda*norm(ShS)/size(ShS,1);
            W(sl,:,:)=pinv(ShS+lamda0*eye(size(ShT,1),size(ShT,1)))*ShT;
        else
            MhM=zeros(kx*ky*ncoil, kx*ky*ncoil,'double');
            for jj=1:cslice
                %generate 5x5 coil sensitivity map for each point and reshape(kx*ky*nc)
                calib0=zeros(nx,ny,kx,ky,ncoil,'double'); 
                for irow=kx2+1:nx+kx2
                    for icol=ky2+1:ny+ky2
                            calib0(irow-kx2,icol-ky2,:,:,1:end)=calib(irow-kx2:irow+kx2,icol-ky2:icol+ky2,1:end,jj);
                    end
                end
                calib0=reshape(calib0,[nx*ny,kx*ky*ncoil]);
                MhM = MhM+calib0'*calib0;   
                if jj==sl
                    Mz=calib0;
                end
            end
            MhT=Mz'*T;
            lamda0=lamda*norm(MhM)/size(MhM,1);
            W(sl,:,:)=pinv(MhM+lamda0*eye(size(MhT,1),size(MhT,1)))*MhT;
        end
    end
    %get each slice from aliasing image
    res=zeros(nx,ny,ncoil,nt,cslice,'double');
    S=zeros(nx,ny,kx,ky,ncoil,nt,'double'); 
    for irow=kx2+1:nx+kx2
        for icol=ky2+1:ny+ky2
            S(irow-kx2,icol-ky2,:,:,1:end,1:end)=kspace(irow-kx2:irow+kx2,icol-ky2:icol+ky2,:,:,1:end,1:end);
         end
    end
    S=reshape(S,[nx*ny,kx*ky*ncoil,nt]);
    for tt=1:nt
        for sl=1:cslice
            tmp=squeeze(S(:,:,tt))*squeeze(W(sl,:,:));
            res(:,:,:,tt,sl)=reshape(tmp,[nx,ny,ncoil]);
        end
    end
    
