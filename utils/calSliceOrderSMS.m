function sliceOrderSMS=calSliceOrderSMS(Nslice,MultiBandFactor)
    if mod(Nslice,2*MultiBandFactor)==0 %if even, then 2,4,6,...,N,1,3,5,...,N-1. if odd, then 1,3,5,...,N,2,4,6,...,N-1.
            if Nslice==2
                sliceOrderSMS = [2,1];
            else
                tt=Nslice/2/MultiBandFactor;
                for k=2:-1:1
                    slice_even=k:2:Nslice;
                    slice_used=slice_even;

                    for i=2:MultiBandFactor
                        slice_tmp=circshift(slice_even,-(i-1)*tt);
                        slice_used=[slice_used;slice_tmp];

                    end
                    slice_used=slice_used(1:MultiBandFactor,1:tt);
                    if k==2
                        sliceOrderSMS=slice_used(:)';
                    else
                        sliceOrderSMS=[sliceOrderSMS,slice_used(:)'];
                    end
                end
            end
        else
            %for lMultiBandFactor=2
            if mod(Nslice,MultiBandFactor)==0 
                %if mod(idivide(prot.OriNslice,int32(2*lMultiBandFactor), 'fix'),2)==0
                if mod(MultiBandFactor,2)==0
                    loopk=2:-1:1;
                else
                    loopk=1:1:2;
                end

                slice_half=Nslice/MultiBandFactor;
                    for k=loopk %1:1:2
                        slice_odd=k:2:Nslice/MultiBandFactor;                    
                        slice_odd_1=zeros(MultiBandFactor-1,size(slice_odd,2),class(slice_odd));
                        slice_odd_1(1,:)=slice_odd+slice_half;
                        slice_used=cat(1,slice_odd,slice_odd_1(1,:));
                        for mk=2:MultiBandFactor-1
                            slice_odd_1(mk,:)=slice_odd_1(mk-1,:)+slice_half;
                            slice_used=cat(1,slice_used,slice_odd_1(mk,:));
                        end

                        if k==loopk(1)
                            sliceOrderSMS=slice_used(:)';
                        else
                            sliceOrderSMS=[sliceOrderSMS,slice_used(:)'];
                        end
                    end
            end
    end
end