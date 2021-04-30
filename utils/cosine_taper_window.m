function tukey_window = cosine_taper_window(k,kc,w,epsilon,k1,option)
% example: tukey_window = cosine_taper_window(128,108,20,2,108,1);1=square;
% example: tukey_window = cosine_taper_window(64,50,20,2,50,2);2=round;
% example: tukey_window = cosine_taper_window(32,25,20,2,25,2);64*64;
% option :
% 1 for 2D separable Tukey window
% 2 for 2D nonseparable(radial) Tukey window
% 3 for 2D nonseparable Tukey window with ellipticallly varying cutoff

D = 2*k;
tukey_window = zeros(D);
tukey_window_1D = zeros(1,D);
if option == 1
    for counter = -(k-1):k
        if abs(counter) < kc
            tukey_window_1D(counter + k) = 1;
        elseif abs(counter) >= (kc + w)
            tukey_window_1D(counter + k) = 0;
        else
            tukey_window_1D(counter + k) = (cos(pi*(abs(counter)-kc)/(2*w))).^2;
        end        
    end
    tukey_window = tukey_window_1D'*tukey_window_1D;
elseif option == 2
    for kx = -(k-1):k
        for ky = -(k-1):k
            counter = sqrt(kx.^2 + ky.^2);
            if abs(counter) < kc
                tukey_window(kx + k,ky + k) = 1;
            elseif abs(counter) >= (kc + w)
                tukey_window(kx + k,ky + k) = 0;
            else
                tukey_window(kx + k,ky + k) = (cos(pi*(abs(counter)-kc)/(2*w))).^2;
            end   
        end
    end
elseif option == 3
    for kx = -(k-1):k 
        for ky = -(k-1):k            
            kc = k1 * sqrt((kx.^2 + ky.^2)/(kx.^2 + ky.^2/epsilon));
            counter = sqrt(kx.^2 + ky.^2);
            if abs(counter) < kc
                tukey_window(kx + k,ky + k) = 1;
            elseif abs(counter) >= (kc + w)
                tukey_window(kx + k,ky + k) = 0;
            else
                tukey_window(kx + k,ky + k) = (cos(pi*(abs(counter)-kc)/(2*w))).^2;
            end   
        end
    end
    tukey_window(k,k) = 1;
end

