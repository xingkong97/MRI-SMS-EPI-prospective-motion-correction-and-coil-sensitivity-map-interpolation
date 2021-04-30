function  [tukey_window,tukey_window_red]=filter_2D(row)
% ************************************************************
%  Matrix
% ************************************************************
kc=row/4; %radius, center_lines/2
w=20; %100---20
epsilon=2;
k1=kc;%center_lines/2  (32-8-15) (24-8-11)  (16-8-7)
option=2;
tukey_window = cosine_taper_window(row/2,kc,w,epsilon,k1,option);
Img_space=ifft2(fftshift(tukey_window));
tukey_window = abs(fftshift(fft2(fftshift(Img_space))));
%figure; mesh(tukey_window);
tukey_window_red=tukey_window;
