function [Recon_BP]=Backprojection(sigMat,n,c,image_width,r_sensor,fs,filter_f,ang_ini,ang_end,ang_step,time_res)
%n = 128;
%c = 1527;  % speed of sound [m/s]
% image_width = 25e-3;                                % width of the image
% r_sensor = 40.5e-3;                                 % distance of the transducer from the centre of rotation [m]
% fs=4e7;
% filter_f = [0.05e6,7e6] ;                           % bandpass filter, lower and upper limits [Hz]
%% 270 degree top
% ang_ini = -0.7853981633;
% ang_end = 3.926990817079;
% ang_step = 0.0184799567858;
% 
% time_res = 2;                                       % time resolution for model-based
% level=4;
% wname='db4';

xsize=[n,n];

dx=image_width/(n-1);
zb=image_width*100/2;

x = (-1)*(n/2-0.5)*dx:dx:(n/2-0.5)*dx;                       
y = (-1)*(n/2-0.5)*dx:dx:(n/2-0.5)*dx; 
[X,Y] = meshgrid(x,y);

len = size(sigMat,1);                               % number of samples in each projection
n_proj = size(sigMat,2);                        % number of projections


ts = 1/fs:1/fs:len/fs;                              % sampling instants

angle_sensor = [ang_ini:ang_step:ang_end];
n_proj=size(angle_sensor,2);
n_angles = 2*n;                                     % number of points for discretizing the curve
limits(1) = r_sensor-(image_width)*sqrt(2)/2;       % limits for the signal
limits(2) = r_sensor+(image_width)*sqrt(2)/2;       % limits for the signal
dx = image_width/n;                                 % increment in x
dt = dx/(time_res*c);                              % increment in t employed to make the model-based reconstruction


%% filter parameters
if( filter_f(2) ~= 0 )
    f_LPF = filter_f(2);
    [b_LPF,a_LPF] = cheby1( 8, .01, 2 * f_LPF/fs * .92 );
end
if( filter_f(1) ~= 0 )
    f_HPF = filter_f(1);
    if( 2 * f_HPF/fs < 0.016 )
        [b_HPF,a_HPF] = cheby1( 2, .01, 2 * f_HPF/fs * 3.3, 'high' );
    else
        [b_HPF,a_HPF] = cheby1( 4, .01, 2 * f_HPF/fs * 1.46, 'high' );
    end
end

% filtering
if( filter_f(2) ~= 0 )
    sigMat = filtfilt(b_LPF, a_LPF, sigMat);
end;
if( filter_f(1) ~= 0 )
    sigMat = filtfilt(b_HPF, a_HPF, sigMat);
end


%% back projection
image_select='full';
xc=0;
yc=0;

tstart=tic;
[Recon_BP,X,Y] = backproject_luis(sigMat,n,r_sensor,angle_sensor,c,image_select,ts,fs,image_width,xc,yc);
save .\wavelet_line\bp Recon_BP
bp_time=toc(tstart)
Recon_BP = -Recon_BP;
