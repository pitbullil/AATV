function [Recon,X,Y] = backproject_luis(sigMat,n,r_sensor,angle_sensor,c0,image_select,t,fs,image_width,xc,yc)

% parameters:
% 
% sigMat - matrix whose columns correspond to the signal measured by each transducer
% n - resolution of the reconstructed image
% r_sensor - distance of the sensor to the centre of rotation [m]
% angle_sensor - vector containing the angular position of the sensor (or sensors) [rad]
% c0 - speed of sound (assumed uniform in tissue and water) [m/s]
% image_select - type of reconstruction - 'direct' the pressure is back-projected, 'derivative' the derivative of the pressure is back-projected, 'full' both terms are back-projected
% t - vector containing the sampling instants [s]
% fs - sampling frequency [Hz]
% image_width - width of the reconstructed image [m]
% xc - horizontal position of the centre of the reconstructed image [m]
% yc - vertical position of the centre of the reconstructed image [m]

if( isempty( strmatch(image_select, {'full' 'direct' 'derivative'}, 'exact') ) )
    warning('MATLAB:valueNotAllowed', 'No valid backprojection image mode selected, using full.');
    image_select = 'full';
end;

sizeT = length(t); % length of the time vector

bp_image = zeros(n,n); % reconstructed image (inizialized)

Dxy = image_width/(n-1); % spatial sampling distance of the reconstructed image
x = [(-1)*(n/2-0.5)*Dxy:Dxy:(n/2-0.5)*Dxy]; x = x+xc; % position of pixels of the reconstructed image in the x direction
y = [(-1)*(n/2-0.5)*Dxy:Dxy:(n/2-0.5)*Dxy]; y = y+yc; % position of pixels of the reconstructed image in the y direction
[X,Y] = meshgrid(x,y); % position of the pixels of the reconstructed image (mesh)

for i = 1:length(angle_sensor)
    
    T_Pt = sigMat(:,i); % term corresponding to the pressure in the back-projection algorithm
    T_dPdt = [T_Pt(2:sizeT)'-T_Pt(1:sizeT-1)' 0]'.*t'*fs; % term corresponding to the derivative of the pressure in the back-projection algorithm
    
    x_sensor = r_sensor*cos(angle_sensor(i)); % horizontal position of the transducer
    y_sensor = r_sensor*sin(angle_sensor(i)); % vertical position of the transducer
    
    xd = X-x_sensor; % horizontal distance of the transducer to the points of the mesh
    yd = Y-y_sensor; % vertical distance of the transducer to the points of the mesh
    d = sqrt(xd.^2+yd.^2); % distance of the transducer to the points of the mesh
    
    tbp = round((d*fs)/c0-t(1)*fs+1); % times required for the acoustic wave to propagate the distance between the sensor and the grid points 
    tbp(tbp>sizeT) = sizeT;
    tbp(tbp<=0) = 1;
    
    switch image_select
        case 'direct'
            bp_image = bp_image + T_Pt(tbp); % update of the reconstructed image for the direct case
        case 'derivative'
            bp_image = bp_image - T_dPdt(tbp); % update of the reconstructed image for the derivative case
        case 'full'
            bp_image = bp_image + T_Pt(tbp) - T_dPdt(tbp); % update of the reconstructed image for the full case
        otherwise
            bp_image = bp_image + T_Pt(tbp) - T_dPdt(tbp); % update of the reconstructed image by default (full case)
    end
    
end

Recon = bp_image; % reconstructed image

