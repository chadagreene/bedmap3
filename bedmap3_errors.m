% This script maps some problem areas in Bedmap3, where the mask indicates
% the presence of ice shelf, yet ice thickness has a zero value. 
% 
% Chad Greene, NASA/JPL, April 3, 2025. 

%% Load data

filename = 'bedmap3.nc'; 
x = ncread(filename,'x'); 
y = ncread(filename,'y'); 
ice_thickness = ncread(filename,'ice_thickness')'; 
mask = ncread(filename,'mask')'; 

%% Plot 
% Uses these three of my: 
% * rgb https://www.mathworks.com/matlabcentral/fileexchange/46872-intuitive-rgb-color-values-from-xkcd?s_tid=srchtitle 
% * maskoverlay https://www.mathworks.com/matlabcentral/fileexchange/98184-maskoverlay?s_tid=srchtitle 
% * labelshelves https://www.mathworks.com/matlabcentral/fileexchange/60246-antarctic-boundaries-grounding-line-and-masks-from-insar 

figure
maskoverlay(x,y,mask==1,'color',0.9*[1 1 1])
hold on
maskoverlay(x,y,mask==2,'color',0.8*[1 1 1])
maskoverlay(x,y,mask==3,'color',0.7*[1 1 1])
maskoverlay(x,y,mask==4,'color',0.6*[1 1 1])

maskoverlay(x,y,mask==3 & ice_thickness==0,'color',rgb('red'),'alpha',0.2)
[~,h] = contour(x,y,ice_thickness==0 & mask==3,[.5 .5],'r'); 
h.Color = rgb('red'); 
axis image
labelshelves('fontsize',6)

%%

% Zoom 1: 
axis([-2088040.00   -1668443.00     812484.00    1232081.00])
exportgraphics(gcf,'bedmap3_errors_zoom_1.jpg','resolution',500)

% Zoom 2: 
axis([277085.43     603554.75   -1817680.25   -1491210.93])
exportgraphics(gcf,'bedmap3_errors_zoom_2.jpg','resolution',500)

% Zoom 3: 
axis([-954356.54    -768308.31   -1367680.41   -1181632.19])
exportgraphics(gcf,'bedmap3_errors_zoom_3.jpg','resolution',500)
