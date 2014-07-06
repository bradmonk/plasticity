function scrollwheel_view
% Example of using a figure WindowScrollWheelFcn callback to rotate
% or zoom a graph, and change font size of an edit text uicontrol.
% Scrolling the uiciontrol text up and down is built in, and the
% myscroll callback does not change or disable that behavior.
% The callback overrides scrolling when right-clicking in the edit box.
% To toggle the text box scrolling behavior, you must click outside the 
% text box first.

% Create a figure with an axes and a edit text box
hf = figure('Visible','off','Units','pixels','Color',[.9 .8 .7]);
set(hf,'Name','Scrollwheel Demo','NumberTitle','off')
ha = axes('Units','normalized','Position',[.025 .2, .95 .75]);
he = uicontrol('Style','edit','Min',1,'Max',100,...
    'Units','normalized','Position',[.025 .025 .95, .15],...
    'FontUnits','points','HorizontalAlignment','left',...
    'BackgroundColor',[.95 .85 .75]);
% Plot a 3-D surfaceplot and customize its appearance
hs = surf(ha,peaks(100),'FaceColor','interp','EdgeColor','none');
colormap('hot')
axis('vis3d','off')
lighting('phong')
camlight('left')
view(25,20)
title('Left-click and scroll to rotate; right-click and scroll to zoom.')
% Add some help text to the edit text box
set(he,'String',help('peaks'))
% Install scroll wheel callback in figure
set(hf,'WindowScrollWheelFcn',{@myscroll,ha,hs,he});
set(hf,'Visible','on')