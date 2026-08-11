function [] = saveoutputfigs(filename)

thisfile='saveoutputfigs';

% Save current directory
pw = what;
cpath = pw.path;
%pd = dir;
%fileparts()

% Get full path to this file (saveoutputfigs.m)
wfiname = which(thisfile);
% Save just the path without the file name
wthisfile = regexprep(wfiname, '(\w*\.m$)', '','once');
w1 = what(wthisfile);

% Simultaniously Save current path and switch path
pathNow = cd(w1.path);

% Load/Save data file
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,filename,'png');
% set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,'STARShiP1','png');
% saveas(gcf, ['outputfigs/STARShiP1.png']);
% printpreview

% Change path back to currently open folder
cd(pathNow)

end