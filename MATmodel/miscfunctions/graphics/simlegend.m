function [varargout] = simlegend(LABS,varargin)
% Useage
% [varargout] = simlegend(LABS,varargin)

add = numel(varargin);
LABS = fliplr(LABS);

%---------------
if add
	
	hP = varargin{1};
	hAc = get(hP,'Children');
	hL = legend(hAc,LABS);
	[hL,hOb,hP,hTxt] = legend;
	hold on
	
%---------------
else

% 	hP = gcf; hA = gca;
% 	hPc = get(hP,'Children');
% 	hAc = get(hA,'Children');
% 	hMOd = findobj('Marker','d');
% 	hMOo = findobj('Marker','o');
% 	hMOx = findobj('Marker','x');
% 	hMOs = findobj('Marker','s');
% 	hMOv = findobj('Marker','v');
% 	hMObs = cat(1,hMOd,hMOo,hMOx,hMOs,hMOv);
% 	hL = legend(hMObs,LABS);
% 	[hL,hOb,hP,hTxt] = legend;
% 	hold on
	
end
%---------------

varargout = {hL,hOb,hP,hTxt};

end

%{

% set(legend,'Location','Best');
% FOT = findobj('Tag','TMxH','Marker','o');
% hAc99 = get(FOT,'Children');
% get(Ph99,'Children')

child_handles = allchild(handle_list)
	child_handles = allchild(handle_list) returns the list of all children 
	(including ones with hidden handles) for each handle. 
	If handle_list is a single element, allchild returns the output in a vector. 
	If handle_list is a vector of handles, the output is a cell array.
	Compare the results these two statements return:
		axes
		get(gca,'Children')
		allchild(gca)

[hL,hOb,hP,LhTxt] = legend
	hL ? Handle of legend axes
	hOb ? Handles of line/patch/text objects used in the legend
	hP ? Handles of lines and other objects used in the plot
	hTxt ? Cell array of the text strings used in the legend


set(legend,'Location','NorthWest')
legend(...,'Location',[left bottom width height])
	Num	Spec|	Axes Location 	|	Current Spec
	-1		| 	Outside right	|	NorthEastOutside
	0		|	Inside axes		|	Best
	1		|	Upper right		|	NorthEast
	2		|	Upper left		|	NorthWest
	3		|	Lower left		|	SouthWest
	4		|	Lower right		|	SouthEast


pch = get(ph,'Children'); 
legend(pch,CMx)
	Each patch/line object in "pch" is labeled using values 
	from each row of the "CMx" matrix (nums) or cell (txt)


legend(h,M,'parameter_name','parameter_value',...)
set(legend_handle, 'Box', 'off')
set(legend_handle, 'Color', 'none')
set(h,'String',{'cos(x)','sin(x)'})

legend
legend('string1','string2',...)
legend(h,'string1','string2',...)
legend(M)
legend(h,M)
legend(M,'parameter_name','parameter_value',...)
legend(h,M,'parameter_name','parameter_value',...)
legend(axes_handle,...)
legend('off'), legend(axes_handle,'off')
legend('toggle'), legend(axes_handle,'toggle')
legend('hide'), legend(axes_handle,'hide')
legend('show'), legend(axes_handle,'show')
legend('boxoff'), legend(axes_handle,'boxoff')
legend('boxon'), legend(axes_handle,'boxon')
legend_handle = legend(...)
legend(...,'Location','location')
legend(...,'Orientation','orientation') 
[legend_h,object_h,plot_h,text_strings] = legend(...)

%}




