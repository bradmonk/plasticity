clc; close all; clear all;


dverts = getfullverts();
dverts(:,1) = [];
dvertsmin = round(min(min(dverts)));
dverts = dverts - dvertsmin;
% dverts = sortrows(dverts,[-3 1]);

dtets = getfulltets() + 1;
dtets(:,1) = [];


%% ---- CREATE DENDRITIC MESH FROM TET & VERT DATA ---
clc; close all

TR = triangulation(dtets,dverts);

[FBtri,FBpoints] = freeBoundary(TR);


f1 = figure(1);
        set(f1,'OuterPosition',[100 100 1200 900],'Color',[1 1 1]);

azview = -32;   elview = 22;

        hax1 = axes('Position',[.05 .05 .9 .9],'Color','none');
            xlim([0 500]); ylim([0 200]); zlim([0 180])
            view([azview elview]); 
            hold on

        hax2 = axes('Position',[.05 .05 .9 .9],'Color','none');
            xlim([0 500]); ylim([0 200]); zlim([0 180])
            view([azview elview]); 
            

    axes(hax1)
trisurf(FBtri,FBpoints(:,1),FBpoints(:,2),FBpoints(:,3), ...
       'FaceColor',[.1 .9 .1],'FaceAlpha', 0.3);

axis off
light('Position',[-80 -15 40]); shading interp;
% set(gca,'CameraPosition',[208 -50 7687]); lighting phong;


%% ---- RUN DIFFUSION ON MESH USING PATH DATA ---

RunDiffusion = 1;
if RunDiffusion

xyz = getpaths2();
xyz = xyz - dvertsmin;

xyz1 = xyz(:,1:3);


set(f1,'CurrentAxes',hax2);     mec = 'none';  mfc = [1 .2 .1];
p2 = scatter3(hax2,xyz(1,1:3:end),xyz(1,2:3:end),xyz(1,3:3:end),90,'filled',...
             'MarkerEdgeColor',mec,'MarkerFaceColor',mfc);
        xlim([0 500]); ylim([0 200]); zlim([0 180])
            view([azview elview]); 


for nn = 1:size(xyz1,1)

    scatter3(hax2,xyz(nn,1:3:end),xyz(nn,2:3:end),xyz(nn,3:3:end),90,'filled',...
             'MarkerEdgeColor',mec,'MarkerFaceColor',mfc);
        xlim([0 500]); ylim([0 200]); zlim([0 180])
            view([azview elview]); 
        drawnow;

end



end

%% ---- VISUAL CHECK ON MESH COMPOSITION ORDER ---

CheckMesh = 0;
if CheckMesh

    f1 = figure(1);
        set(f1,'OuterPosition',[50 150 1200 900]);

        hax1 = axes('Position',[.05 .05 .9 .9],'Color','none');
            view([azview elview]);  xlim([0 500]); ylim([0 200]); zlim([0 180])
            hold on

        hax2 = axes('Position',[.05 .05 .9 .9],'Color','none',...
                    'XTick', [],'YTick',[]);
            view([azview elview]);  xlim([0 500]); ylim([0 200]); zlim([0 180])
            hold on


          axes(hax1)
    ph1 = scatter3(hax1,dverts(:,1), dverts(:,2), dverts(:,3),'.','MarkerEdgeColor',[.8 .8 .8]);
        % camorbit(20,0)

        axes(hax2)
    ph2 = scatter3(hax2,dverts(1,1),dverts(1,2),dverts(1,3),'b.');
        

    for nn = 1:size(dverts,1)

        ph2 = scatter3(hax2,dverts(nn,1),dverts(nn,2),dverts(nn,3),'b.');
        drawnow;


    end

end

%%
