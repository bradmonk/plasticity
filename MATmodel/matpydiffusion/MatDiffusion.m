clc; close all; clear all;
cd(fileparts(which(mfilename)));

%%
clc; close all; clear all;

% %-- GET VERTS AND TETS
% load('serialized_mesh_res_34.mat');
load('serialized_mesh_res_34.mat');
PYdverts = all_vertices;
PYdtets = triangles + 1;

%-- VERT & TET PRE-PROCESSING 
dvertsmin = round(min(min(PYdverts)));
dverts = PYdverts - dvertsmin;
dtets = double(PYdtets);
% TR = triangulation(dtets,dverts);
% [FBtri,FBpoints] = freeBoundary(TR);
FBpoints = dverts;
FBtri = dtets;


min(FBpoints)
min(FBtri)


%% ---- CREATE DENDRITIC MESH FROM TET & VERT DATA ---
clc; close all

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
hts1 = trisurf(FBtri,FBpoints(:,1),FBpoints(:,2),FBpoints(:,3), ...
       'FaceColor',[.1 .9 .1],'FaceAlpha', 0.3);

    %axis off
    title('trisurf of boundaryFacets from alphaShape');
    light('Position',[-193.5 10.8 -17.5]);
    set(hts1,'FaceLighting','flat','EdgeLighting','gouraud');
    shading interp; colormap('hot'); hold on;

%% ---- RUN DIFFUSION ON MESH USING PATH DATA ---

load('vcell3600.mat');

vcell = vcell - dvertsmin;

for ss = 1:size(vcell,1)
    for pp = 1:size(vcell,2)
        for xx = 1:size(vcell,3)

            % k = 6.0;  % Override k.
            % % Run simulation on just a single point.
            % % We need xyz_loc and face_indices for this point.
            % xyz_loc = initial_point;
            % face_indices = [initial_face_index];
            % 
            % [xyz_loc, face_indices] = advance_one_step(...
            %     xyz_loc, face_indices, k, initial_point, initial_face_index, ...
            %     all_vertices, triangles, face_local_bases, neighbor_faces);
            % disp('xyz_loc:'); disp(xyz_loc); disp('face_indices:'); disp(face_indices);

            xyz{ss}(pp,xx) = vcell(ss,pp,xx);

        end
    end
end



    set(f1,'CurrentAxes',hax2);     mec = 'none';  mfc = [.2 .5 .7];
p2 = scatter3(hax2, xyz{1}(:,1), xyz{1}(:,2), xyz{1}(:,3),...
            90,'filled','MarkerEdgeColor',mec,'MarkerFaceColor',mfc);
        xlim([0 500]); ylim([0 200]); zlim([0 180])
            view([azview elview]); 


for nn = 1:size(xyz,2)

    scatter3(hax2, xyz{nn}(:,1), xyz{nn}(:,2), xyz{nn}(:,3),...
            90,'filled','MarkerEdgeColor',mec,'MarkerFaceColor',mfc);
        xlim([0 500]); ylim([0 200]); zlim([0 180])
            view([azview elview]); 
        drawnow;

end


%% -- Analysis
clear xyz1800 zs xx xs

xyzSS = xyz(2601:end);

for ss = 1:size(xyzSS,2)
        xs(:,ss) = xyzSS{ss}(:,1);
        ys(:,ss) = xyzSS{ss}(:,2);
        zs(:,ss) = xyzSS{ss}(:,3);
end

numel(xs)
numel(xs(zs>150))
numel(ys)
numel(ys(zs>150))
numel(zs)
numel(zs(zs>150))


xs = xs(:);
xx = xs( xs>50.1 & xs<449.1 );
hist(xx, 80); hold on


ys = ys(:);
yy = ys( ys>50.1 & ys<149.1 );
hist(ys, 80); hold on

% zzs = reshape(zs,[],3);
% zzs = zzs(:);
% hist(zzs, 80)

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
