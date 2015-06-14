function [varargout] = DiffusionOnMesh(varargin)
%% DiffusionOnMesh.m

% cd(fileparts(which(mfilename)));


if nargin < 1

        lim.x = [-100 400];
        lim.y = [-100 100];
        lim.z = [-100 100];
        lim.v = [-30 20];
        allLims = [lim.x lim.y lim.z lim.v];

elseif nargin == 1

        allLims = varargin{1};

        lim.x = allLims(1:2);
        lim.y = allLims(3:4);
        lim.z = allLims(5:6);
        lim.v = allLims(7:8);

end



doBuildMesh = 0;
doMayaMesh = 0;
doLoadMesh = 1;
doDiffusionOnMesh = 0;

%% ---- LOAD SURFACE MESH


if doLoadMesh

    load('serialized_mesh_res_34.mat');
    PYdverts = all_vertices;
    PYdtets = triangles + 1;

    %-- VERT & TET PRE-PROCESSING 
    % dvertsmin = round(min(min(PYdverts)));
    % dverts = PYdverts - dvertsmin;
    dtets = double(PYdtets);
    FBpoints = PYdverts;
    % FBpoints = dverts;
    FBtri = dtets;

    min(FBpoints)
    min(FBtri)


    f1 = figure(1);
            set(f1,'OuterPosition',[100 100 1200 900],'Color',[1 1 1]);

            hax1 = axes('Position',[.05 .05 .9 .9],'Color','none');
                xlim(lim.x); ylim(lim.y); zlim(lim.z); view(lim.v);
                hold on

            hax2 = axes('Position',[.05 .05 .9 .9],'Color','none');
                xlim(lim.x); ylim(lim.y); zlim(lim.z); view(lim.v);


        axes(hax1)
    hts1 = trisurf(FBtri,FBpoints(:,1),FBpoints(:,2),FBpoints(:,3), ...
           'FaceColor',[.1 .9 .1],'FaceAlpha', 0.3);
            xlabel('µm (x)'); ylabel('µm (y)'); zlabel('µm (z)');

        %axis off
        title('trisurf of boundaryFacets from alphaShape');
        light('Position',[-193.5 10.8 -17.5]);
        set(hts1,'FaceLighting','flat','EdgeLighting','gouraud');
        shading interp; colormap('hot'); hold on;


end; % if doLoadMesh







%%
if doMayaMesh

    % OBJ=read_wobj('intertube.obj');
    % XYZ = OBJ.vertices;
    % AlphShp = alphaShape(XYZ(:,1),XYZ(:,2),XYZ(:,3), 3.2);
    % [AlphShpBF, AlphShpPts] = boundaryFacets(AlphShp);
    % [AlphShpTets,AlphShpVrts] = alphaTriangulation(AlphShp);
    % Tets = AlphShpBF;
    % Vrts = AlphShpPts;

end

%% --- GET MESH OBJECTS


if doBuildMesh

    [XYZ,XYZcell] = TriMeshBuilderFunc();

    TR = delaunayTriangulation(XYZ);

    % TRPoints = TR.Points;
    % TRConLst = TR.ConnectivityList;
    % EDGETRIS = edges(TR);
    % EDGEATTC = edgeAttachments(TR,EDGETRIS);
    % NEIGHBRS = neighbors(TR);
    % VTXATTCH = vertexAttachments(TR);
    % [FBtri,FBpoints] = freeBoundary(TR);
    % bndry = boundary(XYZ,1);

    AlphShp = alphaShape(XYZ(:,1),XYZ(:,2),XYZ(:,3), 1.2);
    % AlphShp.Alpha = AlphShp.Alpha * 1.3;
    [AlphShpBF, AlphShpPts] = boundaryFacets(AlphShp);
    [AlphShpTets,AlphShpVrts] = alphaTriangulation(AlphShp);


    DENshp = alphaShape(XYZcell{1}(:,1),XYZcell{1}(:,2),XYZcell{1}(:,3));
    DENshp.Alpha = DENshp.Alpha * 1.3;
    [DENbf, DENpts] = boundaryFacets(DENshp);

    SPIshp = alphaShape(XYZcell{2}(:,1),XYZcell{2}(:,2),XYZcell{2}(:,3));
    SPIshp.Alpha = SPIshp.Alpha * 1.4;
    [SPIbf, SPIpts] = boundaryFacets(SPIshp);

    SPI2shp = alphaShape(XYZcell{3}(:,1),XYZcell{3}(:,2),XYZcell{3}(:,3));
    SPI2shp.Alpha = SPI2shp.Alpha * 1.4;
    [SPI2bf, SPI2pts] = boundaryFacets(SPI2shp);



    Tets = DENbf;

    Vrts = DENpts;


    % PLOT SURFACE MESH
        f1 = figure(1);
            set(f1,'OuterPosition',[450 100 1100 900],'Color',[1 1 1]);
        hax1 = axes('Position',[.05 .05 .90 .90],'Color','none');
            xlabel('µm (x)'); ylabel('µm (y)'); zlabel('µm (z)');
            xlim(lim.x); ylim(lim.y); zlim(lim.z); view(lim.v);
            hold on;


        axes(hax1);
    scatter3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'.');
        view(lim.v); axis equal; 

        axes(hax1)
    hts1=trisurf(Tets,Vrts(:,1),Vrts(:,2),Vrts(:,3));

        set(hts1,'FaceColor',[.1 .5 .9],'FaceAlpha', 0.6)

    axis equal; title('trisurf of boundaryFacets from alphaShape');
        light('Position',[1000 -80 10]);
        set(hts1,'FaceLighting','flat','EdgeLighting','gouraud');
        %set(hts1,'FaceLighting','gouraud','EdgeLighting','flat');
        shading interp; colormap('hot'); hold on;

        SCALE_FACTOR = 100;
        hax1.XTickLabel = hax1.XTick / SCALE_FACTOR;

end; % if doBuildMesh








%% ---- RUN DIFFUSION ON MESH USING PATH DATA ---
if doDiffusionOnMesh

    load('DiffusionData3600.mat');

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

end; % if doDiffusionOnMesh

%% -- Analysis





varargout = {f1, hax1, hax2, allLims};
%%
end