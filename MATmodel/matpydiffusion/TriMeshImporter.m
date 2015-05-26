%% TriMeshImporter.m

clc; close all; clear all;
cd(fileparts(which(mfilename)));

%% --- GET 3D XYZ POINTS FOR DENDRITE SKELETON

OBJ=read_wobj('dendrite.obj');
XYZ = OBJ.vertices;
XYZcell = {OBJ.vertices};


%% --- PERFORM TRIANGULATION
% DT = DelaunayTri(XYZ);          % Create the tetrahedral mesh
% hullFacets = convexHull(DT);    % Find the facets of the convex hull

TR = delaunayTriangulation(XYZ);


%% --- PROCESS TRIANGULATION

[FBtri,FBpoints] = freeBoundary(TR);

TRPoints = TR.Points;
TRConLst = TR.ConnectivityList;
EDGETRIS = edges(TR);
EDGEATTC = edgeAttachments(TR,EDGETRIS);
NEIGHBRS = neighbors(TR);
VTXATTCH = vertexAttachments(TR);

bndry = boundary(XYZ,1);

%% ---- CREATE ALPHA SHAPES

Alpha_Shape = alphaShape(XYZ(:,1),XYZ(:,2),XYZ(:,3));
Alpha_Shape.Alpha = Alpha_Shape.Alpha * 1.3;
[Alpha_Shape_bf, Alpha_Shape_pts] = boundaryFacets(Alpha_Shape);


%% ---- PLOT SURFACE MESH

    f1 = figure(1);
        set(f1,'OuterPosition',[450 100 1100 900],'Color',[1 1 1]);
    hax1 = axes('Position',[.05 .05 .90 .90],'Color','none');
        xlim([-100 100]); ylim([-100 500]); zlim([-100 180]);
        xlabel('µm (x)'); ylabel('µm (y)'); zlabel('µm (z)');
        vw=[65 20]; view(vw); hold on;


axes(hax1);
scatter3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'.');
view(vw); axis equal; 

    axes(hax1)
hts1=trisurf(Alpha_Shape_bf,Alpha_Shape_pts(:,1),Alpha_Shape_pts(:,2),Alpha_Shape_pts(:,3), 'FaceColor',[.1 .5 .9],'FaceAlpha', 0.6); 

    axis equal; title('trisurf of boundaryFacets from alphaShape');
    light('Position',[1000 -80 10]); %lighting flat;% lighting gouraud; 
    set(hts1,'FaceLighting','flat','EdgeLighting','gouraud');
    set(hts1,'FaceLighting','gouraud','EdgeLighting','flat');
    shading interp; colormap('hot'); hold on;
    % axis off; set(gca,'CameraPosition',[1424 -443 267]);
    % set(gcf,'Renderer','painters','PaperPositionMode','auto');
    % print(gcf,'-dpdf' ,'-painters','/Users/bradleymonk/Desktop/Dendrite')
    % saveas(gcf,'/Users/bradleymonk/Desktop/Dendrite', 'svg')

SCALE_FACTOR = 100;
hax1.XTickLabel = hax1.XTick / SCALE_FACTOR;
hax1.YTickLabel = hax1.YTick / SCALE_FACTOR;
hax1.ZTickLabel = hax1.ZTick / SCALE_FACTOR;


%% ---- CREATE XML FILE FOR DOLFIN

% dolfinXMLalpha(XYZ,'alphaTriangulation',[1.3, 1.4, 1.4])
dolfinXML(XYZ,'delaunayTriangulation')


%%
