%%
clc; close all; clear all;


%% --- GET 3D XYZ POINTS FOR DENDRITE SKELETON
DENxyz = TriMeshFunctions('get_dendrite',1);

SPI1xyz = TriMeshFunctions('get_spine',1);
SPI1xyz(:,2) = SPI1xyz(:,2) + 150;

SPI2xyz = TriMeshFunctions('get_spine',2);
SPI2xyz(:,2) = SPI2xyz(:,2) + 300;


XYZ = [DENxyz; SPI1xyz; SPI2xyz];



%% --- PERFORM TRIANGULATION

DT = DelaunayTri(XYZ);          % Create the tetrahedral mesh
TR = delaunayTriangulation(XYZ);

hullFacets = convexHull(DT);       % Find the facets of the convex hull

[FBtri,FBpoints] = freeBoundary(TR);

bndry = boundary(XYZ,1);

shp = alphaShape(XYZ(:,1),XYZ(:,2),XYZ(:,3));
shp.Alpha = shp.Alpha * 1.3;



%% ---- CREATE ALPHA SHAPES

DENshp = alphaShape(DENxyz(:,1),DENxyz(:,2),DENxyz(:,3));
DENshp.Alpha = DENshp.Alpha * 1.3;
[DENbf, DENpts] = boundaryFacets(DENshp);

SPIshp = alphaShape(SPI1xyz(:,1),SPI1xyz(:,2),SPI1xyz(:,3));
SPIshp.Alpha = SPIshp.Alpha * 1.4;
[SPIbf, SPIpts] = boundaryFacets(SPIshp);

SPI2shp = alphaShape(SPI2xyz(:,1),SPI2xyz(:,2),SPI2xyz(:,3));
SPI2shp.Alpha = SPI2shp.Alpha * 1.4;
[SPI2bf, SPI2pts] = boundaryFacets(SPI2shp);
% ----



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
hts1=trisurf(DENbf,DENpts(:,1),DENpts(:,2),DENpts(:,3), 'FaceColor',[.1 .5 .9],'FaceAlpha', 0.6); 

   axes(hax1)
hts2=trisurf(SPIbf,SPIpts(:,1),SPIpts(:,2),SPIpts(:,3), 'FaceColor',[.1 .5 .9],'FaceAlpha', 0.95);

   axes(hax1)
hts3=trisurf(SPI2bf,SPI2pts(:,1),SPI2pts(:,2),SPI2pts(:,3), 'FaceColor',[.1 .5 .9],'FaceAlpha', 0.95);

%    axes(hax1);
%trisurf(bndry,SPIxyz(:,1),SPIxyz(:,2),SPIxyz(:,3),'FaceColor',[.1 .5 .9],'FaceAlpha', 0.4)

    axis equal; title('trisurf of boundaryFacets from alphaShape');
    light('Position',[1000 -80 10]); %lighting flat;% lighting gouraud; 
    set(hts1,'FaceLighting','flat','EdgeLighting','gouraud');
    set(hts2,'FaceLighting','gouraud','EdgeLighting','flat');
    set(hts3,'FaceLighting','gouraud','EdgeLighting','flat');
    shading interp; colormap('hot'); hold on;
    % axis off; set(gca,'CameraPosition',[1424 -443 267]);
    % set(gcf,'Renderer','painters','PaperPositionMode','auto');
    % print(gcf,'-dpdf' ,'-painters','/Users/bradleymonk/Desktop/Dendrite')
    % saveas(gcf,'/Users/bradleymonk/Desktop/Dendrite', 'svg')

SCALE_FACTOR = 100;
hax1.XTickLabel = hax1.XTick / SCALE_FACTOR;
hax1.YTickLabel = hax1.YTick / SCALE_FACTOR;
hax1.ZTickLabel = hax1.ZTick / SCALE_FACTOR;



%% --- NOTES


%{


% --- PLOT TRIANGULATION

        clc; close all

    f1 = figure(1);
    set(f1,'OuterPosition',[400 100 1000 900],'Color',[1 1 1]);
    hax1 = axes('Position',[.03 .53 .45 .45],'Color','none');
    hax2 = axes('Position',[.53 .53 .45 .45],'Color','none');
    hax3 = axes('Position',[.03 .03 .45 .45],'Color','none');
    hax4 = axes('Position',[.53 .03 .45 .45],'Color','none');
    vw = [45 15];

axes(hax1);
scatter3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'.');
view(vw); axis equal; title('scatter3 xyz points');

% axes(hax2);
% tetramesh(DT,'FaceColor',[0.6875 0.8750 0.8984],'FaceAlpha',0.3);
% view(vw); axis equal; axis off; title('tetramesh DelaunayTri');

axes(hax3);
trisurf(hullFacets,DT.X(:,1),DT.X(:,2),DT.X(:,3),'FaceColor','c')
view(vw); axis equal; axis off; title('trisurf convexHull(DT)');

axes(hax4);
trisurf(bndry,XYZ(:,1),XYZ(:,2),XYZ(:,3),'FaceColor',[.1 .9 .9],'FaceAlpha', 0.3)
set(gca,'CameraPosition',[208 -50 7687]); lighting phong;
view(vw); axis equal; axis off; title('boundary(XYZ)');




    f2 = figure(2);
    set(f2,'OuterPosition',[450 100 1000 900],'Color',[1 1 1]);
    hax1 = axes('Position',[.03 .53 .45 .45],'Color','none');
    hax2 = axes('Position',[.53 .53 .45 .45],'Color','none');
    hax3 = axes('Position',[.03 .03 .45 .45],'Color','none');
    hax4 = axes('Position',[.53 .03 .45 .45],'Color','none');
    vw = [45 15];


axes(hax1);
plot(shp,'FaceColor',[.1 .9 .1],'FaceAlpha', 0.3)
view(vw); axis equal; axis off; title('plot alphaShape(XYZ)');

axes(hax2);
trisurf(FBtri,FBpoints(:,1),FBpoints(:,2),FBpoints(:,3), ...
       'FaceColor',[.1 .9 .1],'FaceAlpha', 0.3);
view(vw); axis equal; axis off; title('trisurf freeBoundary(TR)');

axes(hax3);
trisurf(hullFacets,DT.X(:,1),DT.X(:,2),DT.X(:,3),'FaceColor','c')
view(vw); axis equal; axis off; title('trisurf convexHull(DT)');

axes(hax4);
trisurf(bndry,XYZ(:,1),XYZ(:,2),XYZ(:,3),'FaceColor',[.1 .9 .9],'FaceAlpha', 0.3)
set(gca,'CameraPosition',[208 -50 7687]); lighting phong;
view(vw); axis equal; axis off; title('boundary(XYZ)');




% xyzout = dendrite_xyzout;
% xp = xyzout(:,1);
% yp = xyzout(:,2);
% zp = xyzout(:,3);
% 
% shp1 = alphaShape(xp,yp,zp);
% disp(shp1.Alpha)
% shp1.Alpha = shp1.Alpha * 1.1;
% 
% 
% 
% xyzout = spine_xyzout;
% xp = xyzout(:,1);
% yp = xyzout(:,2);
% zp = xyzout(:,3);
% 
% shp2 = alphaShape(xp,yp,zp);
% disp(shp2.Alpha)
% shp2.Alpha = shp2.Alpha * 1.1;
% 
% 
% shp = shp1;
% shp.Points(end+1,:) = [xp, yp, zp];
% 
% % shp = [shp1 shp2];
% 
% figure(1)
% plot(shp,'FaceColor',[.1 .9 .1],'FaceAlpha', 0.3)
% axis vis3d; title('alphaShape');


% figure(1)
% plot3(xp, yp, zp);
% axis square;
% xlabel('x'); ylabel('y'); zlabel('z')
% 
% figure(1)
% scatter3(xp, yp, zp,'.');
% xlabel('x'); ylabel('y'); zlabel('z')
% axis equal;



% plot(shp,'FaceColor',[.1 .9 .9],'FaceAlpha', 0.3)
% light('Position',[-80 -15 40]); lighting phong; %shading interp;
% % set(gca,'CameraPosition',[208 -50 7687]); lighting phong;
% axis vis3d; axis off; title('alphaShape');

% alphaSpectrum	Alpha values giving distinct alpha shapes
% criticalAlpha	Alpha radius defining a critical transition in the shape
% numRegions	Number of regions in alpha shape
% inShape	Determine if point is inside alpha shape
% alphaTriangulation	Triangulation that fills alpha shape
% boundaryFacets	Boundary facets of alpha shape
% perimeter	Perimeter of 2-D alpha shape
% area	Area of 2-D alpha shape
% surfaceArea	Surface area of 3-D alpha shape
% volume	Volume of 3-D alpha shape
% plot	Plot alpha shape

%}


%%

