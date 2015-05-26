function [] = dolfinXML(XYZ,type)



%% GET TRIANGULATION POINTS AND CONNECTIVITY LIST

switch type

    case 'delaunayTriangulation'

        TR = delaunayTriangulation(XYZ);
        vrts = TR.Points;
        tets = TR.ConnectivityList;


        % DT = DelaunayTri(XYZ);          % Create the tetrahedral mesh
        % hullFacets = convexHull(DT);    % Find the facets of the convex hull

    case 'alphaTriangulation'

        alph = alphaShape(XYZ(:,1),XYZ(:,2),XYZ(:,3));
        alph.Alpha = alph.Alpha * 1.4;
        [tets,vrts] = alphaTriangulation(alph);
        %[tets, vrts] = boundaryFacets(alph);

    otherwise
        disp('error: bad input type')
end



%%

% keyboard

% 2144

%% CREATE XML FILE

dolfinNode = com.mathworks.xml.XMLUtils.createDocument('dolfin');
xml_dolfin = dolfinNode.getDocumentElement;
xml_dolfin.setAttribute('xmlns:dolfin','bradleymonk.com/Dolfin');


xml_mesh = dolfinNode.createElement('mesh');
xml_mesh.setAttribute('celltype','tetrahedron');
xml_mesh.setAttribute('dim','3');
xml_dolfin.appendChild(xml_mesh);


xml_vertices = dolfinNode.createElement('vertices');
xml_vertices.setAttribute('size',num2str(numel(vrts(:,1))));
xml_mesh.appendChild(xml_vertices);

Vtx_indx = 1:numel(vrts(:,1));
for nn = 1:numel(Vtx_indx)

    vertex = dolfinNode.createElement('vertex');
    
    vertex.setAttribute('index',num2str(Vtx_indx(nn)-1));
    vertex.setAttribute('x',sprintf('%0.15e',vrts(nn,1)));
    vertex.setAttribute('y',sprintf('%0.15e',vrts(nn,2)));
    vertex.setAttribute('z',sprintf('%0.15e',vrts(nn,3)));
    
    % vertex.appendChild(meshNode.createTextNode(num2str(nn)));
    xml_vertices.appendChild(vertex);
end



xml_cells = dolfinNode.createElement('cells');
xml_cells.setAttribute('size',num2str(numel(tets(:,1))));
xml_mesh.appendChild(xml_cells);

Tet_indx = 1:numel(tets(:,1));
for mm = 1:numel(Tet_indx)

    tetrahedron = dolfinNode.createElement('tetrahedron');
    
    tetrahedron.setAttribute('index',num2str(Tet_indx(mm)-1));
    tetrahedron.setAttribute('v0',num2str(tets(mm,1)-1));
    tetrahedron.setAttribute('v1',num2str(tets(mm,2)-1));
    tetrahedron.setAttribute('v2',num2str(tets(mm,3)-1));
    tetrahedron.setAttribute('v3',num2str(tets(mm,4)-1));
    
    % vertex.appendChild(meshNode.createTextNode(num2str(mm)));
    xml_cells.appendChild(tetrahedron);
end


%% EXPORT XML FILE
xmlwrite('mesh_res_32_full.xml',dolfinNode);
% type('mesh_res_32_mat.xml');



%% EXAMPLE OUTPUT FROM DOLPHIN

%{

<?xml version="1.0"?>
<dolfin xmlns:dolfin="http://fenicsproject.org">
  <mesh celltype="tetrahedron" dim="3">
    <vertices size="1870">
      <vertex index="0" x="-5.000000000000000e+01" y="-4.619397662556434e+01" z="-3.086582838174552e+01" />
      <vertex index="1" x="-5.000000000000000e+01" y="-4.042197291960878e+01" z="-2.081837216228893e+01" />
      <vertex index="2" x="-5.000000000000000e+01" y="-3.254910668586541e+01" z="-1.234164554963562e+01" />
      <vertex index="3" x="-5.000000000000000e+01" y="-2.297610052168810e+01" z="-5.859577156852141e+00" />
      <vertex index="4" x="-5.000000000000000e+01" y="-1.218627864236704e+01" z="-1.698403082141455e+00" />
      <vertex index="5" x="-5.000000000000000e+01" y="-7.225567482079125e-01" z="-7.116562813618588e-02" />
      <vertex index="1864" x="1.859184868755403e+01" y="4.894117244518976e+01" z="-6.007788068750699e+01" />
      <vertex index="1865" x="3.222183094478699e+02" y="2.761188043682799e+01" z="-2.143059757669406e+01" />
      <vertex index="1866" x="2.438962742823372e+01" y="1.775357536524982e+01" z="-3.703916026306041e+01" />
      <vertex index="1867" x="1.356967176292075e+02" y="3.917781572457769e+01" z="-5.200056445235793e+01" />
      <vertex index="1868" x="2.119462013994452e+02" y="4.140200098750685e+01" z="-2.201253981819792e+01" />
      <vertex index="1869" x="6.981906748912394e+00" y="1.101749587945243e+01" z="-7.316158960612056e+01" />
    </vertices>
    <cells size="7759">
      <tetrahedron index="0" v0="224" v1="325" v2="1334" v3="1576" />
      <tetrahedron index="1" v0="1111" v1="1201" v2="1265" v3="1427" />
      <tetrahedron index="2" v0="92" v1="129" v2="329" v3="1306" />
      <tetrahedron index="3" v0="1091" v1="1314" v2="1529" v3="1646" />
      <tetrahedron index="4" v0="1280" v1="1357" v2="1390" v3="1810" />
      <tetrahedron index="5" v0="199" v1="260" v2="1694" v3="1718" />
      <tetrahedron index="7751" v0="1171" v1="1431" v2="1448" v3="1867" />
      <tetrahedron index="7752" v0="95" v1="613" v2="1701" v3="1868" />
      <tetrahedron index="7753" v0="1147" v1="1331" v2="1446" v3="1869" />
      <tetrahedron index="7754" v0="1132" v1="1200" v2="1331" v3="1869" />
      <tetrahedron index="7755" v0="1132" v1="1200" v2="1501" v3="1869" />
      <tetrahedron index="7756" v0="1132" v1="1501" v2="1640" v3="1869" />
      <tetrahedron index="7757" v0="1132" v1="1230" v2="1640" v3="1869" />
      <tetrahedron index="7758" v0="1230" v1="1501" v2="1640" v3="1869" />
    </cells>
  </mesh>
</dolfin>

%}





%{
%% CREATE XML FILE

dolfinNode = com.mathworks.xml.XMLUtils.createDocument('dolfin');
xml_dolfin = dolfinNode.getDocumentElement;
xml_dolfin.setAttribute('xmlns:dolfin','bradleymonk.com/Dolfin');


xml_mesh = dolfinNode.createElement('mesh');
xml_mesh.setAttribute('celltype','tetrahedron');
xml_mesh.setAttribute('dim','3');
xml_dolfin.appendChild(xml_mesh);


xml_vertices = dolfinNode.createElement('vertices');
xml_vertices.setAttribute('size',num2str(numel(vrts(:,1))));
xml_mesh.appendChild(xml_vertices);

Vtx_indx = 1:numel(vrts(:,1));
for nn = 1:numel(Vtx_indx)

    vertex = dolfinNode.createElement('vertex');
    
    vertex.setAttribute('index',num2str(Vtx_indx(nn)-1));
    vertex.setAttribute('x',sprintf('%0.15e',vrts(nn,1)));
    vertex.setAttribute('y',sprintf('%0.15e',vrts(nn,2)));
    vertex.setAttribute('z',sprintf('%0.15e',vrts(nn,3)));
    
    % vertex.appendChild(meshNode.createTextNode(num2str(nn)));
    xml_vertices.appendChild(vertex);
end



xml_cells = dolfinNode.createElement('cells');
xml_cells.setAttribute('size',num2str(numel(tets(:,1))));
xml_mesh.appendChild(xml_cells);

Tet_indx = 1:numel(tets(:,1));
for mm = 1:numel(Tet_indx)

    tetrahedron = dolfinNode.createElement('tetrahedron');
    
    tetrahedron.setAttribute('index',num2str(Tet_indx(mm)-1));
    tetrahedron.setAttribute('v0',sprintf('%0.15e',tets(mm,1)));
    tetrahedron.setAttribute('v1',sprintf('%0.15e',tets(mm,2)));
    tetrahedron.setAttribute('v2',sprintf('%0.15e',tets(mm,3)));
    tetrahedron.setAttribute('v3',sprintf('%0.15e',tets(mm,4)));
    
    % vertex.appendChild(meshNode.createTextNode(num2str(mm)));
    xml_cells.appendChild(tetrahedron);
end
%}


end