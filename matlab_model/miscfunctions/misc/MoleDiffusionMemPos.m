% Molecular Diffusion with Membrane keepin track of Position

% Determine the grid size and generate the diffusion.
  N = input('Enter Time:');
  nRedMolecules=4;  %Number of Red Particles
  permRed=0.1;  %permeability of Potassium ions
  nBlueMolecules=4;
  permBlue=0.5; %permeability of Sodium ions
  nGreenMolecules=4;
  permGreen=0.1; %permeability of Chloride ions
  membrane=0; %Membrane at x=0
  pos='In';
  [a b]=DiffusionMemPos(N,membrane,permRed,pos);
  xRed=[a];
  xRedLength=length(xRed);
  yRed=[b];
  redIn=2;
  for red=2:nRedMolecules
      [c d] = DiffusionMemPos(N,membrane,permRed,pos);
      xRed=[xRed;c];
      yRed=[yRed;d];
      if red>redIn-0.5;
          pos='Out';
      end
  end
  pos='In';
  [e f]=DiffusionMemPos(N,membrane,permBlue,pos); %Initial array for blue particles
  xBlue=[e]; %Initial X matrix for Blue Particles
  xBlueLength=length(xBlue);
  yBlue=[f]; %Initial y Martix for Blue particles
  blueIn=2;
  for blue=2:nBlueMolecules
      [g h] = DiffusionMemPos(N,membrane,permBlue,pos);
      xBlue=[xBlue;g]; %Concatinate each x row
      yBlue=[yBlue;h]; %Concatinate each y row
      if blue>blueIn-0.5
          pos='Out';
      end
  end
  pos='In';
  [i j]=DiffusionMemPos(N,membrane,permGreen,pos); %Initial array for green particles
  xGreen=[i]; %Initial X matrix for Green Particles
  xGreenLength=length(xGreen);
  yGreen=[j]; %Initial y Martix for Green particles
  greenIn=2;
  for green=2:nGreenMolecules
      [k l] = DiffusionMemPos(N,membrane,permBlue,pos);
      xGreen=[xGreen;k]; %Concatinate each x row
      yGreen=[yGreen;l]; %Concatinate each y row
      if green>greenIn-0.5
          pos='Out';
      end
  end
% Create the figure window
  close all
  figure
  set(gcf,'position',[150 50 600 600])
  hold on
   
% Set the boundary
M = 20;
  
% Set the axes
  axis([-20 20 -20 20])
  axis equal square manual
% Label in and out
  yText=M-1;
  xTextMin=-M+.5;
  xTextMax=M-7.5;
  text(xTextMin,yText,'Inside Cell')
  text(xTextMax,yText,'Outside Cell')
% Set membrane
  xMem=membrane;
  yMin=-20;
  yMax=20;
  blueConOut=0; %Conc of blue ions outside the cell
% Animates the diffusion
mov=avifile('MoleDiffMemPos.avi');
  for col=2:xRedLength %By column
      plot([xMem xMem],[yMin yMax],'-c')
      for row=1:nRedMolecules %By Row
          plot(xRed(row,col-1),yRed(row,col-1),'.w','Markersize',20)
          plot(xRed(row,col),yRed(row,col),'.r','Markersize',20)
      end
      for row=1:nBlueMolecules %By Row
          plot(xBlue(row,col-1),yBlue(row,col-1),'.w','Markersize',20)
          plot(xBlue(row,col),yBlue(row,col),'.b','Markersize',20)
      end
      for row=1:nGreenMolecules %By Row
          plot(xGreen(row,col-1),yGreen(row,col-1),'.w','Markersize',20)
          plot(xGreen(row,col),yGreen(row,col),'.g','Markersize',20)          
      end
      % pause(.1)
      F=getframe(gcf);
      mov=addframe(mov,F);
  end
  mov=close(mov);
  
