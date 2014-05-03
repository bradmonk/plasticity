% Diffusion with membrane

% Determine the grid size and generate the diffusion.
  N = input('Enter Time:');
  nRedMolecules=4;  %Number of K Particles
  permRed=.001;  %permeability of K ions
  nBlueMolecules=4;
  permBlue=.5; %
  membrane=2; %Membrane at x=2
  [a b]=DiffusionMem(N,membrane,permRed);
  xRed=[a];
  xRedLength=length(xRed);
  yRed=[b];
  for red=2:nRedMolecules
      [c d] = DiffusionMem(N,membrane,permRed);
      xRed=[xRed;c];
      yRed=[yRed;d];
      red= red + 1;
  end
  [e f]=DiffusionMem(N,membrane,permBlue); %Initial array for Na particles
  xBlue=[e]; %Initial X matrix for Na Particles
  xBlueLength=length(xBlue);
  yBlue=[f]; %Initial y Martix for Na particles
  for blue=2:nBlueMolecules
      [g h] = DiffusionMem(N,membrane,permBlue);
      xBlue=[xBlue;g]; %Concatinate each x row
      yBlue=[yBlue;h]; %Concatinate each y row
      blue = blue + 1; %Update Index
  end
    
% Create the figure window
  close all
  figure
  set(gcf,'position',[150 50 600 600])
  hold on
  
% Draw the boundary
  M = N+.5;

% Set the axes
  axis([-M M -M M])
  axis equal square manual
% Label in and out
  yText=M-1;
  xTextMin=-M+.5;
  xTextMax=M-3.5;
  text(xTextMin,yText,'Inside Cell')
  text(xTextMax,yText,'Outside Cell')
% Set membrane
  xMem=membrane;
  yMin=-M;
  yMax=M;
  % Animates the diffusion...
  % mov=avifile('DiffMem1s.avi');
  for col=2:xRedLength %By column
      plot([xMem xMem],[yMin yMax],'-c')
	  drawnow
      for row=1:nRedMolecules %By Row
          plot(xRed(row,col-1),yRed(row,col-1),'.w','Markersize',20)
          plot(xRed(row,col),yRed(row,col),'.r','Markersize',20)
      end
      for row=1:nBlueMolecules
          plot(xBlue(row,col-1),yBlue(row,col-1),'.w','Markersize',20)
          plot(xBlue(row,col),yBlue(row,col),'.b','Markersize',20)          
      end
      % pause(.1)
      % F=getframe(gcf);
      % mov=addframe(mov,F);
  end
  % mov=close(mov);
  
  