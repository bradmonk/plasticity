% Diffusion of two molecules.

% Determine the grid size and generate the diffusion.
  N = input('Enter Time:');
  nRedMolecules=4;  %Number of K Particles
  nBlueMolecules=4;
  [a b]=Diffusion(N);
  xRed=[a];
  xRedLength=length(xRed);
  yRed=[b];
  for red=2:nRedMolecules
      [c d] = Diffusion(N);
      xRed=[xRed;c];
      yRed=[yRed;d];
      red= red + 1;
  end
  [e f]=Diffusion(N); %Initial array for Na particles
  xBlue=[e]; %Initial X matrix for Na Particles
  xBlueLength=length(xBlue);
  yBlue=[f]; %Initial y Martix for Na particles
  for blue=2:nBlueMolecules
      [g h] = Diffusion(N);
      xBlue=[xBlue;g]; %Concatinate each x row
      yBlue=[yBlue;h]; %Concatinate each y row
      blue = blue + 1; %Update Index
  end
    
% Create the figure window
  close all
  figure
  set(gcf,'position',[150 50 600 600])
  hold on
 
  
% boundary
  M = 20;
  
% Set the axes
  axis([-M M -M M])
  axis equal square manual
% Animates the diffusion...
mov=avifile('Diff3.avi');
  for col=2:xRedLength %By column
      for row=1:nRedMolecules %By Row
          plot(xRed(row,col-1),yRed(row,col-1),'.w','Markersize',20)
          plot(xRed(row,col),yRed(row,col),'.r','Markersize',20)
      end
      for row=1:nRedMolecules %By Row
          plot(xBlue(row,col-1),yBlue(row,col-1),'.w','Markersize',20)
          plot(xBlue(row,col),yBlue(row,col),'.b','Markersize',20)          
      end
      pause(.1)
      F=getframe(gcf);
      mov=addframe(mov,F);
  end
  mov=close(mov);
