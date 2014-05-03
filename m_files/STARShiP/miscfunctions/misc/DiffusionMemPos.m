 function [x y] = DiffusionMemPos(N,membrane,permeability,position)
% N is a positive integer.
% Simulates a 2D Diffusion with semipermeable membrane
% x and y are row vectors with the property that (x(k),y(k)) is the
% location of the token after N Time.
%Membrane is where the molcules cannot cross
%Permeability is the permeability of the membrane to that molecule
%Position is whether it starts inside or out of the cell

% Initializations...
if (strcmp('In',position))
    k = 0; xc = -10; yc = 0;
else
    k=0; xc=10; yc=0;
end
% In general, (xc,yc) is the location after N time.

% Boundary 
xMinBoundary=-20;
xMaxBoundary=20;
yMinBoundary=-20;
yMaxBoundary=20;

  while k<N 
%    Standing at (xc,yc). 
%    Random determines new location...
     if rand < .5
        if (rand < .5)&&(xc < xMaxBoundary-2)
            if (xc~=membrane)
                xc = xc + 1;
            end
        else if (xc>xMinBoundary+2) && (rand>permeability)   %else if (xc>xMinBoundary+2) || (rand>permeability)
                xc = xc - 1;
            elseif (xc<xMaxBoundary-2)
                xc = xc +1;
            end
        end
     else
        if (rand < .5)&&(yc < yMaxBoundary-1)
           yc = yc + 1;
        else if yc>yMinBoundary+2
                yc = yc - 1;
            else
                yc = yc +1;
            end
        end
     end
%    Save location...
     k = k + 1; x(k) = xc; y(k) = yc; 
  end
