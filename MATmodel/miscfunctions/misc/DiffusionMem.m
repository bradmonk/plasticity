 function [x y] = DiffusionMem(N,membrane,permeability)
% N is a positive integer.
% Simulates a 2D Diffusion with semipermeable membrane
% x and y are row vectors with the property that (x(k),y(k)) is the
% location of the token after N Time.
%Membrane is where the molcules cannot cross
% Initializations...
  k = 0; xc = 0; yc = 0;
  
% In general, (xc,yc) is the location after N time.

  while k<N 
%    Standing at (xc,yc). 
%    Rand to determine new location...
     if rand < .5
        if (rand < .5&& xc~=membrane)
           xc = xc + 1;
        elseif rand<permeability
            xc = xc +1;
        else
           xc = xc - 1;
        end
     else
        if rand < .5
           yc = yc + 1;
        else
           yc = yc - 1;
        end
     end
%    Save location...
     k = k + 1; x(k) = xc; y(k) = yc; 
  end