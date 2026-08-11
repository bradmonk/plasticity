 function [x y] = Diffusion(N)
% N is a positive integer.
% Simulates Diffusion
% x and y are row vectors with the property that (x(k),y(k)) is the
% location of the token after n time

% Initializations...
  k = 0; xc = 0; yc = 0;
  
% In general, (xc,yc) is the location after N Times.

  while k<N 
%    Standing at (xc,yc). 
%    Random determines new location...
     if rand < .5
        if rand < .5
           xc = xc + 1;
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