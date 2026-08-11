function Bounce (Anzahl, Dauer, xMax, yMax)
% draws the path of one diffusing particle.

% default setting
if (nargin < 4)
    xMax = 40;
    yMax = 30;
end
if (nargin < 2)
    Anzahl = 20;
    Dauer = 1000;
end

Position = zeros(2, Anzahl);
Richtung = zeros(1, Anzahl);
% the variable for the path
Kurs = zeros(2, Dauer);
% initialising the particles (random starting points)
for i = 1:Anzahl
    Richtung(1, i) = rand(1)*2*pi;
    Position(1, i) = rand(1)*xMax;
    Position(2, i) = rand(1)*yMax;
end;

clf;
axis square; axis equal;
for i = 1:Dauer
    hold off;
    plot([0, 0, xMax, xMax, 0], [0, yMax, yMax, 0, 0], 'b-');
    hold on;
    plot([-10, -10, xMax+10, xMax+10, -10], [-10, yMax+10, yMax+10, -10, -10], 'b-');
    axis tight;
    % remembering the path
    Kurs(1, i) = Position(1, 1);
    Kurs(2, i) = Position(2, 1);
    for j = 1:Anzahl
       if (j == 1) 
           % plot the special particle
           plot(Position(1, j), Position(2, j), 'r.');
       else
           % plot all other particles
           plot(Position(1, j), Position(2, j), 'b.');
       end;
       for t = 1:Anzahl
           % is there another particle that hits the current one?
           if (t ~= j) & (sqrt((Position(1, j)-Position(1, t))^2 + (Position(2, j)-Position(2, t))^2))<=2
               if (Position(1, t) ~= Position(1, j))
                  Lot = atan((Position(2, t)-Position(2, j))/(Position(1, t)-Position(1, j)));
                else
                  Lot = pi/2;
               end;
               alpha1 = Richtung(1, j);
               alpha2 = Richtung(1, t);
               Richtung(1, j) = mod(2*Lot - alpha2, 2*pi);
               Richtung(1, t) = mod(2*Lot - alpha1, 2*pi);
           end;
       end;
       % calculating the new positions
       Position(1, j) = Position(1, j) + 2*cos(Richtung(1, j));
       Position(2, j) = Position(2, j) + 2*sin(Richtung(1, j));
       
       % is there a wall?
       if (Position(1, j) >= xMax)
           Richtung(1, j) = pi - Richtung(1, j);
           Position(1, j) = xMax;
       end;
       if (Position(1, j) <= 0)
           Richtung(1, j) = pi - Richtung(1, j);
           Position(1, j) = 0;
       end;
       if (Position(2, j) >= yMax)
           Richtung(1, j) = 2*pi - Richtung(1, j);
           Position(2, j) = yMax;
       end;
       if (Position(2, j) <= 0)
           Richtung(1, j) = 2*pi - Richtung(1, j);
           Position(2, j) = 0;
       end;
       Richtung(1, j) = mod(Richtung(1, j), 2*pi);
    end;
    pause (1/24);
end;
pause;
clf;
x = 1:Dauer;
hold on;
% plotting the path
plot([0, 0, xMax, xMax, 0], [0, yMax, yMax, 0, 0], 'b-');
hold on;
plot([-10, -10, xMax+10, xMax+10, -10], [-10, yMax+10, yMax+10, -10, -10], 'b-');
hold on;   
plot(Kurs(1, x), Kurs(2, x), 'r-');
axis tight;
hold off;