function [P] = annulus(r,R,xf,Xf,yf,Yf)
% Creates a annulus patch object and returns the handle.  Input arguments 
% are the inner radius, outer radius, inner x offset, outer x offset, inner
% y offset and outer Y offset.  Changes to the edgecolor and linestyle are
% allowed, and will preserve the correct look of the annulus
t = linspace(0,2*pi,200);
x = xf + r*cos(t);
y = yf + r*sin(t);
X = Xf + R*cos(t);
Y = Yf + R*sin(t);
P = patch([x X],[y Y],[1,.5,.5],'linestyle','non','facealph',.5);
L(1) = line(x,y,'color','k');
L(2) = line(X,Y,'color','k');
axis equal
plistener(P,'edgecolor',@edgec) % listeners for changes to props.
plistener(P,'linestyle',@lnstl)
      function [] = plistener(ax,prp,func)
          % Sets the properties. From proplistener by Yair Altman.
          psetact = 'PropertyPostSet';
          hC = handle(ax);
          hSrc = hC.findprop(prp);
          hl = handle.listener(hC, hSrc, psetact, {func,ax});
          p = findprop(hC, 'Listeners__');
          if isempty(p)
              p = schema.prop(hC, 'Listeners__', 'handle vector');
              set(p,'AccessFlags.Serialize', 'off', ...
                  'AccessFlags.Copy', 'off',...
                  'FactoryValue', [], 'Visible', 'off');
          end
          hC.Listeners__ = hC.Listeners__(ishandle(hC.Listeners__));
          hC.Listeners__ = [hC.Listeners__; hl];
      end
      function [] = edgec(varargin)
          St = get(varargin{3},'edgecolor');
          set(L,'color',St)
      end
      function [] = lnstl(varargin)
          St = get(varargin{3},'linestyle');
          set(varargin{3},'linestyle','none')
          set(L,'linestyle',St)
      end
end