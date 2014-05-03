function [varargout] = CalciumFunction(varargin)


CaE = 4;		% Ca Elimination Rate per 5 steps
	CaDots = 50;	% N Ca dots; must be less than CaRT
	CaRT = 100;		% Ca Release Time, create Ca dots every CaRT steps



%===========================================%
% CALCIUM MATRIX(2,CaNdots)
%-------------------------------------------%
CaD = .01;                    
Cad = 2;                      
CadT = 1;                     
Cak = sqrt(Cad*CaD*CadT); 

CaNdots = 50;
Caxyl = ones(2,CaNdots);

for j = 1:CaNdots  
	Caxyl(1,:)=(col1F+3);
	Caxyl(2,:)=(row2F-3);
end 
%===========================================%
%{
	%-------------------------------%
    %     DO MANUAL STEP SIZES
    %-------------------------------%
	if MANstepsize
		GluR2xyds = stepsize(Ndots, Lx);
	end % xyd = DIRxyd(xyd);
	%}
	
	%{
	%-------------------------------%
    %   CALCIUM FUNCTION
    %-------------------------------%
	if doCalcium
    [Caxyl] = CalciumFun(Cak, CaNdots, Caxyl);
	
	% Every RT steps, make CaDots appear at Caxyl
	if mod(stepN, CaRT) == 0
		CaNdots = CaDots;
		Caxyl = ones(2,CaNdots);
		for j = 1:CaNdots  
			Caxyl(1,:)=(col1F+PSD2CNTR);
			Caxyl(2,:)=(row2F-PSD2CNTR);
		end 
	end
	
	% Delete Ca dots
	if mod(stepN, 5) == 0
		if CaNdots>0
		CaNdots = CaNdots-CaE;
		Caxyl = Caxyl(:,1:end-CaE);
		end
	end
	
	PSD2CaT = size(Caxyl,2); % Get total Ca dots
	PSD1CaT = 0;			 % Get total Ca dots
	
	
	%In Case S-Clustering is Shut Off
	if doUse(3) == 0		
	S1Ds = (S1sumOrig+PSD1CaT)/100;
	S2Ds = (S2sumOrig+PSD2CaT)/100;
	Dr_PSD1 = D/(S1Ds);
	Dr_PSD2 = D/(S2Ds);
	PSD1 = 1/sqrt(Dr_PSD1); 
	PSD2 = 1/sqrt(Dr_PSD2);
	end
	
	end % end doCalcium
	%}



%{
	%-------------------------------%
    %     FRAP Functoin
    %-------------------------------%
	if doFRAP
	if mod(stepN, 1002) == 0
		[Ndots GluR2xyl GluR1xyl GluR1Ndots] = FRAPfun(Ndots, GluR2xyl,... 
		GluR1xyl, GluR1Ndots, row1F, col1F, row1L, col1L, row2L, row2F);
	end
	end
	%}




end