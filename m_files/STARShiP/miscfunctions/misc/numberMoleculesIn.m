 function [chargeIn chargeOut] =numberMoleculesIn(N,membrane,xRed,permRed,nRedMolecules)
% N is a positive integer.
% Indicates directional concentration based change
% x is the new xRed, y the new xBlue and z the new xGreen
% How many positive and negative are in cell
chargeIn=zeros(1,N);
chargeOut=zeros(1,N);

% Number of positive in and out of cell
for cRed=1:length(xRed)
    for rRed=1:nRedMolecules
        if xRed(rRed,cRed)>membrane
            chargeOut(1,cRed)=chargeOut(1,cRed) + 1;
        else
            chargeIn(1,cRed)=chargeIn(1,cRed) + 1;
        end
    end
end
