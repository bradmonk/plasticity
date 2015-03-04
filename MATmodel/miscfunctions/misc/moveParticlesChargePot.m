function [x chargeIn chargeOut]=moveParticlesChargePot(N,membrane,xRed,permRed,nRedMolecules,positiveIn, positiveOut, negativeIn, negativeOut, VK)
%Based upon charge in previous step move away from similiar charge &&Vion
%Change In Charge Concentration
%Each quarter closer to |58| add .1 per quarter to rand to impede movement
%If charge diff add .2
chargeIn=zeros(1,N);
chargeOut=zeros(1,N);
impedanceIon=abs((VK/14.5))*0.1;
for rRed=1:nRedMolecules
    for cRed=2:nRedMolecules
        if ((positiveIn(1,cRed-1)-negativeIn(1,cRed-1))>(positiveOut(1,cRed-1)-negativeOut(1,cRed-1)))
            if (((rand+.20)>0.5)&&(xRed(rRed,cRed)<0))
                if (rand<permRed)&&(rand>impedanceIon)
                    xRed(rRed,cRed)=xRed(rRed,cRed)+1;
                end
                if xRed(rRed,cRed)>membrane
                    chargeIn(1,cRed)=chargeIn(1,cRed) - 1;
                    chargeOut(1,cRed)=chargeOut(1,cRed) + 1;
                end
            end
        elseif(((rand+.20)>0.5)&&(xRed(rRed,cRed)>0))
            if (rand<permRed)&&(rand>impedanceIon)
                xRed(rRed,cRed)=xRed(rRed,cRed)-1;
            end
            if xRed(rRed,cRed)<=membrane
                chargeIn(1,cRed)=chargeIn(1,cRed) + 1;
                chargeOut(1,cRed)=chargeOut(1,cRed) - 1;
            end
        end
    end
end
x=xRed;