function [inbox] = inbox(LB,RT,data)

if LB(1)>RT(1)
	LBt=LB;
	RTt=RT;
	LB(1)=RTt(1);
	RT(1)=LBt(1);
end
if LB(2)>RT(2)
	LBt=LB;
	RTt=RT;
	LB(2)=RTt(2);
	RT(2)=LBt(2);
end


sz = numel(data(1,:));
inbox = zeros(1,sz);

for s = 1:sz
	
	if data(1,s) > LB(1) && data(1,s) < RT(1) &&...
	   data(2,s) > LB(2) && data(2,s) < RT(2)
		
		inbox(s) = 1;
	end
   
end

inbox = inbox>0;

end