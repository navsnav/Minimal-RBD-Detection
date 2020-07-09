function     cost_opt = find_optimise_cost(rf,posterior,Ytrn)

    
firstparam = linspace(1,0,6);     %list of places to search for first parameter
secondparam = linspace(1,0,6);         %list of places to search for second parameter
thirdparam = linspace(1,2,6);     %list of places to search for third parameter
fourthparam = linspace(1,2,6);         %list of places to search fourth second parameter


[F,S,Q,R] = ndgrid(firstparam, secondparam,thirdparam,fourthparam);
fitresult = arrayfun(@(p1,p2,p3,p4) fittingfunction(p1,p2,p3,p4,rf,posterior,Ytrn), F,S,Q,R); %run a fitting on every pair fittingfunction(F(J,K), S(J,K))



[maxval, maxidx] = max(fitresult(:));
bestFirst = F(maxidx);
bestSecond = S(maxidx);
bestThird = Q(maxidx);
bestFourth = R(maxidx);


cost_opt = [0,1,bestFirst;1,0,bestSecond;bestThird,bestFourth,0];
% cost_opt = [0,1,0;1,0,0;1,1,0];

end


