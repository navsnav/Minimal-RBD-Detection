function fitresult = fittingfunction(p1,p2,p3,p4,rf,posterior,Ytrn)

    cost = [0,1,p1;1,0,p2;p3,p4,0];

    expCost = posterior*cost;
    [~,classIndex] = min(expCost,[],2);
    Yhat2 = rf.ClassNames(classIndex);
    Yhat2 = str2num(cell2mat(Yhat2));
    
    metrics  = process_classification_results2(Yhat2==5, Ytrn==5);
    
    ConfMat_Class_Summary = confusionmat(Yhat2, Ytrn, 'order', [0 2 5]);
    kappa = kappa_result(ConfMat_Class_Summary);    
%     fitresult = metrics(6); %F1
%     fitresult = metrics(7); %Kappa
    fitresult = metrics(4); %Precision
%     fitresult = kappa; %Precision
    
    

end