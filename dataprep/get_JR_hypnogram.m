function [hyp, sleep_start] = get_JR_hypnogram(filename)

    s = tdfread(filename,'\t');
    headers = fieldnames(s);
    numHeader = length(headers);
    j = 1;

for k=2:2:numHeader    
    
    if findstr(headers{k-1},'Epoch')
       Epochs = getfield(s,headers{k-1});
    end
    Epochs = (Epochs-1).*30;
    
    if findstr(headers{k},'Stage')
       Stages = getfield(s,headers{k});
    end
    
    if length(Epochs) == length(Stages)
        for i=1:length(Epochs)
            switch strtrim(Stages(i,:))
                case 'Unscored'
                    hyp(j,2) = Epochs(i);
                    hyp(j,1) = 11;
                    j=j+1;
                case 'UNS'
                    hyp(j,2) = Epochs(i);
                    hyp(j,1) = 11;   
                    j=j+1;
                case 'Stage 1'
                    hyp(j,2) = Epochs(i);
                    hyp(j,1) = 1;
                    j=j+1;
                case 'Stage 2'
                    hyp(j,2) = Epochs(i);
                    hyp(j,1) = 2;
                    j=j+1;
                case 'Stage 3'
                    hyp(j,2) = Epochs(i);
                    hyp(j,1) = 3;     
                    j=j+1;
                case 'Stage 4'
                    hyp(j,2) = Epochs(i);
                    hyp(j,1) = 4;    
                    j=j+1;
                case 'REM'
                    hyp(j,2) = Epochs(i);
                    hyp(j,1) = 5;   
                    j=j+1;
                case 'Wake'
                    hyp(j,2) = Epochs(i);
                    hyp(j,1) = 0;    
                    j=j+1;
                case 'MT'
                    hyp(j,2) = Epochs(i);
                    hyp(j,1) = 7;    
                    j=j+1;
                otherwise
                    
            end
            
        end
    end

end

% sleep_start = round(h(1)+ 0.0169*m(1)); % get roughly the time the person went to bed, multiply by 0.0169 to get from 60 min -> 100 for rounding :)
sleep_start = 0;
end % end of get_hypnogram_function