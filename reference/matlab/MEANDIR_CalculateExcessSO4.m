% when able, this script calculates excess SO4 from inversion results.
calculateExcessSO4 = 1; % this is a parameter that can fail if conditions are not met

% calculate excess SO4
if  sum(ismember(ObsList,'SO4')) > 0;  calculateExcessSO4 = 0; disp('did not calculate excess SO4 because SO4 was included in the inversion');                                end
if  isequal(NormalizationType,'SO4');  calculateExcessSO4 = 0; disp('did not calculate excess SO4 because SO4 was included in the inversion as the normalization variable');  end
if  EMHaveSO4 == 0;                    calculateExcessSO4 = 0; disp('did not calculate excess SO4 because the end-members lack SO4 values');                                  end

if  calculateExcessSO4==1
    ExcessSO4  = []; 
    fExcessSO4 = [];
    rivrSO4    = river.model_variable.SO4;

    for i=1:s
        % these should be sum and not nansum
        % if they are nansum, then four 'NaN' will be summed to a 0
        ExcessSO4entries = rivrSO4(i).*(1-sum(results.SO4(:,:,i),1)); 
        fExcessentries   = (1-sum(results.SO4(:,:,i),1));
        ExcessSO4(i,:)   = ExcessSO4entries; 
        fExcessSO4(i,:)  = fExcessentries;  
    end

    % absolute excess SO4
    ExcessSO4_median(1:s,1)         = nanmedian(ExcessSO4',1)';
    ExcessSO4_pct05(1:s,1)          = prctile(ExcessSO4',05)';
    ExcessSO4_pct25(1:s,1)          = prctile(ExcessSO4',25)';
    ExcessSO4_pct75(1:s,1)          = prctile(ExcessSO4',75)';
    ExcessSO4_pct95(1:s,1)          = prctile(ExcessSO4',95)';

    % fractional excess SO4
    fExcessSO4_median(1:s,1)        = nanmedian(fExcessSO4',1)';
    fExcessSO4_pct05(1:s,1)         = prctile(fExcessSO4',05)';
    fExcessSO4_pct25(1:s,1)         = prctile(fExcessSO4',25)';
    fExcessSO4_pct75(1:s,1)         = prctile(fExcessSO4',75)';
    fExcessSO4_pct95(1:s,1)         = prctile(fExcessSO4',95)';

    % save excess SO4
    river.ExcessSO4.absolute.all    = ExcessSO4;
    river.ExcessSO4.absolute.median = ExcessSO4_median;
    river.ExcessSO4.absolute.pct05  = ExcessSO4_pct05;
    river.ExcessSO4.absolute.pct25  = ExcessSO4_pct25;
    river.ExcessSO4.absolute.pct75  = ExcessSO4_pct75;
    river.ExcessSO4.absolute.pct95  = ExcessSO4_pct95;
    
    % save fractional excess SO4
    river.ExcessSO4.fraction.all    = fExcessSO4;
    river.ExcessSO4.fraction.median = fExcessSO4_median;
    river.ExcessSO4.fraction.pct05  = fExcessSO4_pct05;
    river.ExcessSO4.fraction.pct25  = fExcessSO4_pct25;
    river.ExcessSO4.fraction.pct75  = fExcessSO4_pct75;
    river.ExcessSO4.fraction.pct95  = fExcessSO4_pct95;

    % remove clutter from the workspace
    clear ExcessSO4_*;
    clear fExcessSO4_*;

    % tell the use that excess SO4 was calculated
    disp('calculated excess SO4 and saved into field "ExcessSO4"');
    
end % end of primary if statement