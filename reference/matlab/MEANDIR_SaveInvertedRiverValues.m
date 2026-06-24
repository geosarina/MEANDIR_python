% this script is used to save the entries of InvRivDataAgg, which contains information on the river chemistry inverted in each simulation
% this script loops over the fields of InvRivDataAgg and saves the data into the field InvertedRiverValues of Riverdatasource

fn = fieldnames(InvRivDataAgg); 
fn = fn(~ismember(fn,'all_uncut'));
fn = fn(~ismember(fn,'all_errorcut'));
for i=1:length(fn)    
   
    activedattest = eval(sprintf('InvRivDataAgg.%s',fn{i}));
    if sum(~isnan(activedattest(:)))>0
                          
        % pull the active data
        activedat = eval(sprintf('InvRivDataAgg.%s',fn{i}));
        
        % if there is isotopic data and it was converted to isotopic ratios values for inversion, convert it back to delta values for saving results
        if isequal(fn{i},'d7Li')  & sum(ismember(ConvertDelta2RList,'d7Li'))  >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d7Li'));  activedat = (activedat/conv-1)*1000;  end % d7Li
        if isequal(fn{i},'d18O')  & sum(ismember(ConvertDelta2RList,'d18O'))  >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d18O'));  activedat = (activedat/conv-1)*1000;  end % d18O
        if isequal(fn{i},'d26Mg') & sum(ismember(ConvertDelta2RList,'d26Mg')) >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d26Mg')); activedat = (activedat/conv-1)*1000;  end % d26Mg
        if isequal(fn{i},'d30Si') & sum(ismember(ConvertDelta2RList,'d30Si')) >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d30Si')); activedat = (activedat/conv-1)*1000;  end % d30Si
        if isequal(fn{i},'d34S')  & sum(ismember(ConvertDelta2RList,'d34S'))  >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d34S'));  activedat = (activedat/conv-1)*1000;  end % d34S
        if isequal(fn{i},'d42Ca') & sum(ismember(ConvertDelta2RList,'d42Ca')) >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d42Ca')); activedat = (activedat/conv-1)*1000;  end % d42Ca
        if isequal(fn{i},'d44Ca') & sum(ismember(ConvertDelta2RList,'d44Ca')) >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d44Ca')); activedat = (activedat/conv-1)*1000;  end % d44Ca
        if isequal(fn{i},'d56Fe') & sum(ismember(ConvertDelta2RList,'d56Fe')) >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d56Fe')); activedat = (activedat/conv-1)*1000;  end % d56Fe
        if isequal(fn{i},'d98Mo') & sum(ismember(ConvertDelta2RList,'d98Mo')) >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d98Mo')); activedat = (activedat/conv-1)*1000;  end % d98Mo
        if isequal(fn{i},'d13C')  & sum(ismember(ConvertDelta2RList,'d13C'))  >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d13C'));  activedat = (activedat/conv-1)*1000;  end % d13C 
        
        % calculate percentiles of activedat
        inverted_median = []; inverted_median(1:s,1) = nanmedian(activedat',1);
        inverted_pct05  = []; inverted_pct05(1:s,1)  = prctile(activedat',05)';
        inverted_pct25  = []; inverted_pct25(1:s,1)  = prctile(activedat',25)';
        inverted_pct75  = []; inverted_pct75(1:s,1)  = prctile(activedat',75)';
        inverted_pct95  = []; inverted_pct95(1:s,1)  = prctile(activedat',95)';
        
        % save the inverted values
        evalin('base',[sprintf('river.InvertedRiverValues.%s.all    =',fn{i}) 'activedat;'])
        evalin('base',[sprintf('river.InvertedRiverValues.%s.median =',fn{i}) 'inverted_median;'])
        evalin('base',[sprintf('river.InvertedRiverValues.%s.pct05  =',fn{i}) 'inverted_pct05;'])
        evalin('base',[sprintf('river.InvertedRiverValues.%s.pct25  =',fn{i}) 'inverted_pct25;'])
        evalin('base',[sprintf('river.InvertedRiverValues.%s.pct75  =',fn{i}) 'inverted_pct75;'])
        evalin('base',[sprintf('river.InvertedRiverValues.%s.pct95  =',fn{i}) 'inverted_pct95;'])

        % if specified, save the uncut data
        if saveuncutdata==1 & isequal(IterateOver,'End-members') 
            evalin('base',[sprintf('river.InvertedRiverValues.all_uncut.%s =',fn{i}) sprintf('InvRivDataAgg.all_uncut.%s',fn{i}) ';']);
        end
        
    else
        % if all the data are NaN, just save the value of NaN
        evalin('base',[sprintf('river.InvertedRiverValues.%s.all    =',fn{i}) 'NaN;'])
        evalin('base',[sprintf('river.InvertedRiverValues.%s.median =',fn{i}) 'NaN;'])
        evalin('base',[sprintf('river.InvertedRiverValues.%s.pct05  =',fn{i}) 'NaN;'])
        evalin('base',[sprintf('river.InvertedRiverValues.%s.pct25  =',fn{i}) 'NaN;'])
        evalin('base',[sprintf('river.InvertedRiverValues.%s.pct75  =',fn{i}) 'NaN;'])
        evalin('base',[sprintf('river.InvertedRiverValues.%s.pct95  =',fn{i}) 'NaN;'])
    end     
end

% display that the river data has been saved
disp('saved the river data used in the inversion into the field "InvertedRiverValues"');

% remove clutter from the workspace
clear inverted_median;
clear inverted_pct05;
clear inverted_pct25;
clear inverted_pct75;
clear inverted_pct95;
clear InvRivDataAgg;
