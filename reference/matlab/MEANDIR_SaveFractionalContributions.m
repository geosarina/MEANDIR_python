% this script saves the fractional contributions of each end-member to each dissolved variable
% the file is called near the end of MEANDIR, when simulation results are being saved into fields of the river structure

%% (1) iterate over each ion ratio and each end-member
for k=1:nOL
    
    % load the results for this variable
    activeresults = eval(sprintf('results.%s',ObsList{k})); % pull all data for this ion

    % iterate over the number of end-members
    for j = 1:nEM 
         
        % refresh the temporary storage vectors
        medians      = [];
        zerofives    = [];
        twentyfives  = [];
        seventyfives = [];
        ninetyfives  = [];
        
        % pull the data for the active result, reshape,
        % and calculate parameters of interest
        allact = activeresults(j,:,:); 
        allact = reshape(allact,size(allact,2),size(allact,3));
        medians(1:s,1)      = nanmedian(allact,1);
        zerofives(1:s,1)    = prctile(allact,05)';
        twentyfives(1:s,1)  = prctile(allact,25)';
        seventyfives(1:s,1) = prctile(allact,75)';
        ninetyfives(1:s,1)  = prctile(allact,95)';
            
        % if the current variable is an isotope, update the name of the
        % model field to reflect that the f values are the contribution
        % to the product of abundance and isotope information
        if     sum(indxisopos==k)==0;                                           sameobsname =  ObsList{k};
        elseif isequal(ObsList{k},'d7Li');                                      sameobsname = [ObsList{k} '_x_Li'];   % d7Li and Li
        elseif isequal(ObsList{k},'d18O');                                      sameobsname = [ObsList{k} '_x_SO4'];  % d18O and SO4
        elseif isequal(ObsList{k},'d26Mg');                                     sameobsname = [ObsList{k} '_x_Mg'];   % d26Mg and Mg
        elseif isequal(ObsList{k},'d30Si');                                     sameobsname = [ObsList{k} '_x_Si'];   % d30Si and Si
        elseif isequal(ObsList{k},'d34S');                                      sameobsname = [ObsList{k} '_x_SO4'];  % d34S and SO4
        elseif isequal(ObsList{k},'d42Ca');                                     sameobsname = [ObsList{k} '_x_Ca'];   % d42Ca and Ca
        elseif isequal(ObsList{k},'d44Ca');                                     sameobsname = [ObsList{k} '_x_Ca'];   % d44Ca and Ca
        elseif isequal(ObsList{k},'d56Fe');                                     sameobsname = [ObsList{k} '_x_Fe'];   % d56Fe and Fe
        elseif isequal(ObsList{k},'Sr8786');                                    sameobsname = [ObsList{k} '_x_Sr'];   % Sr8786 and Sr
        elseif isequal(ObsList{k},'Os8788');                                    sameobsname = [ObsList{k} '_x_Os'];   % Os8788 and Os
        elseif isequal(ObsList{k},'d98Mo');                                     sameobsname = [ObsList{k} '_x_Mo'];   % d98Mo and Mo
        elseif isequal(ObsList{k},'d13C') & isequal(carbonisotopematch,'HCO3'); sameobsname = [ObsList{k} '_x_HCO3']; % d13C and HCO3
        elseif isequal(ObsList{k},'Fmod') & isequal(carbonisotopematch,'HCO3'); sameobsname = [ObsList{k} '_x_HCO3']; % Fmod and HCO3
        elseif isequal(ObsList{k},'d13C') & isequal(carbonisotopematch,'DIC');  sameobsname = [ObsList{k} '_x_DIC'];  % d13C and DIC
        elseif isequal(ObsList{k},'Fmod') & isequal(carbonisotopematch,'DIC');  sameobsname = [ObsList{k} '_x_DIC'];  % Fmod and DIC
        end

        % save the fractional contributions
        evalin('base',[sprintf('river.fraction.%s.%s.median =',sameobsname,EMList0{j}) 'medians;']);      % save the median fractional contributions
        evalin('base',[sprintf('river.fraction.%s.%s.pct05  =',sameobsname,EMList0{j}) 'zerofives;']);    % save the 05th percentile of fractional contributions
        evalin('base',[sprintf('river.fraction.%s.%s.pct25  =',sameobsname,EMList0{j}) 'twentyfives;']);  % save the 25th percentile of fractional contributions
        evalin('base',[sprintf('river.fraction.%s.%s.pct75  =',sameobsname,EMList0{j}) 'seventyfives;']); % save the 75th percentile of fractional contributions
        evalin('base',[sprintf('river.fraction.%s.%s.pct95  =',sameobsname,EMList0{j}) 'ninetyfives;']);  % save the 95th percentile of fractional contributions
    end
    evalin('base',[sprintf('river.fraction.all.%s  =',sameobsname) 'activeresults;']);
end

% save the misfit data as well
river.misfit.modelobservations = results.misfit;
river.misfit.costfunction      = results.misfit_costfunction;

% consider saving into "fraction" all of the results, including unculled simulations
if  saveuncutdata == 1 & isequal(IterateOver,'End-members') 
    river.misfit.modelobservations_uncut = results.all_uncut.misfit;
    river.misfit.costfunction_uncut      = results.all_uncut.misfit_costfunction;     
    resnames = fieldnames(results.all_uncut);
    resnames = resnames(~ismember(resnames,'misfit'));
    resnames = resnames(~ismember(resnames,'misfit_costfunction'));
    for jj=1:length(resnames)
        adata = eval(sprintf('results.all_uncut.%s',resnames{jj}));
        evalin('base',[sprintf('river.fraction.all_uncut.%s  =',resnames{jj}) sprintf('results.all_uncut.%s',resnames{jj}) ';']);
    end
end

%% (2) SO4
if  sum(ismember(ObsList,'SO4'))==0
    
    % iterate over the number of end-members
    for j = 1:nEM   
        
        % refresh the temporary storage vectors
        medians      = [];
        zerofives    = [];
        twentyfives  = [];
        seventyfives = [];
        ninetyfives  = [];
        
        % pull the data for SO4, reshape,
        % and calculate parameters of interest
        allact = results.SO4(j,:,:); 
        allact = reshape(allact,size(allact,2),size(allact,3));
        medians(1:s,1)      = nanmedian(allact,1);
        zerofives(1:s,1)    = prctile(allact,05)';
        twentyfives(1:s,1)  = prctile(allact,25)';
        seventyfives(1:s,1) = prctile(allact,75)';
        ninetyfives(1:s,1)  = prctile(allact,95)';
                  
        % save the fractional contributions
        evalin('base',[sprintf('river.fraction.SO4.%s.median =',EMList0{j}) 'medians;']);
        evalin('base',[sprintf('river.fraction.SO4.%s.pct05  =',EMList0{j}) 'zerofives;']);
        evalin('base',[sprintf('river.fraction.SO4.%s.pct25  =',EMList0{j}) 'twentyfives;']);
        evalin('base',[sprintf('river.fraction.SO4.%s.pct75  =',EMList0{j}) 'seventyfives;']);
        evalin('base',[sprintf('river.fraction.SO4.%s.pct95  =',EMList0{j}) 'ninetyfives;']);
    end    
    river.fraction.all.SO4  = results.SO4;
end

%% (3) normalization
% iterate over the number of end-members
for j = 1:nEM
    
    % refresh the temporary storage vectors
    medians      = [];
    zerofives    = [];
    twentyfives  = [];
    seventyfives = [];
    ninetyfives  = [];
    
    % pull the data for the normalization, reshape,
    % and calculate parameters of interest    
    allact = results.norm(j,:,:); 
    allact = reshape(allact,size(allact,2),size(allact,3));
    medians(1:s,1)      = nanmedian(allact,1);
    zerofives(1:s,1)    = prctile(allact,05)';
    twentyfives(1:s,1)  = prctile(allact,25)';
    seventyfives(1:s,1) = prctile(allact,75)';
    ninetyfives(1:s,1)  = prctile(allact,95)';
           
    % save the fractional contributions
    evalin('base',[sprintf('river.fraction.norm.%s.median =',EMList0{j}) 'medians;']);
    evalin('base',[sprintf('river.fraction.norm.%s.pct05  =',EMList0{j}) 'zerofives;']);
    evalin('base',[sprintf('river.fraction.norm.%s.pct25  =',EMList0{j}) 'twentyfives;']);
    evalin('base',[sprintf('river.fraction.norm.%s.pct75  =',EMList0{j}) 'seventyfives;']);
    evalin('base',[sprintf('river.fraction.norm.%s.pct95  =',EMList0{j}) 'ninetyfives;']);
end   
river.fraction.all.norm = results.norm;

%% (4) tell the user that the fractional contributions have been saved
disp('saved the fractional contributions of each end-member into the field "fraction"');