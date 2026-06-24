% this script calculates and saves the sum of fractional contributions to the observations in ObsList

clear iterlist;
iterlist = ObsList; 
if sum(ismember(iterlist,'SO4')) >0; else; iterlist{end+1} = 'SO4';  end
if sum(ismember(iterlist,'norm'))>0; else; iterlist{end+1} = 'norm'; end

% a value for double-checking mass balace of inversion results
stillmb = 1; 

% calculate the sum of parameters
for i=1:length(iterlist)

    % get the active results
    activeresults = eval(sprintf('results.%s',iterlist{i}));

    % define some matricies to hold relevant parameters
    workingmedian = [];
    workingmean   = [];
    working05     = [];
    working25     = [];
    working75     = [];
    working95     = [];
    datamatrix    = [];

    % sum the active results for this sample
    for j=1:s        
        workingmedian(j,1) = nanmedian(sum(activeresults(:,:,j),1));
        workingmean(j,1)   = nanmean(sum(activeresults(:,:,j),1));
        working05(j,1)     = prctile(sum(activeresults(:,:,j),1),05);
        working25(j,1)     = prctile(sum(activeresults(:,:,j),1),25);
        working75(j,1)     = prctile(sum(activeresults(:,:,j),1),75);
        working95(j,1)     = prctile(sum(activeresults(:,:,j),1),95);
        datamatrix(j,:)    = sum(activeresults(:,:,j),1);
    end

    % for isotopic information, change the names to reflect the
    % product of isotopic ratios and elemental abundances
    if     isequal(iterlist{i},'d7Li');                                      savename = [iterlist{i} '_x_Li'];   % d7Li and Li
    elseif isequal(iterlist{i},'d18O');                                      savename = [iterlist{i} '_x_SO4'];  % d18O and SO4
    elseif isequal(iterlist{i},'d26Mg');                                     savename = [iterlist{i} '_x_Mg'];   % d26Mg and Mg
    elseif isequal(iterlist{i},'d30Si');                                     savename = [iterlist{i} '_x_Si'];   % d30Si and Si
    elseif isequal(iterlist{i},'d34S');                                      savename = [iterlist{i} '_x_SO4'];  % d34S and SO4
    elseif isequal(iterlist{i},'d42Ca');                                     savename = [iterlist{i} '_x_Ca'];   % d42Ca and Ca
    elseif isequal(iterlist{i},'d44Ca');                                     savename = [iterlist{i} '_x_Ca'];   % d44Ca and Ca
    elseif isequal(iterlist{i},'d56Fe');                                     savename = [iterlist{i} '_x_Fe'];   % d56Fe and Fe
    elseif isequal(iterlist{i},'Sr8786');                                    savename = [iterlist{i} '_x_Sr'];   % Sr8786 and Sr
    elseif isequal(iterlist{i},'Os8788');                                    savename = [iterlist{i} '_x_Os'];   % Os8788 and Os
    elseif isequal(iterlist{i},'d98Mo');                                     savename = [iterlist{i} '_x_Mo'];   % d98Mo and Mo
    elseif isequal(iterlist{i},'d13C') & isequal(carbonisotopematch,'HCO3'); savename = [iterlist{i} '_x_HCO3']; % d13C and HCO3
    elseif isequal(iterlist{i},'Fmod') & isequal(carbonisotopematch,'HCO3'); savename = [iterlist{i} '_x_HCO3']; % Fmod and HCO3
    elseif isequal(iterlist{i},'d13C') & isequal(carbonisotopematch,'DIC');  savename = [iterlist{i} '_x_DIC'];  % d13C and DIC
    elseif isequal(iterlist{i},'Fmod') & isequal(carbonisotopematch,'DIC');  savename = [iterlist{i} '_x_DIC'];  % Fmod and DIC
    else                                 savename = iterlist{i};
    end

    % save results for mass balance
    evalin('base',[sprintf('river.massbalance.median.%s = ', savename)  'workingmedian;']);
    evalin('base',[sprintf('river.massbalance.mean.%s   = ', savename)  'workingmean;']);
    evalin('base',[sprintf('river.massbalance.pct05.%s  = ', savename)  'working05;']);
    evalin('base',[sprintf('river.massbalance.pct25.%s  = ', savename)  'working25;']);
    evalin('base',[sprintf('river.massbalance.pct75.%s  = ', savename)  'working75;']);
    evalin('base',[sprintf('river.massbalance.pct95.%s  = ', savename)  'working95;']);
    evalin('base',[sprintf('river.massbalance.all.%s  = ',   savename)  'datamatrix;']);
end

% tell the user that the mass balance calculations are saved
if stillmb == 1 & isequal(IterateOver,'Samples')
   disp('saved mass balance calculations into the field "massbalance"');
end