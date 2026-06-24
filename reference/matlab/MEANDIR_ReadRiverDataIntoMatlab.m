% this script reads river observations into MATLAB
% the script is mostly just a large "if" statement that tests whether "Riverdatasource" is equal to a particular keyword
% if so, the data corresponding to that keyword is read into MATLAB using the function MEANDIR_StandardInputs

% define the conversion factors for moving between units
% of concentration and charge equivalents 
conc2equi = MEANDIR_DefineConc2Equi; 

% define the fields of a data structure
river = struct;
river = MEANDIR_DefineBasicFields;

% read in the default river data 
% (upon download, this is the Gaillardet and others (1999) data)
if isequal(Riverdatasource,'MyDefaultRiverData')
       [~, ~, all] = xlsread('MyDefaultRiverDataSpreadsheet','MyRiverData');
       MEANDIR_StandardInput;

% Gislason and others (1996)
elseif isequal(Riverdatasource,'SG96RiverData')
       [~, ~, all] = xlsread('RiverDataSpreadsheet_Gislason_etal_1996','MyRiverData');
       MEANDIR_StandardInput

% Gaillardet and others (1999)
elseif isequal(Riverdatasource,'JG99RiverData')
       [~, ~, all] = xlsread('RiverDataSpreadsheet_Gaillardet_etal_1999','MyRiverData');
       MEANDIR_StandardInput

% Torres and others (2016)
elseif isequal(Riverdatasource,'MT16RiverData')
       [~, ~, all] = xlsread('RiverDataSpreadsheet_Torres_etal_2016','MyRiverData');
       MEANDIR_StandardInput

% Burke and others (2018)
elseif isequal(Riverdatasource,'AB18RiverData')
       [~, ~, all] = xlsread('RiverDataSpreadsheet_Burke_etal_2018','MyRiverData');
       MEANDIR_StandardInput;
       river.info.dis = cell2mat(river.info.extrafield1);

% Horan and others (2019)
elseif isequal(Riverdatasource,'KH19Riv')
       [~, ~, all] = xlsread('RiverDataSpreadsheet_Horan_etal_2019','MyRiverData');
       MEANDIR_StandardInput

% Kemeny Alaska observations
elseif isequal(Riverdatasource,'PCKAlaskaData')
       [~, ~, all] = xlsread('RiverDataSpreadsheet_Kemeny_etal_2023','MyRiverData');
       MEANDIR_StandardInput       

else
       disp(sprintf('catastrophic failure! no river data was found for scenario "%s". Check "MEANDIR_ReadRiverDataIntoMatlab" to find your data.',Riverdatasource));
end

%% initial data processing
% construct and save the normalization variable
fn                 = fieldnames(river.observations);                            % get the fieldnames of the river observations
norm_equi          = zeros(size(eval(sprintf('river.observations.%s',fn{1})))); % define a vector of 0 values
norm_conc          = zeros(size(eval(sprintf('river.observations.%s',fn{1})))); % define a vector of 0 values
ObsInNormalization = cat(2,CationsInNormalization,AnionsInNormalization,NeutralInNormalization)';
ObsInNormalization = ObsInNormalization(~ismember(ObsInNormalization,'NaN'));
for i=1:length(ObsInNormalization)    
    norm_equi = norm_equi + eval(sprintf('river.observations.%s_equi',ObsInNormalization{i}));
    norm_conc = norm_conc + eval(sprintf('river.observations.%s_conc',ObsInNormalization{i}));
end
river.observations.norm_conc = norm_conc;
river.observations.norm_equi = norm_equi;
if     isequal(EMUnits,'conc'); river.model_variable.norm = river.observations.norm_conc;
elseif isequal(EMUnits,'equi'); river.model_variable.norm = river.observations.norm_equi;
end

% move the variables for the inversion from "observations" into the field "model_variable"
for i=1:length(ObsList)
    evalin('base',[sprintf('river.model_variable.%s =',    ObsList{i}) sprintf('river.observations.%s_%s',    ObsList{i},EMUnits) ';']);
    evalin('base',[sprintf('river.model_variable.%s_asd =',ObsList{i}) sprintf('river.observations.%s_%s_asd',ObsList{i},EMUnits) ';']);
end

% even if SO4 is not in the inversion, put SO4 observations (even if all NaN values) into the subfield "model_variables"
if sum(ismember(ObsList,'SO4'))==0
    river.model_variable.SO4     = eval(sprintf('river.observations.SO4_%s',EMUnits));      % now assign that data to the variable that the model will use
    river.model_variable.SO4_asd = eval(sprintf('river.observations.SO4_%s_asd',EMUnits));  % find the relevant data
end

% ask how many samples have all of the relevant data
smprep = []; 
for i=1:length(ObsList)
    activedat   = eval(sprintf('river.model_variable.%s',ObsList{i}));
    smprep(:,i) = ~isnan(activedat);
end
allzero = sum(sum(smprep,2)==0);
allwork = sum(sum(smprep,2)==length(ObsList));
numbersamples        = length(sum(smprep,2));
functionalsamplelist = (sum(smprep,2)==length(ObsList)); 
disp(sprintf('imported river data, found %i samples in total and %i samples with all the data required for the requested inversion',numbersamples,allwork));
river.info.numbersampleswithdata = allwork;

% check is ClCritical is possible
if isequal(PrecProcessing,'ClCrit') & ClCriticalValuesGiven==1
    if sum(~isnan(river.observations.ClCr))==0
       disp(sprintf('catastrophic failure! all ClCrit values are NaN, simulation "%s" will only return zeros.',ScenarioList{ScenarioListIndex}));
       river = NaN;
    end
end

% when AdjustRiverObs is set to 1, make sure data exists for analytical uncertainties
if  AdjustRiverObs==1
    EM = MEANDIR_ReadEndMembersIntoMatlab;
    for i=1:length(ObsList)  
        erroronthision = eval(sprintf('river.model_variable.%s_asd',ObsList{i}));
        if  sum(~isnan(erroronthision))==0
            FailCase = 0; 
            disp(sprintf('warning on simulation "%s"! simulation is set to include error on sample measurements ("AdjustRiverObs" is equal to 1), but no error values are found for %s. \nEither include error, remove %s, or set "AdjustRiverObs" to equal zero.',ScenarioList{ScenarioListIndex},ObsList{i},ObsList{i}));
        end  
    end
end

% remove fields that are entirely NaN. keep the fields that are not 100%
% NaN, even if those variables are not included in the inversion
% (this way the data are carried forward for off-line analysis)
fn = fieldnames(river.observations);
for ii=1:length(fn)
    adat = eval(sprintf('river.observations.%s',fn{ii}));
    if     isequal(class(adat),'cell')
           if sum(~ismember(adat,'NaN'))==0
                river.observations = rmfield(river.observations,fn{ii});
           end
    elseif isequal(class(adat),'double')
           if sum(~isnan(adat))==0
                river.observations = rmfield(river.observations,fn{ii});
           end
    end
end

% remove variables to prevent workspace clutter
clear raw*;
clear Unc_raw*;
clear up*;
clear Unc_up*;
clear num;
clear txt;
clear all;
clear tempstruct;
clear smprep;
clear sumions_equi;
clear header;
clear fn;
clear charge;
clear activedat;
clear activeion;