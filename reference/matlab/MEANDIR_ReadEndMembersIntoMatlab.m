function EM = MEANDIR_ReadEndMembersIntoMatlab

% this function reads information on end-member chemistry into MATLAB
% the workbook containing information on end-member chemistry is "MEANDIR_UserEntries" and the spreadsheet name is "MEANDIR_Endmembers"

% open the spreadsheet containing the weathering end-members and identify the header (full names) 
% and secondary names (whatever is written on the second row of the spreadsheet) for each end-member
[num txt all] = xlsread('MEANDIR_UserEntries','MEANDIR_Endmembers');
header     = all(1,:);     for i=1:length(header);      if isnan(header{i});      header{i}      = 'NaN'; end; end
secondname = all(2,1:end); for i=1:length(secondname);  if isnan(secondname{i});  secondname{i}  = 'NaN'; end; end
secondname = secondname(~ismember(secondname,'NaN'));       

% check for the required header entries
if sum(ismember(header,'End-member group (add a custom name)'))~=1
   disp('warning! problem in "MEANDIR_ReadEndMembersIntoMatlab". the value of the A1 cell in the end-members spreadsheet must be "End-member group (add a custom name)"');
   disp('with the current inputs, MEANDIR either cannot find that entry in the first row of the end-member spreadsheet or that value occurs more than once');
end
if sum(ismember(header,'NumeratorVariable'))~=1
   disp('warning! problem in "MEANDIR_ReadEndMembersIntoMatlab". the value of the C1 cell in the end-members spreadsheet must be "NumeratorVariable');
   disp('with the current inputs, MEANDIR either cannot find that entry in the first row of the end-member spreadsheet or that value occurs more than once');
end
if sum(ismember(header,'VariableNameForMatlab'))~=1
   disp('warning! problem in "MEANDIR_ReadEndMembersIntoMatlab". the value of the F1 cell in the end-members spreadsheet must be "VariableNameForMatlab"');
   disp('with the current inputs, MEANDIR either cannot find that entry in the first row of the end-member spreadsheet or that value occurs more than once');
end

% find the names of all the end-member groups
study  = all(3:end,ismember(header,'End-member group (add a custom name)'));
for i=1:length(study)
    if isnan(study{i})
        study{i} = 'NaN'; 
    end
end
allstudynames = unique(study);
namesofEMsets = allstudynames(~ismember(allstudynames,'NaN'))';

% declare a structure EM to hold information on end-member chemistry
clear EM;    % clear the EM structure if it exists
EM = struct; % define a structure

% iterate over the list of unique EM sets
for EMsetindx = 1:length(namesofEMsets)  
   
    % find and the names of all variables
    variablename = all(logical([0; 0; ismember(study,namesofEMsets{EMsetindx})]),ismember(header,'VariableNameForMatlab'));   
    varnumlist   = all(logical([0; 0; ismember(study,namesofEMsets{EMsetindx})]),ismember(header,'NumeratorVariable'));
          
    % find the active data
    activedata     = all(logical([0; 0; ismember(study,namesofEMsets{EMsetindx})]),7:length(secondname)+6);  
    addvaluematrix = NaN(size(activedata));
        
    [xx, yy] = size(activedata);
    % check the data for 'sample' entries
    for xxi=1:xx
    for yyi=1:yy
        if isequal(class(activedata{xxi,yyi}),'char')
            activechar    = activedata{xxi,yyi}; 
            checkhashsign = ismember(activechar,'#'); % Check if there is a # sign in the entry
            if     sum(checkhashsign)==1
                   pre_hashsignchar = activechar(1:find(checkhashsign)-1);   % Get the values prior to the hash mark
                   post_hashsignstr = activechar(find(checkhashsign)+1:end); % Get the values after the hash mark
                   
                   % search for the target word, plus a few variants and possible spelling mistakes
                   if  isequal(pre_hashsignchar,'Sample') | isequal(pre_hashsignchar,'sample') | isequal(pre_hashsignchar,'Samples') | isequal(pre_hashsignchar,'samples') | isequal(pre_hashsignchar,'sampls')  | isequal(pre_hashsignchar,'Sampls') | isequal(pre_hashsignchar,'sampl')  | isequal(pre_hashsignchar,'Sampl')
                       % positions 1:4 are the distribution type, and
                       % positions 5:7 give the sub-type (min, max, men, sig)
                       if     isequal(variablename{xxi}(5:7),'Min')
                              activedata{xxi,yyi} = -123456789;
                       elseif isequal(variablename{xxi}(5:7),'Max')
                              activedata{xxi,yyi} = 123456789;
                       elseif isequal(variablename{xxi}(5:7),'Men')
                              activedata{xxi,yyi} = 123456789i;
                       elseif isequal(variablename{xxi}(5:7),'Sig')
                              disp(sprintf('warning! the standard deviation of an end-member distribution cannot be defined relative to the samples. please fix entry "%s" in end-member "%s", which currently has the value "%s"',variablename{xxi},secondname{yyi},activedata{xxi,yyi}));                          
                       end
                    
                       % add the actual value for offset into the addvalumatrix
                       post_hashsignnum    = str2num(post_hashsignstr);
                       if ~isempty(post_hashsignnum)
                           addvaluematrix(xxi,yyi) = post_hashsignnum;
                       else
                           disp(sprintf('warning! please fix entry "%s" in end-member "%s", which currently has the problematic value "%s".\nbecause MEANDIR did not find a clear number after the # symbol in the end-member spreadsheet, it will use a value of 0.',variablename{xxi},secondname{yyi},post_hashsignstr)) 
                           addvaluematrix(xxi,yyi) = 0;
                       end
                  
                   % search for the target word, plus a few variants and possible spelling mistakes
                   elseif isequal(pre_hashsignchar,'frac') | isequal(pre_hashsignchar,'Frac') | isequal(pre_hashsignchar,'fractionation') | isequal(pre_hashsignchar,'Fractionation') | isequal(pre_hashsignchar,'Fractination') | isequal(pre_hashsignchar,'fracs') 
                       % positions 1:4 are the distribution type, and
                       % positions 5:7 give the sub-type (min, max, men, sig)
                       if     isequal(variablename{xxi}(5:7),'Min')
                              activedata{xxi,yyi} = -314159265;
                       elseif isequal(variablename{xxi}(5:7),'Max')
                              activedata{xxi,yyi} = 314159265;
                       elseif isequal(variablename{xxi}(5:7),'Men')
                              activedata{xxi,yyi} = 314159265i;
                       elseif isequal(variablename{xxi}(5:7),'Sig')
                              disp(sprintf('warning! the standard deviation of an end-member distribution cannot be defined relative to the samples. please fix entry "%s" in end-member "%s", which currently has the value "%s"',variablename{xxi},secondname{yyi},activedata{xxi,yyi}));                          
                       end
                    
                       % add the actual value for offset into the addvalumatrix
                       post_hashsignnum    = str2num(post_hashsignstr);
                       if ~isempty(post_hashsignnum)
                           addvaluematrix(xxi,yyi) = post_hashsignnum;
                       else
                           disp(sprintf('warning! Please fix entry "%s" in end-member "%s", which currently has the problematic value "%s".\nbecause MEANDIR did not find a clear number after the # symbol in the end-member spreadsheet, it will use a value of 0.',variablename{xxi},secondname{yyi},post_hashsignstr)) 
                           addvaluematrix(xxi,yyi) = 0;
                       end

                    else
                        disp(sprintf('warning! Unknown string entry in the end-members spreadsheet. the only acceptable non-numeric entries are "sample" or "fractionation" followed by "#" and a number, such as "fractionation#0" or "sample#3".\nplease fix entry "%s" in end-member "%s", which currently has the problematic value "%s".',variablename{xxi},secondname{yyi},pre_hashsignchar))
                        activedata{xxi,yyi}      = 0;
                        addvaluematrix(xxi,yyi)  = 0;
                   end
            else
                % access this portion of the loop if there is a character entry without a # sign
                if     sum(checkhashsign)>1
                       disp('warning! wach end-member string entry needs only a SINGLE # sign! MEANDIR could not read the end-members sucessfully');
                       activedata{xxi,yyi}     = 0;
                       addvaluematrix(xxi,yyi) = 0;
                elseif sum(checkhashsign)<1
                       disp(sprintf('warning! wnd-member string entries need a # sign! current entry is: "%s". MEANDIR could not read the end-members sucessfully',activechar));
                       activedata{xxi,yyi}     = 0;
                       addvaluematrix(xxi,yyi) = 0;
                end
            end  % end instructions for if there is a # sign
        end      % end instructions for if the value is a string
    end          % end iteration loop over yy
    end          % end iteration loop over xx
    activedata = cell2mat(activedata);  
    
    % for active end-members (where there is at least 1 value 
    % that is not NaN), turn other NaN values to zeros
    numnonnan = sum(~isnan(activedata),1);    
    for nanchecki=1:length(numnonnan)
        if numnonnan(nanchecki)>0
           if sum(isnan(activedata(:,nanchecki)))>0
              activedata(isnan(activedata(:,nanchecki)),nanchecki) = 0;
           end
        end
    end    
    
    [l, w]         = size(activedata);
    rmindx         = (l==sum(isnan(activedata)));
    activedata     = activedata(:,~rmindx);
    EMshortname    = secondname(~rmindx);
    addvaluematrix = addvaluematrix(:,~rmindx);
    
    % from the list of variables, find the name that variable should be saved under (no slashes, sigmas, only one underscore)
    savenamelist = {};
    for i=1:length(variablename)
            aname = variablename{i}; 
            if isequal(class(aname),'logical')
               disp(sprintf('warning! one of the end-member names may be incorrect. check the #%i row of EM group %s in the end-members spreadsheet',i,namesofEMsets{EMsetindx}));
            end            
            cname = aname(~(ismember(aname,'/') | ismember(aname,' ')));
            cname(4) = '_';
            savenamelist{i} = cname;
    end

    % save the data into the end-member structure
    for j = 1:length(savenamelist)
        EM = setfield(EM,namesofEMsets{EMsetindx},'disttype',varnumlist{j},savenamelist{j}(1:3)); 
        for k = 1:length(EMshortname)
            
            % save the end-member values
            EM = setfield(EM,namesofEMsets{EMsetindx},EMshortname{k},savenamelist{j},activedata(j,k));

            % save the end-member offset relative to the sample
            addvalue = addvaluematrix(j,k);
            if ~isnan(addvalue)
                EM = setfield(EM,namesofEMsets{EMsetindx},EMshortname{k},[savenamelist{j} '_addvalue'],addvalue);
            end
        end
    end
         
end % end the loop over the names of EM sets
end % end the function
