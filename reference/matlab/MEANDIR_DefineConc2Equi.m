function conc2equi = MEANDIR_DefineConc2Equi

         % this function defines the "conc2equi" structure that has conversion factors for moving between units
         % of charge equivalents and concentrations. for some dissolved variables, such as Ca or Na, this is a
         % simple stoichiometric factor of +2 or +1, respectively. it is possible to use uncharged species in inversions
         % conducted with charge eqvuialents, such as Si, by defining a conversion factor of +1. in this case having the 
         % normalization be in units of "equi" really means a mixture of charge and concentration. what matters is that each 
         % variable is always treated the same way in defining the river observation vector and the matrix of end-member compositions. 

         % open the spreadsheet containing the conversion factors
         [~, ~, all] = xlsread('MEANDIR_UserEntries','MEANDIR_conc2equi');
         header     = all(1,:);
         for i=1:length(header)
             if isnan(header{i})
                 header{i}      = 'NaN'; 
             end
         end

         % find the list of elements, isotopes, ratios, and standards
         variable   = all(2:end,ismember(header,'Variable'));
         factor     = all(2:end,ismember(header,'ConversionFactor'));

         % check the entries and remove unknowns
         for i=1:length(variable);  if isnan(variable{i});     variable{i}  = 'NaN'; end; end
         for i=1:length(factor);    if ~isnumeric(factor{i});  factor{i}    =  NaN;  end; end

         % convert the factor to numbers
         factor = cell2mat(factor);

         % save the factors for conversion of conc to equi
         clear conc2equi; 
         conc2equi = struct;
         for i=1:length(variable)
             conc2equi = setfield(conc2equi,variable{i},factor(i));
         end
         
end % end of function