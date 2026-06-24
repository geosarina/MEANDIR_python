function RiverMatrix0 = MEANDIR_GenerateRiverMatrix(river, nOL, ObsList, indxisopos, carbonisotopematch, ConvertDelta2RList, ConvertDelta2R)

        % this function generates the matrix "RiverMatrix0", which contains the
        % normalized elemental and isotopic ratios to be used during inversion. 

        RiverMatrix0 = []; 
        for i=1:nOL  % iterate over the number of observations in the inversion (recall that nOL = length(ObsList))

            % if the current variable is not isotopic information
            if    ~ismember(i,indxisopos)     
                   numerator         = eval(sprintf('river.model_variable.%s',ObsList{i}));  % get the river data for this variable
                   rivernorm         = river.model_variable.norm;                            % get the normalization 
                   RiverMatrix0(i,:) = numerator./rivernorm;                                 % save the normalized river data

            % if the current variable is isotopic information    
            elseif ismember(i,indxisopos)
                   deltavalue    = eval(sprintf('river.model_variable.%s',ObsList{i}));
                   activeisoname = ObsList{i}; % find the name of the current isotopic variable

                   if isequal(activeisoname,'d7Li');                                       activeionpos = find(ismember(ObsList,'Li'));    if sum(ismember(ConvertDelta2RList,'d7Li')) >0; conversionfactor = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d7Li'));  end; end % if d7Li, find Li 
                   if isequal(activeisoname,'d18O');                                       activeionpos = find(ismember(ObsList,'SO4'));   if sum(ismember(ConvertDelta2RList,'d18O')) >0; conversionfactor = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d18O'));  end; end % if d18O, find SO4 
                   if isequal(activeisoname,'d26Mg');                                      activeionpos = find(ismember(ObsList,'Mg'));    if sum(ismember(ConvertDelta2RList,'d26Mg'))>0; conversionfactor = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d26Mg')); end; end % if d26Mg, find Mg 
                   if isequal(activeisoname,'d30Si');                                      activeionpos = find(ismember(ObsList,'Si'));    if sum(ismember(ConvertDelta2RList,'d30Si'))>0; conversionfactor = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d30Si')); end; end % if d30Si, find Si
                   if isequal(activeisoname,'d34S');                                       activeionpos = find(ismember(ObsList,'SO4'));   if sum(ismember(ConvertDelta2RList,'d34S')) >0; conversionfactor = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d34S'));  end; end % if d34S, find SO4 
                   if isequal(activeisoname,'d42Ca');                                      activeionpos = find(ismember(ObsList,'Ca'));    if sum(ismember(ConvertDelta2RList,'d42Ca'))>0; conversionfactor = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d42Ca')); end; end % if d42Ca, find Ca 
                   if isequal(activeisoname,'d44Ca');                                      activeionpos = find(ismember(ObsList,'Ca'));    if sum(ismember(ConvertDelta2RList,'d44Ca'))>0; conversionfactor = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d44Ca')); end; end % if d44Ca, find Ca 
                   if isequal(activeisoname,'d56Fe');                                      activeionpos = find(ismember(ObsList,'Fe'));    if sum(ismember(ConvertDelta2RList,'d56Fe'))>0; conversionfactor = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d56Fe')); end; end % if d56Fe, find Fe 
                   if isequal(activeisoname,'Sr8786');                                     activeionpos = find(ismember(ObsList,'Sr'));                                                    conversionfactor = NaN;                                                                 end % if Sr87/86, find Sr 
                   if isequal(activeisoname,'Os8788');                                     activeionpos = find(ismember(ObsList,'Os'));                                                    conversionfactor = NaN;                                                                 end % if Os8788, find Os 
                   if isequal(activeisoname,'d98Mo');                                      activeionpos = find(ismember(ObsList,'Mo'));    if sum(ismember(ConvertDelta2RList,'d98Mo'))>0; conversionfactor = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d98Mo')); end; end % if d98Mo, find Mo 
                   if isequal(activeisoname,'d13C') & isequal(carbonisotopematch,'HCO3');  activeionpos = find(ismember(ObsList,'HCO3'));  if sum(ismember(ConvertDelta2RList,'d13C')) >0; conversionfactor = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d13C'));  end; end % if d13C and the carbon is HCO3, find HCO3 
                   if isequal(activeisoname,'Fmod') & isequal(carbonisotopematch,'HCO3');  activeionpos = find(ismember(ObsList,'HCO3'));                                                  conversionfactor = NaN;                                                                 end % if Fmod and the carbon is HCO3, find HCO3 
                   if isequal(activeisoname,'d13C') & isequal(carbonisotopematch,'DIC');   activeionpos = find(ismember(ObsList,'DIC'));   if sum(ismember(ConvertDelta2RList,'d13C')) >0; conversionfactor = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d13C'));  end; end % if d13C and the carbon is DIC, find HCO3 
                   if isequal(activeisoname,'Fmod') & isequal(carbonisotopematch,'DIC');   activeionpos = find(ismember(ObsList,'DIC'));                                                   conversionfactor = NaN;                                                                 end % if Fmod and the carbon is DIC, find HCO3 

                   % evaluate if the river data should be converted to isotopic ratios:
                   if     sum(ismember(ConvertDelta2RList,activeisoname))==0
                          % in this case, the river data should be the
                          % user-defined value (either a delta value or raw ratio)
                          scaledisotope = deltavalue;
                   elseif sum(ismember(ConvertDelta2RList,activeisoname))==1
                          % in this case, the river data should be the raw
                          % ratio, but where the river data spreadsheet
                          % contained information in delta notation
                          scaledisotope = (deltavalue/1000+1)*conversionfactor;
                   end               

                   numerator         = eval(sprintf('river.model_variable.%s',ObsList{activeionpos}));
                   rivernorm         = river.model_variable.norm;
                   RiverMatrix0(i,:) = scaledisotope.*numerator./rivernorm;
            end  % end of checking isotopic information            
        end      % end of iteration over observations
end              % end of function
