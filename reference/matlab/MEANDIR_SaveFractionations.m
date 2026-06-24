% this script saves the fractionations used in the inversion 
% the script first defines matricies to hold data on fractionation factors, then populates those matricies with values
% this script is run near the end of MEANDIR, when the simulation results are being stored within river

% first, identify where there might be fractionation and make matricies to hold that information
fractionation_dist     = [ismember(distEMdist0,'EPS-NOR') | ismember(distEMdist0,'EPS-UNI') | ismember(distEMdist0,'EPS-LGU')];
fractionation_location = find(fractionation_dist);
fractionation_N        = length(fractionation_location);

% if there are fractionations in this inversion, set "foundfrac" equal to 1
found_fractionations = 0; 
if fractionation_N > 0
   found_fractionations = 1; 
end

% loop over each fractionation and define a matrix for that combination of
% end-member and isotope ratios, with dimension of samples x remaining instances
for jj=1:fractionation_N        
       active_location = fractionation_location(jj);
       count = 0;
       while active_location>0
             active_location = active_location-size(fractionation_dist,1);
             count = count+1; 
       end
       active_location = active_location+size(fractionation_dist,1);
       active_iso       = ObsList{active_location};
       active_EM        = EMList0{count};
       evalin('base',[sprintf('river.fractionation.%s.%s  =',active_EM, active_iso) 'NaN(s,remaininginstances);']);
end % end of iterating over instances of fractionation
         
        
% second, iterate through each sample and each simulation to find and save fractionation information
for i=1:s                   % loop over samples
for j=1:remaininginstances  % loop over remaininginstances
        
    % pull each individual inversion matrix
    inversionmatrix               = EMSucessIterMatrixAll(:,:,j,i);
    inversionmatrix_fractionation = EMSucessIterMatrixAll_frac(:,:,j,i);
    
    % proceed if the inversion matrix is not zero
    if sum(~isnan(inversionmatrix(:))) > 0
       
       % find differences between the two inversion matricies
       difference_matrix    = (inversionmatrix~=inversionmatrix_fractionation);
       difference_locations = find(difference_matrix); % locate those differences
       
       % iterate over each instance of fractionation
       for  k=1:length(difference_locations)
           
            % find the positions of the mismatch (where there is fractionation)
            active_location = difference_locations(k);
            count = 0;
            while active_location>0
                  active_location = active_location-size(difference_matrix,1);
                  count = count+1; 
            end
            active_location = active_location+size(difference_matrix,1);
                      
            active_isopos = active_location;
            active_EMpos  = count; 
            active_iso    = ObsList{active_isopos};
            active_EM     = EMList0{active_EMpos};               
            if isequal(ObsList{active_isopos},'d7Li');                                      active_ion = 'Li';   active_ionpos = find(ismember(ObsList,'Li'));   end % if d7Li, find Li
            if isequal(ObsList{active_isopos},'d18O');                                      active_ion = 'SO4';  active_ionpos = find(ismember(ObsList,'SO4'));  end % if d18O, find SO4
            if isequal(ObsList{active_isopos},'d26Mg');                                     active_ion = 'Mg';   active_ionpos = find(ismember(ObsList,'Mg'));   end % if d26Mg, find Mg
            if isequal(ObsList{active_isopos},'d30Si');                                     active_ion = 'Si';   active_ionpos = find(ismember(ObsList,'Si'));   end % if d30Si, find Si
            if isequal(ObsList{active_isopos},'d34S');                                      active_ion = 'SO4';  active_ionpos = find(ismember(ObsList,'SO4'));  end % if d34S, find SO4
            if isequal(ObsList{active_isopos},'d42Ca');                                     active_ion = 'Ca';   active_ionpos = find(ismember(ObsList,'Ca'));   end % if d42Ca, find Ca
            if isequal(ObsList{active_isopos},'d44Ca');                                     active_ion = 'Ca';   active_ionpos = find(ismember(ObsList,'Ca'));   end % if d44Ca, find Ca
            if isequal(ObsList{active_isopos},'d56Fe');                                     active_ion = 'Fe';   active_ionpos = find(ismember(ObsList,'Fe'));   end % if d56Fe, find Fe
            if isequal(ObsList{active_isopos},'Sr8786');                                    active_ion = 'Sr';   active_ionpos = find(ismember(ObsList,'Sr'));   end % if Sr87/86, find Sr
            if isequal(ObsList{active_isopos},'Os8788');                                    active_ion = 'Os';   active_ionpos = find(ismember(ObsList,'Os'));   end % if Os8788, find Os
            if isequal(ObsList{active_isopos},'d98Mo');                                     active_ion = 'Mo';   active_ionpos = find(ismember(ObsList,'Mo'));   end % if d98Mo, find Mo
            if isequal(ObsList{active_isopos},'d13C') & isequal(carbonisotopematch,'HCO3'); active_ion = 'HCO3'; active_ionpos = find(ismember(ObsList,'HCO3')); end % if d13C and carbon source is HCO3, find HCO3
            if isequal(ObsList{active_isopos},'Fmod') & isequal(carbonisotopematch,'HCO3'); active_ion = 'HCO3'; active_ionpos = find(ismember(ObsList,'HCO3')); end % if Fmod and carbon source is HCO3, find HCO3
            if isequal(ObsList{active_isopos},'d13C') & isequal(carbonisotopematch,'DIC');  active_ion = 'DIC';  active_ionpos = find(ismember(ObsList,'DIC'));  end % if d13C and carbon source is DIC, find HCO3
            if isequal(ObsList{active_isopos},'Fmod') & isequal(carbonisotopematch,'DIC');  active_ion = 'DIC';  active_ionpos = find(ismember(ObsList,'DIC'));  end % if Fmod and carbon source is DIC, find HCO3 
           
            % calculate the fractionation factor for this sample and simulation
            concdata          = inversionmatrix_fractionation(active_ionpos,active_EMpos); % concentration data
            isotopeconcdata   = inversionmatrix_fractionation(active_isopos,active_EMpos); % isotope value
            fractionationdata = isotopeconcdata./concdata;
           
            % if delta values and fractionations were converted to isotopic ratios for the inversion, convert them back to delta values for storing inversion results
            if isequal(ObsList{active_isopos},'d7Li')  & sum(ismember(ConvertDelta2RList,'d7Li'))  >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d7Li'));  fractionationdata = (fractionationdata/conv)*1000;  end % d7Li
            if isequal(ObsList{active_isopos},'d18O')  & sum(ismember(ConvertDelta2RList,'d18O'))  >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d18O'));  fractionationdata = (fractionationdata/conv)*1000;  end % d18O
            if isequal(ObsList{active_isopos},'d26Mg') & sum(ismember(ConvertDelta2RList,'d26Mg')) >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d26Mg')); fractionationdata = (fractionationdata/conv)*1000;  end % d26Mg
            if isequal(ObsList{active_isopos},'d30Si') & sum(ismember(ConvertDelta2RList,'d30Si')) >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d30Si')); fractionationdata = (fractionationdata/conv)*1000;  end % d30Si
            if isequal(ObsList{active_isopos},'d34S')  & sum(ismember(ConvertDelta2RList,'d34S'))  >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d34S'));  fractionationdata = (fractionationdata/conv)*1000;  end % d34S
            if isequal(ObsList{active_isopos},'d42Ca') & sum(ismember(ConvertDelta2RList,'d42Ca')) >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d42Ca')); fractionationdata = (fractionationdata/conv)*1000;  end % d42Ca
            if isequal(ObsList{active_isopos},'d44Ca') & sum(ismember(ConvertDelta2RList,'d44Ca')) >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d44Ca')); fractionationdata = (fractionationdata/conv)*1000;  end % d44Ca
            if isequal(ObsList{active_isopos},'d56Fe') & sum(ismember(ConvertDelta2RList,'d56Fe')) >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d56Fe')); fractionationdata = (fractionationdata/conv)*1000;  end % d56Fe
            if isequal(ObsList{active_isopos},'d98Mo') & sum(ismember(ConvertDelta2RList,'d98Mo')) >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d98Mo')); fractionationdata = (fractionationdata/conv)*1000;  end % d98Mo
            if isequal(ObsList{active_isopos},'d13C')  & sum(ismember(ConvertDelta2RList,'d13C'))  >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d13C'));  fractionationdata = (fractionationdata/conv)*1000;  end % d13C

            % save the fractionation data
            evalin('base',[sprintf('river.fractionation.%s.%s(i,j) =',active_EM, active_iso) 'fractionationdata;']);
           
       end  % end iteration over number of fractionations
    end     % end check of whether inversion matrix is all NaN
end         % end iteration of remaining instances
end         % end iteration over samples

% display information for the user
if found_fractionations==1
    disp('saved fractionation factors into the field "fractionation"');    
end