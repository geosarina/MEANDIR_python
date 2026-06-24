% this script processes the matrix of end-member values from successful simulations to calculate
% a median and range of each ion for each end-member for each sample in the analysis

% find the positions of isotopes
isolist = find(isotopeposition);

for j=1:length(EMList0) % loop over end-members
for k=1:length(ObsList) % loop over observations

    EMratios  = EMSucessIterMatrixAll(k,j,:,:); 
    EMratios  = reshape(EMratios,remaininginstances,s);

    if  sum(ismember(isolist,k))>0
        % if there is isotopic information in the inversion, same the standard isotopic information (in units of per mil, or raw
        % ratios) into the end-member, rather than the product of the isotopic information with the concentration information
        activeisoname = ObsList{k};
        if isequal(activeisoname,'d7Li');                                        activeionpos = find(ismember(ObsList,'Li'));   end % if d7Li, find the Li column
        if isequal(activeisoname,'d18O');                                        activeionpos = find(ismember(ObsList,'SO4'));  end % if d18O, find the SO4 column
        if isequal(activeisoname,'d26Mg');                                       activeionpos = find(ismember(ObsList,'Mg'));   end % if d26Mg, find the Mg column
        if isequal(activeisoname,'d30Si');                                       activeionpos = find(ismember(ObsList,'Si'));   end % if d30Si, find the Si column
        if isequal(activeisoname,'d34S');                                        activeionpos = find(ismember(ObsList,'SO4'));  end % if d34S, find the SO4 column
        if isequal(activeisoname,'d42Ca');                                       activeionpos = find(ismember(ObsList,'Ca'));   end % if d42Ca, find the Ca column
        if isequal(activeisoname,'d44Ca');                                       activeionpos = find(ismember(ObsList,'Ca'));   end % if d44Ca, find the Ca column
        if isequal(activeisoname,'d56Fe');                                       activeionpos = find(ismember(ObsList,'Fe'));   end % if d56Fe, find the Fe column
        if isequal(activeisoname,'Sr8786');                                      activeionpos = find(ismember(ObsList,'Sr'));   end % if Sr87/86, find the Sr column
        if isequal(activeisoname,'d98Mo');                                       activeionpos = find(ismember(ObsList,'Mo'));   end % if d98Mo, find the Mo column
        if isequal(activeisoname,'Os8788');                                      activeionpos = find(ismember(ObsList,'Os'));   end % if Os8788, find the Os column
        if isequal(activeisoname,'d13C') & isequal(carbonisotopematch,'HCO3');   activeionpos = find(ismember(ObsList,'HCO3')); end % if d13C and the carbon source is HCO3, find the HCO3 column
        if isequal(activeisoname,'Fmod') & isequal(carbonisotopematch,'HCO3');   activeionpos = find(ismember(ObsList,'HCO3')); end % if Fmod and the carbon source HCO3, find the HCO3 column
        if isequal(activeisoname,'d13C') & isequal(carbonisotopematch,'DIC');    activeionpos = find(ismember(ObsList,'DIC'));  end % if d13C and the carbon source DIC, find the DIC column
        if isequal(activeisoname,'Fmod') & isequal(carbonisotopematch,'DIC');    activeionpos = find(ismember(ObsList,'DIC'));  end % if Fmod and the carbon source DIC, find the DIC column
        isodata  = EMSucessIterMatrixAll(k,j,:,:);
        isodata  = reshape(isodata,remaininginstances,s);
        concdata = EMSucessIterMatrixAll(activeionpos,j,:,:);
        concdata = reshape(concdata,remaininginstances,s);
        EMratios = isodata./concdata;

        % if delta values and fractionations were converted to isotopic ratios, convert them back
        if isequal(activeisoname,'d7Li')  & sum(ismember(ConvertDelta2RList,'d7Li'))  >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d7Li'));  EMratios = (EMratios/conv-1)*1000; end % d7Li
        if isequal(activeisoname,'d18O')  & sum(ismember(ConvertDelta2RList,'d18O'))  >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d18O'));  EMratios = (EMratios/conv-1)*1000; end % d18O
        if isequal(activeisoname,'d26Mg') & sum(ismember(ConvertDelta2RList,'d26Mg')) >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d26Mg')); EMratios = (EMratios/conv-1)*1000; end % d26Mg
        if isequal(activeisoname,'d30Si') & sum(ismember(ConvertDelta2RList,'d30Si')) >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d30Si')); EMratios = (EMratios/conv-1)*1000; end % d30Si
        if isequal(activeisoname,'d34S')  & sum(ismember(ConvertDelta2RList,'d34S'))  >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d34S'));  EMratios = (EMratios/conv-1)*1000; end % d34S
        if isequal(activeisoname,'d42Ca') & sum(ismember(ConvertDelta2RList,'d42Ca')) >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d42Ca')); EMratios = (EMratios/conv-1)*1000; end % d42Ca
        if isequal(activeisoname,'d44Ca') & sum(ismember(ConvertDelta2RList,'d44Ca')) >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d44Ca')); EMratios = (EMratios/conv-1)*1000; end % d44Ca
        if isequal(activeisoname,'d56Fe') & sum(ismember(ConvertDelta2RList,'d56Fe')) >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d56Fe')); EMratios = (EMratios/conv-1)*1000; end % d56Fe
        if isequal(activeisoname,'d98Mo') & sum(ismember(ConvertDelta2RList,'d98Mo')) >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d98Mo')); EMratios = (EMratios/conv-1)*1000; end % d98Mo
        if isequal(activeisoname,'d13C')  & sum(ismember(ConvertDelta2RList,'d13C'))  >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d13C'));  EMratios = (EMratios/conv-1)*1000; end % d13C
    end

    currentmedianlist = [];
    currentlist05     = [];
    currentlist25     = [];
    currentlist75     = [];
    currentlist95     = [];

    currentmedianlist(1:s,1) = nanmedian(EMratios)';
    currentlist05(1:s,1)     = prctile(EMratios,05)';
    currentlist25(1:s,1)     = prctile(EMratios,25)';
    currentlist75(1:s,1)     = prctile(EMratios,75)';
    currentlist95(1:s,1)     = prctile(EMratios,95)';
    allEMinfo = EMratios';      

    % save the end-member values
    if   sum(ismember(isolist,k))>0
         evalin('base',[sprintf('river.EndMembers.median.%s.%s =',EMList0{j},ObsList{k}) 'currentmedianlist;']);
         evalin('base',[sprintf('river.EndMembers.pct05.%s.%s  =',EMList0{j},ObsList{k}) 'currentlist05;']);
         evalin('base',[sprintf('river.EndMembers.pct25.%s.%s  =',EMList0{j},ObsList{k}) 'currentlist25;']);
         evalin('base',[sprintf('river.EndMembers.pct75.%s.%s  =',EMList0{j},ObsList{k}) 'currentlist75;']);
         evalin('base',[sprintf('river.EndMembers.pct95.%s.%s  =',EMList0{j},ObsList{k}) 'currentlist95;']);
         evalin('base',[sprintf('river.EndMembers.all.%s.%s    =',EMList0{j},ObsList{k}) 'allEMinfo;']);
    else
         evalin('base',[sprintf('river.EndMembers.median.%s.%s_%s =',EMList0{j},ObsList{k},NormalizationType) 'currentmedianlist;']);
         evalin('base',[sprintf('river.EndMembers.pct05.%s.%s_%s  =',EMList0{j},ObsList{k},NormalizationType) 'currentlist05;']);
         evalin('base',[sprintf('river.EndMembers.pct25.%s.%s_%s  =',EMList0{j},ObsList{k},NormalizationType) 'currentlist25;']);
         evalin('base',[sprintf('river.EndMembers.pct75.%s.%s_%s  =',EMList0{j},ObsList{k},NormalizationType) 'currentlist75;']);
         evalin('base',[sprintf('river.EndMembers.pct95.%s.%s_%s  =',EMList0{j},ObsList{k},NormalizationType) 'currentlist95;']);
         evalin('base',[sprintf('river.EndMembers.all.%s.%s_%s    =',EMList0{j},ObsList{k},NormalizationType) 'allEMinfo;']);
    end
end % end iteration on ObsList (k=1:length(ObsList))
end % end iteration on EMList0 (j=1:length(EMList0))

% remove clutter from the workspace
clear currentmedianlist;
clear currentlist05;
clear currentlist25;
clear currentlist75;
clear currentlist95;
clear allEMinfo;

% tell the user that information on the chemistry of inversion-constrained
% end-members is saved in the field "EndMembers"
disp('calculated the ionic and isotopic ratios of inversion-constrained end-members, saved in the field "EndMembers"');