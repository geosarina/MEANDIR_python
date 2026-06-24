function [precratios] = MEANDIR_QuantifyPrecRatios(ObsList, EMIterInst0, InvRivDat, ClCrit, precEM_indx, Cl_indx)

    % this script considers whether accounting for ClCritical values, given the current
    % end-member chemistry, river chemistry, and ClCritical value, would result
    % in a negative concentration for any dissolved constituents. the script
    % produces the variable "precratios", which is the ratio of Cl requested to
    % be derived from precipitation (up to the measured river Cl concentration)
    % relative to the amount of Cl in the river that could be supported by each
    % dissolved constituents. for example, one element of precratios may be: 
    % ClCritical/(river Ca * precipitation Cl/Ca). if this ratio is greater than 
    % 1, removing the requested amount of Cl in ClCritical would result in a 
    % negative Ca concentration. in this case, MEANDIR will re-draw the
    % chemistry of the precipitation end-member to respect the user-given value of 
    % ClCritical while also satisfying the condition that Ca does not become negative.

    % Unpack the InvRivDat data
    if sum(ismember(ObsList,'ALK'))  >0; InvRivDatALK  = InvRivDat.ALK;  end % ALK
    if sum(ismember(ObsList,'DIC'))  >0; InvRivDatDIC  = InvRivDat.DIC;  end % DIC
    if sum(ismember(ObsList,'Ca'))   >0; InvRivDatCa   = InvRivDat.Ca;   end % Ca
    if sum(ismember(ObsList,'Mg'))   >0; InvRivDatMg   = InvRivDat.Mg;   end % Mg
    if sum(ismember(ObsList,'Na'))   >0; InvRivDatNa   = InvRivDat.Na;   end % Na
    if sum(ismember(ObsList,'K'))    >0; InvRivDatK    = InvRivDat.K;    end % K
    if sum(ismember(ObsList,'Sr'))   >0; InvRivDatSr   = InvRivDat.Sr;   end % Sr
    if sum(ismember(ObsList,'Fe'))   >0; InvRivDatFe   = InvRivDat.Fe;   end % Fe
    if sum(ismember(ObsList,'Cl'))   >0; InvRivDatCl   = InvRivDat.Cl;   end % Cl
    if sum(ismember(ObsList,'SO4'))  >0; InvRivDatSO4  = InvRivDat.SO4;  end % SO4
    if sum(ismember(ObsList,'NO3'))  >0; InvRivDatNO3  = InvRivDat.NO3;  end % NO3
    if sum(ismember(ObsList,'PO4'))  >0; InvRivDatPO4  = InvRivDat.PO4;  end % PO4
    if sum(ismember(ObsList,'Si'))   >0; InvRivDatSi   = InvRivDat.Si;   end % Si
    if sum(ismember(ObsList,'Ge'))   >0; InvRivDatGe   = InvRivDat.Ge;   end % Ge
    if sum(ismember(ObsList,'Li'))   >0; InvRivDatLi   = InvRivDat.Li;   end % Li
    if sum(ismember(ObsList,'F'))    >0; InvRivDatF    = InvRivDat.F;    end % F
    if sum(ismember(ObsList,'B'))    >0; InvRivDatB    = InvRivDat.B;    end % B
    if sum(ismember(ObsList,'Re'))   >0; InvRivDatRe   = InvRivDat.Re;   end % Re
    if sum(ismember(ObsList,'Mo'))   >0; InvRivDatMo   = InvRivDat.Mo;   end % Mo
    if sum(ismember(ObsList,'Os'))   >0; InvRivDatOs   = InvRivDat.Os;   end % Os
    if sum(ismember(ObsList,'HCO3')) >0; InvRivDatHCO3 = InvRivDat.HCO3; end % HCO3
     
    % define the Cl value to test. this should be the amount of Cl that the remainder of
    % MEANDIR_ClCriticalCorrection is going to remove from this particular river water sample.
    if     ClCrit >= InvRivDatCl;  ClTest = InvRivDatCl; 
    elseif ClCrit <  InvRivDatCl;  ClTest = ClCrit;      
    else                           ClTest = NaN; 
                                   disp('waning! cannot identify the relative size of ClCritical and river Cl. this error was generated within "MEANDIR_QuantifyPrecRatios"');
    end

    % calculate whether this precipitation end-member would result in negative concentrations for any of the variables in IonList. if any of the
    % entries of precratios are >1, the current combination of river chemistry, precipitation chemistry, and ClCritical value would result in a negative concentration.
    count = 0; precratios = []; 
    if sum(ismember(ObsList,'ALK'))  >0; PrecEM_ClALK  = (EMIterInst0(Cl_indx,precEM_indx)./EMIterInst0(ismember(ObsList,'ALK'),  precEM_indx)); count = count+1; precratios(count) = ClTest/(InvRivDatALK  * PrecEM_ClALK);  end % ALK
    if sum(ismember(ObsList,'DIC'))  >0; PrecEM_ClDIC  = (EMIterInst0(Cl_indx,precEM_indx)./EMIterInst0(ismember(ObsList,'DIC'),  precEM_indx)); count = count+1; precratios(count) = ClTest/(InvRivDatDIC  * PrecEM_ClDIC);  end % DIC
    if sum(ismember(ObsList,'Ca'))   >0; PrecEM_ClCa   = (EMIterInst0(Cl_indx,precEM_indx)./EMIterInst0(ismember(ObsList,'Ca'),   precEM_indx)); count = count+1; precratios(count) = ClTest/(InvRivDatCa   * PrecEM_ClCa);   end % Ca
    if sum(ismember(ObsList,'Mg'))   >0; PrecEM_ClMg   = (EMIterInst0(Cl_indx,precEM_indx)./EMIterInst0(ismember(ObsList,'Mg'),   precEM_indx)); count = count+1; precratios(count) = ClTest/(InvRivDatMg   * PrecEM_ClMg);   end % Mg
    if sum(ismember(ObsList,'Na'))   >0; PrecEM_ClNa   = (EMIterInst0(Cl_indx,precEM_indx)./EMIterInst0(ismember(ObsList,'Na'),   precEM_indx)); count = count+1; precratios(count) = ClTest/(InvRivDatNa   * PrecEM_ClNa);   end % Na
    if sum(ismember(ObsList,'K'))    >0; PrecEM_ClK    = (EMIterInst0(Cl_indx,precEM_indx)./EMIterInst0(ismember(ObsList,'K'),    precEM_indx)); count = count+1; precratios(count) = ClTest/(InvRivDatK    * PrecEM_ClK);    end % K
    if sum(ismember(ObsList,'Sr'))   >0; PrecEM_ClSr   = (EMIterInst0(Cl_indx,precEM_indx)./EMIterInst0(ismember(ObsList,'Sr'),   precEM_indx)); count = count+1; precratios(count) = ClTest/(InvRivDatSr   * PrecEM_ClSr);   end % Sr
    if sum(ismember(ObsList,'Fe'))   >0; PrecEM_ClFe   = (EMIterInst0(Cl_indx,precEM_indx)./EMIterInst0(ismember(ObsList,'Fe'),   precEM_indx)); count = count+1; precratios(count) = ClTest/(InvRivDatFe   * PrecEM_ClFe);   end % Fe
    if sum(ismember(ObsList,'SO4'))  >0; PrecEM_ClSO4  = (EMIterInst0(Cl_indx,precEM_indx)./EMIterInst0(ismember(ObsList,'SO4'),  precEM_indx)); count = count+1; precratios(count) = ClTest/(InvRivDatSO4  * PrecEM_ClSO4);  end % SO4
    if sum(ismember(ObsList,'NO3'))  >0; PrecEM_ClNO3  = (EMIterInst0(Cl_indx,precEM_indx)./EMIterInst0(ismember(ObsList,'NO3'),  precEM_indx)); count = count+1; precratios(count) = ClTest/(InvRivDatNO3  * PrecEM_ClNO3);  end % NO3
    if sum(ismember(ObsList,'PO4'))  >0; PrecEM_ClPO4  = (EMIterInst0(Cl_indx,precEM_indx)./EMIterInst0(ismember(ObsList,'PO4'),  precEM_indx)); count = count+1; precratios(count) = ClTest/(InvRivDatPO4  * PrecEM_ClPO4);  end % PO4
    if sum(ismember(ObsList,'Si'))   >0; PrecEM_ClSi   = (EMIterInst0(Cl_indx,precEM_indx)./EMIterInst0(ismember(ObsList,'Si'),   precEM_indx)); count = count+1; precratios(count) = ClTest/(InvRivDatSi   * PrecEM_ClSi);   end % Si
    if sum(ismember(ObsList,'Ge'))   >0; PrecEM_ClGe   = (EMIterInst0(Cl_indx,precEM_indx)./EMIterInst0(ismember(ObsList,'Ge'),   precEM_indx)); count = count+1; precratios(count) = ClTest/(InvRivDatGe   * PrecEM_ClGe);   end % Ge
    if sum(ismember(ObsList,'Li'))   >0; PrecEM_ClLi   = (EMIterInst0(Cl_indx,precEM_indx)./EMIterInst0(ismember(ObsList,'Li'),   precEM_indx)); count = count+1; precratios(count) = ClTest/(InvRivDatLi   * PrecEM_ClLi);   end % Li
    if sum(ismember(ObsList,'F'))    >0; PrecEM_ClF    = (EMIterInst0(Cl_indx,precEM_indx)./EMIterInst0(ismember(ObsList,'F'),    precEM_indx)); count = count+1; precratios(count) = ClTest/(InvRivDatF    * PrecEM_ClF);    end % F
    if sum(ismember(ObsList,'B'))    >0; PrecEM_ClB    = (EMIterInst0(Cl_indx,precEM_indx)./EMIterInst0(ismember(ObsList,'B'),    precEM_indx)); count = count+1; precratios(count) = ClTest/(InvRivDatB    * PrecEM_ClB);    end % B
    if sum(ismember(ObsList,'Re'))   >0; PrecEM_ClRe   = (EMIterInst0(Cl_indx,precEM_indx)./EMIterInst0(ismember(ObsList,'Re'),   precEM_indx)); count = count+1; precratios(count) = ClTest/(InvRivDatRe   * PrecEM_ClRe);   end % ReK
    if sum(ismember(ObsList,'Mo'))   >0; PrecEM_ClMo   = (EMIterInst0(Cl_indx,precEM_indx)./EMIterInst0(ismember(ObsList,'Mo'),   precEM_indx)); count = count+1; precratios(count) = ClTest/(InvRivDatMo   * PrecEM_ClMo);   end % Mo
    if sum(ismember(ObsList,'Os'))   >0; PrecEM_ClOs   = (EMIterInst0(Cl_indx,precEM_indx)./EMIterInst0(ismember(ObsList,'Os'),   precEM_indx)); count = count+1; precratios(count) = ClTest/(InvRivDatOs   * PrecEM_ClOs);   end % Os
    if sum(ismember(ObsList,'HCO3')) >0; PrecEM_ClHCO3 = (EMIterInst0(Cl_indx,precEM_indx)./EMIterInst0(ismember(ObsList,'HCO3'), precEM_indx)); count = count+1; precratios(count) = ClTest/(InvRivDatHCO3 * PrecEM_ClHCO3); end % HCO3

end % end of function