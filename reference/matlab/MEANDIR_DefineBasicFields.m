function out = MEANDIR_DefineBasicFields

         % this function defines the basic vectors of the "observations" field of the structure "river"
         % the script is called in the function MEANDIR_ReadRiverDataIntoMatlab while reading river observations

         % define a structure
         out = struct; 

         % define the variables for concentration and charge equivalents 
         out = setfield(out,'model_variable','ClCr',          NaN);
         out = setfield(out,'observations','ALK_conc',        NaN);      out = setfield(out,'observations','ALK_equi',        NaN); % ALK
         out = setfield(out,'observations','DIC_conc',        NaN);      out = setfield(out,'observations','DIC_equi',        NaN); % DIC
         out = setfield(out,'observations','Ca_conc',         NaN);      out = setfield(out,'observations','Ca_equi',         NaN); % Ca
         out = setfield(out,'observations','Mg_conc',         NaN);      out = setfield(out,'observations','Mg_equi',         NaN); % Mg
         out = setfield(out,'observations','Na_conc',         NaN);      out = setfield(out,'observations','Na_equi',         NaN); % Na
         out = setfield(out,'observations','K_conc',          NaN);      out = setfield(out,'observations','K_equi',          NaN); % K
         out = setfield(out,'observations','Sr_conc',         NaN);      out = setfield(out,'observations','Sr_equi',         NaN); % Sr
         out = setfield(out,'observations','Fe_conc',         NaN);      out = setfield(out,'observations','Fe_equi',         NaN); % Fe
         out = setfield(out,'observations','Cl_conc',         NaN);      out = setfield(out,'observations','Cl_equi',         NaN); % Cl
         out = setfield(out,'observations','SO4_conc',        NaN);      out = setfield(out,'observations','SO4_equi',        NaN); % SO4
         out = setfield(out,'observations','NO3_conc',        NaN);      out = setfield(out,'observations','NO3_equi',        NaN); % NO3
         out = setfield(out,'observations','PO4_conc',        NaN);      out = setfield(out,'observations','PO4_equi',        NaN); % PO4
         out = setfield(out,'observations','Si_conc',         NaN);      out = setfield(out,'observations','Si_equi',         NaN); % Si
         out = setfield(out,'observations','Ge_conc',         NaN);      out = setfield(out,'observations','Ge_equi',         NaN); % Ge
         out = setfield(out,'observations','Li_conc',         NaN);      out = setfield(out,'observations','Li_equi',         NaN); % Li
         out = setfield(out,'observations','F_conc',          NaN);      out = setfield(out,'observations','F_equi',          NaN); % F
         out = setfield(out,'observations','B_conc',          NaN);      out = setfield(out,'observations','B_equi',          NaN); % B
         out = setfield(out,'observations','Re_conc',         NaN);      out = setfield(out,'observations','Re_equi',         NaN); % Re
         out = setfield(out,'observations','Mo_conc',         NaN);      out = setfield(out,'observations','Mo_equi',         NaN); % Mo
         out = setfield(out,'observations','Os_conc',         NaN);      out = setfield(out,'observations','Os_equi',         NaN); % Os
         out = setfield(out,'observations','HCO3_conc',       NaN);      out = setfield(out,'observations','HCO3_equi',       NaN); % HCO3
         out = setfield(out,'observations','d7Li_conc',       NaN);      out = setfield(out,'observations','d7Li_equi',       NaN); % d7Li
         out = setfield(out,'observations','d13C_conc',       NaN);      out = setfield(out,'observations','d13C_equi',       NaN); % d13C
         out = setfield(out,'observations','d18O_conc',       NaN);      out = setfield(out,'observations','d18O_equi',       NaN); % d18O
         out = setfield(out,'observations','d26Mg_conc',      NaN);      out = setfield(out,'observations','d26Mg_equi',      NaN); % d26Mg
         out = setfield(out,'observations','d30Si_conc',      NaN);      out = setfield(out,'observations','d30Si_equi',      NaN); % d30Si
         out = setfield(out,'observations','d34S_conc',       NaN);      out = setfield(out,'observations','d34S_equi',       NaN); % d34S
         out = setfield(out,'observations','d42Ca_conc',      NaN);      out = setfield(out,'observations','d42Ca_equi',      NaN); % d42Ca
         out = setfield(out,'observations','d44Ca_conc',      NaN);      out = setfield(out,'observations','d44Ca_equi',      NaN); % d44Ca
         out = setfield(out,'observations','d56Fe_conc',      NaN);      out = setfield(out,'observations','d56Fe_equi',      NaN); % d56Fe
         out = setfield(out,'observations','Sr8786_conc',     NaN);      out = setfield(out,'observations','Sr8786_equi',     NaN); % Sr8786
         out = setfield(out,'observations','d98Mo_conc',      NaN);      out = setfield(out,'observations','d98Mo_equi',      NaN); % d98Mo
         out = setfield(out,'observations','Os8788_conc',     NaN);      out = setfield(out,'observations','Os8788_equi',     NaN); % Os8788
         out = setfield(out,'observations','Fmod_conc',       NaN);      out = setfield(out,'observations','Fmod_equi',       NaN); % Fmod

         % define the uncertainty variables for concentration and charge equivalents 
         out = setfield(out,'observations','ALK_conc_asd',    NaN);      out = setfield(out,'observations','ALK_equi_asd',    NaN); % ALK
         out = setfield(out,'observations','DIC_conc_asd',    NaN);      out = setfield(out,'observations','DIC_equi_asd',    NaN); % DIC
         out = setfield(out,'observations','Ca_conc_asd',     NaN);      out = setfield(out,'observations','Ca_equi_asd',     NaN); % Ca
         out = setfield(out,'observations','Mg_conc_asd',     NaN);      out = setfield(out,'observations','Mg_equi_asd',     NaN); % Mg
         out = setfield(out,'observations','Na_conc_asd',     NaN);      out = setfield(out,'observations','Na_equi_asd',     NaN); % Na
         out = setfield(out,'observations','K_conc_asd',      NaN);      out = setfield(out,'observations','K_equi_asd',      NaN); % K
         out = setfield(out,'observations','Sr_conc_asd',     NaN);      out = setfield(out,'observations','Sr_equi_asd',     NaN); % Sr
         out = setfield(out,'observations','Fe_conc_asd',     NaN);      out = setfield(out,'observations','Fe_equi_asd',     NaN); % Fe
         out = setfield(out,'observations','Cl_conc_asd',     NaN);      out = setfield(out,'observations','Cl_equi_asd',     NaN); % Cl
         out = setfield(out,'observations','SO4_conc_asd',    NaN);      out = setfield(out,'observations','SO4_equi_asd',    NaN); % SO4
         out = setfield(out,'observations','NO3_conc_asd',    NaN);      out = setfield(out,'observations','NO3_equi_asd',    NaN); % NO3
         out = setfield(out,'observations','PO4_conc_asd',    NaN);      out = setfield(out,'observations','PO4_equi_asd',    NaN); % PO4
         out = setfield(out,'observations','Si_conc_asd',     NaN);      out = setfield(out,'observations','Si_equi_asd',     NaN); % Si
         out = setfield(out,'observations','Ge_conc_asd',     NaN);      out = setfield(out,'observations','Ge_equi_asd',     NaN); % Ge
         out = setfield(out,'observations','Li_conc_asd',     NaN);      out = setfield(out,'observations','Li_equi_asd',     NaN); % Li
         out = setfield(out,'observations','F_conc_asd',      NaN);      out = setfield(out,'observations','F_equi_asd',      NaN); % F
         out = setfield(out,'observations','B_conc_asd',      NaN);      out = setfield(out,'observations','B_equi_asd',      NaN); % B
         out = setfield(out,'observations','Re_conc_asd',     NaN);      out = setfield(out,'observations','Re_equi_asd',     NaN); % Re
         out = setfield(out,'observations','Mo_conc_asd',     NaN);      out = setfield(out,'observations','Mo_equi_asd',     NaN); % Mo
         out = setfield(out,'observations','Os_conc_asd',     NaN);      out = setfield(out,'observations','Os_equi_asd',     NaN); % Os
         out = setfield(out,'observations','HCO3_conc_asd',   NaN);      out = setfield(out,'observations','HCO3_equi_asd',   NaN); % HCO3
         out = setfield(out,'observations','d7Li_conc_asd',   NaN);      out = setfield(out,'observations','d7Li_equi_asd',   NaN); % d7Li
         out = setfield(out,'observations','d13C_conc_asd',   NaN);      out = setfield(out,'observations','d13C_equi_asd',   NaN); % d13C
         out = setfield(out,'observations','d18O_conc_asd',   NaN);      out = setfield(out,'observations','d18O_equi_asd',   NaN); % d18O
         out = setfield(out,'observations','d26Mg_conc_asd',  NaN);      out = setfield(out,'observations','d26Mg_equi_asd',  NaN); % d26Mg
         out = setfield(out,'observations','d30Si_conc_asd',  NaN);      out = setfield(out,'observations','d30Si_equi_asd',  NaN); % d30Si
         out = setfield(out,'observations','d34S_conc_asd',   NaN);      out = setfield(out,'observations','d34S_equi_asd',   NaN); % d34S
         out = setfield(out,'observations','d42Ca_conc_asd',  NaN);      out = setfield(out,'observations','d42Ca_equi_asd',  NaN); % d42Ca
         out = setfield(out,'observations','d44Ca_conc_asd',  NaN);      out = setfield(out,'observations','d44Ca_equi_asd',  NaN); % d44Ca
         out = setfield(out,'observations','d56Fe_conc_asd',  NaN);      out = setfield(out,'observations','d56Fe_equi_asd',  NaN); % d56Fe
         out = setfield(out,'observations','Sr8786_conc_asd', NaN);      out = setfield(out,'observations','Sr8786_equi_asd', NaN); % Sr8786
         out = setfield(out,'observations','d98Mo_conc_asd',  NaN);      out = setfield(out,'observations','d98Mo_equi_asd',  NaN); % d98Mo
         out = setfield(out,'observations','Os8788_conc_asd', NaN);      out = setfield(out,'observations','Os8788_equi_asd', NaN); % Os8788
         out = setfield(out,'observations','Fmod_conc_asd',   NaN);      out = setfield(out,'observations','Fmod_equi_asd',   NaN); % Fmod

end % end of function