% this script prepares the results of the inversion model to be
% unpacked into 'results' structures, when IterateOver equals "End-members"
    
% the IterEM_allsmpfrac.* have size of (nEM, numberiterations, s)
% the IterEM_allInvRivDat.* have size of (s, numberiteration)
                   IterEM_allsmpfrac = struct;
                   IterEM_allsmpfrac.norm   = IterEM_smpfracStorenorm;   IterEM_allInvRivDat.norm   = IterEM_InvRivDataStore_norm;       % norm
                   IterEM_allsmpfrac.SO4    = IterEM_smpfracStoreSO4;    IterEM_allInvRivDat.SO4    = IterEM_InvRivDataStore_SO4;        % SO4
                   IterEM_allsmpfrac.d34S   = IterEM_smpfracStored34S;   IterEM_allInvRivDat.d34S   = IterEM_InvRivDataStore_d34S;       % d34S
if go.ALK    == 1; IterEM_allsmpfrac.ALK    = IterEM_smpfracStoreALK;    IterEM_allInvRivDat.ALK    = IterEM_InvRivDataStore_ALK;    end % ALK
if go.DIC    == 1; IterEM_allsmpfrac.DIC    = IterEM_smpfracStoreDIC;    IterEM_allInvRivDat.DIC    = IterEM_InvRivDataStore_DIC;    end % DIC
if go.Ca     == 1; IterEM_allsmpfrac.Ca     = IterEM_smpfracStoreCa;     IterEM_allInvRivDat.Ca     = IterEM_InvRivDataStore_Ca;     end % Ca
if go.Mg     == 1; IterEM_allsmpfrac.Mg     = IterEM_smpfracStoreMg;     IterEM_allInvRivDat.Mg     = IterEM_InvRivDataStore_Mg;     end % Mg
if go.Na     == 1; IterEM_allsmpfrac.Na     = IterEM_smpfracStoreNa;     IterEM_allInvRivDat.Na     = IterEM_InvRivDataStore_Na;     end % Na
if go.K      == 1; IterEM_allsmpfrac.K      = IterEM_smpfracStoreK;      IterEM_allInvRivDat.K      = IterEM_InvRivDataStore_K;      end % K
if go.Sr     == 1; IterEM_allsmpfrac.Sr     = IterEM_smpfracStoreSr;     IterEM_allInvRivDat.Sr     = IterEM_InvRivDataStore_Sr;     end % Sr
if go.Fe     == 1; IterEM_allsmpfrac.Fe     = IterEM_smpfracStoreFe;     IterEM_allInvRivDat.Fe     = IterEM_InvRivDataStore_Fe;     end % Fe
if go.Cl     == 1; IterEM_allsmpfrac.Cl     = IterEM_smpfracStoreCl;     IterEM_allInvRivDat.Cl     = IterEM_InvRivDataStore_Cl;     end % Cl
if go.NO3    == 1; IterEM_allsmpfrac.NO3    = IterEM_smpfracStoreNO3;    IterEM_allInvRivDat.NO3    = IterEM_InvRivDataStore_NO3;    end % NO3
if go.PO4    == 1; IterEM_allsmpfrac.PO4    = IterEM_smpfracStorePO4;    IterEM_allInvRivDat.PO4    = IterEM_InvRivDataStore_PO4;    end % PO4
if go.Si     == 1; IterEM_allsmpfrac.Si     = IterEM_smpfracStoreSi;     IterEM_allInvRivDat.Si     = IterEM_InvRivDataStore_Si;     end % Si
if go.Ge     == 1; IterEM_allsmpfrac.Ge     = IterEM_smpfracStoreGe;     IterEM_allInvRivDat.Ge     = IterEM_InvRivDataStore_Ge;     end % Ge
if go.Li     == 1; IterEM_allsmpfrac.Li     = IterEM_smpfracStoreLi;     IterEM_allInvRivDat.Li     = IterEM_InvRivDataStore_Li;     end % Li
if go.F      == 1; IterEM_allsmpfrac.F      = IterEM_smpfracStoreF;      IterEM_allInvRivDat.F      = IterEM_InvRivDataStore_F;      end % F
if go.B      == 1; IterEM_allsmpfrac.B      = IterEM_smpfracStoreB;      IterEM_allInvRivDat.B      = IterEM_InvRivDataStore_B;      end % B
if go.Re     == 1; IterEM_allsmpfrac.Re     = IterEM_smpfracStoreRe;     IterEM_allInvRivDat.Re     = IterEM_InvRivDataStore_Re;     end % Re
if go.Mo     == 1; IterEM_allsmpfrac.Mo     = IterEM_smpfracStoreMo;     IterEM_allInvRivDat.Mo     = IterEM_InvRivDataStore_Mo;     end % Mo
if go.Os     == 1; IterEM_allsmpfrac.Os     = IterEM_smpfracStoreOs;     IterEM_allInvRivDat.Os     = IterEM_InvRivDataStore_Os;     end % Os
if go.HCO3   == 1; IterEM_allsmpfrac.HCO3   = IterEM_smpfracStoreHCO3;   IterEM_allInvRivDat.HCO3   = IterEM_InvRivDataStore_HCO3;   end % HCO3
if go.d7Li   == 1; IterEM_allsmpfrac.d7Li   = IterEM_smpfracStored7Li;   IterEM_allInvRivDat.d7Li   = IterEM_InvRivDataStore_d7Li;   end % d7Li
if go.d13C   == 1; IterEM_allsmpfrac.d13C   = IterEM_smpfracStored13C;   IterEM_allInvRivDat.d13C   = IterEM_InvRivDataStore_d13C;   end % d13C
if go.d18O   == 1; IterEM_allsmpfrac.d18O   = IterEM_smpfracStored18O;   IterEM_allInvRivDat.d18O   = IterEM_InvRivDataStore_d18O;   end % d18O
if go.d26Mg  == 1; IterEM_allsmpfrac.d26Mg  = IterEM_smpfracStored26Mg;  IterEM_allInvRivDat.d26Mg  = IterEM_InvRivDataStore_d26Mg;  end % d26Mg
if go.d30Si  == 1; IterEM_allsmpfrac.d30Si  = IterEM_smpfracStored30Si;  IterEM_allInvRivDat.d30Si  = IterEM_InvRivDataStore_d30Si;  end % d30Si
if go.d42Ca  == 1; IterEM_allsmpfrac.d42Ca  = IterEM_smpfracStored42Ca;  IterEM_allInvRivDat.d42Ca  = IterEM_InvRivDataStore_d42Ca;  end % d42Ca
if go.d44Ca  == 1; IterEM_allsmpfrac.d44Ca  = IterEM_smpfracStored44Ca;  IterEM_allInvRivDat.d44Ca  = IterEM_InvRivDataStore_d44Ca;  end % d44Ca
if go.d56Fe  == 1; IterEM_allsmpfrac.d56Fe  = IterEM_smpfracStored56Fe;  IterEM_allInvRivDat.d56Fe  = IterEM_InvRivDataStore_d56Fe;  end % d56Fe
if go.Sr8786 == 1; IterEM_allsmpfrac.Sr8786 = IterEM_smpfracStoreSr8786; IterEM_allInvRivDat.Sr8786 = IterEM_InvRivDataStore_Sr8786; end % Sr8786
if go.d98Mo  == 1; IterEM_allsmpfrac.d98Mo  = IterEM_smpfracStored98Mo;  IterEM_allInvRivDat.d98Mo  = IterEM_InvRivDataStore_d98Mo;  end % d98Mo
if go.Os8788 == 1; IterEM_allsmpfrac.Os8788 = IterEM_smpfracStoreOs8788; IterEM_allInvRivDat.Os8788 = IterEM_InvRivDataStore_Os8788; end % Os8788
if go.Fmod   == 1; IterEM_allsmpfrac.Fmod   = IterEM_smpfracStoreFmod;   IterEM_allInvRivDat.Fmod   = IterEM_InvRivDataStore_Fmod;   end % Fmod
clear IterEM_smpfracStore*;     % clear variables that are now saved in other structures
clear IterEM_InvRivDataStore_*; % clear variables that are now saved in other structures
               
% organize the results for fractional contributions. this is effectivly just re-packing the structure so that it can be
% operated upon by MEANDIR_UnpackInversionResults analogous to the case when IterateOver equals "Samples"
smprescell = {}; for ii=1:s; smprescell{ii,1} = struct; end    
fn = fieldnames(IterEM_allsmpfrac);
for ii = 1:s
for jj = 1:length(fn)
    alldat = eval(sprintf('IterEM_allsmpfrac.%s',fn{jj}));
    dat = alldat(:,:,ii);        
    evalin('base',[sprintf('smprescell{%i,1}.smpresults%s =',ii,fn{jj}) '[];']);
    evalin('base',[sprintf('smprescell{%i,1}.smpresults%s =',ii,fn{jj}) 'dat;']);
end
end
clear IterEM_allsmpfrac; % remove workspace clutter

% organize the results for inverted data. this is effectivly just re-packing the structure so that it can be
% operated upon by MEANDIR_UnpackInversionResults analogous to the case when IterateOver equals "Samples"
InvRivDatacell = {}; for ii=1:s; InvRivDatacell{ii,1} = struct; end    
fn = fieldnames(IterEM_allInvRivDat);
for ii = 1:s
for jj = 1:length(fn)
    alldat = eval(sprintf('IterEM_allInvRivDat.%s',fn{jj}));
    dat = alldat(ii,:)';        
    evalin('base',[sprintf('InvRivDatacell{%i,1}.%s =',ii,fn{jj}) '[];']);
    evalin('base',[sprintf('InvRivDatacell{%i,1}.%s =',ii,fn{jj}) 'dat;']);
end
end
clear IterEM_InvRivDataStore; % remove workspace clutter