% when IterateOver equals End-members, this script defines the matricies IterEM_smpfracStore* and IterEM_InvRivDataStore_* 
% note the "Store" in the name of these variables. these variables are used to move results out of the outer parfor loop

                    IterEM_smpfracStorenorm   = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_norm   = NaN(s, numberiterations);       %  norm
                    IterEM_smpfracStoreSO4    = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_SO4    = NaN(s, numberiterations);       %  SO4 
                    IterEM_smpfracStored34S   = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_d34S   = NaN(s, numberiterations);       %  d34S
if go.ALK    == 1;  IterEM_smpfracStoreALK    = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_ALK    = NaN(s, numberiterations);  end  %  ALK
if go.DIC    == 1;  IterEM_smpfracStoreDIC    = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_DIC    = NaN(s, numberiterations);  end  %  DIC
if go.Ca     == 1;  IterEM_smpfracStoreCa     = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_Ca     = NaN(s, numberiterations);  end  %  Ca
if go.Mg     == 1;  IterEM_smpfracStoreMg     = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_Mg     = NaN(s, numberiterations);  end  %  Mg
if go.Na     == 1;  IterEM_smpfracStoreNa     = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_Na     = NaN(s, numberiterations);  end  %  Na
if go.K      == 1;  IterEM_smpfracStoreK      = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_K      = NaN(s, numberiterations);  end  %  K
if go.Sr     == 1;  IterEM_smpfracStoreSr     = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_Sr     = NaN(s, numberiterations);  end  %  Sr
if go.Fe     == 1;  IterEM_smpfracStoreFe     = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_Fe     = NaN(s, numberiterations);  end  %  Fe
if go.Cl     == 1;  IterEM_smpfracStoreCl     = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_Cl     = NaN(s, numberiterations);  end  %  Cl
if go.NO3    == 1;  IterEM_smpfracStoreNO3    = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_NO3    = NaN(s, numberiterations);  end  %  NO3
if go.PO4    == 1;  IterEM_smpfracStorePO4    = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_PO4    = NaN(s, numberiterations);  end  %  PO4
if go.Si     == 1;  IterEM_smpfracStoreSi     = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_Si     = NaN(s, numberiterations);  end  %  Si
if go.Ge     == 1;  IterEM_smpfracStoreGe     = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_Ge     = NaN(s, numberiterations);  end  %  Ge
if go.Li     == 1;  IterEM_smpfracStoreLi     = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_Li     = NaN(s, numberiterations);  end  %  Li
if go.F      == 1;  IterEM_smpfracStoreF      = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_F      = NaN(s, numberiterations);  end  %  F
if go.B      == 1;  IterEM_smpfracStoreB      = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_B      = NaN(s, numberiterations);  end  %  B
if go.Re     == 1;  IterEM_smpfracStoreRe     = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_Re     = NaN(s, numberiterations);  end  %  Re
if go.Mo     == 1;  IterEM_smpfracStoreMo     = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_Mo     = NaN(s, numberiterations);  end  %  Mo
if go.Os     == 1;  IterEM_smpfracStoreOs     = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_Os     = NaN(s, numberiterations);  end  %  Os
if go.HCO3   == 1;  IterEM_smpfracStoreHCO3   = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_HCO3   = NaN(s, numberiterations);  end  %  HCO3
if go.d7Li   == 1;  IterEM_smpfracStored7Li   = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_d7Li   = NaN(s, numberiterations);  end  %  d17Li
if go.d13C   == 1;  IterEM_smpfracStored13C   = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_d13C   = NaN(s, numberiterations);  end  %  d13C
if go.d18O   == 1;  IterEM_smpfracStored18O   = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_d18O   = NaN(s, numberiterations);  end  %  d18O
if go.d26Mg  == 1;  IterEM_smpfracStored26Mg  = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_d26Mg  = NaN(s, numberiterations);  end  %  d26Mg
if go.d30Si  == 1;  IterEM_smpfracStored30Si  = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_d30Si  = NaN(s, numberiterations);  end  %  d30Si
if go.d42Ca  == 1;  IterEM_smpfracStored42Ca  = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_d42Ca  = NaN(s, numberiterations);  end  %  d42Ca
if go.d44Ca  == 1;  IterEM_smpfracStored44Ca  = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_d44Ca  = NaN(s, numberiterations);  end  %  d44Ca
if go.d56Fe  == 1;  IterEM_smpfracStored56Fe  = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_d56Fe  = NaN(s, numberiterations);  end  %  d56Fe
if go.Sr8786 == 1;  IterEM_smpfracStoreSr8786 = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_Sr8786 = NaN(s, numberiterations);  end  %  Sr8786
if go.d98Mo  == 1;  IterEM_smpfracStored98Mo  = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_d98Mo  = NaN(s, numberiterations);  end  %  d98Mo
if go.Os8788 == 1;  IterEM_smpfracStoreOs8788 = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_Os8788 = NaN(s, numberiterations);  end  %  Os8788
if go.Fmod   == 1;  IterEM_smpfracStoreFmod   = NaN(nEM, numberiterations, s);  IterEM_InvRivDataStore_Fmod   = NaN(s, numberiterations);  end  %  Fmod