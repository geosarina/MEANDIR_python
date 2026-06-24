function [InvRivDat, RiverColumn0] = MEANDIR_AdjustRiverDataForAnalyticalError(i,river,ObsList,RiverColumn0,AdjustRiverObs,InvRivDat, ObsInNormalization, carbonisotopematch, ConvertDelta2RList, ConvertDelta2R, go)

        % this function updates the river observations stored in RiverColumn0 to account for analytical error during determination of major ion concentration
        % if AdjustRiverObs = 0, the output RiverColumn0 should equal the input RiverColumn0. note that RiverColumn0 contains information on the ratios 
        % of observations, while InvRivDat contains information on the observations themselves

        % (1) unpack the InvRivDat structure   
                           InvRivDatnorm   = InvRivDat.norm;       % normalization
                           InvRivDatSO4    = InvRivDat.SO4;        % SO4
                           InvRivDatd34S   = InvRivDat.d34S;       % d34S
        if go.ALK    == 1; InvRivDatALK    = InvRivDat.ALK;    end % ALK
        if go.DIC    == 1; InvRivDatDIC    = InvRivDat.DIC;    end % DIC
        if go.Ca     == 1; InvRivDatCa     = InvRivDat.Ca;     end % Ca
        if go.Mg     == 1; InvRivDatMg     = InvRivDat.Mg;     end % Mg
        if go.Na     == 1; InvRivDatNa     = InvRivDat.Na;     end % Na
        if go.K      == 1; InvRivDatK      = InvRivDat.K;      end % K
        if go.Sr     == 1; InvRivDatSr     = InvRivDat.Sr;     end % Sr
        if go.Fe     == 1; InvRivDatFe     = InvRivDat.Fe;     end % Fe
        if go.Cl     == 1; InvRivDatCl     = InvRivDat.Cl;     end % Cl
        if go.NO3    == 1; InvRivDatNO3    = InvRivDat.NO3;    end % NO3
        if go.PO4    == 1; InvRivDatPO4    = InvRivDat.PO4;    end % PO4
        if go.Si     == 1; InvRivDatSi     = InvRivDat.Si;     end % Si
        if go.Ge     == 1; InvRivDatGe     = InvRivDat.Ge;     end % Ge
        if go.Li     == 1; InvRivDatLi     = InvRivDat.Li;     end % Li
        if go.F      == 1; InvRivDatF      = InvRivDat.F;      end % F
        if go.B      == 1; InvRivDatB      = InvRivDat.B;      end % B
        if go.Re     == 1; InvRivDatRe     = InvRivDat.Re;     end % Re
        if go.Mo     == 1; InvRivDatMo     = InvRivDat.Mo;     end % Mo
        if go.Os     == 1; InvRivDatOs     = InvRivDat.Os;     end % Os
        if go.HCO3   == 1; InvRivDatHCO3   = InvRivDat.HCO3;   end % HCO3
        if go.d7Li   == 1; InvRivDatd7Li   = InvRivDat.d7Li;   end % d7Li
        if go.d13C   == 1; InvRivDatd13C   = InvRivDat.d13C;   end % d13C
        if go.d18O   == 1; InvRivDatd18O   = InvRivDat.d18O;   end % d18O
        if go.d26Mg  == 1; InvRivDatd26Mg  = InvRivDat.d26Mg;  end % d26Mg
        if go.d30Si  == 1; InvRivDatd30Si  = InvRivDat.d30Si;  end % d30Si
        if go.d42Ca  == 1; InvRivDatd42Ca  = InvRivDat.d42Ca;  end % d42Ca
        if go.d44Ca  == 1; InvRivDatd44Ca  = InvRivDat.d44Ca;  end % d44Ca
        if go.d56Fe  == 1; InvRivDatd56Fe  = InvRivDat.d56Fe;  end % d56Fe
        if go.Sr8786 == 1; InvRivDatSr8786 = InvRivDat.Sr8786; end % Sr8786
        if go.d98Mo  == 1; InvRivDatd98Mo  = InvRivDat.d98Mo;  end % d98Mo
        if go.Os8788 == 1; InvRivDatOs8788 = InvRivDat.Os8788; end % Os8788
        if go.Fmod   == 1; InvRivDatFmod   = InvRivDat.Fmod;   end % Fmod

        % (2) if accounting for measurement error, shift the data
        if     AdjustRiverObs==1
              if sum(ismember(ObsList,'ALK'))   >0; InvRivDatALK    = river.model_variable.ALK(i)    + river.model_variable.ALK_asd(i)    .*randn; end % ALK
              if sum(ismember(ObsList,'DIC'))   >0; InvRivDatDIC    = river.model_variable.DIC(i)    + river.model_variable.DIC_asd(i)    .*randn; end % DIC
              if sum(ismember(ObsList,'Ca'))    >0; InvRivDatCa     = river.model_variable.Ca(i)     + river.model_variable.Ca_asd(i)     .*randn; end % Ca
              if sum(ismember(ObsList,'Mg'))    >0; InvRivDatMg     = river.model_variable.Mg(i)     + river.model_variable.Mg_asd(i)     .*randn; end % Mg
              if sum(ismember(ObsList,'Na'))    >0; InvRivDatNa     = river.model_variable.Na(i)     + river.model_variable.Na_asd(i)     .*randn; end % Na
              if sum(ismember(ObsList,'K'))     >0; InvRivDatK      = river.model_variable.K(i)      + river.model_variable.K_asd(i)      .*randn; end % K
              if sum(ismember(ObsList,'Sr'))    >0; InvRivDatSr     = river.model_variable.Sr(i)     + river.model_variable.Sr_asd(i)     .*randn; end % Sr
              if sum(ismember(ObsList,'Fe'))    >0; InvRivDatFe     = river.model_variable.Fe(i)     + river.model_variable.Fe_asd(i)     .*randn; end % Fe
              if sum(ismember(ObsList,'Cl'))    >0; InvRivDatCl     = river.model_variable.Cl(i)     + river.model_variable.Cl_asd(i)     .*randn; end % Cl
              if sum(ismember(ObsList,'SO4'))   >0; InvRivDatSO4    = river.model_variable.SO4(i)    + river.model_variable.SO4_asd(i)    .*randn; end % SO4
              if sum(ismember(ObsList,'NO3'))   >0; InvRivDatNO3    = river.model_variable.NO3(i)    + river.model_variable.NO3_asd(i)    .*randn; end % NO3
              if sum(ismember(ObsList,'PO4'))   >0; InvRivDatPO4    = river.model_variable.PO4(i)    + river.model_variable.PO4_asd(i)    .*randn; end % PO4
              if sum(ismember(ObsList,'Si'))    >0; InvRivDatSi     = river.model_variable.Si(i)     + river.model_variable.Si_asd(i)     .*randn; end % Si
              if sum(ismember(ObsList,'Ge'))    >0; InvRivDatGe     = river.model_variable.Ge(i)     + river.model_variable.Ge_asd(i)     .*randn; end % Ge
              if sum(ismember(ObsList,'Li'))    >0; InvRivDatLi     = river.model_variable.Li(i)     + river.model_variable.Li_asd(i)     .*randn; end % Li
              if sum(ismember(ObsList,'F'))     >0; InvRivDatF      = river.model_variable.F(i)      + river.model_variable.F_asd(i)      .*randn; end % F
              if sum(ismember(ObsList,'B'))     >0; InvRivDatB      = river.model_variable.B(i)      + river.model_variable.B_asd(i)      .*randn; end % B
              if sum(ismember(ObsList,'Re'))    >0; InvRivDatRe     = river.model_variable.Re(i)     + river.model_variable.Re_asd(i)     .*randn; end % Re
              if sum(ismember(ObsList,'Mo'))    >0; InvRivDatMo     = river.model_variable.Mo(i)     + river.model_variable.Mo_asd(i)     .*randn; end % Mo
              if sum(ismember(ObsList,'Os'))    >0; InvRivDatOs     = river.model_variable.Os(i)     + river.model_variable.Os_asd(i)     .*randn; end % Os
              if sum(ismember(ObsList,'HCO3'))  >0; InvRivDatHCO3   = river.model_variable.HCO3(i)   + river.model_variable.HCO3_asd(i)   .*randn; end % HCO3
              if sum(ismember(ObsList,'d7Li'))  >0; InvRivDatd7Li   = river.model_variable.d7Li(i)   + river.model_variable.d7Li_asd(i)   .*randn; end % d7Li
              if sum(ismember(ObsList,'d13C'))  >0; InvRivDatd13C   = river.model_variable.d13C(i)   + river.model_variable.d13C_asd(i)   .*randn; end % d13C
              if sum(ismember(ObsList,'d18O'))  >0; InvRivDatd18O   = river.model_variable.d18O(i)   + river.model_variable.d18O_asd(i)   .*randn; end % d18O
              if sum(ismember(ObsList,'d26Mg')) >0; InvRivDatd26Mg  = river.model_variable.d26Mg(i)  + river.model_variable.d26Mg_asd(i)  .*randn; end % d26Mg
              if sum(ismember(ObsList,'d30Si')) >0; InvRivDatd30Si  = river.model_variable.d30Si(i)  + river.model_variable.d30Si_asd(i)  .*randn; end % d30Si
              if sum(ismember(ObsList,'d34S'))  >0; InvRivDatd34S   = river.model_variable.d34S(i)   + river.model_variable.d34S_asd(i)   .*randn; end % d34S
              if sum(ismember(ObsList,'d42Ca')) >0; InvRivDatd42Ca  = river.model_variable.d42Ca(i)  + river.model_variable.d42Ca_asd(i)  .*randn; end % d42Ca
              if sum(ismember(ObsList,'d44Ca')) >0; InvRivDatd44Ca  = river.model_variable.d44Ca(i)  + river.model_variable.d44Ca_asd(i)  .*randn; end % d44Ca
              if sum(ismember(ObsList,'d56Fe')) >0; InvRivDatd56Fe  = river.model_variable.d56Fe(i)  + river.model_variable.d56Fe_asd(i)  .*randn; end % d56Fe
              if sum(ismember(ObsList,'Sr8786'))>0; InvRivDatSr8786 = river.model_variable.Sr8786(i) + river.model_variable.Sr8786_asd(i) .*randn; end % Sr8786
              if sum(ismember(ObsList,'d98Mo')) >0; InvRivDatd98Mo  = river.model_variable.d98Mo(i)  + river.model_variable.d98Mo_asd(i)  .*randn; end % d98Mo
              if sum(ismember(ObsList,'Os8788'))>0; InvRivDatOs8788 = river.model_variable.Os8788(i) + river.model_variable.Os8788_asd(i) .*randn; end % Os8788
              if sum(ismember(ObsList,'Fmod'))  >0; InvRivDatFmod   = river.model_variable.Fmod(i)   + river.model_variable.Fmod_asd(i)   .*randn; end % Fmod

              % generate the new, shifted normalization variable
              InvRivDatnorm = 0;
              for i=1:length(ObsInNormalization)  
                  InvRivDatnorm = InvRivDatnorm + eval(sprintf('InvRivDat%s',ObsInNormalization{i}));
              end

        % (3) do not shift the data, if they should not account for measurement error. 
        % rather, in this case just re-define the variables of InvRivDat*
        elseif AdjustRiverObs==0
                                                    InvRivDatnorm       = river.model_variable.norm(i);       % norm
              if sum(ismember(ObsList,'ALK'))   >0; InvRivDatALK        = river.model_variable.ALK(i);    end % ALK
              if sum(ismember(ObsList,'DIC'))   >0; InvRivDatDIC        = river.model_variable.DIC(i);    end % DIC
              if sum(ismember(ObsList,'Ca'))    >0; InvRivDatCa         = river.model_variable.Ca(i);     end % Ca
              if sum(ismember(ObsList,'Mg'))    >0; InvRivDatMg         = river.model_variable.Mg(i);     end % Mg
              if sum(ismember(ObsList,'Na'))    >0; InvRivDatNa         = river.model_variable.Na(i);     end % Na
              if sum(ismember(ObsList,'K'))     >0; InvRivDatK          = river.model_variable.K(i);      end % K
              if sum(ismember(ObsList,'Sr'))    >0; InvRivDatSr         = river.model_variable.Sr(i);     end % Sr
              if sum(ismember(ObsList,'Fe'))    >0; InvRivDatFe         = river.model_variable.Fe(i);     end % Fe
              if sum(ismember(ObsList,'Cl'))    >0; InvRivDatCl         = river.model_variable.Cl(i);     end % Cl
              if sum(ismember(ObsList,'SO4'))   >0; InvRivDatSO4        = river.model_variable.SO4(i);    end % SO4
              if sum(ismember(ObsList,'NO3'))   >0; InvRivDatNO3        = river.model_variable.NO3(i);    end % NO3
              if sum(ismember(ObsList,'PO4'))   >0; InvRivDatPO4        = river.model_variable.PO4(i);    end % PO4
              if sum(ismember(ObsList,'Si'))    >0; InvRivDatSi         = river.model_variable.Si(i);     end % Si
              if sum(ismember(ObsList,'Ge'))    >0; InvRivDatGe         = river.model_variable.Ge(i);     end % Ge
              if sum(ismember(ObsList,'Li'))    >0; InvRivDatLi         = river.model_variable.Li(i);     end % Li
              if sum(ismember(ObsList,'F'))     >0; InvRivDatF          = river.model_variable.F(i);      end % F
              if sum(ismember(ObsList,'B'))     >0; InvRivDatB          = river.model_variable.B(i);      end % B
              if sum(ismember(ObsList,'Re'))    >0; InvRivDatRe         = river.model_variable.Re(i);     end % Re
              if sum(ismember(ObsList,'Mo'))    >0; InvRivDatMo         = river.model_variable.Mo(i);     end % Mo
              if sum(ismember(ObsList,'Os'))    >0; InvRivDatOs         = river.model_variable.Os(i);     end % Os
              if sum(ismember(ObsList,'HCO3'))  >0; InvRivDatHCO3       = river.model_variable.HCO3(i);   end % HCO3
              if sum(ismember(ObsList,'d7Li'))  >0; InvRivDatd7Li       = river.model_variable.d7Li(i);   end % d7Li
              if sum(ismember(ObsList,'d13C'))  >0; InvRivDatd13C       = river.model_variable.d13C(i);   end % d13C
              if sum(ismember(ObsList,'d18O'))  >0; InvRivDatd18O       = river.model_variable.d18O(i);   end % d18O
              if sum(ismember(ObsList,'d26Mg')) >0; InvRivDatd26Mg      = river.model_variable.d26Mg(i);  end % d26Mg
              if sum(ismember(ObsList,'d30Si')) >0; InvRivDatd30Si      = river.model_variable.d30Si(i);  end % d30Si
              if sum(ismember(ObsList,'d34S'))  >0; InvRivDatd34S       = river.model_variable.d34S(i);   end % d34S
              if sum(ismember(ObsList,'d42Ca')) >0; InvRivDatd42Ca      = river.model_variable.d42Ca(i);  end % d42Ca
              if sum(ismember(ObsList,'d44Ca')) >0; InvRivDatd44Ca      = river.model_variable.d44Ca(i);  end % d44Ca
              if sum(ismember(ObsList,'d56Fe')) >0; InvRivDatd56Fe      = river.model_variable.d56Fe(i);  end % d56Fe
              if sum(ismember(ObsList,'Sr8786'))>0; InvRivDatSr8786     = river.model_variable.Sr8786(i); end % Sr8786
              if sum(ismember(ObsList,'d98Mo')) >0; InvRivDatd98Mo      = river.model_variable.d98Mo(i);  end % d98Mo
              if sum(ismember(ObsList,'Os8788'))>0; InvRivDatOs8788     = river.model_variable.Os8788(i); end % Os8788
              if sum(ismember(ObsList,'Fmod'))  >0; InvRivDatFmod       = river.model_variable.Fmod(i);   end % Fmod
        end

        % (4) if delta notation should be converted to isotopic ratios
        if sum(ismember(ConvertDelta2RList,'d7Li')) >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d7Li'));  InvRivDatd7Li  = (InvRivDatd7Li/1000+1)*conv;  end % d7Li
        if sum(ismember(ConvertDelta2RList,'d13C')) >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d13C'));  InvRivDatd13C  = (InvRivDatd13C/1000+1)*conv;  end % d13C
        if sum(ismember(ConvertDelta2RList,'d18O')) >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d18O'));  InvRivDatd18O  = (InvRivDatd18O/1000+1)*conv;  end % d17O
        if sum(ismember(ConvertDelta2RList,'d26Mg'))>0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d26Mg')); InvRivDatd26Mg = (InvRivDatd26Mg/1000+1)*conv; end % d26Mg
        if sum(ismember(ConvertDelta2RList,'d30Si'))>0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d30Si')); InvRivDatd30Si = (InvRivDatd30Si/1000+1)*conv; end % d30Si
        if sum(ismember(ConvertDelta2RList,'d34S')) >0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d34S'));  InvRivDatd34S  = (InvRivDatd34S/1000+1)*conv;  end % d34S
        if sum(ismember(ConvertDelta2RList,'d42Ca'))>0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d42Ca')); InvRivDatd42Ca = (InvRivDatd42Ca/1000+1)*conv; end % d42Ca
        if sum(ismember(ConvertDelta2RList,'d44Ca'))>0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d44Ca')); InvRivDatd44Ca = (InvRivDatd44Ca/1000+1)*conv; end % d44Ca
        if sum(ismember(ConvertDelta2RList,'d56Fe'))>0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d56Fe')); InvRivDatd56Fe = (InvRivDatd56Fe/1000+1)*conv; end % d56Fe
        if sum(ismember(ConvertDelta2RList,'d98Mo'))>0; conv = ConvertDelta2R.ratio(ismember(ConvertDelta2R.isotope,'d98Mo')); InvRivDatd98Mo = (InvRivDatd98Mo/1000+1)*conv; end % d98Mo
    
        % (5) save the data within a structure
                           InvRivDat.norm       = InvRivDatnorm;       % norm
                           InvRivDat.SO4        = InvRivDatSO4;        % SO4
                           InvRivDat.d34S       = InvRivDatd34S;       % d34S
        if go.ALK    == 1; InvRivDat.ALK        = InvRivDatALK;    end % ALK
        if go.DIC    == 1; InvRivDat.DIC        = InvRivDatDIC;    end % DIC
        if go.Ca     == 1; InvRivDat.Ca         = InvRivDatCa;     end % Ca
        if go.Mg     == 1; InvRivDat.Mg         = InvRivDatMg;     end % Mg
        if go.Na     == 1; InvRivDat.Na         = InvRivDatNa;     end % Na
        if go.K      == 1; InvRivDat.K          = InvRivDatK;      end % K
        if go.Sr     == 1; InvRivDat.Sr         = InvRivDatSr;     end % Sr
        if go.Fe     == 1; InvRivDat.Fe         = InvRivDatFe;     end % Fe
        if go.Cl     == 1; InvRivDat.Cl         = InvRivDatCl;     end % Cl
        if go.NO3    == 1; InvRivDat.NO3        = InvRivDatNO3;    end % NO3
        if go.PO4    == 1; InvRivDat.PO4        = InvRivDatPO4;    end % PO4
        if go.Si     == 1; InvRivDat.Si         = InvRivDatSi;     end % Si
        if go.Ge     == 1; InvRivDat.Ge         = InvRivDatGe;     end % Ge
        if go.Li     == 1; InvRivDat.Li         = InvRivDatLi;     end % Li
        if go.F      == 1; InvRivDat.F          = InvRivDatF;      end % F
        if go.B      == 1; InvRivDat.B          = InvRivDatB;      end % B
        if go.Re     == 1; InvRivDat.Re         = InvRivDatRe;     end % Re
        if go.Mo     == 1; InvRivDat.Mo         = InvRivDatMo;     end % Mo
        if go.Os     == 1; InvRivDat.Os         = InvRivDatOs;     end % Os
        if go.HCO3   == 1; InvRivDat.HCO3       = InvRivDatHCO3;   end % HCO3
        if go.d7Li   == 1; InvRivDat.d7Li       = InvRivDatd7Li;   end % d7Li
        if go.d13C   == 1; InvRivDat.d13C       = InvRivDatd13C;   end % d13C
        if go.d18O   == 1; InvRivDat.d18O       = InvRivDatd18O;   end % d18O
        if go.d26Mg  == 1; InvRivDat.d26Mg      = InvRivDatd26Mg;  end % d26Mg
        if go.d30Si  == 1; InvRivDat.d30Si      = InvRivDatd30Si;  end % d30Si
        if go.d42Ca  == 1; InvRivDat.d42Ca      = InvRivDatd42Ca;  end % d42Ca
        if go.d44Ca  == 1; InvRivDat.d44Ca      = InvRivDatd44Ca;  end % d44Ca
        if go.d56Fe  == 1; InvRivDat.d56Fe      = InvRivDatd56Fe;  end % d56Fe
        if go.Sr8786 == 1; InvRivDat.Sr8786     = InvRivDatSr8786; end % Sr8786
        if go.d98Mo  == 1; InvRivDat.d98Mo      = InvRivDatd98Mo;  end % d98Mo
        if go.Os8788 == 1; InvRivDat.Os8788     = InvRivDatOs8788; end % Os8788
        if go.Fmod   == 1; InvRivDat.Fmod       = InvRivDatFmod;   end % Fmod

        % (6) re-define the matrix of normalized river data, having now accounted for analytical error
        if sum(ismember(ObsList,'ALK'))   >0;                                       RiverColumn0(ismember(ObsList,'ALK'),1)    = (InvRivDatALK./  InvRivDatnorm);                    end  %  ALK
        if sum(ismember(ObsList,'DIC'))   >0;                                       RiverColumn0(ismember(ObsList,'DIC'),1)    = (InvRivDatDIC./  InvRivDatnorm);                    end  %  DIC
        if sum(ismember(ObsList,'Ca'))    >0;                                       RiverColumn0(ismember(ObsList,'Ca'),1)     = (InvRivDatCa./   InvRivDatnorm);                    end  %  Ca
        if sum(ismember(ObsList,'Mg'))    >0;                                       RiverColumn0(ismember(ObsList,'Mg'),1)     = (InvRivDatMg./   InvRivDatnorm);                    end  %  Mg
        if sum(ismember(ObsList,'Na'))    >0;                                       RiverColumn0(ismember(ObsList,'Na'),1)     = (InvRivDatNa./   InvRivDatnorm);                    end  %  Na
        if sum(ismember(ObsList,'K'))     >0;                                       RiverColumn0(ismember(ObsList,'K'),1)      = (InvRivDatK./    InvRivDatnorm);                    end  %  K
        if sum(ismember(ObsList,'Sr'))    >0;                                       RiverColumn0(ismember(ObsList,'Sr'),1)     = (InvRivDatSr./   InvRivDatnorm);                    end  %  Sr
        if sum(ismember(ObsList,'Fe'))    >0;                                       RiverColumn0(ismember(ObsList,'Fe'),1)     = (InvRivDatFe./   InvRivDatnorm);                    end  %  Fe
        if sum(ismember(ObsList,'Cl'))    >0;                                       RiverColumn0(ismember(ObsList,'Cl'),1)     = (InvRivDatCl./   InvRivDatnorm);                    end  %  Cl
        if sum(ismember(ObsList,'SO4'))   >0;                                       RiverColumn0(ismember(ObsList,'SO4'),1)    = (InvRivDatSO4./  InvRivDatnorm);                    end  %  SO4
        if sum(ismember(ObsList,'NO3'))   >0;                                       RiverColumn0(ismember(ObsList,'NO3'),1)    = (InvRivDatNO3./  InvRivDatnorm);                    end  %  NO3
        if sum(ismember(ObsList,'PO4'))   >0;                                       RiverColumn0(ismember(ObsList,'PO4'),1)    = (InvRivDatPO4./  InvRivDatnorm);                    end  %  PO4
        if sum(ismember(ObsList,'Si'))    >0;                                       RiverColumn0(ismember(ObsList,'Si'),1)     = (InvRivDatSi./   InvRivDatnorm);                    end  %  Si
        if sum(ismember(ObsList,'Ge'))    >0;                                       RiverColumn0(ismember(ObsList,'Ge'),1)     = (InvRivDatGe./   InvRivDatnorm);                    end  %  Ge
        if sum(ismember(ObsList,'Li'))    >0;                                       RiverColumn0(ismember(ObsList,'Li'),1)     = (InvRivDatLi./   InvRivDatnorm);                    end  %  Li
        if sum(ismember(ObsList,'F'))     >0;                                       RiverColumn0(ismember(ObsList,'F'),1)      = (InvRivDatF./    InvRivDatnorm);                    end  %  F
        if sum(ismember(ObsList,'B'))     >0;                                       RiverColumn0(ismember(ObsList,'B'),1)      = (InvRivDatB./    InvRivDatnorm);                    end  %  B
        if sum(ismember(ObsList,'Re'))    >0;                                       RiverColumn0(ismember(ObsList,'Re'),1)     = (InvRivDatRe./   InvRivDatnorm);                    end  %  Re
        if sum(ismember(ObsList,'Mo'))    >0;                                       RiverColumn0(ismember(ObsList,'Mo'),1)     = (InvRivDatMo./   InvRivDatnorm);                    end  %  Mo
        if sum(ismember(ObsList,'Os'))    >0;                                       RiverColumn0(ismember(ObsList,'Os'),1)     = (InvRivDatOs./   InvRivDatnorm);                    end  %  Os
        if sum(ismember(ObsList,'HCO3'))  >0;                                       RiverColumn0(ismember(ObsList,'HCO3'),1)   = (InvRivDatHCO3./ InvRivDatnorm);                    end  %  HCO3
        if sum(ismember(ObsList,'d7Li'))  >0;                                       RiverColumn0(ismember(ObsList,'d7Li'),1)   = (InvRivDatd7Li)  .*(InvRivDatLi./  InvRivDatnorm);  end  %  d7Li
        if sum(ismember(ObsList,'d18O'))  >0;                                       RiverColumn0(ismember(ObsList,'d18O'),1)   = (InvRivDatd18O)  .*(InvRivDatSO4./ InvRivDatnorm);  end  %  d18O
        if sum(ismember(ObsList,'d26Mg')) >0;                                       RiverColumn0(ismember(ObsList,'d26Mg'),1)  = (InvRivDatd26Mg) .*(InvRivDatMg./  InvRivDatnorm);  end  %  d26Mg
        if sum(ismember(ObsList,'d30Si')) >0;                                       RiverColumn0(ismember(ObsList,'d30Si'),1)  = (InvRivDatd30Si) .*(InvRivDatSi./  InvRivDatnorm);  end  %  d30Si
        if sum(ismember(ObsList,'d34S'))  >0;                                       RiverColumn0(ismember(ObsList,'d34S'),1)   = (InvRivDatd34S)  .*(InvRivDatSO4./ InvRivDatnorm);  end  %  d34S
        if sum(ismember(ObsList,'d42Ca')) >0;                                       RiverColumn0(ismember(ObsList,'d42Ca'),1)  = (InvRivDatd42Ca) .*(InvRivDatCa./  InvRivDatnorm);  end  %  d42Ca
        if sum(ismember(ObsList,'d44Ca')) >0;                                       RiverColumn0(ismember(ObsList,'d44Ca'),1)  = (InvRivDatd44Ca) .*(InvRivDatCa./  InvRivDatnorm);  end  %  d44Ca
        if sum(ismember(ObsList,'d56Fe')) >0;                                       RiverColumn0(ismember(ObsList,'d56Fe'),1)  = (InvRivDatd56Fe) .*(InvRivDatFe./  InvRivDatnorm);  end  %  d56Fe
        if sum(ismember(ObsList,'Sr8786'))>0;                                       RiverColumn0(ismember(ObsList,'Sr8786'),1) = (InvRivDatSr8786).*(InvRivDatSr./  InvRivDatnorm);  end  %  Sr8786
        if sum(ismember(ObsList,'d98Mo')) >0;                                       RiverColumn0(ismember(ObsList,'d98Mo'),1)  = (InvRivDatd98Mo) .*(InvRivDatMo./  InvRivDatnorm);  end  %  d98MO
        if sum(ismember(ObsList,'Os8788'))>0;                                       RiverColumn0(ismember(ObsList,'Os8788'),1) = (InvRivDatOs8788).*(InvRivDatOs./  InvRivDatnorm);  end  %  Os8788
        if sum(ismember(ObsList,'d13C'))  >0 & isequal(carbonisotopematch,'HCO3');  RiverColumn0(ismember(ObsList,'d13C'),1)   = (InvRivDatd13C)  .*(InvRivDatHCO3./InvRivDatnorm);  end  %  d13C for HCO3
        if sum(ismember(ObsList,'Fmod'))  >0 & isequal(carbonisotopematch,'HCO3');  RiverColumn0(ismember(ObsList,'Fmod'),1)   = (InvRivDatFmod)  .*(InvRivDatHCO3./InvRivDatnorm);  end  %  Fmod for HCO3
        if sum(ismember(ObsList,'d13C'))  >0 & isequal(carbonisotopematch,'DIC');   RiverColumn0(ismember(ObsList,'d13C'),1)   = (InvRivDatd13C)  .*(InvRivDatDIC./ InvRivDatnorm);  end  %  d13C for DIC
        if sum(ismember(ObsList,'Fmod'))  >0 & isequal(carbonisotopematch,'DIC');   RiverColumn0(ismember(ObsList,'Fmod'),1)   = (InvRivDatFmod)  .*(InvRivDatDIC./ InvRivDatnorm);  end  %  Fmod for DIC

end % end of the function