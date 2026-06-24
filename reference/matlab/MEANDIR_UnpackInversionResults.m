%% MEANDIR_UnpackInversionResults
% for each sample, unpack cells storing data and organize the results of the inversion.

%% (1) declare results structures
results = struct;
if go.ALK    == 1; results.ALK                 = NaN(nEM,resultssize,s); end % ALK
if go.DIC    == 1; results.DIC                 = NaN(nEM,resultssize,s); end % DIC
if go.Ca     == 1; results.Ca                  = NaN(nEM,resultssize,s); end % Ca
if go.Mg     == 1; results.Mg                  = NaN(nEM,resultssize,s); end % Mg
if go.Na     == 1; results.Na                  = NaN(nEM,resultssize,s); end % Na
if go.K      == 1; results.K                   = NaN(nEM,resultssize,s); end % K
if go.Sr     == 1; results.Sr                  = NaN(nEM,resultssize,s); end % Sr
if go.Fe     == 1; results.Fe                  = NaN(nEM,resultssize,s); end % Fe
if go.Cl     == 1; results.Cl                  = NaN(nEM,resultssize,s); end % Cl
if go.NO3    == 1; results.NO3                 = NaN(nEM,resultssize,s); end % NO3
if go.PO4    == 1; results.PO4                 = NaN(nEM,resultssize,s); end % PO4
if go.Si     == 1; results.Si                  = NaN(nEM,resultssize,s); end % Si
if go.Ge     == 1; results.Ge                  = NaN(nEM,resultssize,s); end % Ge
if go.Li     == 1; results.Li                  = NaN(nEM,resultssize,s); end % Li
if go.F      == 1; results.F                   = NaN(nEM,resultssize,s); end % F
if go.B      == 1; results.B                   = NaN(nEM,resultssize,s); end % B
if go.Re     == 1; results.Re                  = NaN(nEM,resultssize,s); end % Re
if go.Mo     == 1; results.Mo                  = NaN(nEM,resultssize,s); end % Mo
if go.Os     == 1; results.Os                  = NaN(nEM,resultssize,s); end % Os
if go.HCO3   == 1; results.HCO3                = NaN(nEM,resultssize,s); end % HCO3
if go.d7Li   == 1; results.d7Li                = NaN(nEM,resultssize,s); end % d7Li
if go.d13C   == 1; results.d13C                = NaN(nEM,resultssize,s); end % d13C
if go.d18O   == 1; results.d18O                = NaN(nEM,resultssize,s); end % d18O
if go.d26Mg  == 1; results.d26Mg               = NaN(nEM,resultssize,s); end % d26Mg
if go.d30Si  == 1; results.d30Si               = NaN(nEM,resultssize,s); end % d30Si
if go.d42Ca  == 1; results.d42Ca               = NaN(nEM,resultssize,s); end % d42Ca
if go.d44Ca  == 1; results.d44Ca               = NaN(nEM,resultssize,s); end % d44Ca
if go.d56Fe  == 1; results.d56Fe               = NaN(nEM,resultssize,s); end % d56Fe
if go.Sr8786 == 1; results.Sr8786              = NaN(nEM,resultssize,s); end % Sr8786
if go.d98Mo  == 1; results.d98Mo               = NaN(nEM,resultssize,s); end % d98Mo
if go.Os8788 == 1; results.Os8788              = NaN(nEM,resultssize,s); end % Os8786
if go.Fmod   == 1; results.Fmod                = NaN(nEM,resultssize,s); end % Fmod
                   results.SO4                 = NaN(nEM,resultssize,s);     % SO4
                   results.d34S                = NaN(nEM,resultssize,s);     % d34S
                   results.norm                = NaN(nEM,resultssize,s);     % normalization
                   results.misfit              = NaN(s,resultssize);         % misfit from comparing model results and observations
                   results.misfit_costfunction = NaN(s,resultssize);         % misfit from cost function evaluation

%% (2) unpack the samples
for i=1:s  
    smprestemp    = smprescell{i,1};
    InvRivDattemp = InvRivDatacell{i,1};    
    
    if go.ALK    == 1; results.ALK(:,1:resultssize,i)    = smprestemp.smpresultsALK;    InvRivDataAgg.ALK(i,1:resultssize)    = InvRivDattemp.ALK;    end % ALK
    if go.DIC    == 1; results.DIC(:,1:resultssize,i)    = smprestemp.smpresultsDIC;    InvRivDataAgg.DIC(i,1:resultssize)    = InvRivDattemp.DIC;    end % DIC
    if go.Ca     == 1; results.Ca(:,1:resultssize,i)     = smprestemp.smpresultsCa;     InvRivDataAgg.Ca(i,1:resultssize)     = InvRivDattemp.Ca;     end % Ca
    if go.Mg     == 1; results.Mg(:,1:resultssize,i)     = smprestemp.smpresultsMg;     InvRivDataAgg.Mg(i,1:resultssize)     = InvRivDattemp.Mg;     end % Mg
    if go.Na     == 1; results.Na(:,1:resultssize,i)     = smprestemp.smpresultsNa;     InvRivDataAgg.Na(i,1:resultssize)     = InvRivDattemp.Na;     end % Na
    if go.K      == 1; results.K(:,1:resultssize,i)      = smprestemp.smpresultsK;      InvRivDataAgg.K(i,1:resultssize)      = InvRivDattemp.K;      end % K
    if go.Sr     == 1; results.Sr(:,1:resultssize,i)     = smprestemp.smpresultsSr;     InvRivDataAgg.Sr(i,1:resultssize)     = InvRivDattemp.Sr;     end % Sr
    if go.Fe     == 1; results.Fe(:,1:resultssize,i)     = smprestemp.smpresultsFe;     InvRivDataAgg.Fe(i,1:resultssize)     = InvRivDattemp.Fe;     end % Fe
    if go.Cl     == 1; results.Cl(:,1:resultssize,i)     = smprestemp.smpresultsCl;     InvRivDataAgg.Cl(i,1:resultssize)     = InvRivDattemp.Cl;     end % Cl
    if go.NO3    == 1; results.NO3(:,1:resultssize,i)    = smprestemp.smpresultsNO3;    InvRivDataAgg.NO3(i,1:resultssize)    = InvRivDattemp.NO3;    end % NO3
    if go.PO4    == 1; results.PO4(:,1:resultssize,i)    = smprestemp.smpresultsPO4;    InvRivDataAgg.PO4(i,1:resultssize)    = InvRivDattemp.PO4;    end % PO4
    if go.Si     == 1; results.Si(:,1:resultssize,i)     = smprestemp.smpresultsSi;     InvRivDataAgg.Si(i,1:resultssize)     = InvRivDattemp.Si;     end % Si
    if go.Ge     == 1; results.Ge(:,1:resultssize,i)     = smprestemp.smpresultsGe;     InvRivDataAgg.Ge(i,1:resultssize)     = InvRivDattemp.Ge;     end % Ge
    if go.Li     == 1; results.Li(:,1:resultssize,i)     = smprestemp.smpresultsLi;     InvRivDataAgg.Li(i,1:resultssize)     = InvRivDattemp.Li;     end % Li
    if go.F      == 1; results.F(:,1:resultssize,i)      = smprestemp.smpresultsF;      InvRivDataAgg.F(i,1:resultssize)      = InvRivDattemp.F;      end % F
    if go.B      == 1; results.B(:,1:resultssize,i)      = smprestemp.smpresultsB;      InvRivDataAgg.B(i,1:resultssize)      = InvRivDattemp.B;      end % B
    if go.Re     == 1; results.Re(:,1:resultssize,i)     = smprestemp.smpresultsRe;     InvRivDataAgg.Re(i,1:resultssize)     = InvRivDattemp.Re;     end % Re
    if go.Mo     == 1; results.Mo(:,1:resultssize,i)     = smprestemp.smpresultsMo;     InvRivDataAgg.Mo(i,1:resultssize)     = InvRivDattemp.Mo;     end % Mo
    if go.Os     == 1; results.Os(:,1:resultssize,i)     = smprestemp.smpresultsOs;     InvRivDataAgg.Os(i,1:resultssize)     = InvRivDattemp.Os;     end % Os
    if go.HCO3   == 1; results.HCO3(:,1:resultssize,i)   = smprestemp.smpresultsHCO3;   InvRivDataAgg.HCO3(i,1:resultssize)   = InvRivDattemp.HCO3;   end % HCO3
    if go.d7Li   == 1; results.d7Li(:,1:resultssize,i)   = smprestemp.smpresultsd7Li;   InvRivDataAgg.d7Li(i,1:resultssize)   = InvRivDattemp.d7Li;   end % d7Li
    if go.d13C   == 1; results.d13C(:,1:resultssize,i)   = smprestemp.smpresultsd13C;   InvRivDataAgg.d13C(i,1:resultssize)   = InvRivDattemp.d13C;   end % d13C
    if go.d18O   == 1; results.d18O(:,1:resultssize,i)   = smprestemp.smpresultsd18O;   InvRivDataAgg.d18O(i,1:resultssize)   = InvRivDattemp.d18O;   end % d18O
    if go.d26Mg  == 1; results.d26Mg(:,1:resultssize,i)  = smprestemp.smpresultsd26Mg;  InvRivDataAgg.d26Mg(i,1:resultssize)  = InvRivDattemp.d26Mg;  end % d26Mg
    if go.d30Si  == 1; results.d30Si(:,1:resultssize,i)  = smprestemp.smpresultsd30Si;  InvRivDataAgg.d30Si(i,1:resultssize)  = InvRivDattemp.d30Si;  end % d30Si
    if go.d42Ca  == 1; results.d42Ca(:,1:resultssize,i)  = smprestemp.smpresultsd42Ca;  InvRivDataAgg.d42Ca(i,1:resultssize)  = InvRivDattemp.d42Ca;  end % d42Ca
    if go.d44Ca  == 1; results.d44Ca(:,1:resultssize,i)  = smprestemp.smpresultsd44Ca;  InvRivDataAgg.d44Ca(i,1:resultssize)  = InvRivDattemp.d44Ca;  end % d44Ca
    if go.d56Fe  == 1; results.d56Fe(:,1:resultssize,i)  = smprestemp.smpresultsd56Fe;  InvRivDataAgg.d56Fe(i,1:resultssize)  = InvRivDattemp.d56Fe;  end % d56Fe
    if go.Sr8786 == 1; results.Sr8786(:,1:resultssize,i) = smprestemp.smpresultsSr8786; InvRivDataAgg.Sr8786(i,1:resultssize) = InvRivDattemp.Sr8786; end % Sr8786
    if go.d98Mo  == 1; results.d98Mo(:,1:resultssize,i)  = smprestemp.smpresultsd98Mo;  InvRivDataAgg.d98Mo(i,1:resultssize)  = InvRivDattemp.d98Mo;  end % d98Mo
    if go.Os8788 == 1; results.Os8788(:,1:resultssize,i) = smprestemp.smpresultsOs8788; InvRivDataAgg.Os8788(i,1:resultssize) = InvRivDattemp.Os8788; end % Os8788
    if go.Fmod   == 1; results.Fmod(:,1:resultssize,i)   = smprestemp.smpresultsFmod;   InvRivDataAgg.Fmod(i,1:resultssize)   = InvRivDattemp.Fmod;   end % Fmod
                       results.SO4(:,1:resultssize,i)    = smprestemp.smpresultsSO4;    InvRivDataAgg.SO4(i,1:resultssize)    = InvRivDattemp.SO4;        % SO4
                       results.d34S(:,1:resultssize,i)   = smprestemp.smpresultsd34S;   InvRivDataAgg.d34S(i,1:resultssize)   = InvRivDattemp.d34S;       % d34S
                       results.norm(:,1:resultssize,i)   = smprestemp.smpresultsnorm;   InvRivDataAgg.norm(i,1:resultssize)   = InvRivDattemp.norm;       % norm
end


%% (3) save the other information
for iterindx_unpack = 1:numberiterations
    
    % iterindx ranges from 1 to 1 when IterOver is sample, because numberiterations is 1
    % iterindx from 1 to a substantial number when IterOver is end-members, because numberiterations is not 1
    
    if     isequal(IterateOver,'Samples');     primindx = 1:resultssize;    scndindx = 1;
    elseif isequal(IterateOver,'End-members'); primindx = iterindx_unpack;  scndindx = iterindx_unpack;
    end

    for i=1:s  
        successcounter           = successcountercell{i,iterindx_unpack};
        badmbcounter             = badmbcountercell{i,iterindx_unpack};
        iterationcounter         = iterationcountercell{i,iterindx_unpack};
        EMSucessIter0Matrix      = EMSucessIter0Matrixcell{i,iterindx_unpack};
        EMSucessIter0Matrix_frac = EMSucessIter0Matrixcell_frac{i,iterindx_unpack};
        RiverColumnSucessMatrix  = RiverColumnSucessMatrixcell{i,iterindx_unpack};
        smpresmisfit             = smpresmisfitcell{i,iterindx_unpack};
        smpfunmisfit             = smpfunmisfitcell{i,iterindx_unpack};

        % save the end-member matrix
        EMSucessIterMatrixAll(:,:,primindx,i)      = EMSucessIter0Matrix;
        EMSucessIterMatrixAll_frac(:,:,primindx,i) = EMSucessIter0Matrix_frac;

        % save information on misfit between samples and standards
        results.misfit(i,primindx)                 = smpresmisfit;
        results.misfit_costfunction(i,primindx)    = smpfunmisfit;

        % make aggregate vectors describing simulation results
        successfullsimulations(i,scndindx)         = successcounter;
        badmbcount(i,scndindx)                     = badmbcounter;
        iterationcounts(i,scndindx)                = iterationcounter;
    end
end

%% (4) remove workspace clutter
clear smprescell
clear InvRivDatacell; 
clear successcountercell
clear badmbcountercell
clear iterationcountercell
clear EMSucessIter0Matrixcell
clear EMSucessIter0Matrixcell_frac
clear RiverColumnSucessMatrixcell
clear smpresmisfitcell
clear smpfunmisfitcell