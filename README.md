# TcellCOBRA
MATLAB codes for reconstruction of T cell metabolic models based on the multi-omics data

# New Codes
```MATLAB
% initializing COBRA
initCobraToolbox(false) %don't update the toolbox
changeCobraSolver ('gurobi', 'all');

% Reading curated recon 3 model
modelFileName = 'Recon3DModel_301.mat';
modelDirectory = getDistributedModelFolder(modelFileName); %Look up the folder for the distributed Models.
modelFileName= [modelDirectory filesep modelFileName]; % Get the full path. Necessary to be sure, that the right model is loaded
model = readCbModel(modelFileName);

% Model manipulation

% Puniya et al model
modelPuniyaFileName = 'C:\Users\MSI\cobratoolbox\Puniya2021_TNM1055.xml'
modelPuniya = readSBML(modelPuniyaFileName,1000)
modelPuniyaEffFileName = 'C:\Users\MSI\cobratoolbox\T1M1133.xml'
modelPuniyaEff = readSBML(modelPuniyaEffFileName,1000)

% Objective function from macrophage model

% reading marophage model
modelMacroFileName = 'C:\Users\MSI\cobratoolbox\msb201068-s1.xml'
model_macrophage = readSBML(modelMacroFileName,1000)
model = addReaction(model,'biomass_reaction_Mphage','reactionFormula',char(printRxnFormula(modelPuniya,'biomass_reaction_Mphage')))
modelnew = changeObjective(model,'biomass_reaction_Mphage',1.0)
%modelnew = model

% reading Recon 2.2 model
modelRecon2FileName = 'C:\Users\MSI\cobratoolbox\MODEL1603150001_url.xml'
modelRecon2 = readSBML(modelRecon2FileName,1000)

% Adding new reactions
modelnew = addReaction(modelnew,'PRPNCOAHYDm','reactionFormula',char(printRxnFormula(model_macrophage,'PRPNCOAHYDm')))
modelnew = addReaction(modelnew,'PPCOAOm','reactionFormula',char(printRxnFormula(model_macrophage,'PPCOAOm')))
modelnew = addReaction(modelnew,'ARTPLM1','reactionFormula',char(printRxnFormula(model_macrophage,'ARTPLM1')))
modelnew = addReaction(modelnew,'ARTPLM2','reactionFormula',char(printRxnFormula(model_macrophage,'ARTPLM2')))
modelnew = addReaction(modelnew,'PIt2m','reactionFormula',char(printRxnFormula(model_macrophage,'PIt2m')))
modelnew = addReaction(modelnew,'r0941','reactionFormula',char(printRxnFormula(modelRecon2,'r0941')))
modelnew = addReaction(modelnew,'r1464','reactionFormula',char(printRxnFormula(modelRecon2,'r1464')))
modelnew = addReaction(modelnew,'AIRCr_PRASCS','reactionFormula',char(printRxnFormula(modelRecon2,'AIRCr_PRASCS')))
modelnew = addReaction(modelnew,'r0683','reactionFormula',char(printRxnFormula(modelRecon2,'r0683')))
modelnew = addReaction(modelnew,'RTOT_3','reactionFormula',char(printRxnFormula(modelRecon2,'RTOT_3')))
modelnew = addReaction(modelnew,'ARTFR13','reactionFormula',char(printRxnFormula(modelRecon2,'ARTFR13')))
modelnew = addExchangeRxn(modelnew, 'ca2[e]')
modelnew = addExchangeRxn(modelnew, 'cl[e]')

% Removing duplicates
[modelnewfinal, removedRxn, rxnRelationship] = checkDuplicateRxn(modelnew, 'S', 1, 1)

% Macrophage Exchange Constraints
exRxns = printUptakeBound(modelnewfinal)
modeltest = changeRxnBounds(modelnewfinal,modelnewfinal.rxns(exRxns),0,'l')

modeltest = changeRxnBounds(modeltest,'EX_arg_L[e]',-0.02375,'l');
modeltest = changeRxnBounds(modeltest,'EX_but[e]',-0.0038,'l');
modeltest = changeRxnBounds(modeltest,'EX_glc_D[e]',-0.2718,'l');
modeltest = changeRxnBounds(modeltest,'EX_gln_L[e]',-0.0765,'l');
modeltest = changeRxnBounds(modeltest,'EX_hco3[e]',-10,'l');
modeltest = changeRxnBounds(modeltest,'EX_his_L[e]',-0.01,'l');
modeltest = changeRxnBounds(modeltest,'EX_ile_L[e]',-0.01,'l');
modeltest = changeRxnBounds(modeltest,'EX_leu_L[e]',-0.0362216,'l');
modeltest = changeRxnBounds(modeltest,'EX_lys_L[e]',-0.01,'l');
modeltest = changeRxnBounds(modeltest,'EX_met_L[e]',-0.01,'l');
modeltest = changeRxnBounds(modeltest,'EX_o2[e]',-0.3066,'l');
modeltest = changeRxnBounds(modeltest,'EX_ocdca[e]',-0.1,'l');
modeltest = changeRxnBounds(modeltest,'EX_ocdcea[e]',-0.0192,'l');
modeltest = changeRxnBounds(modeltest,'EX_phe_L[e]',-0.01,'l');
modeltest = changeRxnBounds(modeltest,'EX_pi[e]',-10,'l');
modeltest = changeRxnBounds(modeltest,'EX_thr_L[e]',-0.01,'l');
modeltest = changeRxnBounds(modeltest,'EX_trp_L[e]',-0.01,'l');
modeltest = changeRxnBounds(modeltest,'EX_ttdca[e]',-0.1,'l');
modeltest = changeRxnBounds(modeltest,'EX_val_L[e]',-0.01,'l');
modeltest = changeRxnBounds(modeltest,'EX_ca2[e]',-1000,'l');
modeltest = changeRxnBounds(modeltest,'EX_cl[e]',-1000,'l');
modeltest = changeRxnBounds(modeltest,'EX_h[e]',-1000,'l');
modeltest = changeRxnBounds(modeltest,'EX_h2o[e]',-1000,'l');
modeltest = changeRxnBounds(modeltest,'EX_k[e]',-1000,'l');
modeltest = changeRxnBounds(modeltest,'EX_na1[e]',-1000,'l');
modeltest = changeRxnBounds(modeltest,'EX_so4[e]',-1000,'l');
modeltest = changeRxnBounds(modeltest,'EX_ca2[e]',1000,'u');
modeltest = changeRxnBounds(modeltest,'EX_cl[e]',1000,'u');
modeltest = changeRxnBounds(modeltest,'EX_h[e]',1000,'u');
modeltest = changeRxnBounds(modeltest,'EX_h2o[e]',1000,'u');
modeltest = changeRxnBounds(modeltest,'EX_k[e]',1000,'u');
modeltest = changeRxnBounds(modeltest,'EX_na1[e]',1000,'u');
modeltest = changeRxnBounds(modeltest,'EX_so4[e]',1000,'u');
modeltest = changeRxnBounds(modeltest,'EX_tag_hs[e]',-0.01,'l');

modeltest = changeRxnBounds(modeltest,'ALATA_L',-0.20960496,'l');
modeltest = changeRxnBounds(modeltest,'LDH_L',-35.73303552,'l');
modeltest = changeRxnBounds(modeltest,'LDH_Lm',-35.73303552,'l');
modeltest = changeRxnBounds(modeltest,'ACACT1r',0.7437132,'u');
modeltest = changeRxnBounds(modeltest,'ASPTA',5.86970496,'u');
modeltest = changeRxnBounds(modeltest,'ASPTAm',5.86970496,'u');
modeltest = changeRxnBounds(modeltest,'BDHm',0.00872784,'u');
modeltest = changeRxnBounds(modeltest,'GLUDxm',4.4962,'u');
modeltest = changeRxnBounds(modeltest,'GLUDym',4.4962,'u');
modeltest = changeRxnBounds(modeltest,'MDH',21.26331648,'u');
modeltest = changeRxnBounds(modeltest,'MDHm',21.26331648,'u');
modeltest = changeRxnBounds(modeltest,'OCOAT1m',1.71127224,'u');
modeltest = changeRxnBounds(modeltest,'DM_13_cis_retn_n_',0,'l');
modeltest = changeRxnBounds(modeltest,'ACACT1',0,'l');
modeltest = changeRxnBounds(modeltest,'ACACT1x',0,'l');
modeltest = changeRxnBounds(modeltest,'ADRNCPT1',0,'l');
modeltest = changeRxnBounds(modeltest,'AKGDm',0,'l');
modeltest = changeRxnBounds(modeltest,'C160CPT1',0,'l');
modeltest = changeRxnBounds(modeltest,'C161CPT1',0,'l');
modeltest = changeRxnBounds(modeltest,'C161CPT12',0,'l');
modeltest = changeRxnBounds(modeltest,'CLPNDCPT1',0,'l');
modeltest = changeRxnBounds(modeltest,'CSm',0,'l');
modeltest = changeRxnBounds(modeltest,'CYOOm3',0,'l');
modeltest = changeRxnBounds(modeltest,'DCSPTN1CPT1',0,'l');
modeltest = changeRxnBounds(modeltest,'DLNLCGCPT1',0,'l');
modeltest = changeRxnBounds(modeltest,'EICOSTETCPT1',0,'l');
modeltest = changeRxnBounds(modeltest,'ELAIDCPT1',0,'l');
modeltest = changeRxnBounds(modeltest,'G6PDH1rer',0,'l');
modeltest = changeRxnBounds(modeltest,'GLPASE2',0,'l');
modeltest = changeRxnBounds(modeltest,'GLUNm',0,'l');
modeltest = changeRxnBounds(modeltest,'GND',0,'l');
modeltest = changeRxnBounds(modeltest,'GNDer',0,'l');
modeltest = changeRxnBounds(modeltest,'HEX1',0,'l');
modeltest = changeRxnBounds(modeltest,'HEXCCPT1',0,'l');
modeltest = changeRxnBounds(modeltest,'L_LACtcm',0,'l');
modeltest = changeRxnBounds(modeltest,'L_LACtm',0,'l');
modeltest = changeRxnBounds(modeltest,'LGNCCPT1',0,'l');
modeltest = changeRxnBounds(modeltest,'LNELDCCPT1',0,'l');
modeltest = changeRxnBounds(modeltest,'LNLCCPT1',0,'l');
modeltest = changeRxnBounds(modeltest,'LNLNCACPT1',0,'l');
modeltest = changeRxnBounds(modeltest,'LNLNCGCPT1',0,'l');
modeltest = changeRxnBounds(modeltest,'NRVNCCPT1',0,'l');
modeltest = changeRxnBounds(modeltest,'PCm',0,'l');
modeltest = changeRxnBounds(modeltest,'PDHm',0,'l');
modeltest = changeRxnBounds(modeltest,'PEPCK',0,'l');
modeltest = changeRxnBounds(modeltest,'PEPCKm',0,'l');
modeltest = changeRxnBounds(modeltest,'PFK',0,'l');
modeltest = changeRxnBounds(modeltest,'PIt2m',0,'l');
modeltest = changeRxnBounds(modeltest,'PYK',0,'l');
modeltest = changeRxnBounds(modeltest,'STRDNCCPT1',0,'l');
modeltest = changeRxnBounds(modeltest,'TETPENT3CPT1',0,'l');
modeltest = changeRxnBounds(modeltest,'TETPENT6CPT1',0,'l');
modeltest = changeRxnBounds(modeltest,'TETTET6CPT1',0,'l');
modeltest = changeRxnBounds(modeltest,'THD1m',0,'l');
modeltest = changeRxnBounds(modeltest,'TMNDNCCPT1',0,'l');
modeltest = changeRxnBounds(modeltest,'TTDCPT2',0,'l');
modeltest = changeRxnBounds(modeltest,'VACCCPT1',0,'l');

modelnaive = modeltest
modelnaive = changeRxnBounds(modelnaive,'EX_eicostet[e]',-0.01,'l');
modelnaive = changeRxnBounds(modelnaive,'EX_hdca[e]',-0.01,'l');
modelnaive = changeRxnBounds(modelnaive,'EX_hdcea[e]',-0.01,'l');
modelnaive = changeRxnBounds(modelnaive,'EX_lnlc[e]',-0.01,'l');
modelnaive = changeRxnBounds(modelnaive,'EX_lnlnca[e]',-0.01,'l');
modelnaive = changeRxnBounds(modelnaive,'EX_lnlncg[e]',-0.01,'l');
modelnaive = changeRxnBounds(modelnaive,'EX_ocdca[e]',-0.01,'l');
modelnaive = changeRxnBounds(modelnaive,'EX_ocdcea[e]',-0.01,'l');

naiveIrrRxns = {'TTDCPT2','TMNDNCCRNt','TMNDNCCPT2','TMNDNCCPT1','TETTET6CRNt',...
                'TETTET6CPT2','TETTET6CPT1','TETPENT6CRNt','TETPENT6CPT2','TETPENT6CPT1',...
                'TETPENT3CRNt','TETPENT3CPT2','TETPENT3CPT1','STRDNCCRNt','STRDNCCPT2',...
                'STRDNCCPT1','RE0583C','r1401','r1400','r0791','r0735','r0652','r0638','r0636','r0633','r0432',...
                'r0431','r0309','PTDCACRNt','PTDCACRNCPT2','PTDCACRNCPT1','PCRNtm','LNLNCGCRNt',...
                'LNLNCGCPT2','LNLNCGCPT1','LNLNCACRNt','LNLNCACPT2','LNLNCACPT1','LNLCCRNt',...
                'LNLCCPT2','LNLCCPT1','LNELDCCRNt','LNELDCCPT2','LNELDCCPT1','HXCOAx','HXCOAc',...
                'HPDCACRNt','HPDCACRNCPT2','HPDCACRNCPT1','HEXDIACtr','ELAIDCRNt','ELAIDCPT2',...
                'ELAIDCPT1','EICOSTETCRNt','EICOSTETCPT2','EICOSTETCPT1','DLNLCGCRNt',...
                'DLNLCGCPT2','DLNLCGCPT1','DCSPTN1CRNt','DCSPTN1CPT2','DCSPTN1CPT1',...
                'CLPNDCRNt','CLPNDCPT2','CLPNDCPT1','C30CPT1','C226CRNt','C226CPT2','C226CPT1',...
                'C204CRNt','C204CPT2','C204CPT1','C181CRNt','C181CPT2','C181CPT1','C180CPT2',...
                'C180CPT1','C161CRNt','C161CRN2t','C161CPT22','C161CPT2','C161CPT12','C161CPT1',...
                'C160CPT2','C160CPT1','ARACHCPT2','ARACHCPT1','ADRNCRNt','ADRNCPT2','ADRNCPT1'}
[modelnaive_part,matchRevNaive,rev2irrev,irrev2rev] = convertToIrreversible(modelnaive,'sRxns',naiveIrrRxns)

expCol_naive = mapGeneToRxn(modelnaive_part,modelnaive_part.genes,ExpressionNaiveInc,parsedGPR_naive,corrRxn_naive)

NaiveTissueModel = GIMME(modelnaive_part,expCol_naive,0.5)

[NaiveTissueModel_rem,removedMetsNaiveIrr,removedRxnsNaiveIrr] = removeDeadEnds(NaiveTissueModel)
[EffTissueModel_rem,removedMetsEffIrr,removedRxnsEffIrr] = removeDeadEnds(EffTissueModel)

% Glucose-dependence test
for i=0:15
    model = EffTissueModel_rem;
    Glc_lb = -0.3+0.02*i
    model = changeRxnBounds(model,{'EX_glc_D[e]'},Glc_lb,'l');
    solution = optimizeCbModel(model,'max').f
end

% PDHm lactate test
for i=0:10
    model = NaiveTissueModel_rem;
    PDHm_ub = 0.01+0.05*i
    model = changeRxnBounds(model,{'PDHm'},PDHm_ub,'u');
    idx = find(ismember(model.rxns,'LDH_L'));
    solution = optimizeCbModel(model,'max').v(idx)
end


% arg test
model = MemTissueModel_rem;
solution = optimizeCbModel(model,'max').f
model = changeRxnBounds(model,{'EX_arg_L[e]'},0,'l');
solution = optimizeCbModel(model,'max').f

% leu test
model = NaiveTissueModel_rem;
solution = optimizeCbModel(model,'max').f
model = changeRxnBounds(model,{'EX_leu_L[e]'},0,'l');
solution = optimizeCbModel(model,'max').f

% ACC1 test
model = NaiveTissueModel_rem;
model = changeRxnBounds(model,{'ACCOAC'},0,'l');
model = changeRxnBounds(model,{'ACCOAC'},0,'u');
solution = optimizeCbModel(model,'max').f

% Glutamine test
for i=0:15
    model = MemTissueModel_rem;
    Gln_lb = -0.075+0.005*i
    model = changeRxnBounds(model,{'EX_gln_L[e]'},Gln_lb,'l');
    solution1 = optimizeCbModel(model,'max').f
    model = changeRxnBounds(model,{'EX_glc_D[e]'},0,'l');
    solution2 = optimizeCbModel(model,'max').f
end

% General Human Test v5
[TestSolution_naive,TestSolutionName_naive] = Test4HumanFctExtv5(NaiveTissueModel_rem,'all')
[TestSolution_eff,TestSolutionName_eff] = Test4HumanFctExtv5(EffTissueModel_rem,'all')
[TestSolution_mem,TestSolutionName_mem] = Test4HumanFctExtv5(MemTissueModel_rem,'all')




[TestSolution_naive,TestSolutionName_naive] = Test4HumanFctExtv5(NaiveTissueModel_rem,'all')
[TestSolution_eff,TestSolutionName_eff] = Test4HumanFctExtv5(EffTissueModel_rem,'all')

modelnaive = changeRxnBounds(modelnaive,naiveExRxns,-1000,'l')
modelnaive = changeRxnBounds(modelnaive,{'EX_glc_D[e]'},-100,'l')
modelnaive = changeRxnBounds(modelnaive,naiveUpTake,-1,'l')


% media conditions

exRxns = printUptakeBound(modelnewfinal)
modeltest = changeRxnBounds(modelnewfinal,modelnewfinal.rxns(exRxns),0,'l')
naiveExRxns = {'EX_ca2[e]','EX_cl[e]','EX_co[e]','EX_fe2[e]','EX_h[e]','EX_h2o[e]','EX_k[e]','EX_na1[e]','EX_nh4[e]',...
               'EX_no2[e]','EX_o2[e]','EX_pi[e]','EX_so4[e]'}
naiveUpTake = {'EX_ala_L[e]','EX_arg_L[e]',...
               'EX_asn_L[e]','EX_asp_L[e]','EX_cys_L[e]','EX_eicostet[e]','EX_gln_L[e]','EX_glu_L[e]','EX_gly[e]',...
               'EX_hdca[e]','EX_hdcea[e]','EX_his_L[e]','EX_ile_L[e]','EX_leu_L[e]','EX_lnlc[e]','EX_lnlnca[e]',...
               'EX_lnlncg[e]','EX_lys_L[e]','EX_met_L[e]','EX_ocdca[e]','EX_ocdcea[e]','EX_orn_D[e]',...
               'EX_phe_L[e]','EX_pro_L[e]','EX_ribflv[e]','EX_ser_L[e]','EX_thr_L[e]','EX_trp_L[e]','EX_tyr_L[e]',...
               'EX_val_L[e]','EX_vitd3[e]','EX_tag_hs[e]'}
modelnaive = changeRxnBounds(modeltest,naiveExRxns,-1000,'l')
modelnaive = changeRxnBounds(modelnaive,{'EX_glc_D[e]'},-100,'l')
modelnaive = changeRxnBounds(modelnaive,naiveUpTake,-1,'l')


EffExRxns = {'EX_ca2[e]','EX_cl[e]','EX_co2[e]','EX_fe2[e]','EX_h[e]','EX_h2o[e]','EX_k[e]','EX_na1[e]',...
             'EX_nh4[e]','EX_no2[e]','EX_pi[e]','EX_so4[e]'}
EffUptakeRxns = {'EX_ala_L[e]','EX_arg_L[e]','EX_asn_L[e]','EX_asp_L[e]','EX_cys_L[e]','EX_gln_L[e]','EX_glu_L[e]','EX_gly[e]',...
                 'EX_his_L[e]','EX_ile_L[e]','EX_leu_L[e]','EX_lys_L[e]','EX_met_L[e]','EX_orn[e]','EX_orn_D[e]',...
                 'EX_phe_L[e]','EX_pro_L[e]','EX_ribflv[e]','EX_ser_L[e]','EX_thr_L[e]','EX_trp_L[e]','EX_tyr_L[e]',...
                 'EX_val_L[e]','EX_vitd3[e]','EX_tag_hs[e]'}
modeleff = changeRxnBounds(modeltest,EffExRxns,-1000,'l')
modeleff = changeRxnBounds(modeleff,{'EX_glc_D[e]'},-100,'l')
modeleff = changeRxnBounds(modeleff,{'EX_o2[e]'},-20,'l')
modeleff = changeRxnBounds(modeleff,EffUptakeRxns,-1,'l')

MemExRxns =   {'EX_ca2[e]','EX_cl[e]','EX_co[e]','EX_co2[e]','EX_fe2[e]','EX_h[e]','EX_h2o[e]','EX_k[e]','EX_na1[e]','EX_nh4[e]',...
               'EX_no2[e]','EX_o2[e]','EX_pi[e]','EX_so4[e]'}
MemUptakeRxns = {'EX_ala_L[e]','EX_arg_L[e]','EX_asn_L[e]','EX_asp_L[e]','EX_cys_L[e]','EX_eicostet[e]','EX_gln_L[e]','EX_glu_L[e]','EX_gly[e]',...
                 'EX_hdca[e]','EX_hdcea[e]','EX_his_L[e]','EX_ile_L[e]','EX_leu_L[e]','EX_lnlc[e]','EX_lnlnca[e]',...
                 'EX_lnlncg[e]','EX_lys_L[e]','EX_met_L[e]','EX_ocdca[e]','EX_ocdcea[e]','EX_orn_D[e]',...
                 'EX_phe_L[e]','EX_pro_L[e]','EX_ribflv[e]','EX_ser_L[e]','EX_thr_L[e]','EX_trp_L[e]','EX_tyr_L[e]',...
                 'EX_val_L[e]','EX_vitd3[e]','EX_1glyc_hs[e]','EX_ak2lgchol_hs[e]','EX_mag_hs[e]','EX_paf_hs[e]',...
                 'EX_pglyc_hs[e]','EX_tag_hs[e]','EX_HC01444[e]','EX_12dgr120[e]','EX_magarachi_hs[e]',...
                 'EX_maglinl_hs[e]','EX_magole_hs[e]','EX_magpalm_hs[e]','EX_magste_hs[e]','EX_glyc2p[e]','EX_glyc[e]'}
modelmem = changeRxnBounds(modeltest,MemExRxns,-1000,'l')
modelmem = changeRxnBounds(modelmem,{'EX_glc_D[e]'},-100,'l')
modelmem = changeRxnBounds(modelmem,{'EX_o2[e]'},-20,'l')
modelmem = changeRxnBounds(modelmem,MemUptakeRxns,-1,'l')

% Directionality changes 
naiveIrrRxns = {'TTDCPT2','TMNDNCCRNt','TMNDNCCPT2','TMNDNCCPT1','TETTET6CRNt',...
                'TETTET6CPT2','TETTET6CPT1','TETPENT6CRNt','TETPENT6CPT2','TETPENT6CPT1',...
                'TETPENT3CRNt','TETPENT3CPT2','TETPENT3CPT1','STRDNCCRNt','STRDNCCPT2',...
                'STRDNCCPT1','RE0583C','r1401','r1400','r0791','r0735','r0652','r0638','r0636','r0633','r0432',...
                'r0431','r0309','PTDCACRNt','PTDCACRNCPT2','PTDCACRNCPT1','PCRNtm','LNLNCGCRNt',...
                'LNLNCGCPT2','LNLNCGCPT1','LNLNCACRNt','LNLNCACPT2','LNLNCACPT1','LNLCCRNt',...
                'LNLCCPT2','LNLCCPT1','LNELDCCRNt','LNELDCCPT2','LNELDCCPT1','HXCOAx','HXCOAc',...
                'HPDCACRNt','HPDCACRNCPT2','HPDCACRNCPT1','HEXDIACtr','ELAIDCRNt','ELAIDCPT2',...
                'ELAIDCPT1','EICOSTETCRNt','EICOSTETCPT2','EICOSTETCPT1','DLNLCGCRNt',...
                'DLNLCGCPT2','DLNLCGCPT1','DCSPTN1CRNt','DCSPTN1CPT2','DCSPTN1CPT1',...
                'CLPNDCRNt','CLPNDCPT2','CLPNDCPT1','C30CPT1','C226CRNt','C226CPT2','C226CPT1',...
                'C204CRNt','C204CPT2','C204CPT1','C181CRNt','C181CPT2','C181CPT1','C180CPT2',...
                'C180CPT1','C161CRNt','C161CRN2t','C161CPT22','C161CPT2','C161CPT12','C161CPT1',...
                'C160CPT2','C160CPT1','ARACHCPT2','ARACHCPT1','ADRNCRNt','ADRNCPT2','ADRNCPT1'}
[modelnaive_part,matchRevNaive,rev2irrev,irrev2rev] = convertToIrreversible(modelnaive,'sRxns',naiveIrrRxns)

EffIrrRxns = {'ADRNCPT2','ADRNCRNt','C160CPT2','C161CPT2','C161CPT22','C161CRN2t','C161CRNt',...
              'C180CPT2','C181CPT2','C181CRNt','C204CPT2','C204CRNt','C226CPT2','C226CRNt',...
              'CLPNDCPT2','CLPNDCRNt','DCSPTN1CPT2','DCSPTN1CRNt','DLNLCGCPT2',...
              'DLNLCGCRNt','EICOSTETCPT2','EICOSTETCRNt','ELAIDCPT2','ELAIDCRNt',...
              'HPDCACRNCPT2','HPDCACRNt','HXCOAc','HXCOAx','LNELDCCPT2','LNELDCCRNt',...
              'LNLCCPT2','LNLCCRNt','LNLNCGCPT2','LNLNCGCRNt','PCRNtm','PTDCACRNCPT2',...
              'PTDCACRNt','r0309','r0633','r0636','r0638','r0652','r0735','r1400','r1401','RE0583C',...
              'STRDNCCPT2','STRDNCCRNt','TMNDNCCPT2','TMNDNCCRNt','TTDCPT2'}
[modeleff_part,matchRevEff,rev2irrev,irrev2rev] = convertToIrreversible(modeleff,'sRxns',EffIrrRxns)

[modelmem_part,matchRevMem,rev2irrev,irrev2rev] = convertToIrreversible(modelmem,'sRxns',EffIrrRxns)

% Update GPR
[parsedGPR_naive_irr,corrRxn_naive_irr] = extractGPRs(modelnaive_irrev)
[parsedGPR_naive,corrRxn_naive] = extractGPRs(modelnaive)

% Read rules files


% Gene expression to reaction expression, regarding GPR
expCol_naive = mapGeneToRxn(modelnaive_part,modelnaive_part.genes,Expression,parsedGPR_naive,corrRxn_naive)
expCol_eff = mapGeneToRxn(modeleff_part,modeleff_part.genes,Expression,parsedGPR_naive,corrRxn_naive)
expCol_mem = mapGeneToRxn(modelmem_part,modelmem_part.genes,ExpressionMem,parsedGPR_naive,corrRxn_naive)

expCol_eff1 = mapGeneToRxn(modeleff_part,modeleff_part.genes,ExpressionEff,parsedGPR_naive,corrRxn_naive)
expCol_eff2 = mapGeneToRxn(modeleff_part,modeleff_part.genes,ExpressionEffIncNaive,parsedGPR_naive,corrRxn_naive)
expCol_eff3 = mapGeneToRxn(modeleff_part,modeleff_part.genes,ExpressionEffPuniya,parsedGPR_naive,corrRxn_naive)
expCol_eff4 = mapGeneToRxn(modeleff_part,modeleff_part.genes,ExpressionEffIncEff,parsedGPR_naive,corrRxn_naive)

expCol_Puniya = mapGeneToRxn(modelnaive,modelnaive.genes,ExpPuniya,parsedGPR_naive,corrRxn_naive)

% GIMME
NaiveTissueModel = GIMME(modelnaive_part,expCol_naive,0.5)
EffTissueModel = GIMME(modeleff_part,expCol_eff,0.5)
MemTissueModel = GIMME(modelmem_part,expCol_mem,0.5)

EffTissueModel1 = GIMME(modeleff_part,expCol_eff1,0.5)
EffTissueModel2 = GIMME(modeleff_part,expCol_eff2,0.5)
EffTissueModel3 = GIMME(modeleff_part,expCol_eff3,0.5)
EffTissueModel4 = GIMME(modeleff_part,expCol_eff4,0.5)

NaiveTissuePuniya = GIMME(modelnaive,expCol_Puniya,0.5)

% Remove dead end reactions
[NaiveTissueModel_rem,removedMetsNaiveIrr,removedRxnsNaiveIrr] = removeDeadEnds(NaiveTissueModel)
[EffTissueModel_rem,removedMetsEffIrr,removedRxnsEffIrr] = removeDeadEnds(EffTissueModel)
[MemTissueModel_rem,removedMetsMemIrr,removedRxnsMemIrr] = removeDeadEnds(MemTissueModel)

outmodel = writeCbModel(NaiveTissueModel_rem,'format','mat')
outmodel = writeCbModel(EffTissueModel,'format','sbml')
outmodel = writeCbModel(MemTissueModel,'format','sbml')

NaiveTissueModel = readCbModel('C:\Users\MSI\cobratoolbox\NaiveFinal.mat')
EffTissueModel = readCbModel('C:\Users\MSI\cobratoolbox\EffFinal.mat')

NaiveTissueModel = changeRxnBounds(NaiveTissueModel,{'EX_glc_D[e]'},-100,'l')
NaiveTissueModel = changeRxnBounds(NaiveTissueModel,naiveUpTake,-1,'l')

[TestSolution_eff1,TestSolutionName_eff1] = Test4HumanFctExtv5(EffTissueModel1,'all')
[TestSolution_eff2,TestSolutionName_eff2] = Test4HumanFctExtv5(EffTissueModel2,'all')
[TestSolution_eff3,TestSolutionName_eff3] = Test4HumanFctExtv5(EffTissueModel3,'all')
[TestSolution_eff4,TestSolutionName_eff4] = Test4HumanFctExtv5(EffTissueModel4,'all')

% Leak test

modelClosed = NaiveTissueModel_rem;
modelClosed = addDemandReaction(modelClosed,modelClosed.mets);

modelexchanges1 = strmatch('Ex_',modelClosed.rxns);
modelexchanges4 = strmatch('EX_',modelClosed.rxns);
modelexchanges2 = strmatch('DM_',modelClosed.rxns);
modelexchanges3 = strmatch('sink_',modelClosed.rxns);
selExc = (find( full((sum(abs(modelClosed.S)==1,1) ==1) & (sum(modelClosed.S~=0) == 1))))';

modelexchanges = unique([modelexchanges1;modelexchanges2;modelexchanges3;modelexchanges4;selExc]);
modelClosed.lb(ismember(modelClosed.rxns,modelClosed.rxns(modelexchanges)))=0;

[LeakMets,modelClosed] = fastLeakTest(modelClosed, modelClosed.rxns(modelexchanges),0);

% ATP Production aerobic
model = NaiveTissueModel_rem;
model.c(find(model.c)) = 0;
model.lb(ismember(model.rxns,'EX_glc_D[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc_D[e]'))=-1;
model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
model.c(ismember(model.rxns,'DM_atp_c_'))=1;
FBA = optimizeCbModel(model,'max')

% Glucose-dependence test
for i=0:15
    model = NaiveTissueModel_rem;
    Glc_lb = -0.3+0.02*i
    model = changeRxnBounds(model,{'EX_glc_D[e]'},Glc_lb,'l');
    %model = changeRxnBounds(model,{'EX_gln_L[e]'},Glc_lb,'l');
    solution = optimizeCbModel(model,'max').f
end

% PDHm flux effect on Lactate Production test

% 

% Gene essentiality
[grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(NaiveTissueModel_rem, 'FBA')
[grRatio_eff, grRateKO_eff, grRateWT_eff, hasEffect_eff, delRxns_eff, fluxSolution_eff] = singleGeneDeletion(EffTissueModel_rem, 'MOMA')
[grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(MemTissueModel_rem, 'MOMA')

[solutionDel, solutionWT, totalFluxDiff, solStatus] = MOMA(modelWT, modelDel, osenseStr, verbFlag, minNormFlag)



% MTA
changeCobraSolver ('tomlab_cplex', 'MIQP');

solution = optimizeCbModel(NaiveTissueModel_rem,'max')
rxnFBS = diffexprs2rxnFBS(NaiveTissueModel_rem, diffexpEff, solution.v, 'SeparateTranscript', '.', 'logFC', 0, 'pval',0.05)
epsilonOptimal = calculateEPSILON(points, rxnFBS, 'unique_epsilon', false, 'minimum', 1e-3)
[TSscore,deletedGenes,Vout] = MTA(NaiveTissueModel_rem, rxnFBS, solution.v, 'FORCE_CPLEX', 0,'epsilon',epsilonOptimal)

solution = optimizeCbModel(NaiveTissueModel_rem,'max')
rxnFBS = diffexprs2rxnFBS(NaiveTissueModel_rem, diffexpMem, solution.v, 'SeparateTranscript', '.', 'logFC', 0, 'pval',0.1)
epsilonOptimal = calculateEPSILON(points, rxnFBS, 'unique_epsilon', false, 'minimum', 1e-3)
[TSscore_mem,deletedGenes_mem,Vout_mem] = MTA(NaiveTissueModel_rem, rxnFBS, solution.v, 'FORCE_CPLEX', 0)

% FVA
[minFlux_naive, maxFlux_naive] = fluxVariability(NaiveTissueModel_rem, 90, 'printLevel', 1);
NaiveTissueModel_WT = NaiveTissueModel_rem
NaiveTissueModel_WT.lb = minFlux_naive;
NaiveTissueModel_WT.ub = maxFlux_naive;

[minFlux_eff, maxFlux_eff] = fluxVariability(EffTissueModel_rem, 90, 'printLevel', 1);
EffTissueModel_WT = EffTissueModel_rem
EffTissueModel_WT.lb = minFlux_eff
EffTissueModel_WT.ub = maxFlux_eff

[minFlux_mem, maxFlux_mem] = fluxVariability(MemTissueModel_rem, 90, 'printLevel', 1);
MemTissueModel_WT = MemTissueModel_rem
MemTissueModel_WT.lb = minFlux_mem
MemTissueModel_WT.ub = maxFlux_mem  



% ACHR Sampling of flux space
warmupPts_WT = createHRWarmup(NaiveTissueModel_WT, 10000);
ACHRSampler(NaiveTissueModel_WT, warmupPts_WT, 'Naive_samples_WT', 1, 10000, 2)

warmupPts_Eff = createHRWarmup(EffTissueModel_WT, 10000);
ACHRSampler(EffTissueModel_WT, warmupPts_Eff, 'Eff_samples_WT', 1, 10000, 2)

warmupPts_Mem = createHRWarmup(MemTissueModel_WT, 10000);
ACHRSampler(MemTissueModel_WT, warmupPts_Mem, 'Mem_samples_WT', 1, 10000, 2)

% Mutant 1 promotes differentiation to effector
Eff_Dels = {'9380.1','5238.1','51703.1','5243.1','27010.1','5833.1','7390.1','6652.1','2806.1','79896.1'}
[NaiveTissueModel_EffMut, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes(NaiveTissueModel_WT, Eff_Dels)

warmupPts_EffMut = createHRWarmup(NaiveTissueModel_EffMut, 10000);
ACHRSampler(NaiveTissueModel_EffMut, warmupPts_EffMut, 'Naive_samples_EffMut', 1, 10000, 2)

Mem_Dels = {'9563.1','7357.1','549.1','114971.1','51703.1','6526.1','79178.1','3988.1','2937.1','10327.1'}
[NaiveTissueModel_MemMut, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes(NaiveTissueModel_WT, Mem_Dels)
ACHRSampler(NaiveTissueModel_MemMut, warmupPts_Mem, 'Naive_samples_MemMut', 1, 10000, 2)

compareSampleTraj({'LDH_L'}, {samples_wt, samples_effmut1}, {NaiveTissueModel_WT,NaiveTissueModel_EffMut1},400)
[stats, pVals] = compareTwoSamplesStat(samples_wt, samples_effmut1, {'ks','rankSum','tTest'});

histogram(samples_wt(find(ismember(NaiveTissueModel_WT.rxns,"PPAm")),:))
hold on
histogram(samples_effmut1(find(ismember(NaiveTissueModel_EffMut1.rxns,"PPAm")),:))
histogram(samples_eff(find(ismember(EffTissueModel_rem.rxns,"PPAm")),:))
hold off

% Mutant 2 promotes differentiation to memory
[NaiveTissueModel_MemMut1, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes(NaiveTissueModel_WT, {'9563.1'})
[NaiveTissueModel_MemMut2, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes(NaiveTissueModel_WT, {'7357.1'})
[NaiveTissueModel_MemMut3, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes(NaiveTissueModel_WT, {'549.1'})
[NaiveTissueModel_MemMut4, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes(NaiveTissueModel_WT, {'114971.1'})
[NaiveTissueModel_MemMut5, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes(NaiveTissueModel_WT, {'51703.1'})
[NaiveTissueModel_MemMut6, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes(NaiveTissueModel_WT, {'8879.1'})
[NaiveTissueModel_MemMut7, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes(NaiveTissueModel_WT, {'3073.1'})
[NaiveTissueModel_MemMut8, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes(NaiveTissueModel_WT, {'4351.1'})
[NaiveTissueModel_MemMut9, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes(NaiveTissueModel_WT, {'4329.1'})

warmupPts_MemMut_1 = createHRWarmup(NaiveTissueModel_MemMut1, 10000);
ACHRSampler(NaiveTissueModel_MemMut1, warmupPts_MemMut_1, 'Naive_samples_MemMut1', 1, 10000, 2)
ACHRSampler(NaiveTissueModel_MemMut2, warmupPts_MemMut_1, 'Naive_samples_MemMut2', 1, 10000, 2)
ACHRSampler(NaiveTissueModel_MemMut3, warmupPts_MemMut_1, 'Naive_samples_MemMut3', 1, 10000, 2)
ACHRSampler(NaiveTissueModel_MemMut4, warmupPts_MemMut_1, 'Naive_samples_MemMut4', 1, 10000, 2)
ACHRSampler(NaiveTissueModel_MemMut5, warmupPts_MemMut_1, 'Naive_samples_MemMut5', 1, 10000, 2)
ACHRSampler(NaiveTissueModel_MemMut6, warmupPts_MemMut_1, 'Naive_samples_MemMut6', 1, 10000, 2)
ACHRSampler(NaiveTissueModel_MemMut7, warmupPts_MemMut_1, 'Naive_samples_MemMut7', 1, 10000, 2)
ACHRSampler(NaiveTissueModel_MemMut8, warmupPts_MemMut_1, 'Naive_samples_MemMut8', 1, 10000, 2)
ACHRSampler(NaiveTissueModel_MemMut9, warmupPts_MemMut_1, 'Naive_samples_MemMut9', 1, 10000, 2)

[NaiveTissueModel_EffMut1, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes(NaiveTissueModel_WT, {'9380.1'})
ACHRSampler(NaiveTissueModel_EffMut1, warmupPts_EffMut_1, 'Naive_samples_EffMut1', 1, 10000, 2)

[NaiveTissueModel_EffMut2, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes(NaiveTissueModel_WT, {'5238.1'})
ACHRSampler(NaiveTissueModel_EffMut2, warmupPts_EffMut_1, 'Naive_samples_EffMut2', 1, 10000, 2)




% Comparison of OXPHOS and Glycolysis
samples_wt = loadSamples('Naive_samples_WT', 1, 10000);
samples_wt2 = loadSamples('Naive_samples_WT2', 1, 10000);
samples_effmut1 = loadSamples('Naive_samples_EffMut1', 1, 10000);
samples_effmut2 = loadSamples('Naive_samples_EffMut2', 1, 10000);
samples_effmut = loadSamples('Naive_samples_EffMut', 1, 10000);

samples_memmut1 = loadSamples('Naive_samples_MemMut1', 1, 10000);
samples_eff = loadSamples('Eff_samples_WT', 1, 10000);
samples_mem = loadSamples('Mem_samples_WT', 1, 10000);

samples_memmut3 = loadSamples('Naive_samples_MemMut3', 1, 10000);
samples_memmut4 = loadSamples('Naive_samples_MemMut4', 1, 10000);
samples_memmut5 = loadSamples('Naive_samples_MemMut5', 1, 10000);
samples_memmut6 = loadSamples('Naive_samples_MemMut6', 1, 10000);
samples_memmut7 = loadSamples('Naive_samples_MemMut7', 1, 10000);
samples_memmut8 = loadSamples('Naive_samples_MemMut8', 1, 10000);
samples_memmut9 = loadSamples('Naive_samples_MemMut9', 1, 10000);

% OXPHOS
energySubSystems = {'Oxidative phosphorylation'};
energyReactions = NaiveTissueModel_WT.rxns(ismember(NaiveTissueModel_WT.subSystems(1:5511),energySubSystems));
[~,energy_rxnID] = ismember(energyReactions,NaiveTissueModel_WT.rxns(1:5511));
reactionNames = NaiveTissueModel_WT.rxnNames(energy_rxnID);
reactionFormulas = printRxnFormula(NaiveTissueModel_WT,energyReactions,0);
T = table(reactionNames,reactionFormulas,'RowNames',energyReactions)

energySubSystems = {'Oxidative phosphorylation'};
energyReactions = EffTissueModel_rem.rxns(ismember(EffTissueModel_rem.subSystems(1:5511),energySubSystems));
[~,energy_rxnID] = ismember(energyReactions,EffTissueModel_rem.rxns(1:5511));
reactionNames = EffTissueModel_rem.rxnNames(energy_rxnID);
reactionFormulas = printRxnFormula(EffTissueModel_rem,energyReactions,0);
T = table(reactionNames,reactionFormulas,'RowNames',energyReactions)

compareSampleTraj(energyReactions, {samples_wt, samples_eff, samples_mem}, {NaiveTissueModel_WT,EffTissueModel_WT,MemTissueModel_WT},400)

energySubSystems = {'Glycolysis/gluconeogenesis'};
energyReactions = MemTissueModel_WT.rxns(ismember(MemTissueModel_WT.subSystems(1:5511),energySubSystems));
[~,energy_rxnID] = ismember(energyReactions,MemTissueModel_WT.rxns(1:5511));
reactionNames = MemTissueModel_WT.rxnNames(energy_rxnID);
reactionFormulas = printRxnFormula(MemTissueModel_WT,energyReactions,0);
T = table(reactionNames,reactionFormulas,'RowNames',energyReactions)

compareSampleTraj({'r0354','FBA','PGI','PYK','HMR_7748','HMR_7749'}, {samples_wt, samples_eff, samples_mem}, {NaiveTissueModel_WT,EffTissueModel_WT,MemTissueModel_WT},400)

compareSampleTraj({'r0354','FBA','PGI','PYK','HMR_7748','HMR_7749'}, {samples_wt, samples_mem, samples_eff}, {NaiveTissueModel_WT,MemTissueModel_rem,EffTissueModel_rem},400)

histogram(samples_wt(find(ismember(NaiveTissueModel_WT.rxns,"CYOOm2i")),:), linspace(0.12,0.28,100), 'FaceColor','#7F7FFF')
hold on
histogram(samples_memmut3(find(ismember(NaiveTissueModel_WT.rxns,"CYOOm2i")),:), linspace(0.12,0.28,100),'FaceColor','#FF7F7F')
histogram(samples_effmut(find(ismember(NaiveTissueModel_WT.rxns,"CYOOm2i")),:), linspace(0.12,0.28,100),'FaceColor','#7FFFFF')
histogram(samples_mem(find(ismember(MemTissueModel_WT.rxns,"CYOOm2i")),:), linspace(0.12,0.28,100),'FaceColor','#FFFF7F')
histogram(samples_eff(find(ismember(EffTissueModel_WT.rxns,"CYOOm2i")),:), linspace(0.12,0.28,100),'FaceColor','#7FFF7F')
xlabel('Flux through OXPHOS (mmol/gDW/hr)')
ylabel('Count')
title('Histogram of OXPHOS fluxes in samples')
legend('Naive','Mutant memory','Mutant effector','Memory','Effector')
hold off

rnd_idx = randi([1 10000],1,1000);
sampleStats_wt = calcSampleStats(samples_wt(find(ismember(NaiveTissueModel_WT.rxns,"CYOOm2i")),rnd_idx))
sampleStats_mem = calcSampleStats(samples_mem(find(ismember(MemTissueModel_WT.rxns,"CYOOm2i")),rnd_idx))
sampleStats_eff = calcSampleStats(samples_eff(find(ismember(EffTissueModel_WT.rxns,"CYOOm2i")),rnd_idx))
sampleStats_effmut = calcSampleStats(samples_effmut(find(ismember(NaiveTissueModel_WT.rxns,"CYOOm2i")),rnd_idx))
sampleStats_memmut1 = calcSampleStats(samples_memmut1(find(ismember(NaiveTissueModel_WT.rxns,"CYOOm2i")),rnd_idx))

sampleStats_mem = calcSampleStats(samples_mem(find(ismember(MemTissueModel_WT.rxns,"FADH2ETC")),rnd_idx))

rnd_idx = randi([1 10000],1,1000);
[stats, pVals] = compareTwoSamplesStat(samples_memmut1(:,rnd_idx),samples_mem(:,rnd_idx),{'rankSum'})
pVals.ttest(find(ismember(NaiveTissueModel_WT.rxns,"FADH2ETC")))
[sampleDiff_mem, sampleRatio] = calcSampleDifference(samples_wt, samples_mem, 1000);
[sampleDiff_mut1, sampleRatio] = calcSampleDifference(samples_wt, samples_memmut1, 1000);

histogram(sampleDiff(find(ismember(EffTissueModel_rem.rxns,"PYK")),:))
ttest2()

rxnNames = {'PGI','PFK','FBP','FBA','TPI','GAPD','PGK','PGM','ENO','PYK','PPS'};
sampleScatterMatrix(rxnNames,modelSampling,samples);

substrateRxns = {'EX_glc_D[e]','EX_ca2[e]','EX_k[e]','EX_na1[e]','EX_ala_L[e]','EX_arg_L[e]','EX_asn_L[e]','EX_asp_L[e]','EX_cys_L[e]','EX_gln_L[e]','EX_glu_L[e]','EX_his_L[e]','EX_ile_L[e]',....
                 'EX_leu_L[e]','EX_lys_L[e]','EX_met_L[e]','EX_phe_L[e]','EX_pro_L[e]','EX_ser_L[e]','EX_thr_L[e]','EX_trp_L[e]','EX_tyr_L[e]','EX_val_L[e]','EX_btn[e]',....
                 'EX_fol[e]','EX_inost[e]','EX_bz[e]','EX_pnto_R[e]','EX_pydxn[e]','EX_ribflv[e]','EX_but[e]','EX_2hb[e]','EX_2m3hbu[e]','EX_3bcrn[e]','EX_4ohbut[e]','EX_c4crn[e]'}
initConcentrations = [2.7,0.1,0.4,8.8,0.43,0.2,0.05,0.02,0.0652,0.3,0.02,0.015,0.05,0.05,0.04,0.015,0.015,0.02,0.030,0.02,0.005,0.028,0.02,0.002,0.001,0.035,0.001,0.00025,0.001,0.0002,0,0,0,0,0,0]
[concentrationMatrix, excRxnNames, timeVec, biomassVec] = dynamicFBA(NaiveTissueModel_WT, substrateRxns, initConcentrations, 0.005, 12, 8, {'EX_glc_D[e]','EX_gln_L[e]','EX_3bcrn[e]','EX_4ohbut[e]','EX_c4crn[e]'}, {'EX_co2[e]', 'EX_o2[e]', 'EX_h2o[e]', 'EX_h[e]'})

NaiveModel_iMAT = iMAT(modelnaive_part, expCol_naive, 0.5, 0.5)
solution = optimizeCbModel(NaiveModel_iMAT)



modelnaive = changeRxnBounds(modelnaive,naiveExRxns,-1000,'l')
modelnaive = changeRxnBounds(modelnaive,{'EX_glc_D[e]'},-100,'l')
modelnaive = changeRxnBounds(modelnaive,naiveUpTake,-1,'l')

% media conditions

exRxns = printUptakeBound(modelnewfinal)
modeltest = changeRxnBounds(modelnewfinal,modelnewfinal.rxns(exRxns),0,'l')
naiveExRxns = {'EX_ca2[e]','EX_cl[e]','EX_co[e]','EX_fe2[e]','EX_h[e]','EX_h2o[e]','EX_k[e]','EX_na1[e]','EX_nh4[e]',...
               'EX_no2[e]','EX_o2[e]','EX_pi[e]','EX_so4[e]'}
naiveUpTake = {'EX_ala_L[e]','EX_arg_L[e]',...
               'EX_asn_L[e]','EX_asp_L[e]','EX_cys_L[e]','EX_eicostet[e]','EX_gln_L[e]','EX_glu_L[e]','EX_gly[e]',...
               'EX_hdca[e]','EX_hdcea[e]','EX_his_L[e]','EX_ile_L[e]','EX_leu_L[e]','EX_lnlc[e]','EX_lnlnca[e]',...
               'EX_lnlncg[e]','EX_lys_L[e]','EX_met_L[e]','EX_ocdca[e]','EX_ocdcea[e]','EX_orn_D[e]',...
               'EX_phe_L[e]','EX_pro_L[e]','EX_ribflv[e]','EX_ser_L[e]','EX_thr_L[e]','EX_trp_L[e]','EX_tyr_L[e]',...
               'EX_val_L[e]','EX_vitd3[e]','EX_tag_hs[e]'}
modelnaive = changeRxnBounds(modeltest,naiveExRxns,-1000,'l')
modelnaive = changeRxnBounds(modelnaive,{'EX_glc_D[e]'},-100,'l')
modelnaive = changeRxnBounds(modelnaive,naiveUpTake,-1,'l')


EffExRxns = {'EX_ca2[e]','EX_cl[e]','EX_co2[e]','EX_fe2[e]','EX_h[e]','EX_h2o[e]','EX_k[e]','EX_na1[e]',...
             'EX_nh4[e]','EX_no2[e]','EX_pi[e]','EX_so4[e]'}
EffUptakeRxns = {'EX_ala_L[e]','EX_arg_L[e]','EX_asn_L[e]','EX_asp_L[e]','EX_cys_L[e]','EX_gln_L[e]','EX_glu_L[e]','EX_gly[e]',...
                 'EX_his_L[e]','EX_ile_L[e]','EX_leu_L[e]','EX_lys_L[e]','EX_met_L[e]','EX_orn[e]','EX_orn_D[e]',...
                 'EX_phe_L[e]','EX_pro_L[e]','EX_ribflv[e]','EX_ser_L[e]','EX_thr_L[e]','EX_trp_L[e]','EX_tyr_L[e]',...
                 'EX_val_L[e]','EX_vitd3[e]','EX_tag_hs[e]'}
modeleff = changeRxnBounds(modeltest,EffExRxns,-1000,'l')
modeleff = changeRxnBounds(modeleff,{'EX_glc_D[e]'},-100,'l')
modeleff = changeRxnBounds(modeleff,{'EX_o2[e]'},-20,'l')
modeleff = changeRxnBounds(modeleff,EffUptakeRxns,-1,'l')

MemExRxns =   {'EX_ca2[e]','EX_cl[e]','EX_co[e]','EX_co2[e]','EX_fe2[e]','EX_h[e]','EX_h2o[e]','EX_k[e]','EX_na1[e]','EX_nh4[e]',...
               'EX_no2[e]','EX_o2[e]','EX_pi[e]','EX_so4[e]'}
MemUptakeRxns = {'EX_ala_L[e]','EX_arg_L[e]','EX_asn_L[e]','EX_asp_L[e]','EX_cys_L[e]','EX_eicostet[e]','EX_gln_L[e]','EX_glu_L[e]','EX_gly[e]',...
                 'EX_hdca[e]','EX_hdcea[e]','EX_his_L[e]','EX_ile_L[e]','EX_leu_L[e]','EX_lnlc[e]','EX_lnlnca[e]',...
                 'EX_lnlncg[e]','EX_lys_L[e]','EX_met_L[e]','EX_ocdca[e]','EX_ocdcea[e]','EX_orn_D[e]',...
                 'EX_phe_L[e]','EX_pro_L[e]','EX_ribflv[e]','EX_ser_L[e]','EX_thr_L[e]','EX_trp_L[e]','EX_tyr_L[e]',...
                 'EX_val_L[e]','EX_vitd3[e]','EX_1glyc_hs[e]','EX_ak2lgchol_hs[e]','EX_mag_hs[e]','EX_paf_hs[e]',...
                 'EX_pglyc_hs[e]','EX_tag_hs[e]','EX_HC01444[e]','EX_12dgr120[e]','EX_magarachi_hs[e]',...
                 'EX_maglinl_hs[e]','EX_magole_hs[e]','EX_magpalm_hs[e]','EX_magste_hs[e]','EX_glyc2p[e]','EX_glyc[e]'}
modelmem = changeRxnBounds(modeltest,MemExRxns,-1000,'l')
modelmem = changeRxnBounds(modelmem,{'EX_glc_D[e]'},-100,'l')
modelmem = changeRxnBounds(modelmem,{'EX_o2[e]'},-20,'l')
modelmem = changeRxnBounds(modelmem,MemUptakeRxns,-1,'l')

% Directionality changes 
naiveIrrRxns = {'TTDCPT2','TMNDNCCRNt','TMNDNCCPT2','TMNDNCCPT1','TETTET6CRNt',...
                'TETTET6CPT2','TETTET6CPT1','TETPENT6CRNt','TETPENT6CPT2','TETPENT6CPT1',...
                'TETPENT3CRNt','TETPENT3CPT2','TETPENT3CPT1','STRDNCCRNt','STRDNCCPT2',...
                'STRDNCCPT1','RE0583C','r1401','r1400','r0791','r0735','r0652','r0638','r0636','r0633','r0432',...
                'r0431','r0309','PTDCACRNt','PTDCACRNCPT2','PTDCACRNCPT1','PCRNtm','LNLNCGCRNt',...
                'LNLNCGCPT2','LNLNCGCPT1','LNLNCACRNt','LNLNCACPT2','LNLNCACPT1','LNLCCRNt',...
                'LNLCCPT2','LNLCCPT1','LNELDCCRNt','LNELDCCPT2','LNELDCCPT1','HXCOAx','HXCOAc',...
                'HPDCACRNt','HPDCACRNCPT2','HPDCACRNCPT1','HEXDIACtr','ELAIDCRNt','ELAIDCPT2',...
                'ELAIDCPT1','EICOSTETCRNt','EICOSTETCPT2','EICOSTETCPT1','DLNLCGCRNt',...
                'DLNLCGCPT2','DLNLCGCPT1','DCSPTN1CRNt','DCSPTN1CPT2','DCSPTN1CPT1',...
                'CLPNDCRNt','CLPNDCPT2','CLPNDCPT1','C30CPT1','C226CRNt','C226CPT2','C226CPT1',...
                'C204CRNt','C204CPT2','C204CPT1','C181CRNt','C181CPT2','C181CPT1','C180CPT2',...
                'C180CPT1','C161CRNt','C161CRN2t','C161CPT22','C161CPT2','C161CPT12','C161CPT1',...
                'C160CPT2','C160CPT1','ARACHCPT2','ARACHCPT1','ADRNCRNt','ADRNCPT2','ADRNCPT1'}
[modelnaive_part,matchRevNaive,rev2irrev,irrev2rev] = convertToIrreversible(modelnaive,'sRxns',naiveIrrRxns)

EffIrrRxns = {'ADRNCPT2','ADRNCRNt','C160CPT2','C161CPT2','C161CPT22','C161CRN2t','C161CRNt',...
              'C180CPT2','C181CPT2','C181CRNt','C204CPT2','C204CRNt','C226CPT2','C226CRNt',...
              'CLPNDCPT2','CLPNDCRNt','DCSPTN1CPT2','DCSPTN1CRNt','DLNLCGCPT2',...
              'DLNLCGCRNt','EICOSTETCPT2','EICOSTETCRNt','ELAIDCPT2','ELAIDCRNt',...
              'HPDCACRNCPT2','HPDCACRNt','HXCOAc','HXCOAx','LNELDCCPT2','LNELDCCRNt',...
              'LNLCCPT2','LNLCCRNt','LNLNCGCPT2','LNLNCGCRNt','PCRNtm','PTDCACRNCPT2',...
              'PTDCACRNt','r0309','r0633','r0636','r0638','r0652','r0735','r1400','r1401','RE0583C',...
              'STRDNCCPT2','STRDNCCRNt','TMNDNCCPT2','TMNDNCCRNt','TTDCPT2'}
[modeleff_part,matchRevEff,rev2irrev,irrev2rev] = convertToIrreversible(modeleff,'sRxns',EffIrrRxns)

[modelmem_part,matchRevMem,rev2irrev,irrev2rev] = convertToIrreversible(modelmem,'sRxns',EffIrrRxns)

% Update GPR
[parsedGPR_naive_irr,corrRxn_naive_irr] = extractGPRs(modelnaive_irrev)
[parsedGPR_naive,corrRxn_naive] = extractGPRs(modelnaive)

% Read rules files


% Gene expression to reaction expression, regarding GPR
expCol_naive = mapGeneToRxn(modelnaive_part,modelnaive_part.genes,Expression,parsedGPR_naive,corrRxn_naive)
expCol_eff = mapGeneToRxn(modeleff_part,modeleff_part.genes,Expression,parsedGPR_naive,corrRxn_naive)
expCol_mem = mapGeneToRxn(modelmem_part,modelmem_part.genes,ExpressionMem,parsedGPR_naive,corrRxn_naive)

expCol_eff1 = mapGeneToRxn(modeleff_part,modeleff_part.genes,ExpressionEff,parsedGPR_naive,corrRxn_naive)
expCol_eff2 = mapGeneToRxn(modeleff_part,modeleff_part.genes,ExpressionEffIncNaive,parsedGPR_naive,corrRxn_naive)
expCol_eff3 = mapGeneToRxn(modeleff_part,modeleff_part.genes,ExpressionEffPuniya,parsedGPR_naive,corrRxn_naive)
expCol_eff4 = mapGeneToRxn(modeleff_part,modeleff_part.genes,ExpressionEffIncEff,parsedGPR_naive,corrRxn_naive)

expCol_Puniya = mapGeneToRxn(modelnaive,modelnaive.genes,ExpPuniya,parsedGPR_naive,corrRxn_naive)

% GIMME
NaiveTissueModel = GIMME(modelnaive_part,expCol_naive,0.5)
EffTissueModel = GIMME(modeleff_part,expCol_eff,0.5)
MemTissueModel = GIMME(modelmem_part,expCol_mem,0.5)

EffTissueModel1 = GIMME(modeleff_part,expCol_eff1,0.5)
EffTissueModel2 = GIMME(modeleff_part,expCol_eff2,0.5)
EffTissueModel3 = GIMME(modeleff_part,expCol_eff3,0.5)
EffTissueModel4 = GIMME(modeleff_part,expCol_eff4,0.5)

NaiveTissuePuniya = GIMME(modelnaive,expCol_Puniya,0.5)

% Remove dead end reactions
[NaiveTissueModel_rem,removedMetsNaiveIrr,removedRxnsNaiveIrr] = removeDeadEnds(NaiveTissueModel)
[EffTissueModel_rem,removedMetsEffIrr,removedRxnsEffIrr] = removeDeadEnds(EffTissueModel)
[MemTissueModel_rem,removedMetsMemIrr,removedRxnsMemIrr] = removeDeadEnds(MemTissueModel)

outmodel = writeCbModel(NaiveTissueModel_rem,'format','mat')
outmodel = writeCbModel(EffTissueModel,'format','sbml')
outmodel = writeCbModel(MemTissueModel,'format','sbml')

NaiveTissueModel = readCbModel('C:\Users\MSI\cobratoolbox\NaiveFinal.mat')
EffTissueModel = readCbModel('C:\Users\MSI\cobratoolbox\EffFinal.mat')

NaiveTissueModel = changeRxnBounds(NaiveTissueModel,{'EX_glc_D[e]'},-100,'l')
NaiveTissueModel = changeRxnBounds(NaiveTissueModel,naiveUpTake,-1,'l')

[TestSolution_eff1,TestSolutionName_eff1] = Test4HumanFctExtv5(EffTissueModel1,'all')
[TestSolution_eff2,TestSolutionName_eff2] = Test4HumanFctExtv5(EffTissueModel2,'all')
[TestSolution_eff3,TestSolutionName_eff3] = Test4HumanFctExtv5(EffTissueModel3,'all')
[TestSolution_eff4,TestSolutionName_eff4] = Test4HumanFctExtv5(EffTissueModel4,'all')

% Leak test

modelClosed = NaiveTissueModel_rem;
modelClosed = addDemandReaction(modelClosed,modelClosed.mets);

modelexchanges1 = strmatch('Ex_',modelClosed.rxns);
modelexchanges4 = strmatch('EX_',modelClosed.rxns);
modelexchanges2 = strmatch('DM_',modelClosed.rxns);
modelexchanges3 = strmatch('sink_',modelClosed.rxns);
selExc = (find( full((sum(abs(modelClosed.S)==1,1) ==1) & (sum(modelClosed.S~=0) == 1))))';

modelexchanges = unique([modelexchanges1;modelexchanges2;modelexchanges3;modelexchanges4;selExc]);
modelClosed.lb(ismember(modelClosed.rxns,modelClosed.rxns(modelexchanges)))=0;

[LeakMets,modelClosed] = fastLeakTest(modelClosed, modelClosed.rxns(modelexchanges),0);

% ATP Production aerobic
model = NaiveTissueModel_rem;
model.c(find(model.c)) = 0;
model.lb(ismember(model.rxns,'EX_glc_D[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc_D[e]'))=-1;
model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
model.c(ismember(model.rxns,'DM_atp_c_'))=1;
FBA = optimizeCbModel(model,'max')

% Glucose-dependence test
for i=0:15
    model = NaiveTissueModel_rem;
    Glc_lb = -0.3+0.02*i
    model = changeRxnBounds(model,{'EX_glc_D[e]'},Glc_lb,'l');
    %model = changeRxnBounds(model,{'EX_gln_L[e]'},Glc_lb,'l');
    solution = optimizeCbModel(model,'max').f
end

% PDHm flux effect on Lactate Production test

% 





% pre-FBA check
preFBAcheck = checkModelPreFBA(NaiveTissueModel)

% Test Human Cell Functions
[TestSolution_1,TestSolutionName_1] = Test4HumanFctExtv5(NaiveTissueModel_1,'all') 
[TestSolution_2,TestSolutionName_2] = Test4HumanFctExtv5(NaiveTissueModel_2,'all') 
[TestSolution_3,TestSolutionName_3] = Test4HumanFctExtv5(NaiveTissueModel_3,'all') 
[TestSolution_4,TestSolutionName_4] = Test4HumanFctExtv5(NaiveTissueModel_4,'all') 
[TestSolution_5,TestSolutionName_5] = Test4HumanFctExtv5(NaiveTissueModel_5,'all') 
[TestSolution_6,TestSolutionName_6] = Test4HumanFctExtv5(NaiveTissueModel_6,'all') 
[TestSolution_2_rem,TestSolutionName_2] = Test4HumanFctExtv5(NaiveTissueModel_2_rem,'all') 
[TestSolution_4_rem,TestSolutionName_4] = Test4HumanFctExtv5(NaiveTissueModel_4_rem,'all') 

[TestSolutionNaive_Irr,TestSolutionName_Naive_Irr] = Test4HumanFctExtv5(NaiveTissueModelIrr,'all') 
[TestSolutionNaive,TestSolutionName_Naive] = Test4HumanFctExtv5(NaiveTissueModel,'all') 

[TestSolutionEff,TestSolutionName_Eff] = Test4HumanFctExtv5(EffTissueModel,'all') 
[TestSolutionEff_rem,TestSolutionName_Eff_rem] = Test4HumanFctExtv5(EffTissueModel_rem,'all') 

[TestSolutionMem_rem,TestSolutionName_Mem_rem] = Test4HumanFctExtv5(MemTissueModel_rem,'all') 

% Check Model Consistency
load('dataEcoli');

options = 'options_iMAT';
load(['options_methods' filesep options]);
iMAT_model = createTissueSpecificModel(model, options);

options = 'options_GIMME';
load(['options_methods' filesep options]);
GIMME_model = createTissueSpecificModel(model, options);

options = 'options_mCADRE';
load(['options_methods' filesep options]);
mCADRE_model = createTissueSpecificModel(model, options);

modelFileName = 'Recon2.0model.mat';
modelDirectory = getDistributedModelFolder(modelFileName); %Look up the folder for the distributed Models.
modelFileName= [modelDirectory filesep modelFileName]; % Get the full path. Necessary to be sure, that the right model is loaded
model_h = readCbModel(modelFileName);

 




# Old Codes
# Required libraries and models 
## Initializing COBRA Toolbox

```MATLAB
initCobraToolbox(false) %don't update the toolbox
changeCobraSolver ('glpk', 'all');
```

## Reading Recon3D curated model
```MATLAB
% Reading curated recon 3 model
modelFileName = 'Recon3DModel_301.mat';
modelDirectory = getDistributedModelFolder(modelFileName); %Look up the folder for the distributed Models.
modelFileName= [modelDirectory filesep modelFileName]; % Get the full path. Necessary to be sure, that the right model is loaded
model = readCbModel(modelFileName);
```

## Reading Macrophage model
```MATLAB
% Reading marophage model
modelMacroFileName = 'C:\Users\Soroush\cobratoolbox\msb201068-s1.xml'
model_macrophage = readSBML(modelMacroFileName,1000)

%% Biomass function from macrophage model incldues some discrepancies in text (dashed, -, and underscores, _)
% We add macrophage biomass function from the Puniya et al model to avoid discrepancies

model = addReaction(model,'biomass_mac','reactionFormula',char(printRxnFormula(model_macrophage,'biomass_Mphage')))
modelnew = changeObjective(model,'biomass_mac',1.0)
```

## Reading Recon 2 Model
```MATLAB
% reading Recon 2.2 model
modelRecon2FileName = 'C:\Users\Soroush\cobratoolbox\MODEL1603150001_url.xml'
modelRecon2 = readSBML(modelRecon2FileName,1000)
```

## Adding reactions from Recon 2.2 and Macrophage model
```MATLAB
% Adding new reactions
modelnew = addReaction(modelnew,'PRPNCOAHYDm','reactionFormula',char(printRxnFormula(model_macrophage,'PRPNCOAHYDm')))
modelnew = addReaction(modelnew,'PPCOAOm','reactionFormula',char(printRxnFormula(model_macrophage,'PPCOAOm')))
modelnew = addReaction(modelnew,'ARTPLM1','reactionFormula',char(printRxnFormula(model_macrophage,'ARTPLM1')))
modelnew = addReaction(modelnew,'ARTPLM2','reactionFormula',char(printRxnFormula(model_macrophage,'ARTPLM2')))
modelnew = addReaction(modelnew,'PIt2m','reactionFormula',char(printRxnFormula(model_macrophage,'PIt2m')))
modelnew = addReaction(modelnew,'r0941','reactionFormula',char(printRxnFormula(modelRecon2,'r0941')))
modelnew = addReaction(modelnew,'r1464','reactionFormula',char(printRxnFormula(modelRecon2,'r1464')))
modelnew = addReaction(modelnew,'AIRCr_PRASCS','reactionFormula',char(printRxnFormula(modelRecon2,'AIRCr_PRASCS')))
modelnew = addReaction(modelnew,'r0683','reactionFormula',char(printRxnFormula(modelRecon2,'r0683')))
modelnew = addReaction(modelnew,'RTOT_3','reactionFormula',char(printRxnFormula(modelRecon2,'RTOT_3')))
modelnew = addReaction(modelnew,'ARTFR13','reactionFormula',char(printRxnFormula(modelRecon2,'ARTFR13')))
modelnew = addExchangeRxn(modelnew, 'ca2[e]')
modelnew = addExchangeRxn(modelnew, 'cl[e]')
```

## Exchange reactions
```MATLAB
exRxns = printUptakeBound(modelnewfinal)
modeltest = changeRxnBounds(modelnewfinal,modelnewfinal.rxns(exRxns),0,'l')
```
### Naive exchange reactions
```MATLAB
% Following exchange reactions are open to uptake (i.e. lb = -1000)

naiveExRxns = {'EX_ca2[e]','EX_cl[e]','EX_co[e]','EX_co2[e]','EX_fe2[e]','EX_h[e]','EX_h2o[e]','EX_k[e]','EX_na1[e]','EX_nh4[e]',...
               'EX_no2[e]','EX_o2[e]','EX_pi[e]','EX_so4[e]'}
modelnaive = changeRxnBounds(modeltest,naiveExRxns,-1000,'l')

% Following exchange reactions are tightly controlled for uptake (i.e. lb =-1)

naiveUptakeRxns = {'EX_ala_L[e]','EX_arg_L[e]','EX_asn_L[e]','EX_asp_L[e]','EX_cys_L[e]','EX_eicostet[e]','EX_gln_L[e]','EX_glu_L[e]','EX_gly[e]',...
                   'EX_hdca[e]','EX_hdcea[e]','EX_his_L[e]','EX_ile_L[e]','EX_leu_L[e]','EX_lnlc[e]','EX_lnlnca[e]',...
                   'EX_lnlncg[e]','EX_lys_L[e]','EX_met_L[e]','EX_ocdca[e]','EX_ocdcea[e]','EX_orn_D[e]',...
                   'EX_phe_L[e]','EX_pro_L[e]','EX_ribflv[e]','EX_ser_L[e]','EX_thr_L[e]','EX_trp_L[e]','EX_tyr_L[e]',...
                   'EX_val_L[e]','EX_vitd3[e]'}
modelnaive = changeRxnBounds(modelnaive,naiveUptakeRxns,-1,'l')

% Glucose uptake is controlled at lb = -100

modelnaive = changeRxnBounds(modelnaive,{'EX_glc_D[e]'},-100,'l')
```
### Effector exchange reactions
```MATLAB
% Following exchange reactions are open to uptake (i.e. lb = -1000)

EffExRxns = {'EX_ca2[e]','EX_cl[e]','EX_co[e]','EX_co2[e]','EX_fe2[e]','EX_h[e]','EX_h2o[e]','EX_k[e]','EX_na1[e]',...
             'EX_nh4[e]','EX_no2[e]','EX_pi[e]','EX_so4[e]'}
modeleff = changeRxnBounds(modeltest,EffExRxns,-1000,'l')

% Following exchange reactions are tightly controlled for uptake (i.e. lb =-1)

EffUptakeRxns = {'EX_ala_L[e]','EX_arg_L[e]','EX_asn_L[e]','EX_asp_L[e]','EX_cys_L[e]','EX_eicostet[e]','EX_gln_L[e]','EX_glu_L[e]','EX_gly[e]',...
                 'EX_hdca[e]','EX_hdcea[e]','EX_his_L[e]','EX_ile_L[e]','EX_leu_L[e]','EX_lnlc[e]','EX_lnlnca[e]',...
                 'EX_lnlncg[e]','EX_lys_L[e]','EX_met_L[e]','EX_ocdca[e]','EX_ocdcea[e]','EX_orn_D[e]',...
                 'EX_phe_L[e]','EX_pro_L[e]','EX_ribflv[e]','EX_ser_L[e]','EX_thr_L[e]','EX_trp_L[e]','EX_tyr_L[e]',...
                 'EX_val_L[e]','EX_vitd3[e]'}
modeleff = changeRxnBounds(modeleff,EffUptakeRxns,-1,'l')

% Glucose uptake is controlled at lb = -100

modeleff = changeRxnBounds(modeleff,{'EX_glc_D[e]'},-100,'l')

% Since Activated T cells are metbolizing aerobically, O2 uptake is controlled at lb = -20

modeleff = changeRxnBounds(modeleff,{'EX_o2[e]'},-20,'l')
```
### Memory exchange reactions
```MATLAB

% Following exchange reactions are open to uptake (i.e. lb = -1000)

MemExRxns = {'EX_ca2[e]','EX_cl[e]','EX_co[e]','EX_co2[e]','EX_fe2[e]','EX_h[e]','EX_h2o[e]','EX_k[e]','EX_na1[e]','EX_nh4[e]',...
             'EX_no2[e]','EX_o2[e]','EX_pi[e]','EX_so4[e]'}
modelmem = changeRxnBounds(modeltest,MemExRxns,-1000,'l')

% Following exchange reactions are tightly controlled for uptake (i.e. lb =-1)
% Uptake of triglyceride (TG) is also permitted in the memory model regarding its FA oxidation behavior

MemUptakeRxns = {'EX_ala_L[e]','EX_arg_L[e]','EX_asn_L[e]','EX_asp_L[e]','EX_cys_L[e]','EX_eicostet[e]','EX_gln_L[e]','EX_glu_L[e]','EX_gly[e]',...
                 'EX_hdca[e]','EX_hdcea[e]','EX_his_L[e]','EX_ile_L[e]','EX_leu_L[e]','EX_lnlc[e]','EX_lnlnca[e]',...
                 'EX_lnlncg[e]','EX_lys_L[e]','EX_met_L[e]','EX_ocdca[e]','EX_ocdcea[e]','EX_orn_D[e]',...
                 'EX_phe_L[e]','EX_pro_L[e]','EX_ribflv[e]','EX_ser_L[e]','EX_thr_L[e]','EX_trp_L[e]','EX_tyr_L[e]',...
                 'EX_val_L[e]','EX_vitd3[e]','EX_1glyc_hs[e]','EX_ak2lgchol_hs[e]','EX_mag_hs[e]','EX_paf_hs[e]',...
                 'EX_pglyc_hs[e]','EX_tag_hs[e]','EX_HC01444[e]','EX_12dgr120[e]','EX_magarachi_hs[e]',...
                 'EX_maglinl_hs[e]','EX_magole_hs[e]','EX_magpalm_hs[e]','EX_magste_hs[e]','EX_glyc2p[e]','EX_glyc[e]'}

modelmem = changeRxnBounds(modelmem,MemUptakeRxns,-1,'l')

% Glucose uptake is controlled at lb = -100

modelmem = changeRxnBounds(modelmem,{'EX_glc_D[e]'},-100,'l')

% Since Activated T cells are metbolizing aerobically, O2 uptake is controlled at lb = -20

modelmem = changeRxnBounds(modelmem,{'EX_o2[e]'},-20,'l')
```

### Directionality change (Rev to Irrev)
```MATLAB
naiveIrrRxns = {'TTDCPT2','TMNDNCCRNt','TMNDNCCPT2','TMNDNCCPT1','TETTET6CRNt',...
                'TETTET6CPT2','TETTET6CPT1','TETPENT6CRNt','TETPENT6CPT2','TETPENT6CPT1',...
                'TETPENT3CRNt','TETPENT3CPT2','TETPENT3CPT1','STRDNCCRNt','STRDNCCPT2',...
                'STRDNCCPT1','RE0583C','r1401','r1400','r0791','r0735','r0652','r0638','r0636','r0633','r0432',...
                'r0431','r0309','PTDCACRNt','PTDCACRNCPT2','PTDCACRNCPT1','PCRNtm','LNLNCGCRNt',...
                'LNLNCGCPT2','LNLNCGCPT1','LNLNCACRNt','LNLNCACPT2','LNLNCACPT1','LNLCCRNt',...
                'LNLCCPT2','LNLCCPT1','LNELDCCRNt','LNELDCCPT2','LNELDCCPT1','HXCOAx','HXCOAc',...
                'HPDCACRNt','HPDCACRNCPT2','HPDCACRNCPT1','HEXDIACtr','ELAIDCRNt','ELAIDCPT2',...
                'ELAIDCPT1','EICOSTETCRNt','EICOSTETCPT2','EICOSTETCPT1','DLNLCGCRNt',...
                'DLNLCGCPT2','DLNLCGCPT1','DCSPTN1CRNt','DCSPTN1CPT2','DCSPTN1CPT1',...
                'CLPNDCRNt','CLPNDCPT2','CLPNDCPT1','C30CPT1','C226CRNt','C226CPT2','C226CPT1',...
                'C204CRNt','C204CPT2','C204CPT1','C181CRNt','C181CPT2','C181CPT1','C180CPT2',...
                'C180CPT1','C161CRNt','C161CRN2t','C161CPT22','C161CPT2','C161CPT12','C161CPT1',...
                'C160CPT2','C160CPT1','ARACHCPT2','ARACHCPT1','ADRNCRNt','ADRNCPT2','ADRNCPT1'}
% Creating the irreversible model
[modelnaive,matchRevNaive,rev2irrev,irrev2rev] = convertToIrreversible(modelnaive,'sRxns',naiveIrrRxns)
```
