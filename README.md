# TcellCOBRA
MATLAB codes for reconstruction of T cell metabolic models based on the multi-omics data

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
model = addReaction(model,'biomass_mac','reactionFormula',char(printRxnFormula(model_macrophage,'biomass_mac')))
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
