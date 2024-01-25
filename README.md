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
naiveExRxns = {'EX_ca2[e]','EX_cl[e]','EX_co[e]','EX_co2[e]','EX_fe2[e]','EX_h[e]','EX_h2o[e]','EX_k[e]','EX_na1[e]','EX_nh4[e]',...
               'EX_no2[e]','EX_o2[e]','EX_pi[e]','EX_so4[e]','EX_glc_D[e]','EX_ala_L[e]','EX_arg_L[e]',...
               'EX_asn_L[e]','EX_asp_L[e]','EX_cys_L[e]','EX_eicostet[e]','EX_gln_L[e]','EX_glu_L[e]','EX_gly[e]',...
               'EX_hdca[e]','EX_hdcea[e]','EX_his_L[e]','EX_ile_L[e]','EX_leu_L[e]','EX_lnlc[e]','EX_lnlnca[e]',...
               'EX_lnlncg[e]','EX_lys_L[e]','EX_met_L[e]','EX_ocdca[e]','EX_ocdcea[e]','EX_orn_D[e]',...
               'EX_phe_L[e]','EX_pro_L[e]','EX_ribflv[e]','EX_ser_L[e]','EX_thr_L[e]','EX_trp_L[e]','EX_tyr_L[e]',...
               'EX_val_L[e]','EX_vitd3[e]'}
modelnaive = changeRxnBounds(modeltest,naiveExRxns,-1000,'l')
```
### Effector exchange reactions
```MATLAB
EffExRxns = {'EX_ca2[e]','EX_cl[e]','EX_co[e]','EX_co2[e]','EX_fe2[e]','EX_h[e]','EX_h2o[e]','EX_k[e]','EX_na1[e]',...
             'EX_nh4[e]','EX_no2[e]','EX_pi[e]','EX_so4[e]','EX_glc_D[e]','EX_o2[e]','EX_ala_L[e]',...
             'EX_arg_L[e]','EX_asn_L[e]','EX_asp_L[e]','EX_cys_L[e]','EX_gln_L[e]','EX_glu_L[e]','EX_gly[e]',...
             'EX_his_L[e]','EX_ile_L[e]','EX_leu_L[e]','EX_lys_L[e]','EX_met_L[e]','EX_orn[e]','EX_orn_D[e]',...
             'EX_phe_L[e]','EX_pro_L[e]','EX_ribflv[e]','EX_ser_L[e]','EX_thr_L[e]','EX_trp_L[e]','EX_tyr_L[e]',...
             'EX_val_L[e]','EX_vitd3[e]'}
modeleff = changeRxnBounds(modeltest,EffExRxns,-1000,'l')
```
