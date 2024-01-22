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

## Adding Reactions from Recon 2.2 and Macrophage model
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
