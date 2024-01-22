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
