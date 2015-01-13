%makes the xml-ready parameter structures for the ndm plugins

%They all have the fields:
% -programName
% -parameters
% -help

%This scripts goes trhough an array of structures with all the fields for
%each program
programPars=struct('name',{},'parametersList',{},'help',{});
parameterNames= {'name',            'value', 'status'} ;


% PROGRAM PARAMETERS
%%% NDM_CLEAN
np=numel(programPars);
programPars(np+1).name   = 'ndm_clean';
programPars(np+1).parametersList={
                {'wideband',  'true' , 'Optional'}
                {'xml', 'true'    , 'Optional' }
                {'spots', 'true'    , 'Optional' }
                {'pos', 'true'    , 'Optional' }
                {'hipass', 'true'    , 'Optional' }
                };
programPars(np+1).help       = 'Clean intermediate files after pre-processing is complete.' ;     

%%% NDM_CONCATENATE
np=numel(programPars);
programPars(np+1).name   = 'ndm_concatenate';
programPars(np+1).parametersList={
                {'spotsSamplingRate',  '' , 'Mandatory'}
                };
programPars(np+1).help       = 'Concatenate all session files (.dat, .pos and .evt) recorded on the same day.' ; 


%%% NDM_EXTRACTSPIKES
np=numel(programPars);
programPars(np+1).name   = 'ndm_extractspikes';
programPars(np+1).parametersList={
                {'thresholdFactor',  1.725  , 'Mandatory'}
                {'refractoryPeriod', 16    , 'Mandatory'}
                {'peakSearchLength', 32    , 'Mandatory'}
                {'start',            0     , 'Mandatory'}
                {'duration'        , 60    , 'Mandatory'}
                };
programPars(np+1).help       = 'Extract spikes from high-pass filtered .fil file (this creates .res and .spk files).' ;     

%%% NDM_HIPASS
np=numel(programPars);
programPars(np+1).name   = 'ndm_hipass';
programPars(np+1).parametersList={
                {'windowHalfLength', 10  , 'Mandatory'}
                };
programPars(np+1).help       = 'High-pass filter a .dat file (required for spike extraction).' ;     

%%% NDM_LFP
np=numel(programPars);
programPars(np+1).name   = 'ndm_lfp';
programPars(np+1).parametersList={
                {'samplingRate', 1250  , 'Mandatory'}
                };
programPars(np+1).help       = 'Downsample a .dat file to create the corresponding LFP file.' ;    

%%% NDM_PCA
np=numel(programPars);
programPars(np+1).name   = 'ndm_pca';
programPars(np+1).parametersList={
                {'before'          , ''      , 'Mandatory'}
                {'after'           , ''      , 'Mandatory'}
                {'extra'           , 'false' , 'Mandatory'}
                };
programPars(np+1).help       = 'Compute principal component analysis (PCA).' ; 

%%% NDM_START
np=numel(programPars);
programPars(np+1).name   = 'ndm_start';
programPars(np+1).parametersList={
                {'suffixes',  ''        , 'Mandatory'}
                {'wideband', 'false'    , 'Mandatory'}
                {'events'  , 'false'    , 'Mandatory'}
                };
programPars(np+1).help       = 'Perform all processing steps for a multiple sets of multiple-session recordings' ;  


% go through the program array and make and store the par structure for
% each.
fn=file_names;
parametersFolder = fn.ndm_def_par;
for ipr=1:numel(programPars)
    %make parameters subStruct
    parameters=struct;
    parametersList=programPars(ipr).parametersList;
    for ip=1:numel(parametersList)
        for ik=1:numel(parameterNames)
            parameters(ip).parameter.(parameterNames{ik})=parametersList{ip}{ik};
        end
    end
   %make and save the structure
   program.name=programPars(ipr).name;
   program.parameters=parameters;
   program.help=programPars(ipr).help;
   progStrFname=fullfile(parametersFolder,['par_' program.name '.mat'])
   params.program=program;
   save(progStrFname,'program');
end