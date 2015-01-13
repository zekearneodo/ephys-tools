% KlustaSuite V2 default parameters
%Filtering Pars:
classdef klustapar < handle
    properties
        experiment_name = '';
        raw_data_files  = [];
        prb_file        = '';
        nbits           = 16;
        voltage_gain    = 10;
        sample_rate     = 20000;
        
        filter_high         = 0.95;
        filter_butter_order = 3;
        filter_lfp_low      = 0;  % LFP filter low-pass frequency
        filter_lfp_high     = 300;  % LFP filter high-pass frequency
    
        chunk_size    = 1 ;% in seconds
        chunk_overlap = 0.015; %in seconds
        
        %spike detection
        nexcerpts = 50;
        excertp_size = 1.;
        threshold_strong_std_factor = 4.5;
        threshold_weak_std_factor = 2.;
        detect_spikes = 'negative';
        
        %connected_component
        connected_component_join_size = 0.00005;
        % Spike extraction

        extract_s_before = 16;
        extract_s_after = 16;
        waveforms_nsamples = extract_s_before + extract_s_after;
        
        % Features
        nfeatures_per_channel = 3  ;% Number of features per channel.
        pca_nwaveforms_max = 10000 ;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % KlustaKwik parameters %
        %%%%%%%%%%%%%%%%%%%%%%%%%
        MaskStarts = 100;
        %MinClusters = 100
        %MaxClusters = 110
        MaxPossibleClusters =  500;
        FullStepEvery =  10;
        MaxIter = 10000;
        RandomSeed =  654;
        Debug = 0;
        SplitFirst = 20;
        SplitEvery = 100;
        PenaltyK = 0;
        PenaltyKLogN = 1;
        Subset = 1;
        PriorPoint = 1;
        SaveSorted = 0;
        SaveCovarianceMeans = 0;
        UseMaskedInitialConditions = 1;
        AssignToFirstClosestMask = 1;
        UseDistributional = 1;
    end

end