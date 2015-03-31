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
        nchannels       = 32;
        
        filter_high_sec     = 0.95;
        filter_butter_order = 3;
        %filter_lfp_low      = 0;  % LFP filter low-pass frequency
        %filter_lfp_high     = 300;  % LFP filter high-pass frequency
    
        chunk_size_sec    = 1 ;% in seconds
        chunk_overlap_sec = 0.015; %in seconds
        
        %saving raw/filtered data
        save_raw  = true
        save_high = true
        save_low  = true
        
        %spike detection
        nexcerpts = 50;
        excerpt_size_sec = 1.;
        threshold_strong_std_factor = 4.;
        threshold_weak_std_factor = 2.;
        weight_power = 2.;
        detect_spikes = 'negative';
        
        %connected_component
        connected_component_join_size_sec = 0.00005;
        % Spike extraction

        extract_s_before = 10;
        extract_s_after = 10;
        waveforms_nsamples = 20;
        
        % Features
        nfeatures_per_channel = 3  ;% Number of features per channel.
        pca_nwaveforms_max = 10000 ;
        features_contiguous = true
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % KlustaKwik parameters %
        %%%%%%%%%%%%%%%%%%%%%%%%%
        KK_MaskStarts = 500;
        KK_MinClusters = 100
        KK_MaxClusters = 110
        KK_MaxPossibleClusters =  1000;
        KK_FullStepEvery =  10;
        KK_MaxIter = 10000;
        KK_RandomSeed =  654;
        KK_Debug = 0;
        KK_SplitFirst = 20;
        KK_SplitEvery = 40;
        KK_PenaltyK = 0;
        KK_PenaltyKLogN = 1;
        KK_Subset = 1;
        KK_PriorPoint = 1;
        KK_SaveSorted = 0;
        KK_SaveCovarianceMeans = 0;
        KK_UseMaskedInitialConditions = 1;
        KK_AssignToFirstClosestMask = 1;
        KK_UseDistributional = 1;
        KK_RamLimitGB = 120;
    end
    
    methods
        
        function obj = klustapar(quiet)
            % if quiet is 'quiet', don't display the defaults
            if nargin > 0 && strcmpi('quiet',quiet)
                return
            else
                disp(obj);
            end
        end
        
        function result = set_value(obj,parName,parValue)
            if ~isprop(obj,parName)
                error('invalid parameter %s',parName)
                result = 1;
            elseif ~isempty(parValue) %don't let it set empty values
                obj.(parName) = parValue;
                result = 0;
            else
                warning('attempted to set empty value for property %s, left default',parName);
                result = 2;
            end
        end
        
        function fullObj = compute_pars(obj)
            % compute all the values that depend on other parameters
            fullObj = obj;
            sRate = obj.sample_rate;
            fullObj.filter_high = obj.filter_high_sec * 0.5 * sRate;
            fullObj.chunk_size  = round(obj.chunk_size_sec * sRate);
            fullObj.chunk_overlap  = round(obj.chunk_overlap_sec * sRate);
            fullObj.excerpt_size  = round(obj.excerpt_size_sec * sRate);
            fullObj.connected_component_join_size  = round(obj.connected_component_join_size_sec * sRate);
            fullObj.waveforms_nsamples = obj.extract_s_before + obj.extract_s_after;
        end
        
        function result = make_par_file(obj,fileName)
            fullPar = compute_pars(obj);
            pars = fields(fullPar);
            fid = fopen(fileName,'w');
            fprintf(fid,'## parameter files for version 2 of klustasuite generated by klustapar.m\n\n');
            %parameters:
            % char  - > print within quotation
            % num   - > print as num
            % logic - > print the string of the value
            for ip = 1:numel(pars)
                par = pars{ip};
                val = fullPar.(par);
                if isnumeric(val)
                    fprintf(fid,'%s = %d\n',par,val);
                elseif ischar(val)
                    fprintf(fid,'%s = ''%s''\n',par,val);
                elseif islogical(val)
                    if true(val)
                        sOut = 'True';
                    else
                        sOut = 'False';
                    end
                    fprintf(fid,'%s = %s\n',par,sOut);
                end
            end
            result = fclose(fid);
            if ~result
                fprintf('Wrote klustasuite v2 parameters in file %s\n',fileName);
            end
        end
    end

end