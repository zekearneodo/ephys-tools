% sniff_stuff
% created by KP, 2014-05-02
% cript for visualizing various factors of sniff from a sess/rec
% uses data_managment_tools_032

function [pr_good_sniff1,snf] = sniff_stuff(mouse,sess,rec)
global dm sn;

% add the package of common files and folders to the path.
% /baseFolder/ephysDataManagement/current/include
includePath=fullfile(fileparts(pwd),'current','include');
addpath(includePath);
% Functions of include that it uses:
%   - file_names(mouse,sess,rec,stat)

dm = data_management_tools_032();
% sn.good_sniff1 = @good_sniff1;

fn = dm.file_names(mouse, sess, rec);
q = load(fn.trial);
sn.tr = q.trial;
q = load(fn.sniffs);
sn.sniff = q.sniff;

pr_good_sniff1 = good_sniff1(sn);



%% Find times each trial from FV on to first inhalation.
% Calculate probability that the first inhalation arrives after PID signal
% asymptotes and before odor presentation ends.
    function pr_good_sniff1 = good_sniff1(sn)
        tr = sn.tr;
        ig=0; CHECK = []; 
        for ikt = 1:length(tr)
            trial = tr(ikt);
            
            % Only look at odor trials with inhalation during stim
            if trial.odorConc <= 0 || isempty(trial.sniffZeroTimes)
                continue
            end
            
            FV_on = trial.odorTimes(1);
            inhalations = trial.sniffZeroTimes(1,:);
            
            i_first_inh = find(inhalations > FV_on, 1);
            
            if isempty(i_first_inh)
                CHECK = [CHECK ikt];
                sprintf('No positive sniffs recorded in duration of trial %i', ikt)
                continue
            end
            
            ig = ig+1;
            grace_period(ig) = inhalations(i_first_inh) - FV_on;
            
            if grace_period(ig) > 3000
                CHECK = [CHECK ikt];    % trials that probably had bad sniff during odor presentation
            end
            
        end
        
        CHECK = find(grace_period > 1050);  % trials that had bad sniff during odor presentation
        
        inhwin = 1100;
        odor_dur = diff(trial.odorTimes);
        rise_time = 100;
        
        figure;
        cdfplot(grace_period)
        xlim([0 inhwin])
        xlabel('Time (ms) from FV on to first inhalation')
        ylabel('Cumulative probability of occurrence')
        
        %     [h,p] = lillietest(grace_period)   % test for normality
        
        [fct,x] = ecdf(grace_period);
        pr_good_sniff1 = fct(find(x>odor_dur,1)) - fct(find(x>rise_time,1));
        
        sprintf('There are %i "good" odor trials.', ig)
        
        plot_sniff(CHECK);
    end


function plot_sniff(CHECK)
    snf = sn.sniff;
    tr = sn.tr;

    for ikt = [CHECK CHECK-50]

        FV_on = tr(ikt).odorTimes(1);

        nsn = find([sn.sniff.t0]>(tr(ikt).start+FV_on),1);

        trst_sn_t = tr(ikt).start - [sn.sniff(nsn).t0];
        FV_sn_t = trst_sn_t + FV_on;

        figure; hold on
        plot(snf(nsn).waveform)
            t0_diff = [snf(nsn-1).t0]-[snf(nsn).t0];
            x_presnf = t0_diff+[1:length(snf(nsn-1).waveform)]
        plot(x_presnf,snf(nsn-1).waveform)
        plot([snf(nsn-1).t_zer]+t0_diff,[0 0 0],'k*')
        plot([snf(nsn).t_zer],[0 0 0],'k*')
        plot([snf(nsn).t_min],[snf(nsn).y_min],'ro')
        plot([snf(nsn).t_max],[snf(nsn).y_max],'ro')
        plot(FV_sn_t,0,'gd','MarkerSize',25)
        tit = sprintf('sniff number %i, trial number %i', nsn, ikt);
        title(tit)
        hold off
    end    
end
end

    

