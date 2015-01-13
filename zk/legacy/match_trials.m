function [eMatched, inMatched]=match_trials(eSyncRun,in)
% estimation synchronization error for each trial shift value
bSyncRun=[trial(in).runTrialStart];

trial_shift = -10:10; 
nts = numel(trial_shift);
error = zeros(1,nts);
alpha = zeros(1,nts);
err_min = 1e10;

for it = 1:nts
    % defining the shifted trigger times
    %positive shift is when the first n trials of the beh were
    %missed
    %negative is when triggers dont have behavior and have to be skipped
    
    if trial_shift(it) >= 0
        %positive shift, first it trials were not recorded in the ephys
        %file
        bTimes = bSyncRun(:,(1+trial_shift(it)):end);
        inShifted= in(1+trial_shift(it):end);
        eTimes = eSyncRun;
    else
        %negative shift, first it trials were not recorded in the h5 file
        bTimes = bSyncRun;
        eTimes = eSyncRun(:,(1+trial_shift(it)):end);
    end
    
    ne = length(eTimes);
    nb = length(bTimes);
    nn = min([ne, nb]);
    
    %  and alinging triggers by the first trial
    eT0=eTimes(1,1); % the t0 of the aligned ephys, segment, relative to beginning of the run
    bT0=bTimes(1,1);
    eTimes = eTimes(:,1:nn) - eTimes(1,1);
    bTimes = bTimes(:,1:nn) - bTimes(1,1);
    ind    = 1:nn;
    alpha(it) = (eTimes(1,:)/bTimes(1,:));
    error(it) = sqrt(norm(eTimes-alpha(it)*bTimes));
    if strcmpi(figures,'plot')
        figure(sf(irec))
        subplot(nrun,3,3*irun)
        plot(trial_shift(it), error(it), 'o', 'Color', 'r', 'LineWidth', 2)
        xlabel('shift');
        ylabel('error sum');
        hold on
    end
    %
    if error(it)<err_min
        err_min = error(it);
        imin    = it;
        if strcmpi(figures,'plot')
            figure(sf(irec))
            subplot(nrun,3,3*irun-2)
            plot(eTimes,alpha(it)*bTimes,'x')
            xlabel('ephys trig');
            ylabel('behav trig');
            title(trial_shift(it))
            subplot(nrun,3,3*irun-1)
            errorVec=sum(eTimes-alpha(it)*bTimes,1);
            plot(errorVec, 'x')
            xlabel('trial');
            ylabel('error');
            figure(runFig(irun))
            clf
            plot(eTimes,ones(size(eTimes))*double(eThresh),'r*')
            hold on
            plot_if(bTimes*(alpha(it)),ones(size(eTimes))*double(eThresh)*1.2,'mx')
            legend('beh Trigs');
            plot(eTRun(eT0:end))
            title(sprintf('Shift estim + align, rec %s run %02d',info.rec(nrec).name, irun), 'FontSize', 10, 'FontWeight', 'bold')       
        end
        trAlign.eT0Run = eT0; %start of the aligned ephys rel. to the run.
        eMatched = eTimes+eT0; %ephys trig times rel to eT0
        inMatched = inShifted(ind); %beh trials that have a matching ephys trial.
    end
    drawnow
end

if strcmpi(figures,'plot')
    figure(sf(irec))
    title(sprintf('Shift estim + align, rec %s run %02d',info.rec(nrec).name, irun), 'FontSize', 10, 'FontWeight', 'bold')
end

end