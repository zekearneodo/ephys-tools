function [trialsBase] = assemble_baseline_trialStruct(mouse,sess,rec)
fn=file_names(mouse,sess,rec);
% load trial structure
q = load(fn.trial);
trial=q.trial; clear q

% load fields from rsm file
q = load(fn.rsm_data);
Sniff = q.Sniff;
FVpin = q.FVpin;
clear q

trialsBase=struct;
otr=0;
for itr=1:numel(trial)
    tr=trial(itr);
    
    if ~isempty(tr.start) && ~isempty(tr.sniffParabZeroTimes)
        sniffZeros_toUse = round(tr.sniffParabZeroTimes(1:2,:));
        
        % make trial unique id
        trialUId = [fn.basename_an 'trial' num2str(tr.start)];
        
        if diff([tr.odorTimes])>0 && tr.odorConc>0
            
            % get 3000ms of sniff data pre FV onset
            FV_recTime = round(tr.odorTimes(1)*1.009) + tr.start;
            
            t2 = FV_recTime;
            t1 = t2-3000;
            if t1<1
                continue
            end
            %         itr
            %         if itr~=1 && t1<(trial(itr-1).start+trial(itr-1).runTrialDur)
            %             t1=t2-1000;
            %         end
            
            FV_rsmCandidates =  max(1,tr.start-2e4) + find(diff(FVpin(max(1,tr.start-2e4):min(tr.start+2e4,end)))>1e4);
            %                 FV_rsmCandidates =  (max(1,tr.start-2e4) + find(diff(FVpin(max(1,tr.start-2e4):tr.start+2e4))>1e4)+1);
            [~,FVtime_idx] =  min(abs(FV_rsmCandidates-FV_recTime));
            FV_rsmTime = FV_rsmCandidates(FVtime_idx);
            
            ms_tr2rsm = diff([FV_recTime FV_rsmTime]);
            if abs(ms_tr2rsm)>3
%                 error('Couldnt find close enough FVon time in rsm_data')
                continue
            else
                otr=otr+1;
                % Create sniff phase vector
                vectorizedSniff = reshape(sniffZeros_toUse+tr.start,1,2*size(sniffZeros_toUse,2));
                vectorizedSniff = vectorizedSniff(vectorizedSniff>t1&vectorizedSniff<t2)-t1;
                
                inhOnsets = sniffZeros_toUse(1,:)+tr.start > t1 & sniffZeros_toUse(1,:)+tr.start < t2;
                pauseOnsets = sniffZeros_toUse(2,:)+tr.start > t1 & sniffZeros_toUse(2,:)+tr.start < t2;
                
                phase_inRange = [inhOnsets; -1.*pauseOnsets];
                vectorizedPhase = phase_inRange(phase_inRange~=0)';
                
                if numel(vectorizedSniff)<5
                    warning('Skipped a trial with few sniffs detected')
                    trialsBase(otr).trialId = trialUId;
                    trialsBase(otr).start = t1;
                    trialsBase(otr).sniffFlow = Sniff(t1+1+ms_tr2rsm:t2+ms_tr2rsm);
                    trialsBase(otr).sniffPhase = nan;
                else
                    
                    if size(vectorizedSniff)==size(vectorizedPhase)
                        vectorizedPhase = [vectorizedPhase(2) vectorizedPhase];
                        vectorizedSniff = [1 vectorizedSniff t2-t1];
                    else
                        error('sniff vectors not same size')
                    end
                    if ~any(diff(vectorizedSniff)>1200) && vectorizedSniff(1)<1200
                        sniffPhase = nan(1,t2-t1);
                        for ib = 1:length(vectorizedSniff)-1
                            sniffPhase(vectorizedSniff(ib):vectorizedSniff(ib+1)) = vectorizedPhase(ib);
                        end
                        
                        ii = find(tr.sniffZeroTimes(1,:)>tr.odorTimes(1)*1.009,1);
                        figure(1); clf; hold on
%                         plot(1:t2-t1,Sniff(t1+ms_tr2rsm:t2+ms_tr2rsm-1))
%                         plot(vectorizedSniff,zeros(size(vectorizedSniff)),'r*')
                        %         fill([t1;t2]*ones(1,2),ones(2,1)*[-3e4 3e4],'-','LineWidth',100,'Color',[1 0.9 0.95])
                        fill([t1 t1 t2 t2],[-3e4 3e4 3e4 -3e4],[1 0.95 0.97],'EdgeColor','none')
                        plot(t1+ms_tr2rsm:t2+ms_tr2rsm+1e3,Sniff(t1+ms_tr2rsm:t2+ms_tr2rsm+1e3))
                        plot((tr.sniffZeroTimes(1,ii)+tr.start)*ones(1,2),[-3e4 3e4],':g','LineWidth',3)
                        plot((t1+1:t2),-1e4.*sniffPhase,'k')
                        
                        trialsBase(otr).trialId = trialUId;
                        trialsBase(otr).start = t1;
                        trialsBase(otr).sniffFlow = Sniff(t1+1+ms_tr2rsm:t2+ms_tr2rsm);
                        trialsBase(otr).sniffPhase = sniffPhase;
                    else
                        trialsBase(otr).trialId = trialUId;
                        trialsBase(otr).start = t1;
                        trialsBase(otr).sniffFlow = Sniff(t1+1+ms_tr2rsm:t2+ms_tr2rsm);
                        trialsBase(otr).sniffPhase = nan;
                        
                        figure(1); clf; hold on
                        plot(1:t2-t1,Sniff(t1+ms_tr2rsm:t2+ms_tr2rsm-1))
                        plot(vectorizedSniff,zeros(size(vectorizedSniff)),'r*')
                        warning('Sniff_analysis may have missed a cycle; check waveform')
                    end
                end
            end
        end
    end
end
save(fullfile(fn.base_folder,'data_Neil',sprintf('%strialsBase.mat',fn.basename_an)), 'trialsBase')
end