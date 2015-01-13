function [dataPoint] = plot_laser_sniff_kp(mouse,sess)
za=post_acute_04_kp;
% mouse='KPanesthM72';
% sess='006';
recs={'a'};
iClu=2;

props={'delay','avgSpikes'} ;

symbols={'v','^','d','+'}

gs=figure
set(gs,'Nextplot','add')
h(1)=subplot(211)
set(h(1),'Nextplot','add')
h(2)=subplot(212)
set(h(2),'Nextplot','add')


gv=figure
set(gv,'Nextplot','add')
i(1)=subplot(211)
set(i(1),'Nextplot','add')
i(2)=subplot(212)
set(i(2),'Nextplot','add')


for irec=1:numel(recs)
    rec=recs{irec}
    fn=za.file_names(mouse, sess, rec);
    dp=load(fullfile(fn.fold_pr_sess,sprintf('%s_%s_%s_%d_datapoints.mat',mouse,sess,rec,iClu)));
    dataPoint=dp.dataPoint
    vstims=unique([dataPoint.stimV]);
    dstims=unique([dataPoint.stimDur]);
    estims=unique([dataPoint.stimEner]);
    
    makeplots(symbols{irec});
    
end



%plot a property grouped by V
%prop='avgSpikes';

    function makeplots(symbol)
        for iprop=1:numel(props)
            prop=props{iprop};
            toplot=zeros(length(estims),length(vstims));
            for iv=1:length(vstims)
                v=vstims(iv);
                leg{iv}=[num2str(v) 'V'];
                subset=dataPoint([dataPoint.stimV]==v);
                eselect=[subset.stimEner];
                full=find(ismember(estims,eselect));
                empty=find(~ismember(estims,eselect));
                for ifull=1:numel(full)
                    toplot(full(ifull),iv)=subset(ifull).(prop);
                end
                for iempty=1:numel(empty)
                    toplot(empty(iempty),iv)=nan;
                end
            end
            
            figure(gs)
            %subplot(2,ceil(numel(props)/2),iprop)
            hold on
            set(h(iprop),'Nextplot','add')
            plot(h(iprop),estims,toplot,symbol, 'MarkerSize',10,'LineWidth',2,'MarkerFaceColor','k')
            xlabel(h(iprop),'energy (mu J)');
            ylabel(h(iprop),sprintf('%s',prop));
            legend(h(iprop),leg)
            
        end
        suptitle(['Mouse ' mouse ', rec ' sess rec '. Units group ' num2str(iClu)]);
        
%         for iprop=1:numel(props)
%             prop=props{iprop};
%             toplot=zeros(length(dstims),length(vstims));
%             for iv=1:length(vstims)
%                 v=vstims(iv);
%                 leg{iv}=[num2str(v) 'V'];
%                 subset=dataPoint([dataPoint.stimV]==v);
%                 dselect=[subset.stimDur];
%                 full=find(ismember(dstims,dselect));
%                 empty=find(~ismember(dstims,dselect));
%                 for ifull=1:numel(full)
%                     toplot(full(ifull),iv)=subset(ifull).(prop);
%                 end
%                 for iempty=1:numel(empty)
%                     toplot(empty(iempty),iv)=nan;
%                 end
%             end
%             
%             figure(gs)
%             %subplot(2,ceil(numel(props)/2),iprop)
%             hold on
%             set(h(iprop),'Nextplot','add')
%             plot(h(iprop),dstims,toplot,symbol, 'MarkerSize',10,'LineWidth',2,'MarkerFaceColor','k')
%             xlabel(h(iprop),'duration (ms)');
%             ylabel(h(iprop),sprintf('%s',prop));
%             legend(h(iprop),leg)
%         end
%         suptitle(['Mouse ' mouse ', rec ' sess rec '. Units group ' num2str(iClu)]);
        
        %plot a property grouped by duration
        for iprop=1:numel(props)
            prop=props{iprop};
            toplot=transpose(zeros(length(dstims),length(vstims)));
            for id=1:length(dstims)
                d=dstims(id);
                leg2{id}=[num2str(d) 'ms'];
                subset=dataPoint([dataPoint.stimDur]==d);
                vselect=[subset.stimV];
                full=find(ismember(vstims,vselect));
                empty=find(~ismember(vstims,vselect));
                for ifull=1:numel(full)
                    toplot(full(ifull),id)=subset(ifull).(prop);
                end
                for iempty=1:numel(empty)
                    toplot(empty(iempty),id)=nan;
                end
            end
            
            figure(gv)
            %subplot(2,ceil(numel(props)/2),iprop)
            hold on
            set(i(iprop),'Nextplot','add')
            plot(i(iprop),vstims,toplot,symbol, 'MarkerSize',10,'LineWidth',2,'MarkerFaceColor','k')
            xlabel(i(iprop),'V stim (V)');
            ylabel(i(iprop),sprintf('%s',prop));
            legend(i(iprop),leg2)
        end
        suptitle(['Mouse ' mouse ', rec ' sess rec '. Units group ' num2str(iClu)]);
        
     end


end

