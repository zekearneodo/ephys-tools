function rp = result_plotting_tools_kp()
% global rp
rp=struct;

% if this is not the computing server, add the include folder to the path
[~,computerName]=system('hostname');
if ~strcmp(strtrim(computerName),'flipper')
    includePath=fullfile(fileparts(pwd),'current','include');
    addpath(includePath);
end
fn=file_names();
%the parameters common to all the functions
% feat: features to plot:
rp.featList={'latency','nSpikes','jitterISI','latencyISI'};
rp.featvarList={'jitterISI','none'};
rp.stimList={'odor','laser'};
rp.option={};
rp.option=[rp.option struct('parName','vrFun'     ,'parDefault','visualize_responses'                              ,'validation',@ischar)];
rp.option=[rp.option struct('parName','vrFunPath' ,'parDefault',fullfile(get_local_disk,'ephysDataManagement','zk'),'validation',@ischar)];
rp.option=[rp.option struct('parName','odorDbFile','parDefault','odorDb.mat'                                       ,'validation',@ischar)];
rp.option=[rp.option struct('parName','odorDbPath','parDefault',fn.fold_st_meta                                    ,'validation',@ischar)];
rp.option=[rp.option struct('parName','odorDbType','parDefault','struct'                                           ,'validation',@ischar)];
rp.option=[rp.option struct('parName','feature'   ,'parDefault','latencyISI'                                       ,'validation',@(x) ischar(x) && any(strcmpi(x,rp.featList)))];
rp.option=[rp.option struct('parName','featureVar','parDefault','jitterISI'                                        ,'validation',@(x) ischar(x) && any(strcmpi(x,rp.featvarList)))];
rp.option=[rp.option struct('parName','onlyLitral','parDefault',true                                               ,'validation',@islogical)];
rp.option=[rp.option struct('parName','forceVR'   ,'parDefault',false                                              ,'validation',@islogical)];
rp.option=[rp.option struct('parName','quality'   ,'parDefault',[]                                                 ,'validation',@isnumeric)];
rp.option=[rp.option struct('parName','concAxis'  ,'parDefault',[-4.5 0.5]                                         ,'validation',@(x) isnumeric(x) && numel(x)==2)];
rp.option=[rp.option struct('parName','latAxis'   ,'parDefault',[0 120]                                            ,'validation',@(x) isnumeric(x) && numel(x)==2)];
rp.option=[rp.option struct('parName','legends'   ,'parDefault',true                                               ,'validation',@islogical)];
rp.option=[rp.option struct('parName','LatvJit'   ,'parDefault',false                                              ,'validation',@islogical)];
rp.option=[rp.option struct('parName','specConc_uM','parDefault',[]                                                ,'validation',@(x) isnumeric(x) && numel(x)==1)];
rp.opl=0;
% input parameters from functions
% rp.option=[rp.option struct('parName','mouse',  'parDefault','',  'validation',@(x)iscell(x) || ischar(x)) ];
% rp.option=[rp.option struct('parName','rec',  'parDefault','',  'validation',@(x)iscell(x) || ischar(x)) ];
% rp.option=[rp.option struct('parName','allTogether'     ,'parDefault','false','validation',@islogical)];
% rp.option=[rp.option struct('parName','subplotHandle'   ,'parDefault',''     ,'validation',@ishandle)];
% 

%the functions  
    rp.plot_cell          = @plot_cell;
    rp.plot_all_cells     = @plot_all_cells;
    rp.plot_everything    = @plot_everything;
    rp.odor_lookup        = @odor_lookup;
    rp.cell_id            = @cell_id;
    rp.subplot_fig        = @subplot_fig;
    rp.plot_odors         = @plot_odors;
end


function [pev] = plot_everything(varargin)
% Plots response parameter vs stimulus parameter 
%  for all cells that have quality>0 in the unit database.
global gotUnits;

% Define input parameters that are specific to this function
inFunOption={};
inFunOption=[inFunOption struct('parName','mouse',  'parDefault','',  'validation',@(x) iscell(x) || ischar(x)) ];

% Parse input arguments
[iP,varargout] = parseInputs(inFunOption,varargin);
% clear inFunOption

%% Find all units with quality>0
% Now the function opens the cell database (/basefolder/Experiment/units)
% and goes through each mat file and keeps the quality units.

% If a list of mice was entered, it will go only through that list; otherwise all.
mouse=iP.mouse;
if isempty(mouse)
    fprintf('\t no mouse specified, going for all of them...\n');
    mouseList = {''};
else
    if iscell(mouse)
        mouseList =mouse;
    else
        mouseList={mouse};
    end
end

fn=file_names();
nCells=0;
for iM=1:numel(mouseList)
    mouse=mouseList{iM};
    moreCells=[];
    cellBaseName=[mouse '*.mat'];
    cellsList=dir(fullfile(fn.fold_unit_db,cellBaseName));
        %function getCell
    moreCells=cellfun(@(x) getCell(fullfile(fn.fold_unit_db,x)),{cellsList.name},'UniformOutput',false);
    moreCells(cellfun('isempty',moreCells))=[];
    cells(nCells+1:nCells+numel(moreCells))=[moreCells{:}];
    nCells=nCells+numel(moreCells);   
end
gotUnits=unique({cells.Id});

% Now sort the cells by mouse and sess to prepare for plotting.
% For each sess, combining across recs, get all the cells that have the
% same quality. Each of these will get a subplot.
mice=unique({cells.mouse});
nSubPlots=0;
for mouse=mice
    cellsMouse=cells(strcmpi(mouse,{cells.mouse}));
    sessions=unique({cellsMouse.sess});
    for sess=sessions
        cellsSess=cellsMouse(strcmpi(sess,{cellsMouse.sess}));
        nSubPlots=nSubPlots+max([cellsSess.quality]);
    end
end
nSubPlots=double(nSubPlots)+1; %add one subplot for the legend

% Begin creating plot
gev=figure(1);
set(gev,'Position',[100, 100, 1600,1200],...
        'Nextplot','add');
figRows=2;  figCols=ceil(nSubPlots/figRows);
for iSubFig=1:nSubPlots
    hS(iSubFig)=subplot(figRows,figCols,iSubFig);
    set(hS(iSubFig),'Nextplot','add');
end

nPlots=0; %count of subplots as filled in
for mouse=mice
    cellsMouse=cells(strcmpi(mouse,{cells.mouse}));
    sessions=unique({cellsMouse.sess});
    for sess=sessions
        cellsSess=cellsMouse(strcmpi(sess,{cellsMouse.sess}));
        %get the cells of the session
        cellsQual=unique([cellsSess.quality]);
        for qCell=cellsQual
            nPlots=nPlots+1;
            % ~~~~~~~~~~~~~~~~~~~~~~~~~
            pev(nPlots)=plot_all_cells(char(mouse),char(sess),varargout{:},'allTogether',true,'quality',qCell,'subplotHandle',hS(nPlots));
            % ~~~~~~~~~~~~~~~~~~~~~~~~~
        end
    end
end

% Create legend in the extra subplot
if strcmpi(iP.odorDbType,'struct')
    q=load(fullfile(iP.odorDbPath,iP.odorDbFile));
    odorDb=q.odor;
    clear q;
else
    error('I still dont know how to handle a %s type odor database',ip.odorDbType);
end
plottedOId=[];
plottedOId=unique([pev.plottedOId]);
for io=1:numel(plottedOId)
    odorN    = plottedOId(io);
    pLeg{io} = odorDb(io).name{1};
    pHan(io) = plot(hS(end),0,0, ...
                [odorDb(odorN).plotSymbol],'MarkerSize',2,'LineWidth',0.1,'Color',odorDb(odorN).colorCode,...
                'MarkerFaceColor',odorDb(odorN).colorCode);
end
axis off
hLeg=legend(hS(end),pLeg,'location','NorthEast');
set(hLeg,'color','none','box','off','units','pixels')
legText=findobj(hLeg,'type','text');
set(legText,'FontSize',12);

% Save figure
if iP.onlyLitral
    litral='litral';
else
    litral='mitral';
end
fileBase=sprintf('all_%s_%s_%s',litral,iP.stimType,iP.feature);
% print(gev,'-depsc',fullfile(fn.fold_unit_db,sprintf('%s.eps',fileBase)));





    %__________________________________________________%
    function theCell=getCell(unit_filename)
        theCell=load(unit_filename);
        if ~isfield(theCell,'quality') || theCell.quality==0 || theCell.light==0
            theCell='';
        elseif theCell.quality>0
            theCell.resp = get_resp_struct(unit_filename,iP.stimType);
        end
    end
end


function [pac] = plot_all_cells(mouse,sess,varargin)
%plots the feature of all cells in a sessions set of recs over all the stimuli set of a same type

inFunOption={};
inFunOption=[inFunOption struct('parName','allTogether'  , 'parDefault',  'true', 'validation',@islogical)];
inFunOption=[inFunOption struct('parName','subplotHandle', 'parDefault',  ''    , 'validation',@ishandle)];
inFunOption=[inFunOption struct('parName','rec',           'parDefault',  ''    , 'validation',@(x) iscell(x) || ischar(x) )];

% Parse input arguments
[iP,varargout] = parseInputs(inFunOption,varargin);

% Define recs to process for this mouse,sess
fn = file_names(mouse, sess,{},'skipPrint');
q = load(fn.sess_info);
rec=iP.rec;
if isempty(rec)
    fprintf('\t no record specified, going for all of them...\n');
    reclist = {q.info.rec.name};
else
    if iscell(rec)
        reclist =rec;
    else
        reclist={rec};
    end
end
clear q fn

fprintf('** Plotting cells for mouse %s, sess %s \n**', mouse, sess);

pLeg={};
plottedOId   = [];
recStr='';
if iP.allTogether
    %make one figure, one subfigure and send all the plots there
    if isempty(iP.subplotHandle)
        gac=figure;
        set(gac,'Position',[100, 100, 1600,1200])
        set(gac,'Nextplot','add')
        subFig=subplot(1,1,1);
    else
        subFig=iP.subplotHandle;
    end

    set(subFig,'Nextplot','add');
    timesFound=0;
    for iRec=1:numel(reclist)
        recName=reclist{iRec};
        % ~~~~~~~~~~~~~~~~
        plotRec(recName);
        % ~~~~~~~~~~~~~~~~
    end
    if timesFound
        %make the legend:
        %get the names of all the odors used
        plottedOdors = [pac.pc.odor];
        plottedOId   = [plottedOdors.oId];
       if iP.legends
            plottedLeg   = [pac.pc.pLeg];
            plottedHan   = [pac.pc.pHan];
            firstAppear  = arrayfun(@(x) find(plottedOId==x,1),unique(plottedOId));
            hLeg=legend(subFig,plottedHan(firstAppear),plottedLeg(firstAppear));
            set(hLeg,'color','none')
            set(hLeg,'box','off')
            set(hLeg,'units','pixels')
            legText=findobj(hLeg,'type','text');
            set(legText,'FontSize',12);
        end
        % title and save
        recStr(end)='';
        supString=sprintf('%s-%s-recs:%s - q: %d',mouse,sess,recStr,iP.quality);
        set(subFig,'Nextplot','add');
        if isempty(iP.subplotHandle)
            title(supString);
            fn=file_names(mouse,sess,recStr);
            fileBase=[fn.basename_an stimType '_' iP.feature];
            print(gac,'-depsc',fullfile(fn.fold_an_sess,sprintf('%s_obs.eps',fileBase)));
        else
            supString=sprintf('%s-%s- q: %d',mouse,sess,iP.quality);
            title(subFig,supString)
        end
        
    else
        if exist('gac')
        close(gac)
        clear gac
        end
    end
    
else % of if iP.allToghether
    for iRec=1:numel(reclist)
        recName=reclist{iRec};
        plotRec(recName);
    end
end

pac.iP=iP;
pac.plottedOId=plottedOId; %save the list of plotted odors for making legends


   %% 
    function plotRec(thisRecName)
        global gotUnits
        fn=file_names(mouse,sess,thisRecName,'noPrint');
        fprintf('\t rec %s ...',thisRecName);
        cellRecBaseName=cell_id(mouse,sess,thisRecName);
        
%         if strcmp(cellRecBaseName,'KPawakeM72_006_a')
%             aaa=454;
%         end
        
        baseMatch = cellfun(@(x) strfind(x,cellRecBaseName),gotUnits,'UniformOutput',false);
        qCells_i = find(cellfun(@(x) ~isempty(x),baseMatch));
        qCells_thisRec=[];
        if ~isempty(qCells_i)
            for iq=1:numel(qCells_i)
                qCells_thisRec = [gotUnits{qCells_i(iq)}];
                iqFL = dir(fullfile(fn.fold_unit_db,[qCells_thisRec '*.mat']));
                cellFileList(iq) = iqFL;
            end
        else 
            return
        end

        cellIdList={cellFileList.name};
        
%         cellFileList=dir(fullfile(fn.fold_unit_db,[cellRecBaseName '*.mat']));

        %if its only about the litral cells
        if iP.onlyLitral
            cellIdList=cellfun(@(x) validate_id(x(1:end-4)),{cellFileList.name},'UniformOutput',false);
            cellIdList(cellfun('isempty',cellIdList))=[];
        end

        cellIdList
        
        if ~isempty(cellIdList)
            recStr=[recStr thisRecName '-'];
            %print one by one
            fprintf('\t found %d cells \n',numel(cellIdList));
            if ~iP.allTogether
                %there is one figure per rec
                gac=figure;
                set(gac,'Position',[100, 100, 1600,1200])
                set(gac,'Nextplot','add')
                suptitle(['Mouse ' mouse ', rec ' sess '-' thisRecName]);
                if numel(cellIdList)>1
                    figRows=2;
                else
                    figRows=1;
                end
                figCols=ceil(numel(cellIdList)/figRows);
                for iSubFig=1:numel(cellIdList)
                    subFig(iSubFig)=subplot(figRows,figCols,iSubFig);
                    set(subFig(iSubFig),'Nextplot','add');
                end

                %get the cells and plot them
                for iCell=1:numel(cellIdList)
                    cellId=cell_id(mouse,sess,thisRecName,cellIdList{iCell});
                    pac.pc=plot_cell(cellId,stimType,'figHandle',subFig(iCell),varargout{:});
                end
                %save the figure
                fileBase=[fn.basename_an stimType '_' iP.feature];
                print(gac,'-depsc',fullfile(fn.fold_an_sess,sprintf('%s_obs.eps',fileBase)));

            else
                %plot everything in the only subfigure
                for iCell=1:numel(cellIdList)
                    timesFound=timesFound+1;
                    cellId=cell_id(mouse,sess,thisRecName,cellIdList{iCell});
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    pac.pc(timesFound)=plot_cell(cellId,iP.stimType,'figHandle',subFig,varargout{:});
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    pLeg=[pLeg pac.pc.pLeg];
                end
            end

        end

    end %function plotRec
%%

    function cellNumber=validate_id(cellId)
        % double check attributes of cells being plotted
        cellNumber=[];% get the unit Info and the unit response structure.
        cellFilePath=fullfile(fn.fold_unit_db,    sprintf('%s.mat',cellId));
        unitInfo=load(cellFilePath);
        if ~unitInfo.light 
            warning('Unit %s is not light responsive, skipping the plot',unitInfo.Id);
        elseif ~isempty(iP.quality) && ~(unitInfo.quality>0)
            warning('Unit %s does not meet quality (%d/%d), skipping the plot',unitInfo.Id,unitInfo.quality,iP.quality);
        else
            cellNumber=cellId(end-1:end);
        end
        if unitInfo.light && unitInfo.quality>0
            cellNumber=cellId(end-1:end);
        end
    end

end %function

function [pc] = plot_cell(cellId,stimType,varargin)
%plots the feature of a cell over all the stimuli set of a same type
%1-get the cell info from the cell database
%2-get the stimulus and the response from the analyzed data response
%structure
%3- plot
global rp
fn=file_names();

% Parse the input
% Required:
% cellId
% Optional

featList=rp.featList;
% stimType: stimulus to plot
stimList=rp.stimList;

inPar=inputParser;
inPar.addRequired('cellId')
inPar.addRequired('stimType');

%general options with simple validation functions
for iOpt=1:numel(rp.option)
    inPar.addParamValue(rp.option{iOpt}.parName,rp.option{iOpt}.parDefault,rp.option{iOpt}.validation);
end
%optional parameters specific to this function:
option={};
option=[option struct('parName','figHandle'     ,'parDefault',''                  ,'validation',@ishandle)];

for iOpt=1:numel(option)
    inPar.addParamValue(option{iOpt}.parName,option{iOpt}.parDefault,option{iOpt}.validation);
end

inPar.parse(cellId,stimType, varargin{:});
iP=inPar.Results

% get the unit Info and the unit response structure.
cellFilePath=fullfile(fn.fold_unit_db,    sprintf('%s.mat',cellId));
unitInfo=load(cellFilePath)
% triple check attributes of cells being plotted
if iP.onlyLitral && ~unitInfo.light 
    warning('Unit %s is not light responsive, skipping the plot',unitInfo.Id);
    pc=-1;
    return;
elseif ~isempty(iP.quality) && ~(unitInfo.quality>0)
    warning('Unit %s does not meet quality (%d/%d), skipping the plot',unitInfo.Id,unitInfo.quality,iP.quality);
    pc=-1;
    return;
end

fn=file_names(unitInfo.mouse,unitInfo.sess,unitInfo.rec);
cluList='';
for iu=1:numel(unitInfo.clu)
    cluList=[cluList num2str(unitInfo.clu(iu),'%02d')];
end

unitRespFn     = sprintf('%s%s_units%s_resp.mat',fn.basename_an,iP.stimType,cluList);
unitRespFullFn = fullfile(fn.fold_an_sess,unitRespFn);
unitResp       = struct;

%check if the file exists and load the responses(stim) structure
if ~exist(unitRespFullFn,'file') || iP.forceVR
    %in the future, check if visualize_responses function handle and path
    %were provided and do the resp.
    %now, just yell at you and exit
    warning('Response structure for clusters %s is not there for %s (or was requested) \n attempting to create it',cluList,fn.basename_an);
    addpath(iP.vrFunPath);
    visualizeRespFunction=str2func(iP.vrFun);
    vr=visualizeRespFunction(unitInfo.mouse,unitInfo.sess,unitInfo.rec,unitInfo.clu,iP.stimType);
end

if strcmpi(iP.odorDbType,'struct')
    q=load(fullfile(iP.odorDbPath,iP.odorDbFile));
    odorDb=q.odor;
    clear q;
else
    error('I still dont know how to handle a %s type odor database',ip.odorDbType);
end

q=load(unitRespFullFn);
unitResp = q.resp;
clear q;
stim     = [unitResp.stim];

%%%%%%%%%%%%
%%% FOR NOW, exclude stimuli with odor + light presented
withLaser = [stim.laserDur]>0;
stim = stim(~withLaser);
if isempty(stim)
    pc=-1;
    return
end
    
%do concentration series plots of the user selected feature (or the
%default:latency).
odorNames = unique({stim.odorName});
%lookup the odors in the odor database to get all the properties
odorIds   = odor_lookup(odorNames,odorDb);
odorConcs = unique([stim.odorConc]);
respMat   = nan(numel(odorNames),numel(odorConcs));

%get the feature and concentration vector for each odor
if isempty(iP.figHandle)
    gs=figure;
    set(gs,'Nextplot','add');
    figHandle=subplot(1,1,1);
    suptitle(['Mouse ' unitInfo.mouse ', rec ' unitInfo.sess '-' unitInfo.rec '. Units ' cluList]);
else
    figHandle=iP.figHandle;
end


% Go through each odor that the current cell was subjected to, get 
% parameters and plot
pLeg={};
pHan=[];
for io=1:numel(odorNames)
    %indices in the response matrix where the stimulus is this odor
    odor(io).oId     =odorIds(io);
    odor(io).idx     =find(strcmpi(odorNames{io},{stim.odorName}));
    odor(io).stim    =stim(odor(io).idx);
    %reorder by concentration
    [~, iCS]         =sort([odor(io).stim.odorConc]);
    odor(io).idx     =odor(io).idx(iCS);
    odor(io).idx_c   =odor(io).idx;
    odor(io).stim    =odor(io).stim(iCS);
    odor(io).conc    =[odor(io).stim.odorConc];
    odor(io).molar   =odor(io).conc*odorDb(odor(io).oId).molarFactor;
    odor(io).molar_c  =odor(io).molar;
    
    odor(io).concIdx =arrayfun(@(x)find(odorConcs==x,1),[odor(io).stim.odorConc]);
    respMat(io,odor(io).concIdx) = [unitResp(odor(io).idx).(iP.feature)];
    
    % Specify range of concentrations to plot
    if isempty(iP.specConc_uM)
        concTarget=50;
        concRange=100;
    else
        concTarget=iP.specConc_uM;
        concRange=0.5*concTarget;
    end

    molars=odor(io).molar *1e6;
    inRange=find(abs(molars-concTarget) < concRange);
    if isempty(inRange)
        continue
    else
        odor(io).idx_c=odor(io).idx(inRange);
        odor(io).molar_c = odor(io).molar(inRange);
    end

    
    % Add datapoint(s) to figure (figHandle)
    hold on
    set(figHandle,'Nextplot','add')
    pLeg{io}=odorNames{io};
    if ~iP.LatvJit  % this will plot feature as a function of odor conc
        pHan(io)=plot(figHandle,log10(odor(io).molar_c*1E6),[unitResp(odor(io).idx_c).(iP.feature)], ...
            [odorDb(odor(io).oId).plotSymbol],'MarkerSize',10,'LineWidth',1,'Color',odorDb(odor(io).oId).colorCode,...
            'MarkerFaceColor',odorDb(odor(io).oId).colorCode);
        xlabel(figHandle,'log[C (uM)]');
        %         plot(odor(io).conc,[unitResp(odor(io).idx).(iP.feature)], ...
        %         '-^','MarkerSize',10,'LineWidth',2,'MarkerFaceColor','k','Color',odorDb(odor(io).oId).colorCode)
        %     xlabel('log[C (uM)]');
        ylabel(figHandle,sprintf('%s',iP.feature));
        set(figHandle,'XLim',iP.concAxis)
        set(figHandle,'YLim',iP.latAxis)
    
    elseif iP.LatvJit  % this will plot var(latencyISI) as a function of latencyISI
        pHan(io)=plot(figHandle,[unitResp(odor(io).idx_c).latencyISI],[unitResp(odor(io).idx_c).jitterISI], ...
            [odorDb(odor(io).oId).plotSymbol],'MarkerSize',10,'LineStyle','none','Color',odorDb(odor(io).oId).colorCode,...
            'MarkerFaceColor',odorDb(odor(io).oId).colorCode);
        xlabel(figHandle,'latencyISI');
        ylabel(figHandle,'jitterISI');
        set(figHandle,'XLim',[0 150])
        set(figHandle,'YLim',[0 100])
        
    end
        rp.opl=rp.opl+1;
%         odorspl
        % Big structure to send to plot_odors
        odorspl(rp.opl).odorName    = char(odorNames{io});
        odorspl(rp.opl).odorconc    = odor(io).conc;
        odorspl(rp.opl).latency     = [unitResp(odor(io).idx_c).latencyISI];
        odorspl(rp.opl).jitter      = [unitResp(odor(io).idx_c).jitterISI];
        odorspl(rp.opl).cellId      = iP.cellId;
        
%         save('stats_allcells','odorspl')

end

    if iP.legends
        hLeg=legend(figHandle,pLeg);
        set(hLeg,'color','none')
        set(hLeg,'box','off')
        set(hLeg,'units','pixels')
        legText=findobj(hLeg,'type','text');
        set(legText,'FontSize',8);
        set(legText,'FontSize',8);
        %     legAx=findobj(hLeg,'type','axes');
        %     set(legAx,'visible','off');
        %     lp=get(hLeg,'outerposition');
        %     set(hLeg,'outerposition',[lp(1), lp(2),lp(3:4)*0.75]);
    end
    title(figHandle,['Units ' cluList]);


% figure
% plot((odorConcs),respMat','-*')
pc.odor=odor;
pc.resp=unitResp;
pc.concMat=respMat;
pc.odorDb=odorDb;
pc.pLeg=pLeg;
pc.pHan=pHan;

end %function plot_cell



function [iP,varargout] = parseInputs(inFunOption,varargin)
global rp;
% Check input 
varargin=varargin{:};   
if mod(numel(varargin),2)
    warning('Number of input arguments is not even.')
end

inPar=inputParser;

%enter 'global' parameters
for iOpt=1:numel(rp.option)
    option=rp.option{iOpt};
    addParameter(inPar,option.parName,option.parDefault,option.validation);
end
%enter 'local' parameters (only for this function)
for iOpt=1:numel(inFunOption)
    option=inFunOption{iOpt};
    addParameter(inPar,option.parName,option.parDefault,option.validation);
end
% enter generic optional parameters
addParameter(inPar,'stimType',rp.stimList{1},@(x) ischar(x) && any(strcmpi(x,rp.stimList)));

% Parse inputs
parse(inPar,varargin{:});

% varargout replicates varargin, but with function-specific inputs removed
varargout=varargin;
rmOpt=[inFunOption{:}];
rmOpt={rmOpt.parName};
rmPar={rmOpt{:}};
for idxPar=1:numel(rmPar)
    par=rmPar{idxPar};
    iCell=find(cellfun(@(c) strcmpi(par,c),varargout(:)));
    varargout(iCell:iCell+1)=[];
end

iP=inPar.Results;
varargout={varargout};
end


function cellId = cell_id(mouse,sess,rec,unit)
    
if isnumeric(mouse)
    mouse = sprintf('%04d', mouse);
end
if isnumeric(sess)
    sess = sprintf('%03d', sess);
end
if isnumeric(rec)
    rec = sprintf('%02d', rec);
end
cellId=sprintf('%s_%s_%s', mouse, sess, rec);

if nargin>3 && ~isempty(unit)
    if isnumeric(unit)
        unit = sprintf('%02d', unit);
    end
    
    cellId=sprintf('%s_%s_%s_%s', mouse, sess, rec,unit);
end

end

function [oIds] = odor_lookup(odorNames,odorDb)
%lookup index (id) of odors in an odor database (structure)
%the database is a structure array with fields:
%name: cell array of strings
%n: int (odor id)
%the last element of the database is name 'unlisted', n=0;

nOdors=numel(odorNames);
%the default is the zero (last) element of the database (unlisted)
oIds=zeros(nOdors,1);

for io=1:nOdors
    oIdAux=find(arrayfun(@(x) any(strcmpi(odorNames(io),x{:})),{odorDb.name}));
    if isempty(oIdAux)
        warning('Odor %s not found in the odor database, using unlisted ',odorNames{io});
        oIdAux=find(arrayfun(@(x) any(strcmpi('unlisted',x{:})),{odorDb.name}));
    end
    oIds(io)=oIdAux;
end
end

function [gf, gs] = subplot_fig(ny, nx, fig)
    if exist('fig', 'var')
        gf = figure(fig); clf
    else
        gf = figure;
    end
    
    dx = 0.9/nx;    dy = 0.9/ny;
    wx = dx*0.95;   wy = dy*0.95;
    for ix = 1:nx
        for iy = 1:ny
            gs(iy,ix) = subplot('position', [0.08+(ix-1)*dx, 0.95-iy*dy, wx, wy]);
        end
    end
end


























