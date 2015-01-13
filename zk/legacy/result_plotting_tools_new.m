% Functions for plotting results of first analysis:
% 
function rp = result_plotting_tools()
%global rp

% if this is not the computing server, add the include folder to the path
[~,computerName]=system('hostname');
if ~strcmp(strtrim(computerName),'flipper')
    includePath=fullfile(fileparts(pwd),'current','include');
    addpath(includePath);
end
%
    rp.plot_cell          = @plot_cell;
    
    rp.odor_lookup        = @odor_lookup;
    rp.cell_id            = @cell_id;
    rp.subplot_fig        = @subplot_fig;
    
end

function [pc] = plot_cell(cellId,varargin)
%plots the feature of a cell over all the stimuli set of a same type
%1-get the cell info from the cell database
%2-get the stimulus and the response from the analyzed data response
%structure
%3- plot
fn=file_names();

% Parse the input
% Required:
% cellId
% Optional
% feat: features to plot:
featList={'latency','nspk','jitter'};
% stimType: stimulus to plot
stimList={'odor','laser'};
% all the optional variable values

odorDbFileDefault='odorDb.mat';
odorDbPathDefault=fn.st_meta;
odorDbTypeDefault='struct';
% optional parameters

inPar=inputParser;
inPar.addRequired('cellId')
inPar.addRequired('stimType');
%Fuction for visualizing responses
inPar.addParamValue('vrFun','visualize_responses',@(x) ischar(x));
inPar.addParamValue('vrFunPath',fullfile(get_local_disk,'ephysDataManagement','zk'),@ischar)
%Odor database
inPar.addParamValue('odorDbFile',odorDbFileDefault,@(x) ischar(x));
inPar.addParamValue('odorDbPath',odorDbPathDefault,@ischar);
inPar.addParamValue('odorDbType',odorDbTypeDefault,@ischar);

inPar.addOptional('feat', featList{1},@(x) ischar(x) && ~any(strcmpi(x,[inPar.Parameters stimList])) && any(strcmpi(x,featList)));


inPar.parse(cellId, varargin{:});
iP=inPar.Results

% get the unit Info and the unit response structure.
cellFilePath=fullfile(fn.fold_unit_db,    sprintf('%s.mat',cellId));
unitInfo=load(cellFilePath)
fn=file_names(unitInfo.mouse,unitInfo.sess,unitInfo.rec)
cluList='';
for iu=1:numel(unitInfo.clu)
    cluList=[cluList num2str(unitInfo.clu(iu),'%02d')];
end

unitRespFn     = sprintf('%s%s_units%s_resp.mat',fn.basename_an,iP.stimType,cluList);
unitRespFullFn = fullfile(fn.fold_an_sess,unitRespFn);
unitResp       = struct;

%check if the file exists and load the responses(stim) structure
if ~exist(unitRespFullFn,'file')
    %in the future, check if visualize_responses function handle and path
    %were provided and do the resp.
    %now, just yell at you and exit
    warning('Response structure for clusters %s is not there for %s \n attempting to create it',cluList,fn.basename_an);
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

%do concentration series plots of the user selected feature (or the
%default:latency).
odorNames = unique({stim.odorName});
%lookup the odors in the odor database to get all the properties
odorIds   = odor_lookup(odorNames,odorDb);
odorConcs = unique([stim.odorConc]);
respMat   = nan(numel(odorNames),numel(odorConcs));
%get the feature and concentration vector fore ach odor
gs=figure
set(gs,'Nextplot','add')
for io=1:numel(odorNames)
    %indices in the response matrix where the stimulus is this odor
    odor(io).idx     =find(strcmpi(odorNames{io},{stim.odorName}));
    odor(io).oId     =odorIds(io);
    odor(io).stim    =stim(odor(io).idx);
    odor(io).conc    =[odor(io).stim.odorConc];
    odor(io).molar   =odor(io).conc*odorDb(odor(io).oId).molarFactor;
    
    odor(io).concIdx =arrayfun(@(x)find(odorConcs==x,1),[odor(io).stim.odorConc]);
    respMat(io,odor(io).concIdx) = [unitResp(odor(io).idx).(iP.feat)];
    
    %add this plot
    figure(gs)
    hold on
    plot(log10(odor(io).molar*1E6),[unitResp(odor(io).idx).(iP.feat)], ...
        '-^','MarkerSize',10,'LineWidth',2,'MarkerFaceColor','k','Color',odorDb(odor(io).oId).colorCode)
    xlabel('log[C (uM)]');
    ylabel(sprintf('%s',iP.feat));
    
end

figure
plot(log(odorConcs),respMat','-*')
pc.odor=odor;
pc.resp=unitResp;
pc.concMat=respMat;
pc.odorDb=odorDb;

end %function plot_cell

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
    oIds(io)=find(arrayfun(@(x) any(strcmpi(odorNames(io),x{:})),{odorDb.name}));    
end
oIds(oIds==numel(odorDb))=0;
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
    if isnumeric(unit)
        unit = sprintf('%02d', unit);
    end
    
    cellId=sprintf('%s_%s_%s_%s', mouse, sess, rec,unit);
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