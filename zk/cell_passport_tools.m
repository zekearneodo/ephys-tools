function cp = cell_passport_tools(doit)
%functions for making a cell passport
%uses visualize_responses

%add the package of common files and folders to the path.
%here are the functions for file_naming, for instance, and any other
%functions that several script bundles will use.
%that is /baseFolder/ephysDataManagement/current/include
[~,computerName]=system('hostname');
if ~strcmp(strtrim(computerName),'flipper')
    includePath=fullfile(fileparts(pwd),'current','include');
    addpath(includePath);
end

    vr = visualize_responses_new();
           
    cp.plot_response = @plot_response;
    
    
    cp.make_passport       = @make_passport;
    cp.setup_fig           = @setup_fig;
    cp.make_response_frame = @make_response_frame;
    cp.make_frame_geometry = @make_frame_geometry;
    
    cp.plot_light_responses    = @plot_light_responses;
    cp.plot_odor_responses     = @plot_odor_responses;
    cp.plot_sniff_distribution = @plot_sniff_distribution;
    cp.plot_unit               = @plot_unit;
    
    cp.get_no_stim_sniffs = @get_no_stim_sniffs;
    cp.make_sniff_dist    = @make_sniffs_dist;
    
    
    if nargin>0 && doit
        global response;
        response = plot_response();
    end
end

function passport = make_passport(unit_meta)


passport = setup_fig();

passport.unit_meta = unit_meta;
passport.fn        = file_names(unit_meta.mouse, unit_meta.sess, unit_meta.rec);

passport = plot_odor_responses(passport);
passport = plot_light_responses(passport);


%make titles and save the figure
suptitle(sprintf('%s',unit_meta.Id));

end

function passport = setup_fig()

%prepare the figure
fg = figure();
set(fg,'Position',[0 0 1600 1200])
set(0, 'currentfigure', fg);
set(fg,'DefaultAxesFontSize',6)
set(fg,'defaulttextinterpreter','none')

passport.fig.fig_handle = fg;

passport.fig.rows = 4;
passport.fig.cols = 7;
%

end

function passport = plot_odor_responses(passport)
%plots the odor responses in the figure
%for now, it just calls visualize_responses_play

a_unit = passport.unit_meta;

[vr, passport] = visualize_responses_play(passport, a_unit.mouse,a_unit.sess,a_unit.rec,a_unit.clu,'odor');
end

function passport = plot_light_responses(passport)
%plots the odor responses in the figure
%for now, it just calls visualize_responses_play

a_unit = passport.unit_meta;

[vr, passport] = visualize_responses_play(passport, a_unit.mouse,a_unit.sess,a_unit.rec,a_unit.clu,'laser');


end

function passport = plot_sniff_distribution(passport)
%plot the sniff duration distribution

a_unit = passport.unit_meta;

%get the sniff distribution

%make the frame

%plot ths sniff distribution

end

function passport = plot_unit(passport)

a_unit = passport.unit_meta;

%get the unit

%make the correlogram

%make the frame for waveform/correlogram

%plot correlogram

%plot waveform

end


%plots a raster and a psth in response to a stimulus
function response_plot = plot_response(response)

%prepare the figure
fg = response.fg.fig_handle;

% this should be set elsewhere
rows = response.fg.rows;
cols = response.fg.cols;
%
dim_x=1./cols;
dim_y=1./rows;

%make the frame with stim-dependent dimensions
if strcmp(response.sType,'odor')
    % make a title with the odor and the concentration
    title_text = sprintf('%s \n %.2e',response.stim.odorName, response.stim.odorConc);
    frame_set = make_frame_geometry(response.stim);
else
    % make a title with the odor and the concentration
title_text = sprintf('%dmv \n %dms',response.stim.laserAmp, response.stim.laserDur);
    frame_set.i =  0;
    frame_set.j = -1; % will plot immediately above lowes concentration
end
frame     = make_response_frame(frame_set.i*dim_x, (3-frame_set.j-1)*dim_y, dim_x, dim_y, fg);

%plot the raster
t1   = response.vp.t1;
t2   = response.vp.t2;
rmax = response.rmax*1.2;

kt = response.ntrial;
yt = max([1, floor(kt/10)*10]);

plot(frame.ax_rast, response.raster.x, response.raster.y, '.', 'MarkerSize',7 * 2*dim_y);
set(frame.ax_rast, 'NextPlot', 'add');
plot(frame.ax_rast, [0,0], [0, kt+1], '--k')

set(frame.ax_rast, 'XLim', [t1, t2], ...
    'YLim',        [0, kt+1], ...
    'YTick',       [0, yt], ...
    'XTick',       [t1, 0, -t1],...
    'XTickLabel',[]);

%title(stim_str, 'FontSize', 10, 'FontWeight', 'bold') insert some text as
%title
%%%%

%plot the psth
    plot(frame.ax_psth, response.t, response.rate, ...
        'LineWidth', 1, ...
        'Color',     response.stim.lineColor, ...
        'LineStyle', response.stim.lineStyle),
    set(frame.ax_psth, 'NextPlot', 'add');
    %plot(t,ones(size(t))*baseline,'.k')
    %plot(t,ones(size(t))*pkThresh,'b')
    
    %if its odor plot the baseline
    if strcmp(response.sType,'odor')
        plot(frame.ax_psth, response.t,response.rateBaseHist,'LineStyle','--','Color',[.5,0.5,0.5])
    end
    plot(frame.ax_psth, response.latency, response.maxFR,'k*')
    plot(frame.ax_psth, [0,0], [0, rmax*1.2], '--k'), hold off
    set(frame.ax_psth, 'XLim',     [t1, t2], ...
        'XTick',         [t1, 0, -t1], ...
        'YLim',          [0, rmax], ...
        'YTick',         0:20:rmax);
%%%% 

% insert the title
title(frame.ax_rast,title_text,'FontSize', 8);

% 
% plot(frame.ax_rast,x,y)
% set(frame.ax_rast, 'XTickLabel',[]);
% plot(frame.ax_psth,x,z)
% 
response_plot.frame = frame;

end

function frame = make_response_frame(pos_x, pos_y, dim_x, dim_y, fig)

% Makes a frame with two axes for plotting psths. Splits the frame into two
% vertical axes
% 
% pos_x, pos_y: position in the fig of the low-left corner of the frame
% dim_x, dim_y: total dimensions of the frame

% activate or create figure
if exist('fig', 'var')
    gf = figure(fig);
else
    gf = figure;
end
set(0, 'currentfigure', gf);  %# for figures
set(gf,'nextPlot', 'add');

%create axes
shift_x = 0.075;
shift_y = 0.075;

ax_rast = axes('Position',[pos_x + shift_x*dim_x, 0.5*dim_y + pos_y - 0.5*shift_y*dim_y,...
        (1-2*shift_x)*dim_x, (0.45-shift_y)*dim_y]);

    set(ax_rast, 'XTick',[]);
    set(ax_rast, 'YTick',[]);

ax_psth = axes('Position',[pos_x + shift_x*dim_x, pos_y + shift_y*dim_y,...
     (1-2*shift_x)*dim_x, (0.45-shift_y)*dim_y]);
 
    set(ax_psth, 'XTick',[]);
    set(ax_psth, 'YTick',[]);

 
%return the frame
frame.fig  = gf;
frame.ax_rast = ax_rast;
frame.ax_psth = ax_psth;


end

function geometry = make_frame_geometry(one_stim)
%retunrs position and dimension for a stimulus

%there is a lookup table with odors and concentratios, and from there comes
%the position of the frame.
%For now:
 % the grid is fixed to 3 concentrations, 7 odors
 % the dimension of the total grid ar 4x7
 % the lookup table is a struct, made by hand or for the experiment with a
 % script
 
 odors = [];
 
 odor.alias = {'2-hydroxy', '2-hydroxyacetophenone'};
 odor.concs = [0.0051, 2.97e-4, 1.78e-5];
 odors=[odors odor];
 
 odor.alias = {'ethyl-tiglate','ethyl_tiglate','ethyl tiglate'};
 odor.concs = [1.3e-5, 1.3e-4, nan];
 odors=[odors odor];
 
 odor.alias = {'4-methyl_acetophenone','4-methylacetophenone'};
 odor.concs = [0.0018, nan, nan];
 odors=[odors odor];
 
 odor.alias = {'menthone'};
 odor.concs = [0.016, 0.0016, 7e-5];
 odors=[odors odor];
 
 odor.alias = {'acetophenone'};
 odor.concs = [7.0e-4, nan , nan];
 odors=[odors odor];
 
 odor.alias = {'benzaldehyde'};
 odor.concs = [1.9e-4, nan , nan];
 odors=[odors odor];
 
 odor.alias = {'2-4-dimethylacetophenone'};
 odor.concs = [0.0058, nan, nan];
 odors=[odors odor];
 
 
 
 %for testing, just use 2hydroxy and menthone:
 odor_index = find( cellfun(@(x) any(strcmpi(one_stim.odorName,x)),{odors.alias}))
 if isempty(odor_index)
    warning('odor %s not identified as a stimulus for the grid', one_stim.odorName);
    odor_index = 0;
 else
    conc_index = find(odors(odor_index).concs> one_stim.odorConc*0.5 & odors(odor_index).concs<one_stim.odorConc*1.5);
    if isempty(conc_index)
        warning('concentriaton %0.4g of odor %s not identified as a stimulus for the grid', one_stim.odorConc, one_stim.odorName);
        conc_index = 0;
    end
 end
 
 geometry.i = odor_index -1 ;
 geometry.j = conc_index -1;
 
end

function sniffs_dist = make_sniffs_dist(sniffs,bin)
%make the distribuition of sniffs
%sniffs: array of sniffs structures
len_range = [0, 1000];
len_axis  = len_range(1):len_range(2);

inh_lengths = nan(1, numel(sniffs));
exh_lengths = nan(1, numel(sniffs));

for i=1:numel(sniffs)
    inh_lengths(i) = round(sniffs(i).t_zer_fit(2))-sniffs(i).t_zer(1);
    exh_lengths(i) = sniffs(i).t_zer(3)-round(sniffs(i).t_zer_fit(2));
end

inh_count = arrayfun(@(x) sum(inh_lengths==x), len_axis);
exh_count = arrayfun(@(x) sum(exh_lengths==x), len_axis);

sniffs_dist.inh = decimate(inh_count, bin);
sniffs_dist.exh = decimate(exh_count, bin);

end

function sniffs = get_no_stim_sniffs(unit_meta)

% First get all the no_stim sniffs
% Then get the distribution

fn      = file_names(unit_meta.mouse, unit_meta.sess, unit_meta.rec);

load(fn.sniffs) % loads to sniff
load(fn.trial)  % loads to trial


%%clean up trials for sessions < 20 (before serial trial number was used)
badTrials=arrayfun(@(x) isempty(trial(x).start),1:numel(trial));
if any(badTrials)
    warning('Some bad trials found in the rec (start was empty)');
end
trial(badTrials)=[];

badSniffs = 0;
no_stim_sniffs = zeros(1,numel(sniff));
for is=1:numel(sniff)
    sn=sniff(is);
    t_inh = sn.t0+sn.t_zer(1);
    %check if it is within a trial
    prev_trial_on = find([trial.start]< t_inh,1,'last');
    % there is no prev trial
    % or the prev trial ended more than 5 secs ago
    if ( isempty(prev_trial_on) | t_inh > (trial(prev_trial_on).start + trial(prev_trial_on).runTrialDur + 5000));
        no_stim_sniffs(is) = 1;
    end
    
end

if badSniffs>0
    warning('%d problems with sniff zeros',badSniffs);
end

sniffs = sniff(find(no_stim_sniffs));

end