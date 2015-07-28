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
    cp.make_info_frame     = @make_info_frame;
    cp.make_frame_geometry = @make_frame_geometry;
    
    cp.plot_light_responses    = @plot_light_responses;
    cp.plot_odor_responses     = @plot_odor_responses;
    cp.plot_sniff_distribution = @plot_sniff_distribution;
    cp.plot_unit               = @plot_unit;
    cp.plot_info               = @plot_info;
    
    cp.get_no_stim_sniffs = @get_no_stim_sniffs;
    cp.make_sniff_dist    = @make_sniff_dist;
    
    cp.load_unit = @load_unit;
    cp.load_info = @load_info;
    cp.read_info = @read_info;
    
    cp.isi_distribution = @isi_distribution;
    
    
    
    
    if nargin>0 && doit
        global response;
        response = plot_response();
    end
end

function passport = make_passport(unit_meta)

close all
passport = setup_fig();

passport.unit_meta = unit_meta;
passport.fn        = file_names(unit_meta.mouse, unit_meta.sess, unit_meta.rec);

passport = plot_odor_responses(passport);
passport = plot_light_responses(passport);
passport = plot_sniff_distribution(passport);
passport = plot_unit(passport);
passport = plot_info(passport);
%make titles and save the figure
suptitle(sprintf('%s',unit_meta.Id));

%save it to pdf
fn = file_names(unit_meta.mouse, unit_meta.sess, unit_meta.rec);
fig_file_name = fullfile(fn.fold_exp_data,[unit_meta.Id '_fig.pdf']);
set(passport.fig.fig_handle,'PaperOrientation','landscape');
set(passport.fig.fig_handle,'PaperUnits','normalized');
set(passport.fig.fig_handle,'PaperPosition', [0 0 1 1]);
set(passport.fig.fig_handle,'PaperType', 'usletter');
% set(passport.fig.fig_handle,'PaperPositionMode','auto'); 
% set(passport.fig.fig_handle,'PaperPosition', [0 0 1 1]);

print(passport.fig.fig_handle, '-dpdf',fig_file_name);

end

function passport = setup_fig()

%prepare the figure
fg = figure();
set(fg,'Position',[0 0 1600 1200])
%set(fg,'Position',[0 0 792 612])
set(0, 'currentfigure', fg);
set(fg,'DefaultAxesFontSize',6)
set(fg,'defaulttextinterpreter','none')

passport.fig.fig_handle = fg;
passport.fig.default_title_font_size = 6;

passport.fig.rows = 4;
passport.fig.cols = 7;
%
passport.fig.sniff.pos        = [2, 3];
passport.fig.sniff.hist_bin   = 5;
passport.fig.sniff.hist_t_max = 700;

passport.fig.unit.pos        = [3, 3];
passport.fig.unit.hist_bin   = 1;
passport.fig.unit.hist_t_max = 30;

passport.fig.info.pos        = [6, 3];

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

%make the frame
%there should be a unified function for makint the geometry
% for now, just make a frame similar to the raster (2-split vertical), 
% and plot the distributions in the lower, leaving the upper for other


dim_x = 1./passport.fig.cols;
dim_y = 1./passport.fig.rows;

pos_x = passport.fig.sniff.pos(1)*dim_x;
pos_y = passport.fig.sniff.pos(2)*dim_y;

frame = make_response_frame(pos_x, pos_y, dim_x, dim_y, passport.fig.fig_handle);


sniffs = get_no_stim_sniffs(passport.unit_meta);
%get the sniff distribution
sd = make_sniff_dist(sniffs, passport.fig.sniff.hist_bin);


%plot the raster
t1   = 0;
t2   = passport.fig.sniff.hist_t_max;
sd_max = max([sd.inh sd.exh])*1.2;
sd_legend = {'inh' 'exh'};
title_text = 'sniff durations';

% Change to new colors.
set(frame.ax_psth, 'ColorOrder', [0 0.75 0; 0.5 0.5 0.5], 'NextPlot', 'replacechildren');

%plot the distribution in the psth axis
    plot(frame.ax_psth, sd.t, [sd.inh; sd.exh],...
        'LineWidth', 2)
    set(frame.ax_psth, 'NextPlot', 'add');
    set(frame.ax_psth, ...
        'XLim',      [t1, t2], ...
        'XTick',     [t1, 100, t2-100], ...
        'YLim',      [0, sd_max], ...
        'YTick',     0:20:sd_max-10);
    
    frame.leg = legend(frame.ax_psth, sd_legend);
    legend(frame.ax_psth,'boxoff');
    
%%%% 

% insert the title
title(frame.ax_psth,title_text,'FontSize', passport.fig.default_title_font_size);

passport.fig.sniff.frame = frame;

end

function passport = plot_unit(passport)
%plot unit:
% iti distribuition
% waveform (tricky, leave for later)

%make the frame
%there should be a unified function for makint the geometry
% for now, just make a frame similar to the raster (2-split vertical), 
% and plot the distributions in the lower, leaving the upper for other

dim_x = 1./passport.fig.cols;
dim_y = 1./passport.fig.rows;

pos_x = passport.fig.unit.pos(1)*dim_x;
pos_y = passport.fig.unit.pos(2)*dim_y;

frame = make_response_frame(pos_x, pos_y, dim_x, dim_y, passport.fig.fig_handle);

%get the unit
a_unit = passport.unit_meta;
u_s = load_unit(a_unit);

%get the distribution
isi = isi_distribution(u_s.times);
%decimate
isi_dist = decimate(isi.dist, passport.fig.unit.hist_bin);
t        = decimate(isi.t, passport.fig.unit.hist_bin);

%plot the distribuition
t1   = -passport.fig.unit.hist_t_max;
t2   = passport.fig.unit.hist_t_max;

isi_d_max = max(isi_dist)*1.2;
title_text = 'isi';

% Change to new colors.
set(frame.ax_psth, 'ColorOrder', [0 0.75 0; 0.5 0.5 0.5], 'NextPlot', 'replacechildren');

%plot the distribution in the psth axis
    area(frame.ax_psth, [-fliplr(t) t], [fliplr(isi_dist) isi_dist],...
        'FaceColor', [1 1 0],...
        'EdgeColor',  [0.5 0.5, 0.5]);
    
    set(frame.ax_psth, ...
        'XLim',      [t1, t2], ...
        'XTick',     [t1+10, 0, t2-10], ...
        'YLim',      [0, isi_d_max], ...
        'YTick',     0:round(isi_d_max/5):isi_d_max-round(isi_d_max/5));
    
% insert the title
title(frame.ax_psth,title_text,'FontSize', passport.fig.default_title_font_size);

%%%% 
%plot all the spikes in an x axis to spot silent periods
%spike count in intervals of some seconds
%rounded count
bs=5000;
binned_times = floor(u_s.times/bs);
binned_t     = 1:max(binned_times);
spike_dist  = arrayfun(@(x) sum(binned_times==x),binned_t);
spike_times = round(binned_t*bs/1000); %this one is in sec!!

t1   = 0;
t2   = max(spike_times);
rmax = max(spike_dist);
title_text = sprintf('spike events (%d sec bin)', round(bs/1000));

plot(frame.ax_rast, spike_times, spike_dist);
set(frame.ax_rast, 'XLim', [t1, t2], ...
    'YLim',        [0, rmax*1.2], ...
    'YTick',       [0, rmax], ...
    'XTick',       [0, round(t2*0.9)]);

title(frame.ax_rast,title_text,'FontSize', passport.fig.default_title_font_size);


passport.fig.unit.frame = frame;

end

function passport = plot_info(passport)


dim_x = 1./passport.fig.cols;
dim_y = 1./passport.fig.rows;

pos_x = passport.fig.info.pos(1)*dim_x;
pos_y = passport.fig.info.pos(2)*dim_y;



title_text = 'unit meta';
frame = make_info_frame(pos_x, pos_y, dim_x, dim_y, passport.fig.fig_handle);
title(frame.ax_info,title_text,'FontSize', passport.fig.default_title_font_size);

set(frame.ax_info, 'XLim', [0 100], 'YLim', [0 100]);

%get unit info
[s_info] = load_info(passport.unit_meta);

%with s_info buildup the text and positions
%example:
label_props = [];
label_prop.label  = 'Mouse Id: ';
if isfield(s_info, 'mouseId')
    label_prop.text = num2str(s_info.mouseId, '%.1f');
else
    label_prop.text = '';
end
label_prop.pos    = [0, 90];
label_props = [label_props label_prop];

label_prop.label  = 'State: ';
label_prop.text = s_info.rec.awake_anesth;
label_prop.pos    = [0, 80];
label_props = [label_props label_prop];

label_prop.label  = 'Probe: ';
label_prop.text   = s_info.rec.electrode_type;
label_prop.pos    = [0, 70];
label_props = [label_props label_prop];

try
label_prop.label  = 'Depth: ';
label_prop.text   = num2str(s_info.rec.site_depth);
label_prop.pos    = [0, 60];
label_props = [label_props label_prop];
catch me
end
try
label_prop.label  = 'Side: ';
label_prop.text   = (s_info.rec.site_side);
label_prop.pos    = [0, 50];
label_props = [label_props label_prop];
catch me
end
axes(frame.ax_info)
axis off
for i=1:numel(label_props)
    lp = label_props(i);
    text(lp.pos(1),lp.pos(2),[lp.label lp.text],'HorizontalAlignment','left');
end

passport.fig.info.frame = frame;

end

function unit = load_unit(unit_meta)

% clu   = unit_meta.clu;
% 
% if numel(clu)>1
%     warning('Unit %s has more than 1 cluster, still dont know how to handle that', unit_meta.Id)
%     unit = [];
%     return
% end
% 
% group = unit_meta.group;
fn = file_names(unit_meta.mouse, unit_meta.sess, unit_meta.rec);
q=load(fn.exp_spikes);

%fn spikes has all the cells of the rec (the used ones)
%look it up by  sessCell, which is the unique identifier of
%the cell across the session and appears only once in te rec

unit=q.unit([q.unit.sessCell]== unit_meta.sessCell);

end

function [s_info] = load_info(unit_meta)

fn = file_names(unit_meta.mouse, unit_meta.sess, unit_meta.rec);
updated_info = read_info(unit_meta.mouse, unit_meta.sess);

if isempty(updated_info)
    q=load(fn.ss_sess_info);
    s_info = q.info;
else
    s_info = updated_info;
end

s_info.rec(~strcmpi(unit_meta.rec, {s_info.rec.name})) = [];

end

function info = read_info(mouse, sess)
% the program read the following inofmration:
% 1 - log file
% 2 - raw data file directory
% 3 - electrode channel config
% 4 - meta data for individual recordings
% It prints the structure rec for debugging.
% It saves the data into ss_fold, sess_info.mat file

    fn = file_names(mouse, sess);
    if ~exist(fn.log, 'file')
        warning('no log file')
        info = [];
        return
    end
        
    fn = file_names(mouse, sess);
    [info, rec] = read_log_file();
    info.rec = rec;
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function [info, rec] = read_log_file()
        info   = struct();
        rec    = struct();
        irec   = 0;
        if ~exist(fn.log, 'file')
            error('no log file')
        end
        fid = fopen(fn.log, 'r');

        neof   = true;
        plevel = 1;
        while neof
            % read a line from log file
            tl = fgetl(fid);
            % if end of the file, finish the prorgam, 
            % if the line starts with ':' - read the parameter
            if tl== -1
                neof =0;
            elseif ~isempty(tl)&&(tl(1)==':')
                % extract parameter name 'var' and its value 'value'
                [var, value] = strtok(tl(2:end), ':');
                var   = strtrim(var);
                value = strtrim(strtok(value(2:end), '%'));
                nvalue = str2double(value);

                
                switch var
                    case 'rec'
                        plevel = 2;
                        irec   = irec + 1;
                        if irec > 1
                            %inherit the structure but empty all the fields
                            rec(irec) = empty_structure(rec(irec-1));
                        end
                        rec(irec).name = value;
                        % do not inherit run structure fields from prev. fields
                        if isfield(rec(irec),'run') && ~isempty(rec(irec).run)
                            rec(irec).run=[];
                        end
                            
                    case 'run'
                        plevel = 3;
                        irun = nvalue;
                        if isnan(irun)
                            error('wrong run number')
                        end
                        if irun > 1
                            %inherit the structure but empty all the fields
                            rec(irec).run(irun) = empty_structure(rec(irec).run(irun-1));
                        end
                        rec(irec).run(irun).num = irun;
                     otherwise
                        % record parameter value, if it is a numerical convert to numerical
                        if ~isnan(nvalue)
                            value = nvalue;
                        end
                        
                        switch plevel
                            case 1
                                info.(var) = value;
                            case 2
                                rec(irec).(var) = value;
                            case 3
                                rec(irec).run(irun).(var) = value;
                        end
                end
            end
        end
        fclose(fid);
%         save('recdbg.mat','rec')
    end
    
 
    function stru = empty_structure(model)
    %returns and empty model structure
    allFields=fields(model);
    for iF=1:numel(allFields)
            stru.(allFields{iF})=[];
    end
    
    end
end

function isi = isi_distribution(spikes1)

t = [1:100];

isi_ms = round(diff(spikes1));

isi.dist   = arrayfun(@(x) sum(isi_ms==x), t);
isi.t      = t;
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
title(frame.ax_rast,title_text,'FontSize', response.fg.default_title_font_size);

% 
% plot(frame.ax_rast,x,y)
% set(frame.ax_rast, 'XTickLabel',[]);
% plot(frame.ax_psth,x,z)
% 
response_plot.frame = frame;

end

function frame = make_info_frame(pos_x, pos_y, dim_x, dim_y, fig)

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
shift_x = 0.05;
shift_y = 0.05;

ax_info = axes('OuterPosition',[pos_x + shift_x*dim_x, pos_y - 0.5*shift_y*dim_y,...
    (1-2*shift_x)*dim_x, (1-shift_y)*dim_y]);
set(ax_info, 'XTick',[]);
set(ax_info, 'YTick',[]);
set(ax_info, 'NextPlot', 'add');
set(ax_info,'DefaultTextFontSize',8)

%return the frame
frame.fig  = gf;
frame.ax_info = ax_info;


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
shift_x = 0.05;
shift_y = 0.05;

ax_rast = axes('OuterPosition',[pos_x + shift_x*dim_x, 0.47*dim_y + pos_y - 0.5*shift_y*dim_y,...
        (1-2*shift_x)*dim_x, (0.5-shift_y)*dim_y]);

    set(ax_rast, 'XTick',[]);
    set(ax_rast, 'YTick',[]);

ax_psth = axes('OuterPosition',[pos_x + shift_x*dim_x, pos_y + shift_y*dim_y,...
     (1-2*shift_x)*dim_x, (0.47-shift_y)*dim_y]);
 
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
 odor_index = find( cellfun(@(x) any(strcmpi(one_stim.odorName,x)),{odors.alias}));
 if isempty(odor_index)
    warning('odor %s not identified as a stimulus for the grid', one_stim.odorName);
    odor_index = 0;
 else
    %conc_index = find(odors(odor_index).concs> one_stim.odorConc*0.5 & odors(odor_index).concs<one_stim.odorConc*1.5);
    [conc_mismatch, conc_index]=min(abs(odors(odor_index).concs/one_stim.odorConc-1));
    if 1-conc_mismatch < 0.5
        warning('concentriaton %0.4g of odor %s poorly identified as a stimulus for the grid', one_stim.odorConc, one_stim.odorName);
        conc_index = 0;
    end
 end
 
 geometry.i = odor_index -1;
 geometry.j = conc_index -1;
 
end

function sniffs_dist = make_sniff_dist(sniffs,bin)
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
sniffs_dist.t   = decimate(len_axis, bin);

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
    good_sniff = true;
    t_inh = sn.t0+sn.t_zer(1);
    %check if it is within a trial
    prev_trial_on = find([trial.start]< t_inh,1,'last');
    %check if the sniff is fine:
    if sn.t_zer_fit(2) <= sn.t_zer(1) | sn.t_zer(3) <= sn.t_zer_fit(2) | any(imag(sn.t_zer_fit)>0)
        good_sniff = false;
    end
    % there is no prev trial
    % or the prev trial ended more than 5 secs ago
    if ( isempty(prev_trial_on) | t_inh > (trial(prev_trial_on).start + trial(prev_trial_on).runTrialDur + 5000)) & good_sniff
        no_stim_sniffs(is) = 1;
    end
    
end

if badSniffs>0
    warning('%d problems with sniff zeros',badSniffs);
end

sniffs = sniff(find(no_stim_sniffs));

end