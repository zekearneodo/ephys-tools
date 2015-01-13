%% Zk modifications on DR sniff_analysis_03.
% 2013-11-12
% set up for computational server
% add record item menu. 
% Needs debugging: the record menu does not run well when the session is
% changed.

function sniff_analysis_04_zk()
global sa


sa = struct();
sa.disk = '';
sa.mouse = '';
sa.session = [];
sa.record= '';
sa.da = da_tools();

init_gui()
% init_mouse_sess(sa.fig.h)

init_subplots();
init_time_menu();
init_param();
init_fcn();
init_mouse_rec()
end


function init_gui
global sa

    sa.fig.h = figure(24); clf
    screen = get(0, 'ScreenSize');
    sa.fig.title_str = 'Snff analisis     disk: %s    mouse: %s   session: %s rec: %s';
    
    set(sa.fig.h,  'Position', 0.9*screen([3,4,3,4])/2, ...
        'MenuBar',  'none', ...
        'ToolBar',  'none', ...
        'Name',     sprintf(sa.fig.title_str, '', '', '',''),...
        'NumberTitle',  'off')
        
end

function init_mouse_rec()
global sa ch

    fig = sa.fig.h;
    
    ch.h = uimenu(fig, 'Label', 'choose');
    
    ch.disk    = init_menu('disk',    @disk_list,    @disk_choose);
    ch.mouse   = init_menu('mouse',   @mouse_list,   @mouse_choose);
    ch.session = init_menu('session', @session_list, @session_choose);
    ch.record  = init_menu('record',  @record_list,  @record_choose);

    
    [~,computerName]=system('hostname');
    
    local_disk  = '/Volumes/flipper/spikefolder';
    
    if strcmp(cellstr(computerName),'flipper')
        local_disk  = '/experiment';
        ch.disk     = fullfile(local_disk,'pr_data')
        sa.fn = file_names()
    else
        local_disk  = '/Volumes/spikefolder';
        ch.disk     = fullfile(local_disk,'pr_data')
        includePath=fullfile(fileparts(pwd),'current','include');
        addpath(includePath);
        sa.fn = file_names()
    end
    
    set_menu('mouse');
    choose_item(ch.mouse.item_menu(1));
    
    function m = init_menu(name, item_list, item_choose)
        m.h = uimenu(ch.h, 'Label', name);
        m.item_list   = item_list;
        m.item_choose = item_choose;
        m.item_menu   = [];
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function set_menu(name)
        delete_menu(name);
        item = ch.(name).item_list();
        for k = 1:length(item)
            ch.(name).item_menu(k) = uimenu(ch.(name).h, ...
                'Label',    item{k}, ...
                'UserData', struct('name', name, 'num', k), ...
                'Callback', @choose_item);
        end
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function delete_menu(name)
            for ki = 1:length(ch.(name).item_menu)
               if ishandle(ch.(name).item_menu(ki))
                    delete(ch.(name).item_menu(ki));
               end
            end
            ch.(name).item_menu = [];
        end
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function item = disk_list()
        ii = 0;     item = cell(0);
        list = dir('/home/zeke/')
        for il = 1:length(list)
            fold = fullfile('/home/zeke/', list(il).name);
            if exist(fold, 'dir')
                ii = ii + 1;
                item{ii} = list(il).name;
            end
        end
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function item = mouse_list()
        ii = 0; item = cell(0);
        list = dir(sa.fn.data_fold);
        for il = 1:length(list)
            if strcmp(list(il).name(1), 'm')
                fold = fullfile(sa.fn.data_fold, list(il).name);
                m_list = dir(fold);
                if ~isempty(findstr([m_list.name], 'rsm'))
                    ii = ii + 1;
                    item{ii} = strtok(list(il).name(2:end), '_');
                end
            end
        end
        item = unique(item);
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function item = session_list()
        ii = 0; item = cell(0);
        list = dir(sa.fn.data_fold)
        str = sprintf('m%s_', sa.mouse); l_str = length(str);
        for il = 1:length(list)
            if ~isempty(findstr(list(il).name, str))
                fold = fullfile(sa.fn.data_fold, list(il).name);
                m_list = dir(fold);
                 if ~isempty(findstr([m_list.name], 'rsm'))
                    ii = ii + 1;
                    item{ii} = list(il).name(l_str+1:end);
                end
            end
        end
    end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function item = record_list()
        %disp('making record list')
        ii = 0; item = cell(0);
        list = dir(sa.fn.data_fold);
        str = sprintf('m%s_%s', sa.mouse,sa.session); l_str = length(str);
        for il = 1:length(list)
            if ~isempty(findstr(list(il).name, str))
                fold = fullfile(sa.fn.data_fold, list(il).name);
                r_list = dir(fullfile(fold,'*rsm.mat'))
                for irl=1:numel(r_list)
                    ii = ii + 1;
                    base=strsplit(r_list(irl).name,'_rsm.mat');
                    recstr=strsplit(base{1},str(2:end));
                    item{ii} = recstr{end}(2:end);
                end
            end
        end
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function choose_item(hObj,~)
        ud    = get(hObj, 'UserData');
        label = get(hObj, 'Label')
        name  = ud.name
        num   = ud.num;
        
%         if ~strcmp(sa.(name), label)
            for ki = 1:length(ch.(name).item_menu);
                set(ch.(name).item_menu(ki), 'Checked', 'off');
            end
            sa.(name) = label;
            disp(['sa.name=' sa.(name)])
            set(hObj, 'Checked', 'on');
            sa.fn = file_names();
            set(sa.fig.h, 'Name', sprintf(sa.fig.title_str, sa.disk, sa.mouse, sa.session, sa.record));
            
            switch name
                case 'disk'
                    set_menu('mouse')
                     choose_item(ch.mouse.item_menu(1));
                case 'mouse'
                    disp('mouse selected')
                    set_menu('session')
                    disp('going choose session now')
                    choose_item(ch.session.item_menu(1));
                case 'session'
                    set_menu('record')
                    disp('session selected')
                    choose_item(ch.record.item_menu(1));
                case 'record'
                    disp('record selected')
                    init_record()
            end
%         end
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function fn = file_names()
        sess = str2double(sa.session);
        fn.disk          = fileparts(ch.disk);
        fn.data_fold     = fullfile(ch.disk);
        fn.sess_fold     = fullfile(fn.data_fold, sprintf('m%s_%s', sa.mouse, sa.session));
        fn.rsm_data      = fullfile(fn.sess_fold, sprintf('%s_%s_%s_rsm.mat', sa.mouse, sa.session,sa.record));
        
        fn.analysis_fold = fullfile(fn.disk, 'Experiment', 'analysis', 'sniff_analysis');
        fn.sniff_waveforms = fullfile(fn.sess_fold, sprintf('%s_%s_%s_sniff.mat', sa.mouse, sa.session,sa.record));
%         fn.sniff_waveforms = fullfile(fn.analysis_fold, 'data', sprintf('sniff_%s_%s.mat', sa.mouse, sa.session));
        fn.spikes        = fullfile(fn.sess_fold, sprintf('%s_%s_spikes.mat', sa.mouse, sa.session));
        fn.sess          = fullfile(fn.sess_fold, sprintf('sess_%s_%s.mat', sa.mouse, sa.session));
        
       
    end
end

function init_subplots()
global sa
    sa.plot.sniff = subplot('position', [0.20 0.12  0.6  0.83], 'NextPlot', 'add');
    sa.plot.hist  = subplot('position', [0.82 0.12  0.15 0.83], 'Nextplot', 'add');

end

function init_time_menu()
global sa
    pos_fig = get(sa.fig.h, 'Position');
    pos_sp  = get(sa.plot.sniff, 'Position');

    x0 = pos_sp(1)*pos_fig(3);
    xw = (pos_sp(1) + pos_sp(3)/2)*pos_fig(3);
    xe = (pos_sp(1) + pos_sp(3))*pos_fig(3);
    scale = 1./pos_fig([3,4,3,4]);

    sa.start_time = 0;
    sa.window     = 10;

    start_time_h = uicontrol(sa.fig.h, 'Style', 'Edit', ...
        'Units',        'normalized',   ...
        'Position',     [x0-35, 20, 70, 30].*scale, ...
        'String',       sprintf('%d', sa.start_time), ...
        'Callback',     @change_plot,   'UserData',     'start time');

    uicontrol(sa.fig.h,   'Style',  'PushButton', ...
        'Units',        'normalized', ...
        'Position',     [x0-70, 20, 30, 30].*scale, ...
        'String',       '<',    'FontWeight',   'bold', ...
        'CallBack',     @change_plot,   'UserData',     '<');

    uicontrol(sa.fig.h,   'Style',  'PushButton', ...
        'Units',        'normalized', ...
        'Position',     [x0-105, 20, 30, 30].*scale, ...
        'String',       '<<',    'FontWeight',   'bold', ...
        'CallBack',     @change_plot,   'UserData',     '<<');

    uicontrol(sa.fig.h,   'Style',  'PushButton', ...
        'Units',        'normalized', ...
        'Position',     [x0+40, 20, 30, 30].*scale, ...
        'String',       '>',    'FontWeight',   'bold', ...
        'CallBack',     @change_plot,   'UserData',     '>');

    uicontrol(sa.fig.h,   'Style',  'PushButton', ...
        'Units',        'normalized', ...
        'Position',     [x0+75, 20, 30, 30].*scale, ...
        'String',       '>>',    'FontWeight',   'bold', ...
        'CallBack',     @change_plot,   'UserData',     '>>');

    window_h = uicontrol(sa.fig.h, 'Style', 'Edit', ...
        'Units',        'normalized', ...
        'Position',     [xw-30, 20, 60, 30].*scale, ...
        'String',       sprintf('%d', sa.window), ...
        'CallBack',     @change_plot,   'UserData',     'window');


    end_time_h = uicontrol(sa.fig.h, 'Style', 'Text', ...
        'Units',        'normalized', ...
        'Position',     [xe-35, 20, 70, 30].*scale, ...
        'String',       sprintf('%d', sa.start_time + sa.window));

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function change_plot(hObj,~)
        str = get(hObj, 'UserData');
        if isfield(sa.data, 'Y')&&(~isempty(sa.data.Y))
            t_max = length(sa.data.Y)/1000;
        else
            t_max = sa.window;
        end
        t0    = sa.start_time;
        tw    = sa.window;
        switch str
            case 'start time'
                sa.start_time = str2double(get(hObj, 'String'));
            case '>'
                t0  = sa.start_time + sa.window;
                sa.start_time = min([t_max - tw, t0]);
            case '>>'
                sa.start_time = t_max - sa.window;
            case '<'
                t0  = sa.start_time - sa.window;
                sa.start_time = max([0, t0]);
            case '<<'
                sa.start_time = 0;
            case 'window'
                tw = str2double(get(hObj, 'String'));
                if tw == 0, tw = sa.window; end
                sa.window = min([tw, t_max - t0]);
        end
        
        set(start_time_h, 'String', sprintf('%7.3f', sa.start_time));
        set(end_time_h,   'String', sprintf('%7.3f', sa.start_time+sa.window));
        set(window_h,     'String', sprintf('%7.3f', sa.window));
        refresh()
    end
end

function init_record()
global sa
    
%     q = load(sa.fn.rsm_data, 'Sniff_DC');
%     sa.data.Y = q.Sniff_DC';
    sa.q = load(sa.fn.rsm_data, 'Sniff');
    sa.data.Y = cast(sa.q.Sniff,'double');
    
    if isfield(sa, 'fcn')
        fcn = fieldnames(sa.fcn);
        for kf = 1:length(fcn)
            if isfield(sa.fcn.(fcn{kf}), 'comp')
                sa.fcn.(fcn{kf}).computed = 0;
                if isfield(sa.fcn.(fcn{kf}), 'draw_menu')
                    set(sa.fcn.(fcn{kf}).draw_menu, 'Enable', 'off', ...
                        'Checked', 'off')
                end
            end
        end
    end   
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function s = read_rsm_data(fnam)
        if ~exist(fnam, 'file')
            warndlg('file does not exist')
            return
        end

        fid = fopen(fnam, 'r');
        [s, n_read] = fread(fid, 'int16');
        fclose(fid);
    end
    
end

function init_param()
global sa
    pos_fig = get(sa.fig.h, 'Position');
    scale = 1./pos_fig([3,4,3,4]);

    set_param(pos_fig(4)-100, 'ampl', 30000);
    set_param(pos_fig(4)-150, 'upper_thr', 5000);
    set_param(pos_fig(4)-200, 'lower_thr', -5000);
    set_param(pos_fig(4)-250, 'hist_bin',  16);
    set_param(pos_fig(4)-300, 'lp_filt',   20);
    set_param(pos_fig(4)-350, 'hp_filt',   0.2);

        function set_param(y0, name, value)
            uicontrol(sa.fig.h, 'Style',  'Text', ...
                'Units',    'normalized',...
                'Position', [10, y0, 60, 30].*scale, ...
                'String',   sprintf('%s :',name), ...
                'HorizontalAlignment',  'right');
            uicontrol(sa.fig.h, 'Style',  'Edit', ...
                'Units',    'normalized',...
                'Position', [70, y0, 60, 30].*scale, ...
                'String',   sprintf('%d', value), ...
                'UserData', name,   'CallBack', @change_param);
            sa.par.(name) = value;
        end

        function change_param(hObj, ~)
            par_name = get(hObj, 'Userdata');
            sa.par.(par_name) = str2double(get(hObj, 'String'));
        end
end

function init_fcn()
global sa

    fm.comp = uimenu(sa.fig.h, 'Label', 'function');
    fm.draw = uimenu(sa.fig.h, 'Label', 'draw');
    uimenu(sa.fig.h, 'Label', 'refresh plot', 'CallBack', @refresh);  
    
    init_raw_data()
    init_filt_data()
    init_hist()
    init_hist_max()
    init_subtract_dc();
    init_thresh()
    init_invert()
    init_peak()
    init_waveforms()
    init_parabola_fit()
    init_add_spikes()
    init_add_stim_state()
    init_save_waveforms()
    

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function init_raw_data
        create_fcn('raw_data',  'draw', @draw, 'raw_data');
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function draw()
            start = round(sa.start_time*1000);
            wind  = round(sa.window*1000);
            in = start + [1:wind];
            figure(sa.fig.h);   subplot(sa.plot.sniff);
            plot([1:wind]/1000, sa.data.Y(in), 'b')
        end
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function init_filt_data()
        create_fcn('filt_data', 'comp', @comp, 'filtered data',...
                                'draw', @draw, 'filtered data');
                            
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function comp(~,~)
            disp('filtered data')
            lp = sa.par.lp_filt;
            hp = sa.par.hp_filt;
            if (hp==0)&&(lp==0)
                sa.data.Yfilt = sa.data.Y;
            else
                if (lp == 0)&&(hp > 0)
                    [bf, af] = butter(2, hp/500, 'high');
                end
                if (lp > 0)&&(hp == 0)
                    [bf, af] = butter(2, lp/500, 'low');
                end
                if (lp > 0)&&(hp > 0)
                    [bf, af] = butter(2, [hp, lp]/500);
                end
    
                sa.data.Yfilt = filtfilt(bf, af, sa.data.Y);
            end
            enable_draw('filt_data')
%             sa.fcn.filt_data.computed = 1;
%             set(sa.fcn.filt_data.draw_menu, 'Enable', 'on')
        end
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function draw()
            start = round(sa.start_time*1000);
            wind  = round(sa.window*1000);
            in = start + [1:wind];
            figure(sa.fig.h);   subplot(sa.plot.sniff);
            plot([1:wind]/1000, sa.data.Yfilt(in), 'k')
        end
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function init_hist()
        create_fcn('hist', 'comp', @comp, 'histogram', ...
                           'draw', @draw, 'histogram');
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function comp(~,~)
%             if ~sa.fcn.filt_data.computed
%                 sa.fcn.filt_data.comp()
%             end
            disp('histogram')
%             y_axis = [-2^15:2^15];
            sa.data.hist_axis = min(sa.data.Y):max(sa.data.Y);
            sa.data.hist      = hist(sa.data.Y, sa.data.hist_axis);
            sa.data.hist([1:100,end-99:end]) = 0;
            enable_draw('hist')
            
        end
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function draw()
            figure(sa.fig.h), subplot(sa.plot.hist)
            y_axis = sa.data.hist_axis;
            ny  = length(y_axis);
            bin = sa.par.hist_bin;
            ny  = bin*floor(ny/bin);
            y_axis = mean(reshape(y_axis(1:ny),       bin, ny/bin), 1);
            y_hist =  sum(reshape(sa.data.hist(1:ny), bin, ny/bin),1);
            plot(y_hist, y_axis, 'b')
        end
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function init_hist_max()
        create_fcn('hist_max', 'comp', @comp, 'hist_max', ...
                               'draw', @draw, 'hist_max');
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function comp(~,~)
            
            if ~sa.fcn.hist.computed
                sa.fcn.hist.comp()
            end
            disp('histogram maximum position')
            h_axis = sa.data.hist_axis;
            hh     = sa.data.hist;
            % finding maxima of the histogram
            [~, im] = max(hh);
            % parabola fit
            in_fit = im+[-20:20];
            p = polyfit(h_axis(in_fit)-h_axis(im), hh(in_fit),2);
            ym = h_axis(im) - p(2)/p(1)/2;
            
            
            sa.data.Y_zero = ym;
            enable_draw('hist_max')
            
        end
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function draw()
            figure(sa.fig.h), subplot(sa.plot.hist)
            plot([0, max(sa.data.hist)], sa.data.Y_zero*[1,1], 'r', 'LineWidth', 2)
        end
    end
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function init_subtract_dc()
        create_fcn('subtract_dc', 'comp', @comp, 'subrtract_dc')
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function comp(~,~)
            
            if ~sa.fcn.hist_max.computed
                sa.fcn.hist_max.comp()
            end
            disp('subtract dc')
            
            
            sa.data.Y = sa.data.Y - sa.data.Y_zero;
            
        end
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function init_thresh()
        create_fcn('thresh', 'draw', @draw, 'thresholds');
        function draw()
            pl = fieldnames(sa.plot);
            for ip = 1:length(pl)
                x_lim = get(sa.plot.(pl{ip}), 'XLim')';
                subplot(sa.plot.(pl{ip}))
                plot(x_lim*[1,1], [1;1]*[sa.par.lower_thr, sa.par.upper_thr], 'r')
            end
        end
            
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function init_invert()
        create_fcn('invert', 'comp', @comp, 'invert')
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function comp(~,~)
            disp('invert')
            sa.data.Y = - sa.data.Y;
            if sa.fcn.filt_data.computed
                sa.data.Yfilt = -sa.data.Yfilt;
            end
            if sa.fcn.hist.computed
                sa.data.hist = fliplr(sa.data.hist);
            end
        end
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function init_peak()
        create_fcn('peak', 'comp', @comp, 'find peaks', ...
                           'draw', @draw, 'peaks');
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function comp(~,~)
            if ~sa.fcn.filt_data.computed
                sa.fcn.filt_data.comp();
            end
            disp('find peaks')
            y = sa.data.Yfilt;
            in_min = find(diff(sign(diff(y.*(y > sa.par.upper_thr))))==-2) + 1;
            in_max = find(diff(sign(diff(y.*(y < sa.par.lower_thr))))== 2) + 1;
            sa.data.peaks = sort([in_min'; in_max']);
            enable_draw('peak')
        end
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function draw()
            start = round(sa.start_time*1000);
            wind  = round(sa.window*1000);
            peaks = sa.data.peaks;
            in = find((peaks > start).*(peaks < start + wind));
            subplot(sa.plot.sniff)
            plot((peaks(in) - start)/1000, sa.data.Yfilt(peaks(in)), 'xr');
        end
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function init_waveforms()
        create_fcn('waveforms', 'comp', @comp, 'waveforms', ...
                                'draw', @draw, 'waveforms')
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function comp(~,~)
            global t_zer
            if ~sa.fcn.peak.computed
                sa.fcn.peak.comp()
            end
            disp('waveforms')
            y     = sa.data.Yfilt;
            peaks = sa.data.peaks;
            
            max_sniff = 1000;
            margin    = 150;
            
            in_sn = find(diff(sign(y(peaks))) == 2);
            n_sniff = length(in_sn)-1;
            
            t_min = peaks(in_sn);      y_min = y(t_min);
            t_max = peaks(in_sn+1);    y_max = y(t_max);
            
            t_zer   = zeros(1, n_sniff+1);
            for is = 1:n_sniff+1
                t_zer(is) = t_min(is) + find(y(t_min(is):t_max(is))>0, 1, 'first');
%                 t_zer(is) = larger_than_zero(t_min(is), t_max(is), 'first');
            end
            
            t_zer = [1, t_zer, length(y)];

            ks = 0;
            
            gw = waitbar(0, 'waveforms');
            for is = 2:n_sniff+1
                if t_zer(is+1)-t_zer(is) < max_sniff
                    ks = ks + 1;
                    t0   = max([t_zer(is-1), t_zer(is) - margin]);
                    tend = min([t_zer(is+1)+margin, t_zer(is+2)]);
                    t_middle = t_zer(is) + find(y(t_zer(is):t_zer(is+1))<0, 1, 'first') - 1;
                    
                
                    sn(ks).t0       = t0;
                    sn(ks).waveform = y(t0+1:tend);
                    sn(ks).nt       = tend - t0;
                    
                    sn(ks).t_zer      = [t_zer(is), t_middle, t_zer(is+1)]' - t0;
%                     sn(ks).t_zer_down = t_middle - t0;
                
                    in1 = 1:sn(ks).t_zer(1);
                    in2 = (sn(ks).t_zer(1)+1):sn(ks).t_zer(3);
                    in3 = (sn(ks).t_zer(3)+1):(tend-t0);
    
                    [sn(ks).y_min(1,1), tmi1] = min(sn(ks).waveform(in1));
                    [sn(ks).y_max(1,1), tma1] = max(sn(ks).waveform(in2));
                    [sn(ks).y_min(2,1), tmi2] = min(sn(ks).waveform(in2));
                    [sn(ks).y_max(2,1), tma2] = max(sn(ks).waveform(in3));
                    
                    sn(ks).t_min = [tmi1;tmi2+sn(ks).t_zer(1)];
                    sn(ks).t_max = [tma1+sn(ks).t_zer(1); tma2+sn(ks).t_zer(2)];
                elseif t_zer(is+1)-t_zer(is) > max_sniff
                    err_msg = sprintf('sniff number %i comes greater than 1000ms after sniff %i',is,is-1);
                    disp(err_msg)
                end
                if is/100 == round(is/100)
                    waitbar((is-1)/n_sniff, gw);
                end
            end
            
            close(gw)
            sa.data.sniff = sn(1:ks);
            enable_draw('waveforms')
           
%             function t = larger_than_zero(t1,t2, first_last)
%                 t = t1 + find(y(t1:t2)>0, 1, first_last);
%             end
%             
        end
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function draw(~,~)
            start = round(sa.start_time*1000);
            wind  = round(sa.window*1000);
            t_zer = ones(3,1)*[sa.data.sniff.t0] + [sa.data.sniff.t_zer];
            t_min = ones(2,1)*[sa.data.sniff.t0] + [sa.data.sniff.t_min];
            t_max = ones(2,1)*[sa.data.sniff.t0] + [sa.data.sniff.t_max];
            y_zer1 = sa.data.Yfilt(t_zer(1,:));
            y_zer2 = sa.data.Yfilt(t_zer(2,:));
            y_min = [sa.data.sniff.y_min];
            y_max = [sa.data.sniff.y_max];
            
            in = find((t_zer(1,:) > start).*(t_zer(1,:) < start + wind));
            
            subplot(sa.plot.sniff)
            plot((t_zer(1,in) - start)/1000, y_zer1(in), 'om', 'LineWidth', 2);
            plot((t_zer(2,in) - start)/1000, y_zer2(in), 'om', 'LineWidth', 2);
            
            plot((t_min(1,in) - start)/1000, y_min(1,in), 'xg', 'LineWidth', 2);
            plot((t_max(1,in) - start)/1000, y_max(1,in), 'xg', 'LineWidth', 2);
            
        end
        
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function init_parabola_fit()
        create_fcn('parab_fit', 'comp', @comp, 'parabola fit', ...
                                'draw', @draw, 'parabolas')
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function comp(~,~)
            if ~sa.fcn.waveforms.computed
                sa.fcn.waveforms.comp()
            end
            disp('parabola fit')
            sn = sa.data.sniff;
            ns = length(sn);
            gb = waitbar(0, 'fitting parabolas');
            for is = 1:ns
                
                % parabola fit for the first maxima and second minima
                nn   = min([40, round(0.8*(sn(is).t_max(1) - sn(is).t_zer(1)))]);
                [pp1, xx1] = parabola_fit(sn(is).t_max(1), nn);

                nn   = min([40, round(0.8*(sn(is).t_zer(3) - sn(is).t_min(2)))]);
                [pp2, xx2] = parabola_fit(sn(is).t_min(2), nn);
                
                if pp2(1) < 0
                    xx2 = [xx1(2), sn(is).t_zer(3)];
                end
                
                if xx2(1) < xx1(2)
                    xx2(1) = (xx1(2)+xx2(1))/2;
                    xx1(2) = xx2(1);
                end
                sn(is).pp_fit = {pp1, pp2};
                sn(is).t_zer_fit = [xx1 xx2];
                
                    
%                 
%                 [sn(is).pp_fit{1}, xx] = parabola_fit(sn(is).t_max(1), nn);
%                 sn(is).t_zer_fit([1,2]) = xx;
%         
%                 nn   = min([40, round(0.8*(sn(is).t_zer(3) - sn(is).t_min(2)))]);
%                 [pp, xx] = parabola_fit(sn(is).t_min(2), nn);
%                     
%                     
%                 [sn(is).pp_fit{2}, xx] = parabola_fit(sn(is).t_min(2), nn);
%                 sn(is).t_zer_fit([3,4]) = xx;
%                 
%                 if sn(is).t_zer_fit(3) < sn(is).t_zer_fit(2)
%                     sn(is).t_zer_fit([2,3]) = mean(sn(is).t_zer_fit([2,3]))+ [-1,1];
%                 end
                if is/100 == round(is/100)
                    waitbar(is/ns, gb);
                end
            end
            
            
            sa.data.sniff = sn;

            close(gb)
            enable_draw('parab_fit')

            function [pp, xx]= parabola_fit(x0, np)
                x = -np:np;
%                 size(x)
%                 size(sn(is).waveform(x0+x))
                pp = polyfit(x, sn(is).waveform(x0+x),2);
                d  = pp(2)/2/pp(1) - pp(3)/pp(1);
                xx = x0 - pp(2)/2/pp(1) + sqrt(d)*[-1,1];
            end
    
        end
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function draw(~,~)
            start = round(sa.start_time*1000);
            wind  = round(sa.window*1000);
            
            t_zer = ones(3,1)*[sa.data.sniff.t0] + [sa.data.sniff.t_zer];
            in = find((t_zer(1,:) > start).*(t_zer(1,:) < start + wind));
            
            for is = in
                sn = sa.data.sniff(is);
                x     = -50:50;
                x0(1) = sn.t_max(1);
                x0(2) = sn.t_min(2);
                for ip =1:2
                    pp  = sn.pp_fit{ip};
                    yy  = pp(1).*x.^2 + pp(2)*x + pp(3);
                    xx0 = x0(ip) + sn.t0 - start;
                    plot((xx0+x)/1000, yy, 'r')
                end
                plot((sn.t_zer_fit + sn.t0 - start)/1000, zeros(1,4), 'xr', 'LineWidth', 2)
                    
            end
        end
    end

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function init_add_spikes()
        create_fcn('add_spikes', 'comp', @comp, 'add spikes')

        function comp(~,~)
            disp('add spikes')
            if ~exist(sa.fn.spikes, 'file')
                disp('no spike file')
                return
            end
            
            % loading spikes
            q = load(sa.fn.spikes);
            nu   = length(q.unit);
            for iu = 1:nu
                spikes{iu} = q.unit(iu).spikeTimes;
            end
            sn = sa.data.sniff;

            % sorting spike for sniffs
            for is = 1:length(sn)
                t0   = sn(is).t0;
                tend = t0 + sn(is).nt;
    
                for iu = 1:nu
                    in = find((spikes{iu} > t0).*(spikes{iu} < tend));
                    sn(is).spikes{iu} = spikes{iu}(in) - t0;
                end
    
                if is/100 == round(is/100)
                    disp(is)
                end
            end
            
            sa.data.sniff = sn;
            set(sa.fcn.add_spikes.comp_menu, 'Checked', 'on')
        end
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function init_add_stim_state()
        create_fcn('add_stim_state', 'comp', @comp, 'add stimulus info', ...
                                     'draw', @draw, 'stimulus state')
        
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function comp(~,~)
            disp('add stim state')
            sn = sa.data.sniff;
            
            if ~exist(sa.fn.sess, 'file')
                disp('no trial info')
                return
            end
            
            q  = load(sa.fn.sess);
            tr = q.trial;

            for is = 1:length(sn)
                sn(is).state = 0;
            end
            % defining state for each sniff 0 - no stimulus, 1 - stimulus
            post_stim = 1000;

            for it = 1:length(tr)
                stim = [10000, 0];
                if isfield(tr, 'odorTime')&&(diff(tr(it).odorTime) > 0)
                    stim(1) = min([stim(1), tr(it).odorTime(1)]);
                    stim(2) = max([stim(2), tr(it).odorTime(2)]);
                end
    
                if isfield(tr, 'laserTime')&&(diff(tr(it).laserTime) > 0);
                    stim(1) = min([stim(1), tr(it).laserTime(1)]);
                    stim(2) = max([stim(2), tr(it).laserTime(2)]);
                end
    
                if stim(2) > stim(1);
                    t0   = tr(it).start + stim(1);
                    tend = tr(it).start + stim(2) + post_stim;
        
                    for is = find(([sn.t0]+[sn.nt] > t0).*([sn.t0] < tend))
                        sn(is).state = 1;
                    end
                    disp(it)
                end
            end
            sa.data.sniff = sn;
            enable_draw('add_stim_state')
        end
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function draw(~,~)
            start = round(sa.start_time*1000);
            wind  = round(sa.window*1000);
            t_zer = ones(3,1)*[sa.data.sniff.t0] + [sa.data.sniff.t_zer];
            
            in = find((t_zer(1,:) > start).*(t_zer(1,:) < start + wind));
            y0 = sa.par.ampl*0.95;
            subplot(sa.plot.sniff)
            for is = in
                if sa.data.sniff(is).state
                    plot((t_zer(:,is)-start)/1000, y0*ones(3,1), 'r', 'LineWidth', 2)
                end
            end
            
        end
        
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function init_save_waveforms()
        create_fcn('save_waveforms', 'comp', @comp, 'save waveforms')
        set(sa.fcn.save_waveforms.comp_menu, 'Separator', 'on')
        
        function comp(~,~)
            if ~sa.fcn.waveforms.computed
                sa.fcn.waveforms.comp();
            end
            disp(['save waveforms to ' sa.fn.sniff_waveforms])
            sniff = sa.data.sniff;
            save(sa.fn.sniff_waveforms, 'sniff')
        end
    end



    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function create_fcn(name, varargin)
        ii = 0;
        while ii < length(varargin)
            ii = ii+1;  str = varargin{ii};
            ii = ii+1;  fcn = varargin{ii};
            ii = ii+1;  label = varargin{ii};
            switch str
                case 'comp'
                    sa.fcn.(name).comp_menu = uimenu(fm.comp, ...
                        'Label',    label, ...
                        'CallBack', fcn);
                    sa.fcn.(name).computed = 0;
                case 'draw'
                    sa.fcn.(name).draw_menu = uimenu(fm.draw, ...
                        'Label',    label, 'CallBack', @check_draw, ...
                        'Enable',   'off', 'UserData',  name);
            end
            sa.fcn.(name).(str) = fcn;
        end
        if ~isfield(sa.fcn.(name), 'comp')
            set(sa.fcn.(name).draw_menu, 'Enable', 'on');
        end
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function check_draw(hObj,~)
        if strcmp(get(hObj, 'Checked'), 'on')
            set(hObj, 'Checked', 'off')
        else
            set(hObj, 'Checked', 'on')
        end
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function enable_draw(name)
        sa.fcn.(name).computed = 1;
        set(sa.fcn.(name).draw_menu, 'Enable', 'on')
        set(sa.fcn.(name).comp_menu, 'Checked', 'on')

    end



end

function refresh(~,~)
global sa
    % clear current plot
    pl = fieldnames(sa.plot);
    for ip = 1:length(pl)
        child = get(sa.plot.(pl{ip}), 'Children');
        for ic = 1:length(child)
            delete(child(ic))
        end
    end
    
    % plot new graphs
    fcn = fieldnames(sa.fcn);
    for kf = 1:length(fcn)
        if isfield(sa.fcn.(fcn{kf}), 'draw')&&strcmp(get(sa.fcn.(fcn{kf}).draw_menu, 'Checked'), 'on')
            sa.fcn.(fcn{kf}).draw()
        end
    end
    
    % set limits
    set(sa.plot.sniff, 'YLim', sa.par.ampl*[-1,1], ...
        'XLim',         [0, sa.window], ...
        'Xgrid',        'on', ...
        'YGrid',        'on');
    set(sa.plot.hist, 'YLim', sa.par.ampl*[-1,1], ...
        'YTickLabel',   [], ...
        'YGrid',        'on');
end