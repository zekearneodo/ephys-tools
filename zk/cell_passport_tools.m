function cp = cell_passport_tools(doit)
%functions for making a cell passport
%uses visualize_responses
    vr = visualize_responses_new();
           
    cp.plot_response = @plot_responses;
    cp.pa = @plot_aids;
    
    if nargin>0 && doit
        global response;
        response = plot_response();
    end
end

function canvas = setup_fig()

end

%plots a raster and a psth in response to a stimulus
function response_plot = plot_response(resp)
pa = plot_aids();

pos_x = 1./4.;
pos_y = 2./7.;

dim_x=1./4.;
dim_y=1./6.;

%prepare the figure
fg = resp.fg;
set(0, 'currentfigure, fg');
set(fg,'DefaultAxesFontSize',8)

%make the frame with stim-dependent dimensions
%frame_set = get_frame_geometry(resp.stim)
frame = pa.make_response_frame(pos_x, pos_y, dim_x, dim_y, h(1));

%plot the raster
kt = resp.ntrial;
yt = max([1, floor(kt/10)*10]);

plot(frame.ax_rast, resp.x, resp.y, '.', 'MarkerSize',7);
set(frame.ax_rast, 'NextPlot', add);
plot([0,0], [0, kt+1], '--k'), hold off
set(frame.ax_rast, 'XLim', [t1, t2], ...
    'YLim',      [0, kt+1], ...
    'YTick',     [0, yt], ...
    'XTick',     [-200,0,200], ...
    'FontSize',  10);
%title(stim_str, 'FontSize', 10, 'FontWeight', 'bold') insert some text as
%title

%plot the psth
    plot(frame.ax_psth, resp.t, resp.rate, ...
        'LineWidth', 1, ...
        'Color',     resp.stim.lineColor, ...
        'LineStyle', resp.stim.lineStyle),
    set(frame.ax_rast, 'NextPlot', add);
    %plot(t,ones(size(t))*baseline,'.k')
    %plot(t,ones(size(t))*pkThresh,'b')
    if strcmp(resp.sType,'odor')
        plot(frame.ax_rast, resp.t,resp.rateBaseHist,'LineStyle','--','Color',[.5,0.5,0.5])
    end
    plot(frame.ax_rast, resp.pkTime, resp.pkRate,'k*')
    plot(frame.ax_rast, [0,0], [0, rmax], '--k'), hold off
    set(gs, 'XLim',     [t1, t2], ...
        'XTick',         [-200,0,200], ...
        'YLim',          [0, rmax], ...
        'YTick',         0:20:rmax, ...
        'FontSize',      10);


plot(frame.ax_rast,x,y)
set(frame.ax_rast, 'XTickLabel',[]);
plot(frame.ax_psth,x,z)

response_plot.frame = frame;

end