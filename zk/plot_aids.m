function pa = plot_aids()

pa.make_response_frame = @make_response_frame;

end




function frame = make_response_frame(pos_x, pos_y, dim_x, dim_y, fig)

% Makes a frame with two axes for plotting psths. Splits the frame into two
% vertical axes
% 
% pos_x, pos_y: position in the fig of the low-left corner of the frame
% dim_x, dim_y: total dimensions of the frame

% activate or create figure
if exist('fig', 'var')
    gf = figure(fig); clf
else
    gf = figure;
end
set(0, 'currentfigure', gf);  %# for figures

%create axes
ax_rast = axes('Position',[0.05+pos_x, 0.5*dim_y + pos_y, 1*dim_x-0.075, 0.45*dim_y-0.05]);
ax_psth = axes('Position',[0.05+pos_x, 0.05 + pos_y, 1*dim_x-0.075, 0.45*dim_y-0.05]);

%return the frame
frame.fig  = gf;
frame.ax_rast = ax_rast;
frame.ax_psth = ax_psth;

end


