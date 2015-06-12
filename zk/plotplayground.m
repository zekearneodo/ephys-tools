function pa = plot_aids()

x=1:1000;
y=2*x;
z=3*x.^2;

pos_x = 0.3;
pos_y = 0.5;

dim_x=0.3;
dim_y=0.5;

close all
h(1)=figure(1);
set(h(1),'DefaultAxesFontSize',8)


plot(ax_rast,x,y)
plot(ax_psth,x,z)

fn.make_response_frame = @make_response_frame;

end



function frame = make_response_frame(pos_x, pos_y, dim_x, dim_y, fig)

% Makes a frame with two axes for plotting psths. Splits the frame into two
% vertical axes
% 
% pos_x, pos_y: position in the fig of the low-left corner of the frame
% dim_x, dim_y: total dimensions of the frame

if exist('fig', 'var')
    gf = figure(fig); clf
else
    gf = figure;
end

ax_rast = axes('Position',[0.05+pos_x, 0.5*dim_y + pos_y, 1*dim_x-0.075, 0.45*dim_y-0.05]);
ax_psth = axes('Position',[0.05+pos_x, 0.05 + pos_y, 1*dim_x-0.075, 0.45*dim_y-0.05]);

frame.fig  = gf;
frame.axes = [ax_rast, ax_psth];

end


