% for plotting box in zoom in 
%
function [x_fig, y_fig] = ds2nfu(x, y)
    ax = gca;
    ax_units = ax.Units;
    ax.Units = 'normalized';
    ax_pos = ax.Position;
    ax.Units = ax_units;

    xlim = get(ax, 'XLim');
    ylim = get(ax, 'YLim');

    x_fig = ax_pos(1) + (x - xlim(1)) / diff(xlim) * ax_pos(3);
    y_fig = ax_pos(2) + (y - ylim(1)) / diff(ylim) * ax_pos(4);
end