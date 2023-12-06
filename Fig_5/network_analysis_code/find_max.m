% detect maximum fsp at every event
function [ind_max] = find_max(nn, ind_event, max_det_win, bins, f_sp, flag_plot)
ind_max = [];
for jj=1:length(ind_event{nn})
    [max1, ind] = max(f_sp{nn}(ind_event{nn}(jj):1:ind_event{nn}(jj)+max_det_win));
    ind_max(jj) = ind_event{nn}(jj) + ind - 1;
    % do you want to see the maxima?
    if flag_plot == 1
        figure(441+nn); plot(bins{nn}, f_sp{nn}, '-ob'); hold on;
        plot(bins{nn}(ind_max), f_sp{nn}(ind_max), 'or', 'linewidth', 4, 'markersize', 15);
    end
end
