function [window_left, window_right] = bin_win_fit(nn, win, bins, ind_max)
window_left = [];
window_right = [];
for jj=1:length(ind_max)
    % left limit...
    if ind_max(jj) - win < 1
        window_left(jj) = 1;
    else
        window_left(jj) = ind_max(jj) - win;
    end
    % right limit...
    if ind_max(jj) + win > length(bins{nn})
        window_right(jj) = length(bins{nn});
    else
        window_right(jj) = ind_max(jj) + win;
    end
end
