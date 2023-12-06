function [ind_event] = detect_threshold_crossing(nn, thresh, dt_event_detection, bin_dt, bins, f_sp, flag_plot)

ind_event = [];
count1 = 0;
dummy_num = 0;
tend = 1;

while tend < length(f_sp{nn})
    count1 = count1 + 1;
    clear dummy_num
    dummy_num = find(f_sp{nn}(tend:length(f_sp{nn})) > thresh(nn), 1);
    if isempty(dummy_num)~=1
        ind_event(count1) = dummy_num + tend - 1;
        tend = ind_event(count1) + round(dt_event_detection/bin_dt);
        % -------------------
        % do you want to see some plots of the event detection...
        if flag_plot == 1
            figure(435); plot(bins{nn}, f_sp{nn}, '-ob'); hold on;
            plot(bins{nn}(ind_event(count1)), f_sp{nn}(ind_event(count1)), 'or', 'linewidth', 4, 'markersize', 15);
            if tend<length(f_sp{nn})
                plot(bins{nn}(tend), f_sp{nn}(tend), 'og', 'linewidth', 4, 'markersize', 15);
            end
            plot([min(bins{nn}) max(bins{nn})], [thresh(nn) thresh(nn)], '--r'); hold on;
            pause;
            close(435);
        end
        % -------------------
    else
        % print('im in here...');
        tend = length(f_sp{nn}) + 1;
    end
end
