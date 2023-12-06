function master_spike_metrics(composition_v)

%% analyze network simulation results focusing around interictal generation in WG1 vs. 4

if nargin<1
    composition_v = cell(1001,1);
    composition_v = {'0' '100' '200'};
end

% extract the composition values
compo_v = zeros(size(composition_v,2),1);

% define a few paths and some parameters
path_s = cell(length(composition_v),1);
for kk = 1:length(path_s)
    path_s{kk} = ['./' num2str(composition_v{kk}) '/'];
    compo_v(kk) = str2num(composition_v{kk});
end
spike_f = 'spikes.csv';

% some more parameters...
bin_dt = 10;                    % binning time step for assessing spike frequency
pre = 0;
sigma_times = 1.5;              % sigma for defining the event detection threshold
dt_event_detection = 100;       % time interval to start detection of next event
flag_plot = 0;                  % flag to plot event detection
flag_plot2 = 0;                 % flag to plot event maxima

max_detect_win = 6;             % number of points after threshold crossing to look for max
max_win = 12;                    % number of bins around the max to calculate stats

% number of cell types to be considered
CT = cell(2,1);
CT{1} = [1 499];
CT{2} = [500 505];

%% Extract some data for analysis...
thresh = zeros(length(path_s), 1);
bins = cell(length(composition_v),1);
f_sp = cell(length(composition_v),1);
f_sp_ct = cell(length(composition_v),1);
ind_event_start = cell(length(composition_v),1);
ind_event_max = cell(length(composition_v),1);
win_max_left = cell(length(composition_v),1);
win_max_right = cell(length(composition_v),1);
fit_param = cell(length(composition_v),1);

for kk = 1:length(path_s)
    [post, spikeTimes, spiked_v, ct_ind, bins{kk}, f_sp_ct{kk}, f_sp{kk}] = spike_metrics_industrial(path_s{kk}, spike_f, bin_dt, pre, 1, CT);  
    thresh(kk) = mean(f_sp{kk}) + sigma_times.*std(f_sp{kk});
    % detect threshold crossing...
    [ind] = detect_threshold_crossing(kk, thresh, dt_event_detection, bin_dt, bins, f_sp, flag_plot);
    ind_event_start{kk}(1:length(ind)) = ind;
    % detect maximum fsp at every event
    [ind_m] = find_max(kk, ind_event_start, max_detect_win, bins, f_sp, flag_plot2);
    ind_event_max{kk}(1:length(ind_m)) = ind_m;
    % determine window (bins) around maximum to fit function (gaussian)
    [win_left, win_right] = bin_win_fit(kk, max_win, bins, ind_event_max{kk});
    win_max_left{kk}(1:length(win_left)) = win_left;
    win_max_right{kk}(1:length(win_right)) = win_right;
    % fit function at every event and extract some numbers...
    [fit_params] = fit_data_func(kk, win_max_left, win_max_right, bins, f_sp, 'gauss1');
    fit_param{kk}(1:size(fit_params,1), 1:size(fit_params,2)) = fit_params;
end

%% Plot and figure parameters...
linew = 2.0;
ax_linew = 1.5;
msize = 12.0;
ax_msize = 15.0;
color_v = ['b' 'r'];

%% Attempt at plotting...
figure(3033); set(gcf,'color','w');
for kk = 1:length(path_s)
    
    [post, spikeTimes, spiked_v, ct_ind] = spike_metrics_industrial(path_s{kk}, spike_f, bin_dt, pre, 1, CT);
        
    % raster plot...
    subplot(length(path_s), 2, 2*(kk-1)+1);
    for j=1:size(ct_ind,1)
        for i = 1:1:length(ct_ind{j})
            plot(round(spikeTimes{ct_ind{j}(i)}) + pre, spiked_v(ct_ind{j}(i)),'.', 'color', color_v(j), 'MarkerSize',4); hold on;
        end
    end
    
    xlabel('time / ms','fontsize', 14);
    ylabel('Cell ID','fontsize', 14);
    box off; set(gca,'linewidth',ax_linew);
    axis([-pre-25 post+25 -25 spiked_v(length(spiked_v))+25]);
    
    subplot(length(path_s), 2, 2*(kk-1)+2);
    for j=1:size(ct_ind,1)
        plot(bins{kk}, f_sp_ct{kk}(:,j), '-', 'color', color_v(j), 'linewidth', linew); hold on;
    end
    plot(bins{kk}, f_sp{kk}, '-k', 'linewidth', linew-1.5);
    plot([-pre-25 post+25], [thresh(kk) thresh(kk)], '--', 'color', [0.3 0.9 0.3], 'linewidth', linew-0.5);
    ylabel('Spike frequency / Hz','fontsize', 14);
    xlabel('time / ms','fontsize', 14);
    box off; set(gca,'linewidth',ax_linew);
    xlim([-pre-25 post+25]);
end

% plot some event statistics... 
fit_a = zeros(length(fit_param),1); fit_o = zeros(length(fit_param),1); fit_w = zeros(length(fit_param),1);
for j = 1:length(fit_param), fit_a(j) = fit_param{j}(1); fit_o(j) = fit_param{j}(2);  fit_w(j) = fit_param{j}(3); end;

figure(3043); set(gcf,'color','w');
subplot(1,3,1);
plot(compo_v, fit_a, '-ok', 'linewidth', linew); hold on;
box off;
axis([-0.05 1.05 -5 105]);
ylabel('event spike frequency amplitude (Hz)','fontsize', 14);
xlabel('WG1 vs. WG4 composition','fontsize', 14);

subplot(1,3,2);
plot(compo_v, fit_o, '-ok', 'linewidth', linew); hold on;
box off;
ylabel('event centroid (timing) / ms','fontsize', 14);
xlabel('WG1 vs. WG4 composition','fontsize', 14);

subplot(1,3,3);
plot(compo_v, fit_w, '-ok', 'linewidth', linew); hold on;
box off;
ylabel('event duration / ms','fontsize', 14);
xlabel('WG1 vs. WG4 composition','fontsize', 14);

%% =======================================================================
% detect threshold crossing...
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

%% =======================================================================
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

%% =======================================================================
% determine window (bins) around maximum to fit function (gaussian)
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

%% =======================================================================
% fit data with function
function [fit_params] = fit_data_func(nn, win_max_left, win_max_right, bins, f_sp, fit_func)

fit_params = zeros(length(win_max_left{nn}), 3);

for j=1:length(win_max_left{nn})
    clear x y f
     x = [win_max_left{nn}(j):1:win_max_right{nn}(j)]';
    %x = bins{j}(win_max_left{nn}(j):1:win_max_right{nn}(j))';
    y = [f_sp{nn}(win_max_left{nn}(j):1:win_max_right{nn}(j))];
    
    f = fit(x, y, fit_func);
    fit_params(j, 1) = f.a1;
    fit_params(j, 2) = f.b1;
    fit_params(j, 3) = f.c1;
    
        figure; set(gcf,'color','w');
        plot(x,y); hold on;
        plot(f,x,y);
        box off;
end




