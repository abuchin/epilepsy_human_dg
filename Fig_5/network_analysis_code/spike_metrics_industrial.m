function [post, spikeTimes, spiked_v, ct_ind, bins, f_sp_ct, f_sp] = spike_metrics_industrial(path2, spike_f, bin_dt, pre, SaveFile_Flag, CT)

%% Load data...
if nargin < 2
    path2 = './100/';
    spike_f = 'spikes.csv';
elseif nargin == 2
    bin_dt = 25;        % binning time step
elseif nargin == 3
    pre = 0;
elseif nargin == 4
    % Do you want to save the output data file?
    SaveFile_Flag = 1;
elseif nargin == 5
    CT_index = 2;
    if CT_index == 1
        ['CTs accounted for are 1']
        CT = cell(1,1);
        CT{1} = [1  505];
        color_v = ['k'];
    elseif CT_index == 2
        ['CTs accounted for are 2']
        CT = cell(2,1);
        CT{1} = [1 499];
        CT{2} = [500 505];
        color_v = ['b' 'r'];
    elseif CT_index == 3
        ['CTs accounted for are 3']
        CT = cell(3,1);
        CT{1} = [1 160];
        CT{2} = [161 500];
        CT{3} = [500 505];
        color_v = ['b' 'g' 'r'];
    end
end

%% start loading data...
spike_data = load([path2 spike_f]);
spike_times = spike_data(:,1);
cell_ID = spike_data(:,2);

% Define spikeTimes cell array...
uni_v = unique(cell_ID);
spikeTimes_all = cell(length(uni_v),1);
spikeTimes = cell(length(uni_v),1);
['Number of cells with at least one spike... ' num2str(length(uni_v))]

for ii=1:1:length(uni_v)
    clear ind
    ind = find(cell_ID == uni_v(ii));
    % Cell output definition in ms!
    spikeTimes_all{uni_v(ii) + 1} = spike_times(ind)';
end

% Exclude cells that did not spike...
spiked_v = find(~cellfun('isempty', spikeTimes_all)); 
spikeTimes = spikeTimes_all(spiked_v,:);
ct_ind = cell(length(CT),1);
if length(CT)>1
    for j=1:length(CT)
        ct_ind{j} = find(spiked_v >= CT{j}(1) & spiked_v <= CT{j}(2));
    end
elseif length(CT)==1
    ct_ind = spiked_v;
end

%% calculate a few spiking metrics...
post = ceil(max([spikeTimes{:}]));
bins = linspace(pre, post, round((post-pre)/bin_dt));

f_sp = zeros(length(bins),1);
f_sp_ct = zeros(length(bins),length(CT));

for i = 1:length(bins)-1
    clear sp_ind
    sp_ind = find([spikeTimes{:}] >= bins(i) & [spikeTimes{:}] < bins(i+1));
    f_sp(i) = length(sp_ind)/bin_dt/1e-3/length(spiked_v);    % Spike frequency in Hz!
    %  Cell type specific metrics...
    if length(CT)>1
        for j=1:length(CT)
            clear sp_ind_ct
            sp_ind_ct = find([spikeTimes{ct_ind{j}}] >= bins(i) & [spikeTimes{ct_ind{j}}] < bins(i+1));
            f_sp_ct(i,j) = length(sp_ind_ct)/bin_dt/1e-3/length(ct_ind{j});    % Spike frequency in Hz!
        end
    end
end

%% saving output files...
if SaveFile_Flag == 1
    out_name = char([path2 'raster.mat']);
    save(out_name, 'spikeTimes');
end








