% Bin 1 = House Upright
% Bin 2 = Face Upright
% Bin 3 = Face Inverted
%cd('/Volumes/T7/Decoding analysis')
clc;
clear all;
addpath('~/Documents/MATLAB/eeglab2024.2.1/');
eeglab

% Generalize Functiong Across Time points
tps = [{'T1'}];

%faces vs houses = 72 v 72 = cross_blocks = 5 72/5 = 14 trials per
%upright v inverted = 72 v 72 = cross_blocks = 5 72/5 = 14 trials per
%upright v inverted = 24 v 24 v 24; blocks=2
%all bins = 24 v 24 v 24 v 24 v 24 v 24 v 24 v 24 v 24 blocks =2

%what type of run is this?

%general parameters
decode_int = 5;
decode_time = [-200 499];  
niter = 100;
bin_list = [{1:2} {2:3} {4:6} {1:9}];
bin_min = [1,2,4,1];
bin_max = [2,3,6,9];

min_trials = [20,20,10,10];
%trial_thresholds = {[20 30 40 50 60 70] [20 30 40 50 60 70] [6 9 12 15 18 21 24] [6 9 12 15 18 21 24]};
%cross_blocks     = {[ 2  3  4  5  6  7] [ 2  3  4  5  6  7] [2 3  4  5  6  7  8] [2 3  4  5  6  7  8]};
trial_thresholds = {[20 30 40 50 60 70] [20 30 40 50 60 70] [10 15 20] [10 15 20]};
cross_blocks     = {[ 2  3  4  5  6  7] [ 2  3  4  5  6  7] [2   3  4] [ 2  3  4]};

%% create folder names
npermin1 = trial_thresholds{1}(1)/cross_blocks{1}(1);
npermin2 = trial_thresholds{3}(1)/cross_blocks{3}(1);

run_label = sprintf('/mvpc_dint-%d_min-%d.%d.%d.%d_nper-%d.%d_cbmin-%d_iter-%d', ...
    decode_int, min_trials(1), min_trials(2), min_trials(3), min_trials(4), ...
    npermin1, npermin2, cross_blocks{1}(1), niter);

save_dir = fullfile(pwd, tps{1}, run_label);

% Create the directory (MATLAB automatically creates intermediate folders)
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

cat_folders  = {'faces vs houses', 'upright vs inverted'};
img_folders  = {'face_identity', 'all_bins'};

% Base MVPC folder types
base_folders = {'mvpc_files', 'mvpc_output'};

% === Create the erp_trials folder at the same level ===
erp_path = fullfile(save_dir, 'trials_per_erp');
if ~exist(erp_path, 'dir'), mkdir(erp_path); end

% === Create the MVPC folders and subfolders ===
for b = 1:length(base_folders)

    % Create category_bins subfolders
    for c = 1:length(cat_folders)
        folder_to_make = fullfile(save_dir, base_folders{b}, 'category_bins', cat_folders{c});
        if ~exist(folder_to_make, 'dir'), mkdir(folder_to_make); end
    end

    % Create image_bins subfolders
    for i = 1:length(img_folders)
        folder_to_make = fullfile(save_dir, base_folders{b}, 'image_bins', img_folders{i});
        if ~exist(folder_to_make, 'dir'), mkdir(folder_to_make); end
    end

end




subs = [{'faces vs houses'} {'upright vs inverted'} {'face_identity'} {'all_bins'}];
load_paths = [{'/best_files/category_bins/'} {'/best_files/category_bins/'} {'/best_files/image_bins/'}  {'/best_files/image_bins/'}];
save_paths = [{[run_label '/mvpc_files/category_bins/faces vs houses/']} {[run_label '/mvpc_files/category_bins/upright vs inverted/']} {[run_label '/mvpc_files/image_bins/face_identity/']} {[run_label '/mvpc_files/image_bins/all_bins/']}];
output_paths = [{[run_label '/mvpc_output/category_bins/faces vs houses/']} {[run_label '/mvpc_output/category_bins/upright vs inverted/']} {[run_label '/mvpc_output/image_bins/face_identity/']} {[run_label '/mvpc_output/image_bins/all_bins/']}];

%% start decoding

for timepoints = 1:length(tps)
    
    for blocks = 1:length(subs)
        
        %run decoding
        files = dir([tps{timepoints} load_paths{blocks} '*.best']);
        %files = files(1);
        runMVPA(files, tps{timepoints}, load_paths{blocks},save_paths{blocks},output_paths{blocks}, bin_list{blocks},bin_min(blocks), bin_max(blocks), min_trials(blocks), trial_thresholds{blocks}, cross_blocks{blocks}, decode_int, decode_time,niter)
        %read in MVPC and save stuff
        files = dir([tps{timepoints} save_paths{blocks} '*.mvpc']);
        %files = files(1);
        exportMVPC(files,tps{timepoints}, save_dir, save_paths{blocks},output_paths{blocks}, subs{blocks})
    
    end

end