function runDecoding(files, timepoint, load_path, save_path, output_path, bins, bmin, bmax, min_trial, thresholds, n_blocks, decode_int, decode_time, niter)

%timepoint='T1';
%load_path=load_paths{3};
%min_trial=5;
%bmin=bin_min(3);
%bmax=bin_max(3);
%thresholds=trial_thresholds{3};
%n_blocks = cross_blocks{3}

for i = 1:length(files)
    eegFile = [files(i).folder '/' files(i).name];
    [eegPath, eegHead, eegExt] = fileparts(eegFile);
        
    % Load Data and check number of trials in bins
    BEST = pop_loadbest( 'filename', [eegHead '.best'], 'filepath', [timepoint load_path] );
        
    not_enough = any(BEST.n_trials_per_bin(bmin:bmax) < min_trial);

    if not_enough
        continue;
    end
    
    % evaluate the lowest number of trials and set the highest cross block
    % whilst maintaining equal trials amongst participants
    low = min(BEST.n_trials_per_bin(bmin:bmax));
    idx = find(low >= thresholds, 1, 'last');
    n_cross_blocks = n_blocks(idx);

    %Run Decoding
    MVPC = pop_decoding( BEST , 'Channels', [ 1:120], 'classcoding', 'OneVsAll', 'Classes', bins, 'Decode_Every_Npoint', decode_int, 'DecodeTimes', decode_time, 'EqualizeTrials', 'classes', 'Method', 'SVM', 'nCrossblocks', n_cross_blocks, 'nIter', niter, 'ParCompute', 'on' );
    MVPC = pop_savemymvpc(MVPC, 'mvpcname', eegHead, 'filename', [eegHead '.mvpc'], 'filepath', [timepoint save_path], 'Warning', 'off');
    MVPC.history = [];
    pop_mvpc2text( MVPC, [timepoint output_path eegHead '_mvpc.txt'],  'time', 'on', 'timeunit',  0.001, 'transpose', 'off');
end

end

