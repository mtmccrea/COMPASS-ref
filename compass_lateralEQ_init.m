function lateralEQ_struct = compass_lateralEQ_init(synthesis_struct, initEQ_struct, debug, verbose)


    % Preprocess the EQ table and build the lateralEQ_struct for synthesis

    % Unpack struct with info about the precomputed EQs
    fwBiasSteps     = initEQ_struct.fwBiasSteps;
    binFreq_latEQ   = initEQ_struct.binFreq_latEQ;
    listeningRadii  = initEQ_struct.listeningRadii;
    lsSpreads       = initEQ_struct.lsSpreads;
    panDirs         = initEQ_struct.panDirs;
    Fc_erb          = initEQ_struct.Fc_erb;
    earWeight       = initEQ_struct.earWeight;
    listeningRadius = initEQ_struct.listeningRadius;
    interpType      = initEQ_struct.interpType;
    correctionTarget= initEQ_struct.correctionTarget;
    eq_table        = initEQ_struct.eq_table;
    k_freq          = initEQ_struct.k_freq;
    k_dir           = initEQ_struct.k_dir;
    k_clip          = initEQ_struct.k_clip;
    lsDesign        = initEQ_struct.lsDesign;
    NpanStepsOut    = numel(fwBiasSteps);
    if initEQ_struct.dim == 2
        lsSpread2D      = initEQ_struct.lsSpread2D; % for 2D regular arrays
        onPoint2D       = initEQ_struct.onPoint2D;
    end
    
    % Only store the sampled bins that are unique, i.e. no repeated low/high 
    % freqs outside ERB band range (they can be queried on demand)
    lateralEQ_struct.binFreq_latEQ = binFreq_latEQ;
    binFreq_latEQ_precomp = unique(binFreq_latEQ);
    lateralEQ_struct.binFreq_latEQ_precomp = binFreq_latEQ_precomp;
    lateralEQ_struct.Nfrq = numel(lateralEQ_struct.binFreq_latEQ);

    lateralEQ_struct.fwBiasSteps = fwBiasSteps;
    lateralEQ_struct.earWeight = earWeight;
    lateralEQ_struct.listeningRadius = listeningRadius;
    lateralEQ_struct.interpType = interpType;
    lateralEQ_struct.correctionTarget = correctionTarget;
    lateralEQ_struct.Fc_erb = Fc_erb;
    lateralEQ_struct.Nerb = numel(Fc_erb);
    lateralEQ_struct.maxSpread = max(lsSpreads);
    lateralEQ_struct.minSpread = min(lsSpreads);
    % lateralEQ_struct.lsSpread = lsSpread; % for 2D

    if nargin < 4, debug = false; end
    if nargin < 5, verbose = false; end
    lateralEQ_struct.debug = debug;
    lateralEQ_struct.verbose = verbose;

    disp('Precomputing EQ for DoA directions...')
    % Preprocess
    % Size: eq_table (NlisteningRadii, NlsSpread, NpanStepsOut, Ndir, Nerb, 2)
    % values correspond to those of the precomputed eq_table
    [d1,d2,d3,d4,d5,d6] = ndgrid( ...
        listeningRadii, lsSpreads, fwBiasSteps, panDirs, Fc_erb, [0 1]);

    switch initEQ_struct.dim
        case 2
            % Table reduction one: fix listening radius, lsSpread, ear weight
            [dq1,dq2,dq3,dq4,dq5,dq6] = ndgrid(   ...
                listeningRadius, ...   % listeningRadii, matching original precalc'd vals or interpolated
                lsSpread2D,      ...   % LS spread, matching original precalc'd vals or interpolated
                fwBiasSteps,     ...   % forward bias: pan "position" within LS pair, formulated as [0,1] from left/backmost speaker
                panDirs,         ...   % positive lateral pan direction ([0,90] deg)
                Fc_erb,          ...   % all ERB band freqs
                earWeight        ...   % evenly weighted between lateral and contralateral ears
                );
            % Interp to >> Size: eq_table_culled (NpanStepsOut, NpanDirs, Fc_erb)
            eq_table_culled = interpn( ...
                d1, d2, d3, d4, d5, d6, eq_table, dq1, dq2, dq3, dq4, dq5, dq6, interpType);
            eq_table_culled = squeeze(eq_table_culled);
            checkForNaNs(eq_table_culled, 'EQ table reduction stage 1');

            % Pre-processing: smoothing (over freq, direction), soft clipping
            for pani = 1:NpanStepsOut
                % eq_table, k_freq, k_dir, k_clip, Ndirs, Nfrq
                eq_table_culled(pani, :, :) = smoothEQTable(    ...
                    squeeze(eq_table_culled(pani,:,:))',        ...
                    k_freq, k_dir, k_clip,                      ...
                    numel(panDirs), numel(Fc_erb)               ...
                    )';
            end
            checkForNaNs(eq_table_culled, 'EQ table reduction stage 1: post-smoothing');
            clear d1 d2 d3 d4 d5 d6 dq1 dq2 dq3 dq4 dq5 dq6 pani 

            % Table reduction two: fix the spectrum to all bins within ERB bands
            [d1, d2, d3]    = ndgrid(fwBiasSteps, panDirs, Fc_erb);
            [dq1, dq2, dq3] = ndgrid(fwBiasSteps, panDirs, binFreq_latEQ_precomp);
            % To size: eq_table_culled (NpanStepsOut, NpanDirs, Fc_erb)
            eq_table_culled = interpn(  ...
                d1, d2, d3, eq_table_culled, dq1, dq2, dq3, interpType);
            eq_table_culled = squeeze(eq_table_culled);
            checkForNaNs(eq_table_culled, 'EQ table reduction stage 2');

            lateralEQ_struct.eq_table_culled = eq_table_culled;
            lateralEQ_struct.interpGrid{1} = dq1;
            lateralEQ_struct.interpGrid{2} = dq2;
            lateralEQ_struct.interpGrid{3} = dq3;
            lateralEQ_struct.lsArrayAz2D = getRing(360/lsSpread2D, onPoint2D);

        case 3
            % Table reduction one: fix listening radius, lsSpread, ear weight
            [dq1,dq2,dq3,dq4,dq5,dq6] = ndgrid(   ...
                listeningRadius, ...   % listeningRadii, matching original precalc'd vals or interpolated
                lsSpreads,       ...   % LS spreads, matching original precalc'd vals or interpolated
                fwBiasSteps,     ...   % forward bias: pan "position" within LS pair, formulated as [0,1] from left/backmost speaker
                panDirs,         ...   % positive lateral pan direction ([0,90] deg)
                Fc_erb,          ...   % all ERB band freqs
                earWeight        ...   % evenly weighted between lateral and contralateral ears
                );
            % Interp to >> Size: eq_table_culled (NlsSpreads, NpanStepsOut, NpanDirs, Fc_erb)
            eq_table_culled = interpn( ...
                d1, d2, d3, d4, d5, d6, eq_table, dq1, dq2, dq3, dq4, dq5, dq6, interpType);
            eq_table_culled = squeeze(eq_table_culled);
            checkForNaNs(eq_table_culled, 'reduction stage 1');
            NlsSpreads = numel(lsSpreads);
            % Pre-processing: smoothing (over freq, direction), soft clipping
            for sprdi = 1:NlsSpreads
                for pani = 1:NpanStepsOut
                    % eq_table, k_freq, k_dir, k_clip, Ndirs, Nfrq
                    eq_table_culled(sprdi,pani, :, :) = smoothEQTable(  ...
                        squeeze(eq_table_culled(sprdi,pani,:,:))',            ...
                        k_freq, k_dir, k_clip,                          ...
                        numel(panDirs), numel(Fc_erb)                   ...
                        )';
                end
            end
            checkForNaNs(eq_table_culled, 'post-smoothing');
            clear d1 d2 d3 d4 d5 d6 dq1 dq2 dq3 dq4 dq5 dq6 pani 

            % Table reduction two: fix the spectrum to all bins within ERB bands
            [d1, d2, d3, d4]     = ndgrid(lsSpreads, fwBiasSteps, panDirs, Fc_erb);
            [dq1, dq2, dq3, dq4] = ndgrid(lsSpreads, fwBiasSteps, panDirs, binFreq_latEQ_precomp);
            % To size: eq_table_culled (NpanStepsOut, NpanDirs, Fc_erb)
            eq_table_culled = interpn(  ...
                d1, d2, d3, d4, eq_table_culled, dq1, dq2, dq3, dq4, interpType);
            eq_table_culled = squeeze(eq_table_culled);
            checkForNaNs(eq_table_culled, 'reduction stage 2');

            lateralEQ_struct.eq_table_culled = eq_table_culled;
            lateralEQ_struct.interpGrid{1} = dq1;
            lateralEQ_struct.interpGrid{2} = dq2;
            lateralEQ_struct.interpGrid{3} = dq3;
            lateralEQ_struct.interpGrid{4} = dq4;

            lateralEQ_struct.lsLatAngs_rad = azel2interauralRad(lsDesign.Directions(:,1), lsDesign.Directions(:,2));
    end
    lateralEQ_struct.k_freq = k_freq;
    lateralEQ_struct.k_dir = k_dir;
    lateralEQ_struct.k_clip = k_clip;

    clear d1 d2 d3 d4 d5 d6 dq1 dq2 dq3 dq4 dq5 dq6 pani 

    
    % Precompute the EQs for each direction of the DoA grid (takes a bit of time)

    Ngrid = size(synthesis_struct.DOAgrid,1);
    doas_aziElevRad = synthesis_struct.DOAgrid_aziElevRad;  % (NdoaPnts, 2)
    gains_nerb = synthesis_struct.vbap_gtable;              % (Nls, NdoaPnts)
    gains_nerb = gains_nerb ./ sum(gains_nerb, 1);          % return base gains to amplitude-normalized

    srcLatAng_rad = zeros(Ngrid,1);
    fwBias        = ones(Ngrid,2) * -1;
    lsSpreads_rad = ones(Ngrid,2) * -1;
    weights       = ones(Ngrid,2) * -1;

    for gi = 1:Ngrid
        [srcLatAng_rad(gi), fwb, sprd_ls, wgt] = ...
            getForwardBias3DArray(doas_aziElevRad(gi,:), lateralEQ_struct.lsLatAngs_rad, gains_nerb(:,gi), verbose);
        if any(fwb > 1)
            if verbose
            warning('fwBias is above 1:\n\tfwBias: %s\n\tsrcLatAng_deg: %.1f\n\tDoA:%s\n\tlsSpreads_deg: %s\n\tgains_nerb(:,doai): %s\n', ...
                join(string(fwb), ' '), ...
                rad2deg(srcLatAng_rad(gi)), ...
                join(string(rad2deg(doas_aziElevRad(gi,:))), ' '), ...
                join(string(rad2deg(sprd_ls)), ' '), ...
                join(string(gains_nerb(:,gi)), ' '));
            end
            fwb = min(1,fwb);
        end
        pairIdx = 1:numel(fwb);
        fwBias(gi,pairIdx) = fwb;
        lsSpreads_rad(gi,pairIdx) = sprd_ls;
        weights(gi,pairIdx) = wgt;
    end

    % Calculate EQs for every direction, given the user-defined smoothing params
    srcLatAng_deg = abs(rad2deg(srcLatAng_rad));
    lsSpreads_deg = rad2deg(lsSpreads_rad);
    lsSpreads_deg = min(lsSpreads_deg, lateralEQ_struct.maxSpread); % clip spread to precomputed max (e.g. 90)

    clear gains_nerb verbose pairIdx gi lsSpreads_rad doas_aziElevRad

    latEQs = zeros(numel(lateralEQ_struct.binFreq_latEQ),Ngrid);
    igrid = lateralEQ_struct.interpGrid;
    for gi = 1:Ngrid
        ws = weights(gi,:);
        pairmask = ws > 0;
        Npair = sum(pairmask);

        [dq1, dq2, dq3, dq4] = ndgrid( ...
            lsSpreads_deg(gi,pairmask), fwBias(gi,pairmask), srcLatAng_deg(gi), lateralEQ_struct.binFreq_latEQ);

        % Size: (Npair, Npair, NbinFrq)
        pairEQ = squeeze( interpn(                          ...
            igrid{1}, igrid{2}, igrid{3}, igrid{4},         ...
            lateralEQ_struct.eq_table_culled,               ...
            dq1, dq2, dq3, dq4, lateralEQ_struct.interpType ...
            ));

        for pairi = 1:Npair
            latEQs(:,gi) = latEQs(:,gi) + ...
                squeeze(pairEQ(pairi,pairi,:)) .* ws(pairi);
        end
    end
    checkForNaNs(latEQs);
    lateralEQ_struct.latEQ_doa = latEQs;
    if debug
        lateralEQ_struct.fwBias_doa = fwBias;
        lateralEQ_struct.lsSpreads_deg_doa = lsSpreads_deg;
        lateralEQ_struct.weights_doa = weights;
    end
    clear dq1 dq2 dq3 dq4 pairEQ pairi igrid 
    disp('EQ precomputed.')
end