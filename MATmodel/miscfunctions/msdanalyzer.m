classdef msdanalyzer
    %%MSDANALYZER a class for simple mean square displacement analysis.
    %
    % msdanalyzer is a MATLAB per-value class dedicated to the mean square
    % displacement analysis of single particle trajectories. After the
    % instantiation, the user can feed the object with the trajectories of
    % several particles. The class allows him to then derive the MSD
    % curves, the velocity autocorrelation functions, to plot them and fit
    % them, and possible correct for drift if required.
    %
    % Jean-Yves Tinevez - Institut Pasteur, 2013 <tinevez at pasteur dot fr>
    
    properties (Constant)
        % Tolerance for binning delays together. Two delays will be binned
        % together if they differ in absolute value by less than
        % 10^-TOLERANCE.
        TOLERANCE = 12;
    end
    
    properties (SetAccess = private)
        % The trajectories stored in a cell array, one T x n_dim per particle
        tracks = {};
        % The dimensionality of the problem
        n_dim
        % Spatial units
        space_units;
        % Time units
        time_units
        % Stores the MSD
        msd
        % Stores the velocty correlation
        vcorr
        % Stores the linear fit of MSD vs t
        lfit
        % Stores the linear fit of the log log plot of MSD vs t
        loglogfit
        % Drift movement
        drift
    end
    
    properties (SetAccess = private, GetAccess = private, Hidden = true)
        % If false, msd needs to be recomputed
        msd_valid = false;
        % If false, vcorr needs to be recomputed
        vcorr_valid = false
        
    end
    
    
    %% Constructor
    methods
        
        function obj = msdanalyzer(n_dim, space_units, time_units)
            %%MSDANALYZER Builds a new MSD analyzer object.
            % obj = MSDANALYZER(dimensionality, space__units, time_units)
            % builds a new empty msdanalyzer object with the specified
            % dimensionality (2 for 2D, 3 for 3D, etc...), spatial units
            % and time units. Units are strings, used only for display and
            % plotting.
            
            if ~isreal(n_dim) || ~isscalar(n_dim) || n_dim < 1 || ~(n_dim == floor(double(n_dim)))
                error('msdanalyzer:BadDimensionality', ...
                    'Dimensionality must be a positive integer, got %f.', ...
                    n_dim);
            end
            obj.n_dim = n_dim;
            obj.space_units = space_units;
            obj.time_units = time_units;
        end
    end
    
    %% Public methods
    methods
        
        function obj = addAll(obj, tracks)
            %%ADDALL Add specified trajectories to msd analyzer.
            %
            % obj = obj.addAll(tracks) adds the given tracks the
            % msdanalyzer object obj.
            %
            % Tracks must be specified as a cell array, one array per
            % track. Each track array must be of size N x (Ndim+1) where N
            % is the number of individual measurements of the particle
            % position and Ndim is the dimensionality specified during
            % object construction. Track array must be arranged with time
            % and space as follow: [ Ti Xi Yi ... ] etc..
            %
            % Adding new tracks to an existing object invalidates its
            % stored MSD values, and will cause them to be recalculated if
            % needed.
            
            tracks = tracks(:);
            n_tracks = numel(tracks);
            
            % Check dimensionality
            for i = 1 : n_tracks
                track = tracks{i};
                track_dim = size(track, 2);
                if track_dim ~= obj.n_dim + 1
                    error('msdanalyzer:addAll:BadDimensionality', ...
                        'Tracks must be of size N x (nDim+1) with [ T0 X0 Y0 ... ] etc...');
                end
            end
            
            % Add to track collection
            obj.tracks = [
                obj.tracks;
                tracks
                ];
            obj.msd_valid = false;
            obj.vcorr_valid = false;
            
        end
        
        
        function varargout = plotTracks(obj, ha, indices, corrected)
            %%PLOTTRACKS Plot the tracks stored in this object.
            %
            % obj.plotTracks plots the particle trajectories stored in the
            % msdanalyzer object obj in the current axes. This method only
            % works for 2D or 3D problems.
            %
            % obj.plotTracks(ha) plots the trajectories in the axes with
            % handle ha.
            %
            % obj.plotTracks(ha, indices) where indices is a vector, allows
            % to specify the track to be plotted, using their indices.
            % Leave the vector empty to plot all trajectories.
            %
            % obj.plotTracks(ha, indices, corrected) where corrected is a a
            % boolean flag, allows to specify whether the plot should
            % display the trajectories corrected for drift (true) or
            % uncorrected (false). A proper drift vector must be computed
            % prior to setting this flag to true. See
            % msdanalyzer.computeDrift.
            %
            % hps = obj.plotTracks(...) returns the handles to the
            % individual line objects created.
            %
            % [hps, ha] = obj.plotTracks(...) returns also the handle to
            % the axes handle the trajectories are plot in.
            
            if obj.n_dim < 2 || obj.n_dim > 3
                error('msdanalyzer:plotTracks:UnsupportedDimensionality', ...
                    'Can only plot tracks for 2D or 3D problems, got %dD.', obj.n_dim);
            end
            
            if nargin < 2
                ha = gca;
            end
            if nargin < 3 || isempty(indices)
                indices = 1 : numel(obj.tracks);
            end
            if nargin < 4
                corrected = false;
            end
            
            n_tracks = numel(indices);
            colors = jet(n_tracks);
            
            hold(ha, 'on');
            hps = NaN(n_tracks, 1);
            
            if obj.n_dim == 2
                % 2D case
                for i = 1 : n_tracks
                    
                    index = indices(i);
                    track = obj.tracks{index};
                    
                    x = track(:,2);
                    y = track(:,3);
                    
                    if corrected && ~isempty(obj.drift)
                        tdrift = obj.drift(:,1);
                        xdrift = obj.drift(:, 2);
                        ydrift = obj.drift(:, 3);
                        t = track(:,1);
                        [~, index_in_drift_time, ~] = intersect(tdrift, t);
                        % Subtract drift position to track position
                        x = x - xdrift(index_in_drift_time);
                        y = y - ydrift(index_in_drift_time);
                    end
                    
                    hps(i) =  plot(ha, x, y, ...
                        'Color', colors(i,:));
                    
                end
                
            else
                % 3D case
                for i = 1 : n_tracks
                    
                    index = indices(i);
                    track = obj.tracks{index};
                    
                    x = track(:,2);
                    y = track(:,3);
                    z = track(:,4);
                    
                    if corrected && ~isempty(obj.drift)
                        tdrift = obj.drift(:,1);
                        xdrift = obj.drift(:, 2);
                        ydrift = obj.drift(:, 3);
                        zdrift = obj.drift(:, 4);
                        t = track(:,1);
                        [~, index_in_drift_time, ~] = intersect(tdrift, t);
                        % Subtract drift position to track position
                        x = x - xdrift(index_in_drift_time);
                        y = y - ydrift(index_in_drift_time);
                        z = z - zdrift(index_in_drift_time);
                    end
                    
                    hps(i) =  plot3(ha, x, y, z, ...
                        'Color', colors(i,:));
                    
                end
                
            end
            
            % Output
            if nargout > 0
                varargout{1} = hps;
                if nargout > 1
                    varargout{2} = ha;
                end
            end
            
        end
        
        
        function varargout = plotDrift(obj, ha)
            %%PLOTDRIFT Plot drift stored in this object.
            %
            % obj.plotDrift plots the calculated drift position in the
            % current axes.
            %
            % obj.plotDrift(ha) plots the drift position in the axes
            % specified by the handle ha.
            %
            % hp = obj.plotDrift(...) returns the plot handle for the
            % created line.
            %
            % [hp, ha] = obj.plotDrift(...) also returns the axes handle.
            
            if obj.n_dim < 2 || obj.n_dim > 3
                error('msdanalyzer:plotDrift:UnsupportedDimensionality', ...
                    'Can only plot drift for 2D or 3D problems, got %dD.', obj.n_dim);
            end
            
            if isempty(obj.drift)
                if nargout > 0
                    varargout{1} = [];
                    if nargout > 1
                        varargout{2} = [];
                    end
                end
                return
            end
            
            if nargin < 2
                ha = gca;
            end
            
            if obj.n_dim == 2
                hp = plot(ha, obj.drift(:,2), obj.drift(:,3), 'k-', ...
                    'LineWidth', 2);
            else
                hp = plot3(ha, obj.drift(:,2), obj.drift(:,3), obj.drift(:,4), 'k-', ...
                    'LineWidth', 2);
            end
            
            if nargout > 0
                varargout{1} = hp;
                if nargout > 1
                    varargout{2} = ha;
                end
            end
            
        end
        
        
        function varargout = labelPlotTracks(obj, ha)
            %%LABELPLOTTRACKS A convenience method to set the axes labels.
            %
            % obj.labelPlotTracks(ha) sets the axis label of the axes with
            % the specified handle ha. It is meant for axes containing the
            % plot of the particles trajectories and their drift.
            %
            % hl = obj.plotTracks(...) returns the handle to the generated
            % labels.
            
            if obj.n_dim < 2 || obj.n_dim > 3
                error('msdanalyzer:labelPlotTracks:UnsupportedDimensionality', ...
                    'Can only label axis for 2D or 3D problems, got %dD.', obj.n_dim);
            end
            
            if nargin < 2
                ha = gca;
            end
            
            hl = NaN(obj.n_dim, 1);
            hl(1) = xlabel(ha, ['X (' obj.space_units ')'] );
            hl(2) = ylabel(ha, ['Y (' obj.space_units ')'] );
            if obj.n_dim ==3
                hl(3) = zlabel(ha, ['Z (' obj.space_units ')'] );
            end
            axis equal
            if nargout > 0
                varargout{1} = hl;
            end
        end
        
        
        function varargout = labelPlotMSD(obj, ha)
            %%LABELPLOTMSD A convenience method to set the axes labels.
            %
            % obj.labelPlotMSD(ha) sets the axis label of the axes with
            % the specified handle ha. It is meant for axes containing the
            % plot of the mean-square-displacement.
            %
            % hl = obj.plotMSD(...) returns the handle to the generated
            % labels.
            
            if nargin < 2
                ha = gca;
            end
            
            hl = NaN(2, 1);
            hl(1) = xlabel(ha, ['Delay (' obj.time_units ')']);
            hl(2) = ylabel(ha, ['MSD (' obj.space_units '^2)']);
            
            xl = xlim(ha);
            xlim(ha, [0 xl(2)]);
            yl = ylim(ha);
            ylim(ha, [0 yl(2)]);
            box(ha, 'off')
            
            if nargout > 0
                varargout{1} = hl;
            end
        end
        
        
        function varargout = labelPlotVCorr(obj, ha)
            %%LABELPLOTVCORR A convenience method to set the axes labels.
            %
            % obj.labelPlotVCorr(ha) sets the axis label of the axes with
            % the specified handle ha. It is meant for axes containing the
            % plot of the velocity autocorrelation.
            %
            % hl = obj.labelPlotVCorr(...) returns the handle to the generated
            % labels.
            %
            % [hl, hline] = obj.labelPlotVCorr(...) also returns the handle
            % to generate y=0 dashed line.
            
            if nargin < 2
                ha = gca;
            end
            
            hl = NaN(2, 1);
            hl(1) = xlabel(ha, ['Delay (' obj.time_units ')']);
            hl(2) = ylabel(ha, 'Normalized velocity autocorrelation');
            
            
            xl = xlim(ha);
            xlim(ha, [0 xl(2)]);
            box(ha, 'off')
            
            hline = line([0 xl(2)], [0 0], ...
                'Color', 'k', ...
                'LineStyle', '--');
            uistack(hline, 'bottom')
            
            if nargout > 0
                varargout{1} = hl;
                if nargout > 1
                    varargout{2} = hline;
                end
            end
        end
        
        
        function obj = computeMSD(obj, indices)
            %%COMPUTEMSD Compute the mean-squared-displacement for this object.
            %
            % obj = obj.computeMSD computes the MSD for all the tracks stored
            % in this object. If a drift correction was computed prior to this
            % method call, it is used to correct positions before MSD
            % calculation.
            %
            % Results are stored in the msd field of this object as a cell
            % array, one cell per particle. The array is a double array of size
            % N x 4, and is arranged as follow: [dt mean std N ; ...] where dt
            % is the delay for the MSD, mean is the mean MSD value for this
            % delay, std the standard deviation and N the number of points in
            % the average.
            %
            % obj = obj.computeMSD(indices) computes the MSD only for the
            % particles with the specified indices. Use an empty array to take
            % all particles.
            
            if nargin < 2 || isempty(indices)
                indices = 1 : numel(obj.tracks);
            end
            
            n_tracks = numel(indices);
            fprintf('Computing MSD of %d tracks... ', n_tracks);
            
            % First, find all possible delays in time vectors.
            % Time can be arbitrary spaced, with frames missings,
            % non-uniform sampling, etc... so we have to do this clean.
            % We use a certain tolerance to bin delays together
            delays = obj.getAllDelays;
            n_delays = numel(delays);
            
            obj.msd = cell(n_tracks, 1);
            if ~isempty(obj.drift)
                tdrift = obj.drift(:,1);
                xdrift = obj.drift(:, 2:end);
            end
            
            fprintf('%4d/%4d', 0, n_tracks);
            
            for i = 1 : n_tracks
                
                fprintf('\b\b\b\b\b\b\b\b\b%4d/%4d', i, n_tracks);
                
                mean_msd    = zeros(n_delays, 1);
                M2_msd2     = zeros(n_delays, 1);
                n_msd       = zeros(n_delays, 1);
                
                index = indices(i);
                track = obj.tracks{index};
                t = track(:,1);
                t = msdanalyzer.roundn(t, msdanalyzer.TOLERANCE);
                X = track(:, 2:end);
                
                % Determine drift correction
                if ~isempty(obj.drift)
                    % Determine target delay index in bulk
                    [~, index_in_drift_time, index_in_track_time] = intersect(tdrift, t);
                    % Keep only track times that can be corrected.
                    X = X(index_in_track_time, :);
                    t = t(index_in_track_time);
                    % Subtract drift position to track position
                    X = X - xdrift(index_in_drift_time, :);
                    
                end
                
                
                n_detections = size(X, 1);
                
                for j = 1 : n_detections - 1
                    
                    % Delay in physical units
                    dt = t(j+1:end) - t(j);
                    dt = msdanalyzer.roundn(dt, msdanalyzer.TOLERANCE);
                    
                    % Determine target delay index in bulk
                    [~, index_in_all_delays, ~] = intersect(delays, dt);
                    
                    % Square displacement in bulk
                    dX = X(j+1:end,:) - repmat(X(j,:), [(n_detections-j) 1] );
                    dr2 = sum( dX .* dX, 2);
                    
                    % Store for mean computation / Knuth
                    n_msd(index_in_all_delays)     = n_msd(index_in_all_delays) + 1;
                    delta = dr2 - mean_msd(index_in_all_delays);
                    mean_msd(index_in_all_delays) = mean_msd(index_in_all_delays) + delta ./ n_msd(index_in_all_delays);
                    M2_msd2(index_in_all_delays)  = M2_msd2(index_in_all_delays) + delta .* (dr2 - mean_msd(index_in_all_delays));
                end
                
                n_msd(1) = n_detections;
                std_msd = sqrt( M2_msd2 ./ n_msd ) ;
                
                % We replace points for which N=0 by Nan, to later treat
                % then as missing data. Indeed, for each msd cell, all the
                % delays are present. But some tracks might not have all
                % delays
                delay_not_present = n_msd == 0;
                mean_msd( delay_not_present ) = NaN;
                
                obj.msd{index} = [ delays mean_msd std_msd n_msd ];
                
            end
            fprintf('\b\b\b\b\b\b\b\b\bDone.\n')
            
            obj.msd_valid = true;
            
        end
        
        
        function varargout = plotMSD(obj, ha, indices, errorbar)
            %% PLOTMSD Plot the mean square displacement curves.
            %
            % obj.plotMSD plots the MSD curves in the current axes.
            %
            % obj.plotMSD(ha) plots the MSD curves in the axes
            % specified by the handle ha.
            %
            % obj.plotMSD(ha, indices) plots the MSD curves for the
            % particles with the specified indices only. Leave empty to
            % plot MSD for all particles.
            %
            % obj.plotMSD(ha, indices, errorbar), where errorbar is a
            % boolean flag, allows to specify whether the curves should be
            % plotted with error bars (equal to standard deviation). It is
            % false by default.
            %
            % hps =  obj.plotMSD(...) returns the handle array for the
            % lines generated.
            %
            % [hps, ha] =  obj.plotMSD(...) also return the axes handle in
            % which the lines were plotted.
            
            if ~obj.msd_valid
                obj = obj.computeMSD;
            end
            
            if nargin < 2
                ha = gca;
            end
            if nargin < 3 || isempty(indices)
                indices = 1 : numel(obj.msd);
            end
            if nargin < 4
                errorbar = false;
            end
            
            n_spots = numel(indices);
            colors = jet(n_spots);
            
            hold(ha, 'on');
            if errorbar
            else
                hps = NaN(n_spots, 1);
            end
            
            for i = 1 : n_spots
                
                index = indices(i);
                
                msd_spot = obj.msd{index};
                t = msd_spot(:,1);
                m = msd_spot(:,2);
                if errorbar
                    s = msd_spot(:,3);
                    hps(i) = msdanalyzer.errorShade(ha, t, m, s, colors(i,:), true);
                else
                    hps(i) = plot(ha, t, m, 'Color', colors(i,:));
                end
                
            end
            
            obj.labelPlotMSD(ha);
            
            if nargout > 0
                varargout{1} = hps;
                if nargout > 1
                    varargout{2} = ha;
                end
            end
            
            
        end
        
        
        function obj = fitLogLogMSD(obj, clip_factor)
            %%FITLOGLOGMSD Fit the log-log MSD to determine behavior.
            %
            % obj = obj.fitLogLogMSD fits each MSD curve stored in this object
            % in a log-log fashion. If x = log(delays) and y = log(msd) where
            % 'delays' are the delays at which the msd is calculated, then this
            % method fits y = f(x) by a straight line y = alpha * x + gamma, so
            % that we approximate the MSD curves by MSD = gamma * delay^alpha.
            % By default, only the first 25% of each MSD curve is considered
            % for the fit,
            %
            % Results are stored in the 'loglogfit' field of the returned
            % object. It is a structure with 3 fields:
            % - alpha: all the values for the slope of the log-log fit.
            % - gamma: all the values for the value at origin of the log-log fit.
            % - r2fit: the adjusted R2 value as a indicator of the goodness of
            % the fit.
            %
            % obj = obj.fitLogLogMSD(clip_factor) does the fit, taking into
            % account only the first potion of each MSD curve specified by
            % 'clip_factor' (a double between 0 and 1). If the value
            % exceeds 1, then the clip factor is understood to be the
            % maximal number of point to take into account in the fit. By
            % default, it is set to 0.25.
            
            if nargin < 2
                clip_factor = 0.25;
            end
            
            if ~obj.msd_valid
                obj = obj.computeMSD;
            end
            n_spots = numel(obj.msd);
            
            if clip_factor < 1
                fprintf('Fitting %d curves of log(MSD) = f(log(t)), taking only the first %d%% of each curve... ',...
                    n_spots, ceil(100 * clip_factor) )
            else
                fprintf('Fitting %d curves of log(MSD) = f(log(t)), taking only the first %d points of each curve... ',...
                    n_spots, round(clip_factor) )
            end
            
            alpha = NaN(n_spots, 1);
            gamma = NaN(n_spots, 1);
            r2fit = NaN(n_spots, 1);
            ft = fittype('poly1');
            
            fprintf('%4d/%4d', 0, n_spots);
            for i_spot = 1 : n_spots
                
                fprintf('\b\b\b\b\b\b\b\b\b%4d/%4d', i_spot, n_spots);
                
                msd_spot = obj.msd{i_spot};
                
                t = msd_spot(:,1);
                y = msd_spot(:,2);
                w = msd_spot(:,4);
                
                % Clip data
                if clip_factor < 1
                    t_limit = 2 : round(numel(t) * clip_factor);
                else
                    t_limit = 2 : min(1+round(clip_factor), numel(t));
                end
                t = t(t_limit);
                y = y(t_limit);
                w = w(t_limit);
                
                % Thrash bad data
                nonnan = ~isnan(y);
                
                t = t(nonnan);
                y = y(nonnan);
                w = w(nonnan);
                
                if numel(y) < 2
                    continue
                end
                
                xl = log(t);
                yl = log(y);
                
                bad_log =  isinf(xl) | isinf(yl);
                xl(bad_log) = [];
                yl(bad_log) = [];
                w(bad_log) = [];
                
                if numel(xl) < 2
                    continue
                end
                
                [fo, gof] = fit(xl, yl, ft, 'Weights', w);
                
                alpha(i_spot) = fo.p1;
                gamma(i_spot) = exp(fo.p2);
                r2fit(i_spot) = gof.adjrsquare;
                
            end
            fprintf('\b\b\b\b\b\b\b\b\bDone.\n')
            
            obj.loglogfit = struct(...
                'alpha', alpha, ...
                'gamma', gamma, ...
                'r2fit', r2fit);
            
        end
        
        
        function obj = fitMSD(obj, clip_factor)
            %%FITMSD Fit all MSD curves by a linear function.
            %
            % obj = obj.fitMSD fits all MSD curves by a straight line 
            %                      y = a * x + b.
            % The fit is therefore rigorously valid only for purely
            % diffusive behavior.
            %
            % Results are stored in the 'fit' field of the returned
            % object. It is a structure with 2 fields:
            % - a: all the values of the slope of the linear fit.
            % - b: all the values for the intersect of the linear fit.
            % - r2fit: the adjusted R2 value as a indicator of the goodness 
            % of the fit.
            %
            % obj = obj.fitMSD(clip_factor) does the fit, taking into
            % account only the first potion of the average MSD curve
            % specified by 'clip_factor' (a double between 0 and 1). If the
            % value exceeds 1, then the clip factor is understood to be the
            % maximal number of point to take into account in the fit. By
            % default, it is set to 0.25.
            
            
            if nargin < 2
                clip_factor = 0.25;
            end
            
            if ~obj.msd_valid
                obj = obj.computeMSD;
            end
            n_spots = numel(obj.msd);
            
            if clip_factor < 1
                fprintf('Fitting %d curves of MSD = f(t), taking only the first %d%% of each curve... ',...
                    n_spots, ceil(100 * clip_factor) )
            else
                fprintf('Fitting %d curves of MSD = f(t), taking only the first %d points of each curve... ',...
                    n_spots, round(clip_factor) )
            end
            
            a = NaN(n_spots, 1);
            b = NaN(n_spots, 1);
            r2fit = NaN(n_spots, 1);
            ft = fittype('poly1');
            
            fprintf('%4d/%4d', 0, n_spots);
            for i_spot = 1 : n_spots
                
                fprintf('\b\b\b\b\b\b\b\b\b%4d/%4d', i_spot, n_spots);
                
                msd_spot = obj.msd{i_spot};
                
                t = msd_spot(:,1);
                y = msd_spot(:,2);
                w = msd_spot(:,4);
                
                % Clip data, never take the first one dt = 0
                if clip_factor < 1
                    t_limit = 2 : round(numel(t) * clip_factor);
                else
                    t_limit = 2 : min(1+round(clip_factor), numel(t));
                end
                t = t(t_limit);
                y = y(t_limit);
                w = w(t_limit);
                
                % Thrash bad data
                nonnan = ~isnan(y);
                x = t(nonnan);
                y = y(nonnan);
                w = w(nonnan);
                
                if numel(y) < 2
                    continue
                end
                
                [fo, gof] = fit(x, y, ft, 'Weights', w);
                
                a(i_spot) = fo.p1;
                b(i_spot) = fo.p2;
                r2fit(i_spot) = gof.adjrsquare;
                
            end
            fprintf('\b\b\b\b\b\b\b\b\bDone.\n')
            
            obj.lfit = struct(...
                'a', a, ...
                'b', b, ...
                'r2fit', r2fit);
            
        end
        
        
        function  varargout = fitMeanMSD(obj, clip_factor)
            %%FITMEANMSD Fit the weighted averaged MSD by a linear function.
            %
            % obj.fitMeanMSD computes and fits the weighted mean MSD by a
            % straight line y = a * x. The fit is therefore valid only for
            % purely diffusive behavior. Fit results are displayed in the
            % command window.
            %
            % obj.fitMeanMSD(clip_factor) does the fit, taking into account
            % only the first potion of the average MSD curve specified by
            % 'clip_factor' (a double between 0 and 1). If the value
            % exceeds 1, then the clip factor is understood to be the
            % maximal number of point to take into account in the fit. By
            % default, it is set to 0.25.
            %
            % [fo, gof] = obj.fitMeanMSD(...) returns the fit object and the
            % goodness of fit.
            
            if nargin < 2
                clip_factor = 0.25;
            end
            
            if ~obj.msd_valid
                obj = obj.computeMSD;
            end
            
            ft = fittype('poly1');
            mmsd = obj.getMeanMSD;
            
            t = mmsd(:,1);
            y = mmsd(:,2);
            w = 1./mmsd(:,3);
            
            % Clip data, never take the first one dt = 0
            if clip_factor < 1
                t_limit = 2 : round(numel(t) * clip_factor);
            else
                t_limit = 2 : min(1+round(clip_factor), numel(t));
            end
            t = t(t_limit);
            y = y(t_limit);
            w = w(t_limit);
            
            [fo, gof] = fit(t, y, ft, 'Weights', w);
            
            ci = confint(fo);
            str = sprintf([
                'Estimating D through linear weighted fit of the mean MSD curve.\n', ...
                'D = %.3e with 95%% confidence interval [ %.3e - %.3e ].\n', ...
                'Goodness of fit: R² = %.3f.' ], ...
                fo.p1/2/obj.n_dim, ci(1)/2/obj.n_dim, ci(2)/2/obj.n_dim, gof.adjrsquare);
            disp(str)
            
            if nargout > 0
                varargout{1} = fo;
                if nargout > 1
                    varargout{2} = gof;
                end
            end
            
        end
        
        
        function msmsd = getMeanMSD(obj, indices)
            %%GETMEANMSD Compute the weighted mean of all MSD curves.
            %
            % msd = obj.getMeanMSD computes and return the weighted mean of all
            % MSD curves stored in this object. All possible delays are first
            % derived, and for each delay, a weighted mean is computed from all
            % the MSD curves stored in this object. Weights are set to be the
            % number of points averaged to generate the mean square
            % displacement value at the given delay. Thus, we give more weight
            % to MSD curves with greater certainty (larger number of elements
            % averaged).
            %
            % Results are returned as a N x 4 double array, and ordered as
            % following: [ dT M STD N ] with:
            % - dT the delay vector
            % - M the weighted mean of MSD for each delay
            % - STD the weighted standard deviation
            % - N the number of degrees of freedom in the weighted mean
            % (see http://en.wikipedia.org/wiki/Weighted_mean)
            %
            % msd = obj.getMeanMSD(indices) only takes into account the MSD
            % curves with the specified indices.
            
            if ~obj.msd_valid
                obj = obj.computeMSD(indices);
            end
            
            if nargin < 2 || isempty(indices)
                indices = 1 : numel(obj.msd);
            end
            
            n_tracks = numel(indices);
            
            % First, collect all possible delays
            all_delays = cell(n_tracks, 1);
            for i = 1 : n_tracks
                index = indices(i);
                all_delays{i} = obj.msd{index}(:,1);
            end
            delays = unique( vertcat( all_delays{:} ) );
            n_delays = numel(delays);
            
            % Collect
            sum_weight          = zeros(n_delays, 1);
            sum_weighted_mean   = zeros(n_delays, 1);
            
            % 1st pass
            for i = 1 : n_tracks
                
                index = indices(i);
                
                t = obj.msd{index}(:,1);
                m = obj.msd{index}(:,2);
                n = obj.msd{index}(:,4);
                
                % Do not tak NaNs
                valid = ~isnan(m);
                t = t(valid);
                m = m(valid);
                n = n(valid);
                
                % Find common indices
                [~, index_in_all_delays, ~] = intersect(delays, t);
                
                % Accumulate
                sum_weight(index_in_all_delays)           = sum_weight(index_in_all_delays)         + n;
                sum_weighted_mean(index_in_all_delays)    = sum_weighted_mean(index_in_all_delays)  + m .* n;
            end
            
            % Compute weighted mean
            mmean = sum_weighted_mean ./ sum_weight;
            
            % 2nd pass: unbiased variance estimator
            sum_weighted_variance = zeros(n_delays, 1);
            sum_square_weight     = zeros(n_delays, 1);
            
            for i = 1 : n_tracks
                
                index = indices(i);
                
                t = obj.msd{index}(:,1);
                m = obj.msd{index}(:,2);
                n = obj.msd{index}(:,4);
                
                % Do not tak NaNs
                valid = ~isnan(m);
                t = t(valid);
                m = m(valid);
                n = n(valid);
                
                % Find common indices
                [~, index_in_all_delays, ~] = intersect(delays, t);
                
                % Accumulate
                sum_weighted_variance(index_in_all_delays)    = sum_weighted_variance(index_in_all_delays)  + n .* (m - mmean(index_in_all_delays)).^2 ;
                sum_square_weight(index_in_all_delays)        = sum_square_weight(index_in_all_delays)      + n.^2;
            end
            
            % Standard deviation
            mstd = sqrt( sum_weight ./ (sum_weight.^2 - sum_square_weight) .* sum_weighted_variance );
            
            % Output [ T mean std Nfreedom ]
            msmsd = [ delays mmean mstd (sum_weight.^2 ./ sum_square_weight) ];
            
        end
        
        
        function varargout = plotMeanMSD(obj, ha, errorbar, indices)
            %%PLOTMEANMSD Plot the weighted mean of the MSD curves.
            %
            % obj,plotMeanMSD computes and plots the weighted of all MSD
            % curves. See msdanalyzer.getMeanMSD.
            %
            % obj,plotMeanMSD(ha) plots the curve in the axes with the
            % specified handle.
            %
            % obj,plotMeanMSD(ha, errorbar) where 'errorbar' is a boolean allow
            % to specify whether to plot the curve with error bars indicating
            % the weighted standard deviation. Default is false.
            %
            % obj,plotMeanMSD(ha, errorbar, indices) computes and plots the
            % mean only fothe MSD curves whose indices are given on the
            % 'indices' array.
            %
            % h = obj,plotMeanMSD(...) returns the handle to the line plotted.
            %
            % [h, ha] = obj,plotMeanMSD(...) also returns the handle of the
            % axes in which the curve was plotted.
            
            if nargin < 4
                indices = [];
            end
            
            msmsd = obj.getMeanMSD(indices);
            
            if nargin < 3
                errorbar = false;
                if nargin < 2
                    ha = gca;
                end
            end
            
            if errorbar
                h = msdanalyzer.errorShade(ha, msmsd(:,1), msmsd(:,2), msmsd(:,3), [0 0 0], false);
                set(h.mainLine, 'LineWidth', 2);
                
            else
                h = plot(ha, msmsd(:,1), msmsd(:,2), 'k', ...
                    'LineWidth', 2);
            end
            
            obj.labelPlotMSD(ha);
            
            if nargout > 0
                varargout{1} = h;
                if nargout > 1
                    varargout{2} = ha;
                end
            end
            
        end
        
        
        function obj = computeDrift(obj, method, extra, interpmethod)
            %%COMPUTEDRIFT Compute and store drift correction.
            %
            % obj = obj.computeDrift(method) computes and stores the drift
            % using one of the 4 following methods:
            %
            % 'clear' does not compute drift and remove any prior drift
            % computation results.
            %
            % 'manual' allow to specify manually the drift vector:
            % obj = obj.computeDrift('manual', dv); where dv is a double array
            % of size N x (nDim+1) (nDim being the problem dimensionality), and
            % must be arranged as following: [ Ti Xi Yi ... ] etc...
            % On top of this, the drift vector must cover all the possible time
            % points specified in the tracks field of this object: It must
            % start before the first point and end after the last one,
            % otherwise an error is thrown.
            %
            % Missing values within these extremities are interpolated using a
            % linear interpolation scheme by default. To specify another
            % interpolation scheme, use the following syntax:
            % obj.computeDrift('manual', dv, interpmethod), with interpmethod
            % being any value accepted by interp1 ('linear', 'nearest',
            % 'spline', 'pchip', 'cubic').
            %
            % 'centroid' derives the drift by computing the center of mass of
            % all particles at each time point. This method best work for a
            % large number of particle and when the same number of particles is
            % found at every time points. It fails silently otherwise.
            %
            % 'velocity' derives drift by computing instantaneous velocities
            % and averaging them all together, at each time point. If the
            % particles are in sufficient number, and if their real movement is
            % indeed uncorrelated, the uncorrelated part will cancel when
            % averaging, leaving only the correlated part. We assume this part
            % is due to the drift. This method is more robust than the
            % 'centroid' method against particle disappearance and appearance.
            %
            % Results are stored in the 'drift' field of the returned object.
            % It is a double array of size N x (nDim+1) (nDim being the problem
            % dimensionality), and must be arranged as following: [ Ti Xi Yi ... ]
            % etc. If present, it will by used for any call to computeMSD and
            % computeVCorr methods.
            
            
            % First compute common time points
            time = obj.getCommonTimes();
            n_times = numel(time);
            n_tracks = numel(obj.tracks);
            
            switch lower(method)
                
                case 'manual'
                    
                    if nargin < 4
                        interpmethod = 'linear';
                    end
                    
                    drift_dim = size(extra, 2);
                    if drift_dim ~= obj.n_dim + 1
                        error('msdanalyzer:computeDrift:BadDimensionality', ...
                            'Drift must be of size N x (nDim+1) with [ T0 X0 Y0 ... ] etc...');
                    end
                    
                    uninterpolated_time = extra(:, 1);
                    uninterpolated_drift = extra(:, 2:end);
                    if min(uninterpolated_time) > min(time) || max(uninterpolated_time) < max(time)
                        error('msdanalyzer:computeDrift:BadTimeVector', ...
                            'For manual drift correction, time vector must cover all time vector from all tracks.');
                    end
                    
                    ldrift = interp1(...
                        uninterpolated_time, ...
                        uninterpolated_drift, ...
                        time, interpmethod);
                    
                    obj.drift = [time ldrift];
                    
                case 'centroid'
                    
                    ldrift = zeros(n_times, obj.n_dim);
                    n_drift = zeros(n_times, 1);
                    for i = 1 : n_tracks
                        
                        t = obj.tracks{i}(:,1);
                        t = msdanalyzer.roundn(t, msdanalyzer.TOLERANCE);
                        
                        % Determine target time index in bulk
                        [~, index_in_all_tracks_time, ~] = intersect(time, t);
                        
                        % Add to mean accum for these indexes
                        n_drift(index_in_all_tracks_time) = n_drift(index_in_all_tracks_time) + 1;
                        ldrift(index_in_all_tracks_time, :) = ldrift(index_in_all_tracks_time, :) + obj.tracks{i}(:, 2:end);
                        
                    end
                    
                    ldrift = ldrift ./ repmat(n_drift, [1 obj.n_dim]);
                    obj.drift = [time ldrift];
                    
                    
                case 'velocity'
                    
                    sum_V = zeros(n_times, obj.n_dim);
                    n_V = zeros(n_times, 1);
                    
                    for i = 1 : n_tracks
                        
                        t = obj.tracks{i}(:,1);
                        t = msdanalyzer.roundn(t, msdanalyzer.TOLERANCE);
                        
                        % Determine target time index in bulk
                        [~, index_in_all_tracks_time, ~] = intersect(time, t);
                        
                        % Remove first element
                        index_in_all_tracks_time(1) = [];
                        
                        % Compute speed
                        V = diff( obj.tracks{i}(:, 2:end) ) ./ repmat(diff(t), [ 1 obj.n_dim]);
                        
                        % Add to mean accum for these indexes
                        n_V(index_in_all_tracks_time) = n_V(index_in_all_tracks_time) + 1;
                        sum_V(index_in_all_tracks_time, :) = sum_V(index_in_all_tracks_time, :) + V;
                        
                    end
                    
                    % Build accumulated drift
                    sum_V(1, :) = 0;
                    n_V(1, :) = 1;
                    % Integrate
                    d_time = [0; diff(time) ];
                    ldrift = cumsum( sum_V ./ repmat(n_V, [1 obj.n_dim]) .* repmat(d_time, [1 obj.n_dim]), 1);
                    obj.drift = [time ldrift];
                    
                case 'clear'
                    
                    obj.drift = [];
                    
                    
                otherwise
                    error('msdanalyzer:computeDriftCorrection:UnknownCorrectionMethod', ...
                        'Unknown correction method %s. Must be ''clear'', ''manual'', ''centroid'' or ''velocity''.', ...
                        method);
            end
            
            obj.msd_valid = false;
            obj.vcorr_valid = false;
            
        end
        
        
        function velocities = getVelocities(obj, indices)
            %%GETVELOCITIES Generate and return the instantaneous velocities.
            %
            % v = obj.getVelocities returns in v the instantaneous velocities
            % calculated over all the particles tracjectories stored in this
            % object, using dX/dt. The velocities are corrected for drift if
            % this object holds a proper drift field.
            %
            % This method returns a cell array, one cell per particle. Arrays
            % are N x (Ndim+1) double arrays, with Ndim the dimensionality set
            % at object creation. Data is organized as follow:  [ Ti Vxi Vyi ... ].
            %
            % v = obj.getVelocities(indices) restrict the calculation over only
            % the particles with specified indices. Use an empty array to use
            % take all.
            
            if nargin < 2 || isempty(indices)
                indices = 1 : numel(obj.tracks);
            end
            
            n_tracks = numel(indices);
            velocities = cell(n_tracks, 1);
            
            for i = 1 : n_tracks
                
                index = indices(i);
                
                t = obj.tracks{index}(:, 1);
                X = obj.tracks{index}(:, 2:end);
                
                % Determine drift correction
                if ~isempty(obj.drift)
                    tdrift = obj.drift(:, 1);
                    xdrift = obj.drift(:, 2:end);
                    % Determine target delay index in bulk
                    [~, index_in_drift_time, index_in_track_time] = intersect(tdrift, t);
                    % Keep only track times that can be corrected.
                    X = X(index_in_track_time, :);
                    t = t(index_in_track_time);
                    % Subtract drift position to track position
                    X = X - xdrift(index_in_drift_time, :);
                    
                end
                
                dX = diff(X, 1) ./ repmat(diff(t), [1 obj.n_dim]);
                velocities{i} = [ t(1:end-1) dX];
            end
            
        end
        
        
        function obj = computeVCorr(obj, indices)
            %%COMPUTEVCORR Compute velocity autocorrelation.
            %
            % obj = obj.computeVCorr computes the velocity autocorrelation for all
            % the particles trajectories stored in this object. Velocity
            % autocorrelation is defined as vc(t) = < v(i+t) x v(i) >, the mean
            % being taken over all possible pairs inside a trajectories.
            %
            % Results are stored in the 'vcorr' field of the returned object.
            % The velocity autocorrelation is stored for each particles in a
            % cell array, one cell per particle. The array is a double array of
            % size N x 4, and is arranged as follow: [dt mean std N ; ...]
            % where dt is the delay for the autocorrelation, mean is the mean
            % autocorrelation value for this delay, std the standard deviation
            % and N the number of points in the average.
            %
            % obj = obj.computeVCorr(indices) computes the velocity
            % autocorrelation only for the particles with the specified
            % indices. Use an empty array to take all particles.
            
            obj.vcorr = cell(numel(obj.tracks), 1);
            
            if nargin < 2 || isempty(indices)
                indices = 1 : numel(obj.tracks);
            end
            
            % Get instantaneous velocities
            velocities = obj.getVelocities(indices);
            delays = obj.getAllDelays(indices);
            n_delays = numel(delays);
            n_tracks = numel(velocities);
            
            fprintf('Computing velocity autocorrelation of %d tracks... ', n_tracks);
            fprintf('%4d/%4d', 0, n_tracks);
            for i = 1 : n_tracks
                fprintf('\b\b\b\b\b\b\b\b\b%4d/%4d', i, n_tracks);
                
                % Holder for mean, std calculations
                sum_vcorr     = zeros(n_delays-1, 1);
                sum_vcorr2    = zeros(n_delays-1, 1);
                n_vcorr       = zeros(n_delays-1, 1);
                
                % Unwrap data
                vc = velocities{i};
                t = vc(:, 1);
                V = vc(:, 2:end);
                n_detections = size(V, 1);
                
                % First compute velocity correleation at dt = 0 over all tracks
                Vc0     = mean( sum( V.^2, 2) );
                
                % Other dts
                for j = 1 : n_detections - 1
                    
                    % Delay in physical units
                    dt = t(j+1:end) - t(j);
                    dt = msdanalyzer.roundn(dt, msdanalyzer.TOLERANCE);
                    
                    % Determine target delay index in bulk
                    [~, index_in_all_delays, ~] = intersect(delays, dt);
                    
                    % Velocity correlation in bulk
                    lvcorr = sum( repmat(V(j, :), [ (n_detections-j) 1]) .* V(j+1:end, :), 2 );
                    
                    % Normalize
                    lvcorr = lvcorr ./ Vc0;
                    
                    % Store for mean computation
                    sum_vcorr(index_in_all_delays)   = sum_vcorr(index_in_all_delays) + lvcorr;
                    sum_vcorr2(index_in_all_delays)  = sum_vcorr2(index_in_all_delays) + (lvcorr .^ 2);
                    n_vcorr(index_in_all_delays)     = n_vcorr(index_in_all_delays) + 1;
                    
                end
                
                mean_vcorr = sum_vcorr ./ n_vcorr;
                std_vcorr = sqrt( (sum_vcorr2 ./ n_vcorr) - (mean_vcorr .* mean_vcorr)) ;
                vcorrelation = [ delays(1:end-1) mean_vcorr std_vcorr n_vcorr ];
                vcorrelation(1,:) = [0 1 0 n_detections];
                
                % Store in object field
                index = indices(i);
                obj.vcorr{index} = vcorrelation;
                
            end
            fprintf('\b\b\b\b\b\b\b\b\bDone.\n')
            
            obj.vcorr_valid = true;
            
        end
        
        
        function msvcorr = getMeanVCorr(obj, indices)
            %%GETMEANVCORR Compute the weighted mean of velocity autocorrelation.
            %
            % msd = obj.getMeanVCorr computes and return the weigthed mean of
            % all velocity autocorrelation curves stored in this object. All
            % possible delays are first derived, and for each delay, a weighted
            % mean is computed from all the velocity autocorrelation curves
            % stored in this object. Weights are set to be the number of points
            % averaged to generate the mean square displacement value at the
            % given delay. Thus, we give more weight to velocity
            % autocorrelation curves with greater certainty (larger number of
            % elements averaged).
            %
            % Results are returned as a N x 4 double array, and ordered as
            % following: [ dT M STD N ] with:
            % - dT the delay vector
            % - M the weighted mean of velocity autocorrelation for each delay
            % - STD the weighted standard deviation
            % - N the number of degrees of freedom in the weighted mean
            % (see http://en.wikipedia.org/wiki/Weighted_mean)
            %
            % msd = obj.getMeanVCorr(indices) only takes into account the
            % velocity autocorrelation curves with the specified indices,
            
            if ~obj.vcorr_valid
                obj = obj.computeVCorr(indices);
            end
            
            if nargin < 2 || isempty(indices)
                indices = 1 : numel(obj.vcorr);
            end
            
            n_tracks = numel(indices);
            
            % First, collect all possible delays
            all_delays = cell(n_tracks, 1);
            for i = 1 : n_tracks
                index = indices(i);
                all_delays{i} = obj.vcorr{index}(:,1);
            end
            delays = unique( vertcat( all_delays{:} ) );
            n_delays = numel(delays);
            
            
            % Collect
            sum_weight          = zeros(n_delays, 1);
            sum_weighted_mean   = zeros(n_delays, 1);
            
            % 1st pass
            for i = 1 : n_tracks
                
                index = indices(i);
                
                t = obj.vcorr{index}(:,1);
                m = obj.vcorr{index}(:,2);
                n = obj.vcorr{index}(:,4);
                
                % Do not tak NaNs
                valid = ~isnan(m);
                t = t(valid);
                m = m(valid);
                n = n(valid);
                
                % Find common indices
                [~, index_in_all_delays, ~] = intersect(delays, t);
                
                % Accumulate
                sum_weight(index_in_all_delays)           = sum_weight(index_in_all_delays)         + n;
                sum_weighted_mean(index_in_all_delays)    = sum_weighted_mean(index_in_all_delays)  + m .* n;
            end
            
            % Compute weighted mean
            mmean = sum_weighted_mean ./ sum_weight;
            
            % 2nd pass: unbiased variance estimator
            sum_weighted_variance = zeros(n_delays, 1);
            sum_square_weight     = zeros(n_delays, 1);
            
            for i = 1 : n_tracks
                
                index = indices(i);
                
                t = obj.vcorr{index}(:,1);
                m = obj.vcorr{index}(:,2);
                n = obj.vcorr{index}(:,4);
                
                % Do not tak NaNs
                valid = ~isnan(m);
                t = t(valid);
                m = m(valid);
                n = n(valid);
                
                % Find common indices
                [~, index_in_all_delays, ~] = intersect(delays, t);
                
                % Accumulate
                sum_weighted_variance(index_in_all_delays)    = sum_weighted_variance(index_in_all_delays)  + n .* (m - mmean(index_in_all_delays)).^2 ;
                sum_square_weight(index_in_all_delays)        = sum_square_weight(index_in_all_delays)      + n.^2;
            end
            
            % Standard deviation
            mstd = sqrt( sum_weight ./ (sum_weight.^2 - sum_square_weight) .* sum_weighted_variance );
            
            % Output [ T mean std Nfreedom ]
            msvcorr = [ delays mmean mstd (sum_weight.^2 ./ sum_square_weight) ];
            
        end
        
        
        function varargout = plotMeanVCorr(obj, ha, errorbar, indices)
            %%PLOTMEANVCORR Plot the weighted mean of the velocity autocorrelation curves.
            %
            % obj,plotMeanVCorr computes and plots the weighted of all velocity
            % autocorrelation curves. See msdanalyzer.getMeanVCorr.
            %
            % obj,plotMeanVCorr(ha) plots the curve in the axes with the
            % specified handle.
            %
            % obj,plotMeanVCorr(ha, errorbar) where 'errorbar' is a boolean
            % allow to specify whether to plot the curve with error bars
            % indicating the weighted standard deviation. Default is false.
            %
            % obj,plotMeanVCorr(ha, errorbar, indices) computes and plots the
            % mean only fothe velocity autocorrelation curves whose indices are
            % given on the 'indices' array.
            %
            % h = obj,plotMeanVCorr(...) returns the handle to the line plotted.
            %
            % [h, ha] = obj,plotMeanVCorr(...) also returns the handle of the
            % axes in which the curve was plotted.
            
            if nargin < 4
                indices = [];
                
                if nargin < 3
                    errorbar = false;
                    if nargin < 2
                        ha = gca;
                    end
                end
            end
            mvc = obj.getMeanVCorr(indices);
            
            if errorbar
                h = msdanalyzer.errorShade(ha, mvc(:,1), mvc(:,2), mvc(:,3), [0 0 0], false);
                set(h.mainLine, 'LineWidth', 2);
                
            else
                h = plot(ha, mvc(:,1), mvc(:,2), 'k', ...
                    'LineWidth', 2);
            end
            
            obj.labelPlotVCorr(ha);
            
            if nargout > 0
                varargout{1} = h;
                if nargout > 1
                    varargout{2} = ha;
                end
            end
            
        end
        
    end
    
    %% Private methods
    
    methods (Access = private)
        
        function time = getCommonTimes(obj)
            
            n_tracks = numel(obj.tracks);
            times = cell(n_tracks, 1);
            for i = 1 : n_tracks
                times{i} = obj.tracks{i}(:,1);
            end
            time = unique( vertcat(times{:}) );
            time = msdanalyzer.roundn(time, msdanalyzer.TOLERANCE);
            
        end
        
        function delays = getAllDelays(obj, indices)
            % First, find all possible delays in time vectors.
            % Time can be arbitrary spaced, with frames missings,
            % non-uniform sampling, etc... so we have to do this clean.
            
            if nargin < 2 || isempty(indices)
                indices = 1 : numel(obj.tracks);
            end
            
            n_tracks = numel(indices);
            all_delays = cell(n_tracks, 1);
            for i = 1 : n_tracks
                index = indices(i);
                track = obj.tracks{index};
                t = track(:,1);
                [T1, T2] = meshgrid(t, t);
                dT = msdanalyzer.roundn(abs(T1(:)-T2(:)), msdanalyzer.TOLERANCE);
                all_delays{i} = unique(dT);
            end
            delays = unique( vertcat(all_delays{:}) );
        end
        
    end
    
    %% Static methods
    
    methods (Static, Access = private)
        
        function wm = weightedmean(x, w)
            wm = sum( x .* w) ./ sum(w);
        end
        
        function sewm = standarderrorweightedmean(x, w)
            n = numel(w);
            wbar = mean(w);
            xbar = sum( x .* w) ./ sum(w);
            sewm = n /((n-1) * sum(w)^2) * (sum( (w.*x - wbar*xbar).^2) ...
                - 2 * xbar * sum( (w-wbar).*(w.*x - wbar*xbar)) ...
                + xbar^2 * sum((w-wbar).^2));
        end
        
        
        
        function H = errorShade(ha, x, y, errBar, col, transparent)
            % Adapted from Rob Campbell code, at:
            % http://www.mathworks.com/matlabcentral/fileexchange/26311-shadederrorbar/content/shadedErrorBar.m
            hold on
            H.mainLine = plot(ha, x, y, 'Color', col);
            
            edgeColor = col + (1-col) * 0.55;
            patchSaturation = 0.15; %How de-saturated or transparent to make the patch
            if transparent
                faceAlpha=patchSaturation;
                patchColor=col;
                set(gcf,'renderer','openGL')
            else
                faceAlpha=1;
                patchColor=col+(1-col)*(1-patchSaturation);
                set(gcf,'renderer','painters')
            end
            
            %Calculate the y values at which we will place the error bars
            uE = y + errBar;
            lE = y - errBar;
            
            %Make the cordinats for the patch
            yP = [ lE ; flipud(uE) ];
            xP = [ x ; flipud(x) ];
            
            invalid = isnan(xP) | isnan(yP) | isinf(xP) | isinf(yP);
            yP(invalid) = [];
            xP(invalid) = [];
            
            
            H.patch = patch(xP, yP, 1, ...
                'Facecolor', patchColor,...
                'Edgecolor', 'none',...
                'Facealpha', faceAlpha, ...
                'Parent', ha);
            
            %Make nice edges around the patch.
            H.edge(1) = plot(ha, x, lE, '-', 'Color', edgeColor);
            H.edge(2) = plot(ha, x, uE, '-', 'Color', edgeColor);
            
            %The main line is now covered by the patch object and was plotted first to
            %extract the RGB value of the main plot line. I am not aware of an easy way
            %to change the order of plot elements on the graph so we'll just remove it
            %and put it back (yuk!)
            delete(H.mainLine)
            H.mainLine = plot(ha, x, y, 'Color', col);
        end
        
    end
    
    %% Static and public functions
    
    methods (Static)
        
        %ROUNDN  Round towards nearest number with Nth decimals.
        %   ROUNDN(X,N) rounds the elements of X to the nearest numbers with the
        %   precision given by N.
        %
        %   Examples:   roundn(8.73,0) = 9
        %               roundn(8.73,1) = 8.7
        %               roundn(8.73,2) = 8.73
        %               roundn(8.73,-1) = 10
        %
        %   See also ROUND
        % Jean-Yves Tinevez - MPI-CBG - August 2007
        
        function Y = roundn(X,N)
            Y = 10^(-N) .* round(X.*10^(N));
        end
        
        function newobj = pool(msda_arr)
            %%POOL Pool the data of several masanalyzer objects in a new one.
            
            if ~isa(msda_arr, 'msdanalyzer')
                error('msdanalyzer:pool:BadArgument', ...
                    'Expect arguments to be of class ''masanalyzer'', got ''%s''.',...
                    class(msda_arr))
            end
            
            n_obj = numel(msda_arr);
            
            % Check dimensionality and units consistency
            lspace_units = msda_arr(1).space_units;
            ltime_units = msda_arr(1).time_units;
            ln_dim = msda_arr(1).n_dim;
            
            for i = 2 : n_obj
                
                obj = msda_arr(i);
                if ~strcmp(obj.space_units, lspace_units)
                    error('msdanalyzer:pool:InconsistentArray', ...
                        'Element %d is inconsistent. Expected space units to be %s, got %s.', ...
                        i, lspace_units, obj.space_units)
                end
                if ~strcmp(obj.time_units, ltime_units)
                    error('msdanalyzer:pool:InconsistentArray', ...
                        'Element %d is inconsistent. Expected time units to be %s, got %s.', ...
                        i, ltime_units, obj.time_units)
                end
                if obj.n_dim ~= ln_dim
                    error('msdanalyzer:pool:InconsistentArray', ...
                        'Element %d is inconsistent. Expected dimensionality to be %d, got %d.', ...
                        i, ln_dim, obj.n_dim)
                end
                
            end
            
            all_msd = cell(n_obj,1 );
            all_vcorr = cell(n_obj,1 );
            
            for i = 1 : n_obj
                obj = msda_arr(i);
                
                obj = obj.computeDrift('velocity');
                obj = obj.computeMSD;
                obj = obj.computeVCorr;
                
                all_msd{i} = obj.msd;
                all_vcorr{i} = obj.vcorr;
            end
            
            newobj = msdanalyzer(ln_dim, lspace_units, ltime_units);
            newobj.msd = vertcat(all_msd{:});
            newobj.vcorr = vertcat(all_vcorr{:});
            newobj.msd_valid = true;
            newobj.vcorr_valid = true;
            
        end
        
    end
    
end

