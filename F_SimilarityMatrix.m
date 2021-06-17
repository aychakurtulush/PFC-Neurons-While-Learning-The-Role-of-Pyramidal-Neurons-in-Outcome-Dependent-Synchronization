
function [varargout] = F_SimilarityMatrix(spkdata,Didxs,T,bin,varargin)

[r, ~] = size(Didxs);
if r==1; Didxs = Didxs'; end   % make column vector

% defaults  
% bin = 0.001;    % 1ms bins
opts.BLmeth = 'Gaussian';   % use Gaussian window around spikes 
opts.BLpars = 0.01;        % std dev is 10 ms
opts.Dmeth = 'corrcoef';     % use correlation coefficient
opts.nlimit = 6;       % minimum number of nodes  
opts.modopts = '';   % set no options for the all eigenvector method...
opts.blnS = 0;       % run clustering algorithm by default

if nargin >= 5
    if isstruct(opts) 
        tempopts = varargin{1}; 
        fnames = fieldnames(tempopts);
        for i = 1:length(fnames)
            opts = setfield(opts,fnames{i},getfield(tempopts,fnames{i}));
        end
    end
end

% storage
% Rn = cell(numel(opts.BLpars),1);        % retained neurons...
% Gmax = struct('grps',[],'grpsizes',[],'ngrps',[],'Q',[]);
% Gcon = struct('grps',[],'grpsizes',[],'ngrps',[],'Q',[]);

spkfcn = cell(numel(opts.BLpars),1);

bins = T(1)+bin/2:bin:T(2);


%% analyse data
for loop = 1:numel(opts.BLpars)
    
    % set up convolution window if using...
    sig = opts.BLpars(loop)/bin;
    switch opts.BLmeth
        case 'Gaussian'
            x = -5*sig:1:5*sig';  % x-axis values of the discretely sampled Gaussian, out to 5xSD
            h = (1/(sqrt(2*pi*sig^2)))*exp(-((x.^2*(1/(2*sig^2))))); % y-axis values of the Gaussian
            shiftbase = floor(length(x)/2);
            % keyboard
        case 'exponential'
            x = 0:1:10*sig';   % spread to 10 times the time constant
            h = exp(-x/sig);
            shiftbase = length(x);
    end      
    h = h ./sum(h); % make sure kernel has unit area, then can use for rate functions

    % convolve window with data-spike trains
    [spkfcn{loop},~] = convolve_spiketrains(spkdata,h,shiftbase,Didxs,bins,bin,T,opts);
%     Nidxs = numel(idxs);
    
    %% now compute selected distance between those functions
    [Sxy{loop}] = constructS(spkfcn{loop},opts);
%     Rn{loop} = 1:Nidxs;
% 
%     % keyboard
%     % do modularity on that graph
%     edges(loop) = sum(sum(Sxy{loop}~=0)); 
%     nodes(loop) = numel(Rn{loop});
%     
%     if ~opts.blnS & nodes(loop) >= opts.nlimit & edges(loop) > log(nodes(loop))   
%         % only group if: (a) asked to do so; 
%         % (b) there are enough nodes and
%         % (c) the number of edges ensures a likely fully-connected graph
%         
%         if numel(opts.modopts) == 1
%             [grps,Gmax(loop).Q,grpscon,Gcon(loop).Q,ctr] = allevsplitConTransitive(Sxy{loop},opts.modopts{1}); 
%         elseif numel(opts.modopts) == 2
%             [grps,Gmax(loop).Q,grpscon,Gcon(loop).Q,ctr]  = allevsplitConTransitive(Sxy{loop},opts.modopts{1},opts.modopts{2}); 
%         else
%             [grps,Gmax(loop).Q,grpscon,Gcon(loop).Q,ctr]  = allevsplitConTransitive(Sxy{loop}); 
%         end
%         
%         % Qmax answers
%         Gmax(loop).ngrps = max(grps);
%         % get sizes 
%         siz = [];
%         for i = 1:Gmax(loop).ngrps
%             siz = [siz numel(find(grps == i))];   
%         end
%         Gmax(loop).grpsizes = siz;
%         
%         % Qcon answers
%         Gcon(loop).ngrps = max(grpscon);
%         % get sizes 
%         siz = [];
%         for i = 1:Gcon(loop).ngrps
%             siz = [siz numel(find(grpscon == i))];   
%         end
%         Gcon(loop).grpsizes = siz;
% 
%     else
%         grps = zeros(numel(Rn{loop}),1);  % nothing to group!
%         Gmax(loop).ngrps = 0;
%         Gmax(loop).grpsizes = [];
%         Gmax(loop).Q = 0;
%         grpscon = zeros(numel(Rn{loop}),1);  % nothing to group!
%         Gcon(loop).ngrps = 0;
%         Gcon(loop).grpsizes = [];
%         Gcon(loop).Q = 0;
%     end
% 
%     % remap from retained index count to ID stamps...
%     Gmax(loop).grps = [idxs(Rn{loop}) grps];  % [ID stamp; Group membership]
%     Gcon(loop).grps = [idxs(Rn{loop}) grpscon];  % [ID stamp; Group membership]
      
end


varargout{1} = Sxy;
varargout{2} = spkfcn;
varargout{3} = shiftbase;


function [spkfcn,idxs] = convolve_spiketrains(spkdata,h,shiftbase,Didxs,bins,bin,T,opts)
    Nidxs = numel(Didxs);
    
    %% go round and compute spike-train binless functions
    spkfcn = zeros(numel(bins),Nidxs);
    nspikes = NaN; % just in case there are no spikes....
    
    for j = 1:Nidxs
        currix = find(spkdata(:,1) == Didxs(j));
        nspikes(j) = numel(currix);
        clear spk
        [spk,~] = spike_train_from_times(spkdata(currix,2),bin,T);
        if nspikes(j) > 0   % only bother doing convolution if there's something to convolve!!
            switch opts.BLmeth
                case 'Gaussian'
                    try
                        y = conv(h,spk);
                        [~, c] = size(y); if c>1; y = y'; end  % for reasons best known to Matlab, certain convolutions will return this as a row vector rather than a column vector
                        shifty = y(shiftbase+1:end-shiftbase);   % shift convolved signal to line up with spike-times
                        if numel(shifty) < numel(bins)
                            % pad with zeros
                            diffbins = numel(bins) - numel(shifty);
                            shifty = [zeros(diffbins,1); shifty]; 
                        end % can occasionally happen with width pars that are not integer multiples of step-size 
                        spkfcn(:,j) = shifty;
                    catch
                        disp('I ran into a problem convolving the Gaussian')
                        keyboard
                    end
                    
                case 'exponential'
                    try
                        y = conv(h,spk);
                        [~, c] = size(y); if c>1; y = y'; end  % for reasons best known to Matlab, certain convolutions will return this as a row vector rather than a column vector
                        % keyboard
                        spkfcn(:,j) = y(1:numel(bins));   % truncate convolved signal to line up with recording time 
                    catch
                        disp('I ran into a problem convolving the exponential')
                        keyboard
                    end
            end
        end
        % keyboard
    end
   
    % keyboard
    
    idxs = Didxs; 
    try
        if any(nspikes == 0)
            % if not firing, then strip out cells
            spkfcn(:,nspikes==0) = []; 
            idxs(nspikes==0) = [];
        end
    catch
        disp('I ran into a problem removing non-firing cells')
        keyboard
    end
    
    % convert to firing rate (spikes/s)
    spkfcn = spkfcn ./ bin; % sums to 1 if spike in every bin

function [Sxy] = constructS(spkfcn,opts)

    switch opts.Dmeth            
        case {'corr','corrcoef'}
            Sxy = eval([opts.Dmeth '(spkfcn);']);     % compute correlation
            if any(isnan(Sxy))
                keyboard
            end
%             Sxy(Sxy < 0) = 0;   % rectify correlation coefficient = similarity....
            % place zeros on the diagonal: no self-connections allowed...
%             Sxy(eye(Nidxs)==1) = 0;
%          case 'cosine'
%             try 
%                 % how many entries in Sxy to compute? 
%                 inds = 1:Nidxs^2;
%                 uniqueinds = triu(reshape(inds,Nidxs,Nidxs)); % indexes of all unique pairs
%                 uniqueinds = uniqueinds - diag(diag(uniqueinds)); % remove self-indices
%                 [ip,jp] = ind2sub([Nidxs,Nidxs],uniqueinds(uniqueinds>0)); % get all [i,j] pairs...
%                 Sxy = zeros(Nidxs);
%                 for prs = 1:numel(ip)
%                     A = spkfcn(:,ip(prs)); B = spkfcn(:,jp(prs));   % current pair
%                     Sxy(ip(prs),jp(prs)) = (A'*B)/(norm(A)*norm(B));
%                     Sxy(jp(prs),ip(prs)) = Sxy(ip(prs),jp(prs));
%                 end
% 
%             catch
%                 keyboard
%             end

        otherwise
            error('Unknown pairwise similarity metric specified') 
    end

function [y,bins] = spike_train_from_times(times,bin_size,T)

% SPIKE_TRAIN_FROM_TIMES create spike train array
%
%   SPIKE_TRAIN_FROM_TIMES(A,BINSIZE,T) converts the spike train represented by event times A (in seconds) into
%   a binary spike train array using bins BINSIZE seconds wide, and where T is a two-element array specifying the start
%   and end times in seconds. 
%   
%   [Y,BINS] = SPIKE_TRAIN_FROM_TIMES(...) returns the spike train Y and the bin centers (times) BINS
%
%   Mark Humphries 22/4/2004

% time_seconds = T(2) - T(1);
bins = T(1)+bin_size/2:bin_size:T(2)-bin_size/2;
y = hist(times(times<T(2)),bins);

