% gac.m
% Growth and Aggregation of Clusters model
%

function [V_arr,total_pop_arr] = gac(lr,la,lc,ls,dt,Tmax,n0,nu_A,nu_F)

%% default values for input parameters

% rate of cluster growth
if ~exist('lr','var')||isempty(lr)
    lr = 1;
end

% rate of cluster aggregation
if ~exist('la','var')||isempty(la)
    la = .1;
end

% rate of cluster collapse
if ~exist('lc','var')||isempty(lc)
    lc = .005;
end

% rate of cluster sprouting
if ~exist('ls','var')||isempty(ls)
    ls = .1;
end

% numerical time step
if ~exist('dt','var')||isempty(dt)
    dt = .1;
end

% total simulation time
if ~exist('Tmax','var')||isempty(Tmax)
    Tmax = 24;
end

% initial starting population of single cells
if ~exist('n0','var')||isempty(n0)
    n0 = 10;
end

% aggregation exponent
if ~exist('nu_A','var')||isempty(nu_A)
    nu_A = 2/3;
end

% fragmentation exponent
if ~exist('nu_F','var')||isempty(nu_F)
    nu_F = 2/3;
end
% shuffle random number generator
rng('shuffle');


%% intialize some arrays

% main object of the simulation:  array of cluster volumes
if numel(n0)==1
    V_arr = ones(1,n0);
elseif size(n0,1) > 1
    disp('gac error:  initial cluster size array n0 must have dim 1xM');
    return
else
    V_arr = n0;
end

% time vector
tvec = 0:dt:Tmax;

% number of simulation steps
numsteps = numel(tvec);

% array for keeping track of total population size
total_pop_arr = zeros(1,numsteps);
total_pop_arr(1) = sum(n0);

% pregenerate lots of random numbers for speed.
number_of_random_numbers_to_pre_generate = round(1000.*numsteps);
lots_of_random_numbers = rand(1,number_of_random_numbers_to_pre_generate);

% counter for how many of the pre-generated random numbers have been used.
n = 0;

%% fixed parameters
% carrying capacity of whole system.  regulates cluster growth in 'global' mode.  
max_total_pop = 10^(4);                       

% logical for printing progress to screen
l_print_progress = false;

%% main simulation.  loop over time.
for s = 2:numsteps
    
    % print update to console
    if l_print_progress
        disp(['step number ' num2str(s) ' number of clusters = ' num2str(numel(V_arr))])
    end   
    
    %% growth
    
    % grow all clusters in one deterministic vectorized step
    V_arr = V_arr + dt.*lr.*V_arr.*(1-sum(V_arr)./max_total_pop);
    
    % if clusters die out, remove them
    V_arr = V_arr(V_arr>=1);
                
    %% collapse 
    
    % collect random numbers.  If we've run out of pregenerated ones,
    % generate some new ones.
    if n + numel(V_arr) < number_of_random_numbers_to_pre_generate
        random_numbers = lots_of_random_numbers(n+1:n+numel(V_arr));
        n = n + numel(V_arr);
    else
        random_numbers = rand(1,numel(V_arr));
    end
    
    % compute probability of collapse as a function of cluster volume.
    % here, the relationship is linear in linear length dimension.
    Pc = lc.*(V_arr).^(1/3).*dt;

    % remove clusters that have collapsed from the cluster array.
    V_arr = V_arr(random_numbers >= Pc);
        
    %% aggregation
    
    % collect random numbers.  If we've run out of pregenerated ones,
    % generate some new ones.
    if n + numel(V_arr) < number_of_random_numbers_to_pre_generate
        random_numbers = lots_of_random_numbers(n+1:n+numel(V_arr));
        n = n + numel(V_arr);
    else
        random_numbers = rand(1,numel(V_arr));
    end
    
    % compute the probability of aggregation in this timestep. 
    % made product kernel 5/1/18
    Pa = la.*dt.*V_arr.^nu_A;
    
    % find ids of clusters that will aggregate in this timestep
    agg_ids = find(random_numbers < Pa);
    
    % loop over aggregates
    for c = 1:numel(agg_ids)
              
        % choose a cluster at random 
        %random_cluster_id = randi([1,numel(V_arr)]);
        % made product kernel 5/1/18: pick a cluster based on its size.
        % construct a discrete probability distribution for clusters by
        % normalizing the cluster size array and inverse project uniform
        % random numbers via the cumulative distribution function.

        rn = rand();
        random_cluster_id = find(rn <= cumsum(V_arr.^nu_A./sum(V_arr.^nu_A)),1);
        
        if random_cluster_id ~= agg_ids(c)
            
            % add cluster sizes
            V_arr(agg_ids(c)) = V_arr(agg_ids(c)) + V_arr(random_cluster_id);
            
            % mark randomly chosen cluster for removal with a 0.  Using the
            % product kernel code above, this guarantees that this cluster
            % will not be chosen again.
            V_arr(random_cluster_id) = 0;
            
            % if the randomly chosen cluster was part of the original
            % set of aggregators, remove it from agg_ids.
            if ismember(agg_ids,random_cluster_id)
                agg_ids = agg_ids(~ismember(agg_ids,random_cluster_id));
            end
            
        end
            
        
    end
    
    % removed clusters that have been absorbed
    V_arr = V_arr(V_arr ~= 0);
    
    %% sprout
    
    % collect random numbers.  If we've run out of pregenerated ones,
    % generate some new ones.
    if n + numel(V_arr) < number_of_random_numbers_to_pre_generate
        random_numbers = lots_of_random_numbers(n+1:n+numel(V_arr));
        n = n + numel(V_arr);
    else
        random_numbers = rand(1,numel(V_arr));
    end
 
    % compute probability of sprouting.  require clusters to be of at least
    % size 2 to sprout.
    Ps = ls.*dt.*(V_arr >= 2).*(V_arr).^(nu_F);
    
    % find ids of clusters that sprouted
    sprout_ids = find(random_numbers < Ps);
    
    % reduce sprouted clusters by 1
    V_arr(sprout_ids) = V_arr(sprout_ids) - 1;
    
    % compute number of newly sprouted clusters
    num_sprouts = sum(random_numbers < Ps);
    
    % add this collection of new single cells to cluster array
    V_arr = [V_arr, ones(1,num_sprouts)];
    
    %% other updates
    
    % update total population array
    total_pop_arr(s) = sum(V_arr);
        
        

end







end