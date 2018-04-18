% gac.m
% Growth and Aggregation of Clusters model
%

function [V_arr,total_pop_arr] = gac(lr,la,lc,ls,dt,Tmax,n0)

%% default values for input parameters

% rate of cluster growth
if ~exist('lr','var')||isempty(lr)
    lr = 1.4;
end

% rate of cluster aggregation
if ~exist('la','var')||isempty(la)
    la = .1;
end

% rate of cluster collapse
if ~exist('lc','var')||isempty(lc)
    lc = .1;
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

% shuffle random number generator
rng('shuffle');



%% intialize some arrays

% main object of the simulation:  array of cluster volumes
V_arr = ones(1,n0);

% time vector
tvec = 0:dt:Tmax;

% number of simulation steps
numsteps = numel(tvec);

% array for keeping track of total population size
total_pop_arr = zeros(1,numsteps);
total_pop_arr(1) = n0;

% pregenerate lots of random numbers for speed.
number_of_random_numbers_to_pre_generate = 1000.*numsteps;
lots_of_random_numbers = rand(1,number_of_random_numbers_to_pre_generate);

% counter for how many of the pre-generated random numbers have been used.
n = 0;

%% fixed parameters
% method for regulating growth of clusters.
%	"global' =  cluster growth is regulated by the total
%               population size, max_tot_pop.  
%   'local' =   cluster growth is regulated by local
%               constraints, i.e. space or local nutrient depletion, via a maximum
%               volume for a cluster, Vmax.
%
cluster_growth_regulation_method = 'global';

% carrying capacity of whole system.  regulates cluster growth in 'global' mode.  
max_total_pop = 10^(4);                       

% maximum volume of individual clump.  regulates cluster growth in 'local' model, and
% also sets scale for how probability of collapse depends on cluster size.
Vmax = max_total_pop;%1e4;

% NOT USED YET.exponent for determining aggregation probability.  Pa ~ V^theta.
% theta = 0:    all clusters have equal prob. of aggregating.
% theta = 2/3:  Pa ~ surface area of cluster.
% theta = 1:    Pa ~ volume of cluster.
theta = 2/3;

% logical for stochastic vs deterministic growth
l_stochastic_growth = false;

% logical for printing progress to screen
l_print_progress = false;

%% main simulation.  loop over time.
for s = 2:numsteps
    
    % print update to console
    if l_print_progress
        disp(['step number ' num2str(s) ' number of clusters = ' num2str(numel(V_arr))])
    end   
    
    %% growth
    
    switch l_stochastic_growth
        case true
            
            % collect random numbers.  If we've run out of pregenerated ones,
            % generate some new ones.
            if n + numel(V_arr) < number_of_random_numbers_to_pre_generate
                random_numbers = lots_of_random_numbers(n+1:n+numel(V_arr));
                n = n + numel(V_arr);
            else
                random_numbers = rand(1,numel(V_arr));
            end
            
            % compute probability of clusters growing in this timestep.  there are
            % currently two methods.  'global' = cluster growth is regulated by the total
            % population size, max_tot_pop.  'local' = cluster growth is regulated by local
            % constraints, i.e. space or local nutrient depletion, via a maximum
            % volume for a cluster, Vmax.
            switch cluster_growth_regulation_method
                
                case 'global'
                    
                    Pr = lr.*(1-sum(V_arr)./max_total_pop).*dt.*ones(size(V_arr));
                    
                case 'local'
                    
                    Pr = lr.*(1-V_arr./Vmax).*dt;
                    
                case 'both'
                    
                    Pr = lr.*(1-V_arr./Vmax).*dt.*(1-sum(V_arr)./max_total_pop);
                    
                    
            end
            
            % grow clusters in one vectorized step
            V_arr(random_numbers <= Pr) = 2.*V_arr(random_numbers <= Pr);
    
    case false
            
        switch cluster_growth_regulation_method
            
            case 'global'
                % grow clusters in one vectorized step
                V_arr = V_arr + dt.*lr./log2(exp(1)).*V_arr.*(1-sum(V_arr)./max_total_pop);
            
            case 'local'
                
                % grow clusters in one vectorized step
                V_arr = V_arr + dt.*lr./log2(exp(1)).*V_arr.*(1-V_arr./Vmax);
                
            case 'both'
                
                % grow clusters in one vectorized step
                V_arr = V_arr + dt.*lr./log2(exp(1)).*V_arr.*(1-sum(V_arr)./max_total_pop).*(1-V_arr./Vmax);
        end
        
        V_arr = V_arr(V_arr>=1);
                
    end
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
    % here, the relationship is linear and scaled by Vmax.  
    %Pc = lc.*V_arr./Vmax.*dt;
    
    % here, the relationship is linear in linear length dimension and scaled by Vmax.
    Pc = lc.*(V_arr./Vmax).^(1/3).*dt;

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
    Pa = la.*dt;
    
    % find ids of clusters that will aggregate in this timestep
    agg_ids = find(random_numbers < Pa);
    
    % loop over aggregates
    for c = 1:numel(agg_ids)
        
        % agg_ids can change size, need this check
        if c > numel(agg_ids)
            break
        end
        
        % choose a cluster at random
        random_cluster_id = randi([1,numel(V_arr)]);
        
        %if sum(ismember([agg_ids,0],random_cluster_id)) == 0
        if random_cluster_id ~= agg_ids(c)
            
            % add cluster sizes
            V_arr(agg_ids(c)) = V_arr(agg_ids(c)) + V_arr(random_cluster_id);
            
            % mark randomly chosen cluster for removal with a 0
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
    Ps = ls.*dt.*(V_arr >= 2).*(V_arr).^(2/3);
    
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