% Program:  multigac_tau_constant_kernels.m
%
% cluster_sizes array now becomes num_species x num_clusters array,
% elements define number of cells of each species in each cluster

function [cluster_sizes,abundance_arr,tvec,total_num_clumps_arr] = multigac_tau_constant_kernels(num_species,reaction_rates,Tmax,n0,carrying_capacities,interaction_matrix,tau)

%% default values for input parameters

if ~exist('num_species','var')||isempty(num_species)
    num_species = 2;
end

if ~exist('reaction_rates','var')||isempty(reaction_rates)
    growth_rate = 1;
    aggregation_rate = 0.1;
    expulsion_rate = 0.1;
    fragmentation_rate = 10;
    reaction_rates = repmat([growth_rate, aggregation_rate, expulsion_rate, fragmentation_rate],num_species,1);
elseif size(reaction_rates,2) ~= 4
    disp('error: reaction_rates needs 4 columns');
    return
elseif size(reaction_rates,1) == 1
    reaction_rates = repmat(reaction_rates,num_species,1);
elseif size(reaction_rates,1) ~=num_species
    disp('error: reactions rates needs either one row or num_species rows');
    return 
end

growth_rate = reaction_rates(:,1);
aggregation_rate = reaction_rates(:,2);
expulsion_rate = reaction_rates(:,3);
fragmentation_rate = reaction_rates(:,4);

% total simulation time
if ~exist('Tmax','var')||isempty(Tmax)
    Tmax = 24;
end

% initial starting population of single cells
if ~exist('n0','var')||isempty(n0)
    n0 = 10;
end

% tau 
if ~exist('tau','var')||isempty(tau)
    tau = .005;
end

% carrying capacities
if ~exist('carrying_capacities','var')||isempty(carrying_capacities)
    carrying_capacities = 1e4.*ones(num_species,1);
elseif numel(carrying_capacities) == 1
    carrying_capacities = repmat(carrying_capacities,num_species,1);
end

% interaction_matrix 
if ~exist('interaction_matrix','var')||isempty(interaction_matrix)
    interaction_matrix = -1.*ones(num_species);
end

% shuffle random number generator
rng('shuffle');


%% intialize some arrays

% main object of the simulation:  array of cluster volumes
if numel(n0)==1
    cluster_sizes = zeros(num_species,round(num_species.*n0));
    for k = 1:num_species
        cluster_sizes(k,[(k-1).*round(n0) + 1]:[k.*round(n0)]) = 1;
    end
    
    % array for keeping track of number of clumps over time
    total_num_clumps_arr = round(num_species.*n0);

elseif size(n0,1) ~= num_species
    disp('gac error:  initial cluster size array n0 must have num_species rows');
    return
else
    cluster_sizes = n0;
    
     % array for keeping track of number of clumps over time
    total_num_clumps_arr = size(n0,2);
    
end


% array for keeping track of total population size
abundance_arr = sum(cluster_sizes,2);

% pregenerate lots of random numbers for speed.
number_of_random_numbers_to_pre_generate = round(1e6);
lots_of_random_numbers = rand(1,number_of_random_numbers_to_pre_generate);
max_allowed_number_of_random_numbers = 1e8;

% counter for how many of the pre-generated random numbers have been used.
n = 0;

% time
t = 0;
tvec = t;

% hour marker
disp_time_marker = 1;
disp_time_increment = 1;


%% fixed parameters

% growth timestep
fraction_of_delta_t = 1;
baseline_dt = .1;

% logical for printing progress to screen
l_print_progress = true;


%% main simulation.  loop over time.
while t < Tmax 
            
    % print update to console
    if l_print_progress && t >= disp_time_marker
        disp(['time = ' num2str(t,2) ' total number of clusters = ' num2str(size(cluster_sizes,2))])
        disp_time_marker = disp_time_marker + disp_time_increment;
    end   
          
   
        
   %% growth (deterministic) - do this first
   
   % construct a timestep for numerical integration.  can change
   % resolution here for better or worse accuracy.  Base timestep is
   % .01, if time between reactions gets small, use a smaller step.
   %dt = min(delta_t./10,.01);
   dt = min(tau.*fraction_of_delta_t,baseline_dt);
   
   % number of steps taken between reactions
   %numsteps = round((delta_t - dt)/dt);
   numsteps = ceil(tau/dt);
   
   
   % update time array.  new variable name is needed to avoid
   % referencing issues with the nested function.
   tvec = [tvec, linspace(t+dt, t+tau, numsteps)];
   
   % temporary arrays to be updated during numerical integration
   tmp_tot_pop = zeros(num_species,numsteps);
   tmp_num_clumps = zeros(1,numsteps);
   
   % growth rate matrix
   growth_rate_mat = repmat(growth_rate,1,size(cluster_sizes,2));
   
   % carrying capacity matrix
   carrying_capacities_mat = repmat(carrying_capacities,1,size(cluster_sizes,2));
   
   % loop over time and update according to growth equation
   for s = 1:numsteps
       
       % grow all clusters in one deterministic vectorized step
       cluster_sizes = cluster_sizes + dt.*growth_rate_mat./carrying_capacities_mat.*cluster_sizes.*(carrying_capacities_mat+interaction_matrix*sum(cluster_sizes,2));
       
       % if clusters die out, remove them
       cluster_sizes(cluster_sizes<=1) = 0;
       
       % update tmp_totpop
       tmp_tot_pop(:,s) = sum(cluster_sizes,2);
       
       % update tmp_num_clumps
       tmp_num_clumps(s) = size(cluster_sizes,2);
       
   end
   
   % append total population array
   abundance_arr = [abundance_arr, tmp_tot_pop];
   
   % append total number of clumps array
   total_num_clumps_arr = [total_num_clumps_arr, tmp_num_clumps];
    
    
   % update time
   t = t + tau;
    %% fragmentation
    % do by each cluster
    num_frag_events_for_each_cluster = poissrnd(mean(fragmentation_rate)*tau,1,size(cluster_sizes,2));
    
    if ~isempty(num_frag_events_for_each_cluster)
        for j = 1:size(cluster_sizes,2)
            this_num_frag_events = num_frag_events_for_each_cluster(j);
            
            if this_num_frag_events > 0
                this_random_species_partition = randi(num_species,this_num_frag_events,1);
                
                [frag_events_by_species,~] = hist(this_random_species_partition,1:num_species);
         
                frag_events_by_species = frag_events_by_species';
                
                frag_events_by_species([cluster_sizes(:,j) - frag_events_by_species] < 1) = 0;
                
                cluster_sizes(:,j) = cluster_sizes(:,j) - frag_events_by_species;
                
                for k = 1:num_species
                    new_singletons = zeros(num_species,frag_events_by_species(k));
                    new_singletons(k,:) = 1;
                    cluster_sizes = [cluster_sizes,new_singletons];
                end
                
                total_num_clumps_arr(end) = size(cluster_sizes,2);
            end
        end
    end
    
    %% aggregation
    % do all in one
    number_of_possible_agg_reactions = round(.5.*(size(cluster_sizes,2).^2 - size(cluster_sizes,2)));
    total_prob_rate_of_an_agg_event_happening  = mean(aggregation_rate).*number_of_possible_agg_reactions;
    
    num_agg_events_to_happen = poissrnd(total_prob_rate_of_an_agg_event_happening*tau);
    if num_agg_events_to_happen > 0
        [random_numbers, n, lots_of_random_numbers] = select_x_random_numbers(2*num_agg_events_to_happen);
        
        agg_ids = reshape(ceil(random_numbers.*size(cluster_sizes,2)),num_agg_events_to_happen,2);
        agg_ids(diff(agg_ids,[],2)==0,:) = [];
        [~,unique_rows,~] = unique(agg_ids(:,2));
        agg_ids = agg_ids(unique_rows,:);
        agg_ids(ismember(agg_ids(:,2),agg_ids(:,1)),:) = [];
        
        for a = 1:size(agg_ids,1)
            cluster_sizes(:,agg_ids(a,1)) = cluster_sizes(:,agg_ids(a,1)) + cluster_sizes(:,agg_ids(a,2));
        end
        
        cluster_sizes(:,agg_ids(:,2)) = [];
        
        
        total_num_clumps_arr(end) = size(cluster_sizes,2);
    end
    
    %% expulsion
    % 
    total_prob_rate_of_explusion_event_happening = mean(expulsion_rate)*size(cluster_sizes,2);
    num_expulsion_events_to_happen = poissrnd(total_prob_rate_of_explusion_event_happening*tau);
    
    if num_expulsion_events_to_happen > 0
        [random_numbers, n, lots_of_random_numbers] = select_x_random_numbers(num_expulsion_events_to_happen);
        
        expulsion_ids = unique(ceil(random_numbers.*size(cluster_sizes,2)));
        
        cluster_sizes(:,expulsion_ids) = [];
        
        
        total_num_clumps_arr(end) = size(cluster_sizes,2);
        abundance_arr(:,end) = sum(cluster_sizes,2);
    end
    
end
    
    
     function [random_numbers, n_out, lots_of_random_numbers_out] = select_x_random_numbers(x)
         
       
         
        if n + x < number_of_random_numbers_to_pre_generate
            random_numbers = lots_of_random_numbers(n+1:n+x);
            n_out = n + x;
            lots_of_random_numbers_out = lots_of_random_numbers;
        else
            disp('gac: ran out of pre-generated rns, generating new ones')

            % generate new ones
            if x < number_of_random_numbers_to_pre_generate
                lots_of_random_numbers_out = rand(1,number_of_random_numbers_to_pre_generate);
            elseif x < max_allowed_number_of_random_numbers
                lots_of_random_numbers_out = rand(1,x);
            else
                disp(['error: number of random numbers exceeds memory limit '])
                return
            end
            
            % collect the one we need
            random_numbers = lots_of_random_numbers_out(1:x);
            
            % reset counter
            n_out = x;
            
        end
        
    end
    
    
    
  

   
        
end





    




