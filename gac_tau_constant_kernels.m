% Program:  gac_tau_constant_kernels.m

function [cluster_sizes,total_pop_arr,tvec,num_clumps_arr] = gac_tau_constant_kernels(growth_rate,aggregation_rate,expulsion_rate,fragmentation_rate,Tmax,n0,max_total_pop)

%% default values for input parameters

% rate of cluster growth
if ~exist('growth_rate','var')||isempty(growth_rate)
    growth_rate = 1;
end

% rate of cluster aggregation
if ~exist('aggregation_rate','var')||isempty(aggregation_rate)
    aggregation_rate = 1;
end

% rate of cluster explusion
if ~exist('expulsion_rate','var')||isempty(expulsion_rate)
    expulsion_rate = .1;%.005;
end

% rate of cluster fragmentation_rate
if ~exist('fragmentation_rate','var')||isempty(fragmentation_rate)
    fragmentation_rate = 20;
end

% total simulation time
if ~exist('Tmax','var')||isempty(Tmax)
    Tmax = 72;
end

% initial starting population of single cells
if ~exist('n0','var')||isempty(n0)
    n0 = 10;
end


% carrying capacity
if ~exist('max_total_pop','var')||isempty(max_total_pop)
    max_total_pop = 1e5;
end


% shuffle random number generator
rng('shuffle');


%% intialize some arrays

% main object of the simulation:  array of cluster volumes
if numel(n0)==1
    cluster_sizes = ones(1,round(n0));
    
    % array for keeping track of number of clumps over time
    num_clumps_arr = n0;

elseif size(n0,1) > 1
    disp('gac error:  initial cluster size array n0 must have dim 1xM');
    return
else
    cluster_sizes = n0;
    
     % array for keeping track of number of clumps over time
    num_clumps_arr = numel(n0);
    
end


% array for keeping track of total population size
total_pop_arr = sum(n0);

% pregenerate lots of random numbers for speed.
number_of_random_numbers_to_pre_generate = round(1e7);
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
% tau -- fixed for now
tau = .05;
                   
% growth timestep
fraction_of_delta_t = 1;
baseline_dt = .1;

% logical for printing progress to screen
l_print_progress = false;


%% main simulation.  loop over time.
while t < Tmax 
            
    % print update to console
    if l_print_progress && t >= disp_time_marker
        disp(['time = ' num2str(t,2) ' number of clusters = ' num2str(numel(cluster_sizes))])
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
   tmp_tot_pop = zeros(1,numsteps);
   tmp_num_clumps = zeros(1,numsteps);
   
   % loop over time and update according to growth equation
   for s = 1:numsteps
       
       % grow all clusters in one deterministic vectorized step
       cluster_sizes = cluster_sizes + dt.*growth_rate.*cluster_sizes.*(1-sum(cluster_sizes)./max_total_pop);
       
       % if clusters die out, remove them
       %cluster_sizes_out = cluster_sizes_out(cluster_sizes_out>=1);
       
       % update tmp_totpop
       tmp_tot_pop(s) = sum(cluster_sizes);
       
       % update tmp_num_clumps
       tmp_num_clumps(s) = numel(cluster_sizes);
       
   end
   
   % append total population array
   total_pop_arr = [total_pop_arr, tmp_tot_pop];
   
   % append total number of clumps array
   num_clumps_arr = [num_clumps_arr, tmp_num_clumps];
    
    
   % update time
   t = t + tau;
    %% fragmentation
    % do by each cluster
    num_frag_events_for_each_cluster = poissrnd(fragmentation_rate*tau,1,numel(cluster_sizes));
    
    if sum(num_frag_events_for_each_cluster) > 0
        num_frag_events_for_each_cluster([cluster_sizes - num_frag_events_for_each_cluster] < 1 ) = 0;
        
        cluster_sizes = cluster_sizes - num_frag_events_for_each_cluster;
        cluster_sizes = [cluster_sizes, ones(1,sum(num_frag_events_for_each_cluster))];
        
        num_clumps_arr(end) = numel(cluster_sizes);
    end
    
    %% aggregation
    % do all in one
    number_of_possible_agg_reactions = round(.5.*(numel(cluster_sizes).^2 - numel(cluster_sizes)));
    total_prob_rate_of_an_agg_event_happening  = aggregation_rate.*number_of_possible_agg_reactions;
    
    num_agg_events_to_happen = poissrnd(total_prob_rate_of_an_agg_event_happening*tau);
    if num_agg_events_to_happen > 0
        [random_numbers, n, lots_of_random_numbers] = select_x_random_numbers(2*num_agg_events_to_happen);
        
        agg_ids = reshape(ceil(random_numbers.*numel(cluster_sizes)),num_agg_events_to_happen,2);
        agg_ids(diff(agg_ids,[],2)==0,:) = [];
        [~,unique_rows,~] = unique(agg_ids(:,2));
        agg_ids = agg_ids(unique_rows,:);
        agg_ids(ismember(agg_ids(:,2),agg_ids(:,1)),:) = [];
        
        for a = 1:size(agg_ids,1)
            cluster_sizes(agg_ids(a,1)) = cluster_sizes(agg_ids(a,1)) + cluster_sizes(agg_ids(a,2));
        end
        
        cluster_sizes(agg_ids(:,2)) = [];
        
        
        num_clumps_arr(end) = numel(cluster_sizes);
    end
    
    %% expulsion
    % 
    total_prob_rate_of_explusion_event_happening = expulsion_rate*numel(cluster_sizes);
    num_expulsion_events_to_happen = poissrnd(total_prob_rate_of_explusion_event_happening*tau);
    
    if num_expulsion_events_to_happen > 0
        [random_numbers, n, lots_of_random_numbers] = select_x_random_numbers(num_expulsion_events_to_happen);
        
        expulsion_ids = unique(ceil(random_numbers.*numel(cluster_sizes)));
        
        cluster_sizes(expulsion_ids) = [];
        
        
        num_clumps_arr(end) = numel(cluster_sizes);
        total_pop_arr(end) = sum(cluster_sizes);
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




%% scrap
%  % compute probability of expulsion as a function of cluster volume.
%     % here, the relationship is linear in linear length dimension.
%     lambda_E = expulsion_rate.*ones(1,numel(cluster_sizes));
%     
%     % create an array of ids labelling every possible explusion reaction
%     ids_for_expulsion = 1:numel(lambda_E);
%     
%     % create an accesory array that labels (1) these reactions as "expulsion"
%     label_for_lambda_E = 1.*ones(1,numel(lambda_E));
%     
%     %aggregation rates
%     number_of_possible_agg_reactions = round(.5.*(numel(cluster_sizes).^2 - numel(cluster_sizes)));
%     lambda_A = aggregation_rate.*ones(1,number_of_possible_agg_reactions);
%    
%      % assign ids to each of the possible aggregation reactions.
%     ids_for_agg = 1:numel(lambda_A);
%     
%     % create labels (2) that keep track of these possible reactions as
%     % "aggregation"
%     label_for_lambda_A = 2.*ones(1,numel(lambda_A));
%         
%     % compute probability of fragmentation.  require clusters to be of at least
%     % size 2 to fragment. contains zeros.
%     lambda_F = fragmentation_rate.*(cluster_sizes >= 2).*ones(1,numel(cluster_sizes));
%     
%      % create labels (3) that keep track of these possible reactions as
%     % "fragmentation" 
%     label_for_lambda_F = 3.*ones(1,numel(lambda_F));
%     
%     % assign ids to each of the possible fragmentation reactions.
%     ids_for_frag = 1:numel(lambda_F);
%     
%     % assemble all of the reaction ids into a single linear array
%     ids_arr = [ids_for_expulsion, ids_for_agg, ids_for_frag];
%     
%     % assemble all reaction probability rates into a single linear array
%     lambda_arr = [lambda_E, lambda_A, lambda_F];
%     
%     % compute total proabability rate of a reaction happening
%     lambda_total = sum(lambda_arr);
%     
%     % assemble all labels into a single linear array
%     all_labels = [label_for_lambda_E, label_for_lambda_A, label_for_lambda_F]; 
%     
%     number_of_aggregation_events = poissrnd(sum(lambda_A)*tau);
%     
%     number_of_fragmentation_events = min(floor(sum(cluster_sizes(cluster_sizes > 2))),poissrnd(sum(lambda_F)*tau));
% 
%     
%     
%     
%     
%     
    
    
    
    
    
    




