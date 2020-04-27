% Program:  gac_tau_constant_kernels.m

function [cluster_sizes,total_pop_arr,tvec,num_clumps_arr] = gac_tau_constant_kernels(growth_rate,aggregation_rate,expulsion_rate,fragmentation_rate,Tmax,n0,max_total_pop,tau,fragmentation_exponent)

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

% tau
if ~exist('tau','var')||isempty(tau)
    tau = 0.001;
end

% fragmentation_exponent
if ~exist('fragmentation_exponent','var')||isempty(fragmentation_exponent)
    fragmentation_exponent = 0;
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

% tau for timestep--fix for now
%tau = .001;

% growth timestep
fraction_of_delta_t = 1;
baseline_dt = .1;

% logical for printing progress to screen
l_print_progress = true;

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
   
   
   new_tvec = linspace(t+dt, t+tau, numsteps);
   
   tmp_tot_pop = exp(growth_rate.*(new_tvec(end)-t))./((1/total_pop_arr(end)) + (1/max_total_pop).*(exp(growth_rate.*(new_tvec(end)-t)) - 1));
   tmp_num_clumps = num_clumps_arr(end).*ones(1,numsteps);
   
   cluster_sizes = exp(growth_rate.*(new_tvec(end)-t))./((1./cluster_sizes) + (total_pop_arr(end)./max_total_pop./cluster_sizes).*(exp(growth_rate.*(new_tvec(end)-t)) - 1));

   % append time array
   tvec = [tvec, new_tvec];
   
   % append total population array
   total_pop_arr = [total_pop_arr, tmp_tot_pop];
   
   % append total number of clumps array
   num_clumps_arr = [num_clumps_arr, tmp_num_clumps];
     
   % update time
   t = t + tau;
   
   %% fragmentation
   % do by each cluster
   total_fragmentation_rate = fragmentation_rate.*(cluster_sizes.*(1-sum(cluster_sizes)./max_total_pop)).^fragmentation_exponent;
   %%%%%%%%%%%%%%% old way
   num_frag_events_for_each_cluster = poissrnd(total_fragmentation_rate.*tau,1,numel(cluster_sizes));
   
   %%%%%%%%%%%%% new way
%    number_of_rns_per_cluster = max(ceil(10*max(total_fragmentation_rate.*tau)./0.7),2);
%    [random_numbers, n, lots_of_random_numbers] = select_x_random_numbers(numel(cluster_sizes)*number_of_rns_per_cluster);
%    random_numbers = reshape(random_numbers,numel(cluster_sizes),number_of_rns_per_cluster);
%    [num_frag_events_for_each_cluster,n,lots_of_random_numbers] = convert_uniform_rns_to_poisson(random_numbers,[total_fragmentation_rate.*tau]');
%     num_frag_events_for_each_cluster = num_frag_events_for_each_cluster';
    %%%%%%%%%%%%%%%%
    
    num_frag_events_for_each_cluster([cluster_sizes - num_frag_events_for_each_cluster] < 1 ) = 0;

    if sum(num_frag_events_for_each_cluster) > 0
        
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
    
%%%%%%%%%%%%% not in use yet
    function [r,n_out,lots_of_random_numbers_out] = convert_uniform_rns_to_poisson(u,lambda)
        l_success = 0;
        counter = zeros(size(u,1),1);
        number_of_rns_to_select_if_needed = size(u,2);%max(ceil(10.*max(lambda)./0.7),2);
        
        if numel(lambda)==1 && size(u,1)>1
            lambda = repmat(lambda,size(u,1),1);
        end
        
        while l_success==0
            tmp = cumprod(u,2) < exp(-lambda);

            [~,r] = max(tmp,[],2);
            r = r-1;
            r(r==0 & tmp(:,1)==0)= NaN;
            
            if sum(isnan(r)) > 0
                [random_numbers, n_out, lots_of_random_numbers_out] = select_x_random_numbers(number_of_rns_to_select_if_needed*size(u,1));
                random_numbers = reshape(random_numbers,size(u,1),number_of_rns_to_select_if_needed);
                u =  repmat(~isnan(r),1,size(u,2)).*u + repmat(isnan(r),1,size(u,2)).*random_numbers;
                counter = counter + isnan(r).*size(u,2);
            else
                l_success=1;
                r = r+counter;
                n_out = n;
                lots_of_random_numbers_out = lots_of_random_numbers;
            end
            
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




 
    
    
    
    
    
    




