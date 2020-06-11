%Copyright (c) 2020 
%3D Systems Packaging Research Center (PRC), Georgia Tech.

% Two-Stage Bayesian Optimization (TSBO) for maximizing a black-box
% function. If you want minimization, multiply objective function by "-1".

%This material is based on work supported by NSF I/UCRC Center for Advanced
%Electronics Through Machine Learning (CAEML)
%For questions and queries, please contact: htorun3@gatech.edu

%Please cite our paper if you use the code:
%H. M. Torun, M. Swaminathan, A. Kavungal Davis and M. L. F. Bellaredj,
%"A Global Bayesian Optimization Algorithm and Its Application to Integrated System Design,"
%in IEEE Transactions on Very Large Scale Integration (VLSI) Systems, vol. 26, no. 4, pp. 792-802, April 2018.
%%

clear all
close all
clc

addpath('HFSS/')
addpath('test_function/')
Rosen2 = @(value)  -1*( ((1-value(:,1)).^2)+(100*((value(:,2)-(value(:,1).^2)).^2)) );
Branin = @(value) -1*braninf(value);
shekel5 = @(value) -1*shekelMf(value,5);
Hart6 = @(value) -1*hart6f(value);
Hart3 = @(value) -1*hart3f(value);
HFSS_inductor = @(value) EMBEDDED_INDUCTOR(value(1,1),value(1,2),value(1,3),value(1,4),value(1,5),value(1,6),value(1,7),value(1,8),value(1,9),value(1,10));

%% Select which function to be optimized. Set Dimensions and Sample Space
f = Branin; dimension=2; sample_space = [-5 10; 0 15]; optima =  -0.397887357752662; 
% f = Rosen2; dimension=2; sample_space = [-5 10; -5 10];optima =  0;
% f = shekel5; dimension=4; sample_space = [0 10; 0 10; 0 10; 0 10]; optima= 10.153195850979039;
% f = Hart6; dimension=6; sample_space = [0 1; 0 1; 0 1; 0 1; 0 1; 0 1]; optima = 3.322368011391339;
% f = Hart3; dimension=3; sample_space = [0 1; 0 1; 0 1]; optima= 3.862782147819745;
% f = HFSS_inductor; dimension = 10; sample_space = [53, 170;   50, 650; 35 170  ;50, 350  ; 2, 20  ; 3, 13  ; 2, 30  ;  0, 100 ;   0.1 1;  35 170];

%% Initilization - Set Parameters of TSBO
%Switching count: "n" in Eq. 16 in TVLSI paper. This defines the amount of
%time you want to spend before switching to 2nd stage. We used this as "10"
%in the paper for challenge functions. But it should not be less than 20 for
%actual design optimization problems just to ensure proper exploration,
% Higher the number, less chance of stucking in local optima, but more time
% required to fine tune.
N_switching_count = 25; 

count_max = 100; %maximum number of function evaluations(simulations)
N_acq_train = 50; %Number of simulations to be used for Learning Acquisition functions block

%Hyperparameters of EI and PI can be used in their default. 
sigma1 = 0.1; %hyperparameter zeta for Expected improvement in Eq. 9 in TVLSI paper
sigma2 = 0.01; %hyperparameter zeta for Probability of Improvement in Eq.8 in TVLSI paper
alpha_tightness = 0.15; %Tightness of region selected for fine tuning, Eq.17 in TVLSI paper
%% Starting Point

% Center of sample space is used for starting point (xs)
% You can provide additional points here by
% adding more sampling points to the variable "total_samples"

regions_total = splitregion(sample_space,dimension);
xs = ((sample_space(:,1)+sample_space(:,2))/2)';

total_samples = xs;
for i = 1:2^dimension
    test_x(i,:) = (regions_total(:,1,i) + regions_total(:,2,i))/2;
end

total_targets = [];
for i = 1:size(total_samples,1)
    total_targets = [total_targets;f(total_samples(i,:))];
end
% initialize partitioning tree
max_depth = count_max;
t = cell(max_depth,1);
for a = 1:max_depth
    t{a}.regions = [];
    t{a}.sum_of_variances = [];
    t{a}.test_x = [];
    t{a}.variances = [];
end
%%

gains = zeros(3,1);
best_samples = [];
best_targets =[];
regions_total = splitregion(sample_space,dimension);
new_region = regions_total;

t{1,1}.regions = new_region;
switching_count = 0;
stuckk = false;

M = 1;
nu = 0.1;
GP_varsigma = @(M) sqrt(2*log(pi^2*(M)^2/(12*nu)));
new_region_index_array = [];
previous_best_sample = 0;
previous_best_target = Inf;
ucb_count = 0; ei_count = 0; pi_count = 0;

tic
count = 0;
while (count <= count_max)
    
    count = count + 1;
    % Number of new regions to be generated
    if(count == 1)
        n_regions = 2^dimension;
    else
        n_regions = n_future_region;
    end
    n_future_region = n_regions - 1 + 2^dimension;
    
    % Train a GP Model using all the samples
    % Here, Matlab's own version is used for robustness to numerical errors
    % In BALDO code (additionally provided), GPML toolbox is used
    % for full-customization capability.
    
    % You can change kernel function to be better suited to the problem at
    % hand. Here, "Matern 5/2 with Automatic Relevance Determination (ARD)"
    % is used. All the results in the paper is with this configuration.
    gprMdl = fitrgp(total_samples,total_targets,'Basis','none','Sigma',0.00001,'ConstantSigma',true,'KernelFunction','ardmatern52','FitMethod','exact','PredictMethod','exact');
    
    % Learning Acquisition Functions: Calculate gains of each acquisition function.
    if(count > 1)
        gains(1) = gains(1) + abs(predict(gprMdl,test_UCB));
        gains(2) = gains(2) + abs(predict(gprMdl,test_EI));
        gains(3) = gains(3) + abs(predict(gprMdl,test_PI));
    end
    
    % Get previous best observed value so far. This is used to determine if
    % a switching to 2nd stage should occur.
    if(count>1)
        previous_best_target = max_of_targets;
        previous_best_sample = max_sample;
    end
    %Get current best value and corresponding input parameters.
    [max_of_targets,max_target_index] = max(total_targets);
    max_sample = total_samples(max_target_index,:);
    
    % predict test points using GP
    [gp_output,sample_std]  = predict(gprMdl,test_x);
    
    %calculate EI, UCB and PI
    Z1 = (gp_output - max_of_targets-sigma1)./(sample_std);
    EI = (-max_of_targets+gp_output-sigma1).*normcdf(Z1) + (sample_std).*normpdf(Z1);
    
    UCB = gp_output + GP_varsigma(M)*(sample_std);
    M = M+1;
    
    PI = normcdf((gp_output - max_of_targets-sigma2)./sample_std);
    
    % Select UCB, EI or PI sequentally for selecting next input parameters
    % to be simulated.
    if(rem(count,3) == 2)
        eta = UCB;
        if (count < N_acq_train)
            ucb_count = ucb_count + 1;
        end
    elseif(rem(count,3) == 0)
        eta = EI;
        if (count < N_acq_train)
            ei_count = ei_count + 1;
        end
    elseif(rem(count,3) == 1)
        eta = PI;
        if (count < N_acq_train)
            pi_count = pi_count +1;
        end
        % else
    end
    
    % Get suggested samples from each acquisiton function for updating
    % gains in Learning Acquisition Functions
    [u_x,u_i] = max(UCB);
    test_UCB = test_x(u_i,:);
    [e_x,e_i] = max(EI);
    test_EI = test_x(e_i,:);
    [p_x,p_i] = max(PI);
    test_PI = test_x(p_i,:);

    % Select best acquisition function after N_acq_train simulations
    if(count > N_acq_train)
        [a,sel] = max(gains);
        if(sel == 1)
            eta = UCB;
            ucb_count = ucb_count + 1;
        elseif sel == 2
            eta = EI;
            ei_count = ei_count + 1;
        else
            eta = PI;
            pi_count = pi_count +1;
        end
    end
    
    % Get new sample to be simulated, "sample_new", and the region it
    % belongs to.
    [~,new_region_index] = max(eta);
    sample_new = test_x(new_region_index,:);
    
    target_new = f(sample_new);
    
    % If switching occurs, get the "small enough region" to be further
    % exploited.
    if(stuckk)
        stuckk = false;
        new_region_temp = t{count-max_target_index+1}.regions(:,:,new_region_index_array(count-max_target_index+1));
    else
        new_region_temp = t{count,1}.regions(:,:,new_region_index);
    end
    % Divide selected region by 2^D new regions and generate candidate
    % points.
    new_region = splitregion(new_region_temp,dimension);
    
    test_temp = squeeze((new_region(:,1,:)+new_region(:,2,:))/2)';
 
    temp_x = test_x;
    test_x = zeros(n_future_region,dimension);
    if(new_region_index == 1)
        test_x(1:2^dimension,:) = test_temp;
        test_x(2^dimension + 1:n_future_region,:) = temp_x(2:n_regions,:);
    else
        test_x(1:new_region_index-1,:) = temp_x(1:new_region_index-1,:);
        test_x(new_region_index:new_region_index+2^dimension-1,:) = test_temp;
        test_x(new_region_index+2^dimension:n_future_region,:) = temp_x(new_region_index+1:n_regions,:);
    end
    
    t{count+1,1}.regions = zeros(dimension,2,n_future_region);    
    if(new_region_index == 1)
        t{count+1,1}.regions(:,:,1:2^dimension) = new_region;
        t{count+1,1}.regions(:,:,2^dimension + 1:n_future_region) = t{count,1}.regions(:,:,2:n_regions);
    else
        t{count+1,1}.regions(:,:,1:new_region_index-1) =  t{count,1}.regions(:,:,1:new_region_index-1);
        t{count+1,1}.regions(:,:,new_region_index:new_region_index+2^dimension-1) = new_region;
        t{count+1,1}.regions(:,:,new_region_index+2^dimension :n_future_region) = t{count,1}.regions(:,:,new_region_index+1:n_regions);
    end
    % Clear previously stored regions for memory saving, can comment out
    % for debugging purposes.
    t{count,1}.regions = [];
   
    % Report outputs
    clc
    target_new
    max_of_targets
    switching_count
    count
    
    % Check if the algorithm is "stuck", i.e. switching count in paper.
    % If stuck for "N_switching_count" times, switch to 2nd stage.
    if ((norm(previous_best_sample-max_sample) <= 0.01*norm(max_sample)))
        switching_count = switching_count +1;
        if(switching_count == N_switching_count)
            count_first = count;
            break
        end
    else
        switching_count = 0;
    end
    
    best_targets = [max_of_targets;best_targets];
    best_samples = [max_sample;best_samples];
    
    if(count == count_max)
        toc        
        break;
    else
        total_samples = [sample_new;total_samples];
        total_targets = [target_new;total_targets];
    end

end
%% 2nd Stage: Pure Exploitation
% This might not be necessary if you dont want a fine tuning of your
% parameters. If you want to skip fine tuning, change "N_switch_count = count_max".

total_samples = [sample_new;total_samples];
total_targets = [target_new;total_targets];
[max_of_targets,max_target_index] = max(total_targets);
max_sample = total_samples(max_target_index,:);

%Get the small enough region to be explotied. See Eq. 17 in TVLSI paper
converged = sample_space;
for i = 1:dimension
    converged(i,:) = [max([sample_space(i,1),(1-alpha_tightness)*max_sample(:,i)]) min([sample_space(i,2),(1+alpha_tightness)*max_sample(:,i)])];
end


converged_region = splitregion33(converged,dimension);
t{count+1,1}.regions = converged_region;
total_samples2 = max_sample;
total_targets2 = max_of_targets;
n_future_region = 3;

M = size(total_targets2,1);

% Uncomment below if you want to see how many times each acquisition function is
% called in the first stage.


%Below is very similar to 1st stage, except instead of 2^D divisions, 3
%divisions occur at each iteration. For better fine tuning, set "size =
%1000" as this is the size of candidate points generated in each region.
%Recommended is "sizee = 500" for a balanced time trade-off.

selection = [zeros(1,dimension)];
sizee = 500;
test_x = zeros(sizee+1,dimension,n_future_region);

for i = 1:n_future_region
    regions = converged_region;
    for a = 1:size(selection,1)
        test_x(a,:,i) = 0.5*(regions(:,1,i)'+regions(:,2,i)') + 0.5*(regions(:,2,i)'-regions(:,1,i)').*selection(a,:);
    end
    test_x(2:sizee+1,:,i) = regions(:,1,i)' + (regions(:,2,i)'-regions(:,1,i)').*rand(sizee,dimension);
end


while(count <= count_max)
    count = count + 1;
    n_regions = n_future_region;
    n_future_region = n_regions - 1 + 3;
    gprMdl = fitrgp(total_samples2,total_targets2,'Basis','none','Sigma',0.00001,'ConstantSigma',true,'KernelFunction','ardmatern52','FitMethod','exact','PredictMethod','exact');
    if(count > 1)
        gains(1) = gains(1) + abs(predict(gprMdl,test_UCB));
        gains(2) = gains(2) + abs(predict(gprMdl,test_EI));
        gains(3) = gains(3) + abs(predict(gprMdl,test_PI));
    end
    if(count > 1)
        [max_of_targets,max_target_index] = max(total_targets2);
        max_sample = total_samples2(max_target_index,:);
    end

    gp_output = zeros(sizee+1,1,n_regions); sample_std = gp_output;
    for i = 1:n_regions
        [gp_output(:,:,i),sample_std(:,:,i)]  = predict(gprMdl,test_x(:,:,i));
    end

    Z1 = zeros(sizee+1,1,n_regions);
    Z2 = Z1; acquisition1 = Z1;
    for i = 1:n_regions
        Z1(:,:,i) = (gp_output(:,:,i) - max_of_targets-sigma1)./(sample_std(:,:,i));
        Z2(:,:,i) = (gp_output(:,:,i) - max_of_targets-sigma2)./(sample_std(:,:,i));
        acquisition1(:,:,i) = (-max_of_targets+gp_output(:,:,i)-sigma1).*normcdf(Z1(:,:,i)) + (sample_std(:,:,i)).*normpdf(Z1(:,:,i));
    end
    
    EI = acquisition1;
    PI = Z2;
    UCB = gp_output + GP_varsigma(M)*(sample_std);
    M = M + 1;
    
    if(count <= N_acq_train)
        if(rem(count,3) == 0)
            eta = UCB;
            ucb_count = ucb_count + 1;
        elseif(rem(count,3) == 1)
            eta = EI;
            ei_count = ei_count + 1;
        elseif(rem(count,3) == 2)
            eta = PI;
            pi_count = pi_count +1;
        end
    end
    [~,u_i] = max(UCB); [~,u_i2] = max(max(UCB));
    test_UCB = test_x(u_i(u_i2),:,u_i2);
    [~,e_i] = max(EI); [~,e_i2] = max(max(EI));
    test_EI = test_x(e_i(e_i2),:,e_i2);
    [~,p_i] = max(PI); [~,p_i2] = max(max(PI));
    test_PI = test_x(p_i(p_i2),:,p_i2);
    
    if(count > 50)
        [~,sel] = max(gains);
        if(sel == 1)
            eta = UCB;
            ucb_count = ucb_count + 1;
        elseif (sel == 2)
            eta = EI;
            ei_count = ei_count + 1;
            
        else
            eta = PI;
            pi_count = pi_count +1;
        end
    end
    
    % You can just use UCB for 2nd stage as we are just exploiting,
    % Comment out "eta=UCB" below if you selected alpha_tightness > 0.2 as there still
    % can be benefit from Learning Acquisition Functions 
    eta = UCB;
    [in_region_max,in_region_index] = max(eta);
    [max_UCB,new_region_index] = max(max(eta));
    sample_new = test_x(in_region_index(new_region_index),:,new_region_index);
    target_new = f(sample_new);
    
    new_region_temp = t{count,1}.regions(:,:,new_region_index);
    new_region = splitregion33(new_region_temp,dimension);
    

    t{count+1,1}.regions = zeros(dimension,2,n_future_region);
    if(new_region_index == 1)
        t{count+1,1}.regions(:,:,1:3) = new_region;
        t{count+1,1}.regions(:,:,3 + 1:n_future_region) = t{count,1}.regions(:,:,2:n_regions);
    else
        t{count+1,1}.regions(:,:,1:new_region_index-1) =  t{count,1}.regions(:,:,1:new_region_index-1);
        t{count+1,1}.regions(:,:,new_region_index:new_region_index+3-1) = new_region;
        t{count+1,1}.regions(:,:,new_region_index+3 :n_future_region) = t{count,1}.regions(:,:,new_region_index+1:n_regions);
    end
    t{count,1}.regions = [];
    
    test_x = zeros(sizee+1,dimension,n_future_region);
    
    selection = [zeros(1,dimension)];
    for i = 1:n_future_region
        regions = t{count+1,1}.regions;
        for a = 1:size(selection,1)
            test_x(a,:,i) = 0.5*(regions(:,1,i)'+regions(:,2,i)') + 0.5*(regions(:,2,i)'-regions(:,1,i)').*selection(a,:);
        end
        test_x(2:sizee+1,:,i) = regions(:,1,i)' + (regions(:,2,i)'-regions(:,1,i)').*rand(sizee,dimension);
    end
    best_targets = [max_of_targets;best_targets];
    best_samples = [max_sample;best_samples];
    
    clc
    target_new
    max_of_targets
    max_sample
    count
    
    
    total_samples2 = [sample_new;total_samples2];
    total_targets2 = [target_new;total_targets2];
end


toc
%%
fid = figure;
h=gca;
plot(1:1:length(best_targets),flip(best_targets),'-b','LineWidth',3)
% h.FontWeight = 'bold';
h.FontSize = 24;
h.XGrid = 'on'; %h.XMinorGrid = 'on';
h.YGrid = 'on'; %h.YMinorGrid = 'on';
xlabel('Number of Function Evaluations');
ylabel('Best Observed Value');

set(fid, 'Units', 'inches'); % set units
PaperWidth = 9.5; % paper width
PaperHeight = 7.5; % paper height
set(fid, 'PaperSize', [PaperWidth PaperHeight]); % set paper size
afFigurePosition = [1 1 PaperWidth PaperHeight]; % set figure on screen [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition); % set figure position on paper [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto'); %
set(gcf, 'Renderer', 'painters');
set(gca, 'Units','normalized','Position',[0.125 0.13 0.845 0.85]); % fit axes within figure