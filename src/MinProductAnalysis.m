%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Senior Design - Min and Product Processing Analysis
% Louisiana Tech University
% Pablo Johnson, Daniel Sartori, Tyler Trosclair, John Willis
% Sponsored by Dr. Kaushallya Adhikari
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set flag to one if you want to plot the peak side lobe peak 
% and pause the execution. You shouldn't really need to set this to one
% unless something is changed with the PSL part. 
tic
%% Set variables
calc_min = 1; % if 1, plot the results from minimum processing
calc_product = 0; % if 1, plot the results from product processing
close_graphs = 0; % close the graphs after they are generated. 

for spacing = 1:1
flag = 0;
diff_min_prod = [];

% Here are the coprime pairs we will test
% Coprimes.Pairs = {[2,3], [3,4], [4,5], [5,6], [6, 7], [8, 9], [9, 10]};
Coprimes.Pairs = GenerateCoprimePairs(3,7,spacing);

% Here is our preallocation for the Min and Prod data
Coprimes.MinData = cell(1,length(Coprimes.Pairs));
Coprimes.ProdData = cell(1,length(Coprimes.Pairs)); 



% max number of periods added. 
max_multiples = 5;

multipliers = 0:max_multiples;


% For each coprime pair...  
for i = 1:length(Coprimes.Pairs)
    % ... set the M and N. 
    M = Coprimes.Pairs{i}(1);
    N = Coprimes.Pairs{i}(2);
    % Initialize variables and preallocate space 
    U1 = N; 
    U2 = M;
    total_Bmin_max = zeros(1,max_multiples-1);
    total_Bprod_max = zeros(1,max_multiples -1);
    % Getting data from increasing each subarray by one period
    for add = multipliers
        % set add1, add2 to however many more periods 
        %   we want to add for this iteration
        add1 = (U1*add);
        add2 = (U2*add);

        % We need Bmin, Bprod, and u to look in the correct range for the peak
        % side lobe
        [Bmin, Bprod, u] = ProductMinBeampattern(M, N, U1, U2, add1, add2);
        % Our search range is defined as where the main lobe goes from the 
        % bottom of the main lobe up to 1 
        %(only need to search one side since it's symetric for our case)
        u_search_range = u > (2/(U1*(M+add1)));

        % Search the actual range for the PSL for Product/Min Processing
        [Bmin_max, Bmin_max_loc] = max(Bmin(u_search_range));
        [Bprod_max, Bprod_max_loc] = max(Bprod(u_search_range));

        % This is to find set the index of B*_max_loc to be in relation to
        % the entirity of Bmin, not just the indices searched. 
        start_index = find(u_search_range == 1,1)-1;
        Bmin_max_loc = start_index + Bmin_max_loc;
        Bprod_max_loc = start_index + Bprod_max_loc;

        % Append our found total_B*_max to our total vector
        total_Bmin_max(add+1) = Bmin_max;
        total_Bprod_max(add+1) = Bprod_max;
        
        % If flag is 1, then we make sure our PSL finding is accurate. 
        if flag == 1
            gcf;
            hold on;
            plot(u(Bmin_max_loc), Bmin_max, 'r*');
            plot(u(Bprod_max_loc), Bprod_max,'b*');
            pause;
            close gcf;
        else
            % else, close the figure that is generated from
            % ProductMinBeampattern
            close gcf;
        end

    end
    % Put our new Bmin and Bprod data into the corresponding cell and field
    Coprimes.MinData{i} = total_Bmin_max;
    Coprimes.ProdData{i} = total_Bprod_max;
end
%{
%% Plot the found B*_maxes
% % % % This section is useful, but the mesh display below is a better
% representation of the data
figure;
% Just a vector to show on the graph which multiplier corresponds to what
% data point. 
multiplier = 2:max_multiples;
hold on;
% We generate a 3D plot 
for i = 1:length(Coprimes.Pairs)
    plot3(i*ones(size(Coprimes.MinData{i})), multiplier , Coprimes.MinData{i}, ':');
    plot3(i*ones(size(Coprimes.ProdData{i})), multiplier, Coprimes.ProdData{i}, '--');
    
end
view(3);
set(gca,'Ydir','reverse')
%}
    %% Plot Min Processing mesh
    if calc_min == 1
        fig = figure;
        hold on;
        % Initialize our three matrices for our data points from the various
        % coprime pairs. 
        X = [];
        Y = [];
        Z = [];
        % Iterage over the cell array, combining data into rows. 
        for i = 1:length(Coprimes.Pairs)
            X = [X;i*ones(size(Coprimes.MinData{i}))];
            Y = [Y;multipliers];
            Z = [Z;Coprimes.MinData{i}];
        end
        % This is the target dB. The following lines of code find the point in the
        % data closest to the target dB and plot a red star to denote it. 
        target_dB = 13;
        Min_target = min(min(abs(Z+target_dB)));
        [row, col] = find(abs(Z + target_dB) == Min_target);

        % Plot the mesh
        mesh(X,Y,Z);
        title(['Min Processing PSLs with coprime difference = ', num2str(spacing)]);
        xlabel('Coprime Pair');
        ylabel('Additional Periods');
        zlabel('PSL Power (dB)');
        view(3); % Set the view to isometric
        set(gca,'Ydir','reverse'); % I reversed the y axis for readability 
            %(this essentially just rotates the graph in a certain way)
            
        % Plot our target_dB point
        plot3(X(row,col), Y(row,col), Z(row, col), 'r*');
        
        % Custom data cursor for readability
        dcm_obj = datacursormode(fig);
        set(dcm_obj, 'UpdateFcn', {@myupdatefcn,Coprimes.Pairs});
        
        % Save the power from the two period power for each coprime pair
        Coprimes.Two_Period_Power{spacing} = Z(:,2);
         
    end
    
    Min_Z_mesh = Z;
    %% Plot Product processing mesh
    if calc_product == 1
        % I leave out the comments for this section since it's just a 
        % repeat of the product processing mesh section.
        figure;
        hold on;
        X = [];
        Y = [];
        Z = [];

        for i = 1:length(Coprimes.Pairs)
            X = [X;i*ones(size(Coprimes.ProdData{i}))];
            Y = [Y;multipliers];
            Z = [Z;Coprimes.ProdData{i}];
        end
        target_dB = 13;
        Prod_target = min(min(abs(Z+target_dB)));
        [row, col] = find(abs(Z + target_dB) == Prod_target);
        mesh(X,Y,Z);
        title(['Product Processing PSLs with coprime difference = ', num2str(spacing)]);
        xlabel('Coprime Pair');
        ylabel('Additional Periods');
        zlabel('PSL Power (dB)');
        view(3);
        set(gca,'Ydir','reverse')
        hold on;
        plot3(X(row,col), Y(row,col), Z(row, col), 'r*');
        
        dcm_obj = datacursormode(fig);
        set(dcm_obj, 'UpdateFcn', {@myupdatefcn,Coprimes.Pairs});
    end
    if calc_min == 1 && calc_product == 1
        Prod_Z_mesh = Z;
        diff_min_prod = [diff_min_prod , min(min(Prod_Z_mesh - Min_Z_mesh))];
        Coprimes.Diff_min_prod{spacing} = diff_min_prod;
    end
    if close_graphs == 1
        close all
    end
end
toc;