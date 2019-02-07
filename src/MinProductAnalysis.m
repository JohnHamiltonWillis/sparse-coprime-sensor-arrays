%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Senior Design - Min and Product Processing Analysis
% Louisiana Tech University
% Pablo Johnson, Daniel Sartori, Tyler Trosclair, John Willis
% Sponsored by Dr. Kaushallya Adhikari
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datetime
% file_loc = 'C:\Users\ELEN 479\Documents\Generated_Tables\00001_res\';
% mkdir(file_loc);
% mkdir([file_loc,'Product']); 
% mkdir([file_loc,'Minimum']);
% mkdir([file_loc,'Difference']);
% mkdir([file_loc,'Figures']);
% mkdir([file_loc,'Figures\Minimum']);
% mkdir([file_loc,'Figures\Product']);
tic
%% Set variables
calc_min = 1; % if 1, calculate the results from minimum processing
plot_min = 1; % if 1, plot the results from minimum processing
calc_product = 1; % if 1, plot the results from product processing
plot_prod = 1; % if 1, calculate the results from product processing
close_graphs = 0; % close the graphs after they are generated. 
save_min = 0;
save_min_fig = 0;
save_prod = 0;
save_prod_fig = 0;

spacing_min = 1;
spacing_max = 1;

for spacing = spacing_min:spacing_max
flag = 0; % Used to verify PSL finding is accurate
diff_min_prod = [];

% Here are the coprime pairs we will test
% Coprimes.Pairs = {[2,3], [3,4], [4,5], [5,6], [6, 7], [8, 9], [9, 10]};
min_int = 3;
max_int = 5;

Coprimes.Pairs = GenerateCoprimePairs(min_int,max_int,spacing);

% Here is our preallocation for the Min and Prod data
Coprimes.MinData = cell(1,length(Coprimes.Pairs));
Coprimes.ProdData = cell(1,length(Coprimes.Pairs)); 



% max number of periods added. 
max_multiples = 10;

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
            Bmin_max
            Bprod_max
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
        Y = Y+1;
        % This is the target dB. The following lines of code find the point in the
        % data closest to the target dB and plot a red star to denote it. 
        target_dB = 13;
        Min_target = min(min(abs(Z+target_dB)));
        [row, col] = find(abs(Z + target_dB) == Min_target);
        
        if plot_min == 1
            % Plot the mesh
            fig = figure;
            hold on;
            mesh(X,Y,Z);
            title(['Min Processing PSLs with coprime difference = ', num2str(spacing)]);
            xlabel('Coprime Pair');
            ylabel('Periods');
            zlabel('PSL Power (dB)');
            view(3); % Set the view to isometric
            set(gca,'Ydir','reverse'); % I reversed the y axis for readability 
                %(this essentially just rotates the graph in a certain way)

            % Plot our target_dB point
            plot3(X(row,col), Y(row,col), Z(row, col), 'r*');

            % Custom data cursor for readability
            dcm_obj = datacursormode(fig);
            set(dcm_obj, 'UpdateFcn', {@MinProductAnalysis_updatefcn,Coprimes.Pairs});
            if save_min_fig
                filename = [file_loc,'Figures\Minimum\', num2str(min_int), '_', num2str(max_int), '_', num2str(spacing)];
                savefig(gcf,filename, 'compact');
            end
        end
        
        pair_names = cell(1,length(Coprimes.Pairs));
        for i = 1:length(Coprimes.Pairs)
            pair_names{i} = ['(',num2str(Coprimes.Pairs{i}(1)), ',',... 
                num2str(Coprimes.Pairs{i}(2)),')'];
        end
        
        col_names = cell(1,length(Z(1,:)));
        for i = 1:length(Z(1,:))
            col_names{i} = ['P',num2str(i)];
        end
        
        z_table = array2table(Z,'RowNames', pair_names, 'VariableNames', col_names);
        if save_min
            filename = [file_loc, 'Minimum\', num2str(min_int), '_', num2str(max_int), '_', num2str(spacing)];
            save(filename, 'z_table','X','Y','Z');
        end
        % Save the power from the two period power for each coprime pair
        Coprimes.Two_Period_Power{spacing} = Z(:,2);
        Coprimes.Min_PSL_table{spacing} = z_table;
         
    end
    
    Min_Z_mesh = Z;
    %% Plot Product processing mesh
    if calc_product == 1
        % I leave out the comments for this section since it's just a 
        % repeat of the product processing mesh section.

        X = [];
        Y = [];
        Z = [];

        for i = 1:length(Coprimes.Pairs)
            X = [X;i*ones(size(Coprimes.ProdData{i}))];
            Y = [Y;multipliers];
            Z = [Z;Coprimes.ProdData{i}];
        end
        Y = Y+1;
        target_dB = 13;
        Prod_target = min(min(abs(Z+target_dB)));
        [row, col] = find(abs(Z + target_dB) == Prod_target);
        
        if plot_prod == 1
            fig = figure;
            hold on;    
            mesh(X,Y,Z);
            title(['Product Processing PSLs with coprime difference = ', num2str(spacing)]);
            xlabel('Coprime Pair');
            ylabel('Periods');
            zlabel('PSL Power (dB)');
            view(3);
            set(gca,'Ydir','reverse')
            hold on;
            plot3(X(row,col), Y(row,col), Z(row, col), 'r*');

            dcm_obj = datacursormode(fig);
            set(dcm_obj, 'UpdateFcn', {@MinProductAnalysis_updatefcn,Coprimes.Pairs});
            if save_prod_fig 
                filename = [file_loc, 'Figures\Product\', num2str(min_int), '_', num2str(max_int), '_', num2str(spacing)];
                savefig(gcf,filename, 'compact');
            end
        end
                pair_names = cell(1,length(Coprimes.Pairs));
        for i = 1:length(Coprimes.Pairs)
            pair_names{i} = ['(',num2str(Coprimes.Pairs{i}(1)), ',',... 
                num2str(Coprimes.Pairs{i}(2)),')'];
        end
        
        col_names = cell(1,length(Z(1,:)));
        for i = 1:length(Z(1,:))
            col_names{i} = ['P',num2str(i)];
        end
        
        z_table = array2table(Z,'RowNames', pair_names, 'VariableNames', col_names);
        if save_prod
            filename = [file_loc, 'Product\', num2str(min_int), '_', num2str(max_int), '_', num2str(spacing)];
            save(filename, 'z_table','X','Y','Z');
        end
        Coprimes.Prod_PSL_table{spacing} = z_table;
    end
    
%% Difference Calculations
    if calc_min == 1 && calc_product == 1
        Prod_Z_mesh = Z;
        diff_min_prod = [diff_min_prod , min(min(Prod_Z_mesh - Min_Z_mesh))];
        Coprimes.Diff_min_prod{spacing} = diff_min_prod;
    end
    if close_graphs == 1
        close all
    end
end

if save_min && save_prod
    filename = [file_loc,'Difference\', num2str(min_int), '_', num2str(max_int), '_', num2str(spacing)];
    save(filename, 'diff_min_prod');
end
toc;