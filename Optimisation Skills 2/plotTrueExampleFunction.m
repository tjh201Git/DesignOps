function plotTrueExampleFunction(func)
% PLOTTRUEEXAMPLEFUNCTION Plots a 3D cube with a 4th dimension heatmap.
%
%   Args:
%       func: A function handle, e.g., @(coords) sum(coords, 1), 
%             that accepts ONE 3xN array (coords) and returns a 1xN array (W).

    %% 1. Grid Setup
    nodes = 20;
    lin = linspace(0, 1, nodes);

    % Generate the 3D grid (70x70x70 arrays)
    [XX, YY, ZZ] = meshgrid(lin, lin, lin);

    %% 2. Data Preparation for Input
    % Flatten the 3D arrays into column vectors (N x 1)
    x_flat = XX(:);
    y_flat = YY(:);
    z_flat = ZZ(:);

    % Create the single 3xN coordinate array as requested: [xvals; yvals; zvals]
    % Transposing the flat vectors (x_flat') makes them 1xN, 
    % resulting in a 3xN matrix when concatenated vertically.
    coords_array = [x_flat, y_flat, z_flat];

    %% 3. Calculate the 4th Dimension (W)
    % Call the user-provided function handle with the single coordinate array
    try
        W_flat = func(coords_array);
    catch ME
        error('Plotting:FunctionError', ...
              ['Error calling the provided function handle. ' ...
               'Ensure your function accepts ONE 3xN coordinate array ' ...
               'and returns a 1xN array of values. Details: %s'], ...
               ME.message);
    end
    
    % Ensure W_flat is a column vector (N x 1) for scatter3 plotting
    if isrow(W_flat)
        W_flat = W_flat';
    end
    
    % Basic validation
    if length(W_flat) ~= length(x_flat)
         error('Plotting:SizeMismatch', ...
              'The function output size (%d) does not match the input point count (%d).', ...
              length(W_flat), length(x_flat));
    end

    %% 4. Plotting
    figure; % Open a new figure window

    % Create the 3D scatter plot
    % S=1 is the marker size; W_flat is the color data
    scatter3(x_flat, y_flat, z_flat, 1, W_flat, 'filled');

    % Set axes and title
    xlabel('X position');
    ylabel('Y position');
    zlabel('Z position');
    title_str = sprintf('3D Cube Heatmap: Function %s', func2str(func));
    title(title_str);

    % Customization
    axis tight;
    grid on;
    view(3); % Set to 3D view

    % Add a color bar
    h = colorbar;
    ylabel(h, 'Function Output (W)');
end