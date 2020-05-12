function quiverC2D(x,y,u,v,maxNumArrows)
%quiverC2D creates a 2D quiver plot and adds a color coding. The color coding is
%given by the absolut values of the component vectors. Large values result in colors 
%from the upper end of the used colormap. Plotting parameters have to be changed within 
%the function in this version. Values equal to NaN or inf are set to zero.
%In case of complex arrays the absolute value is used.
% 
%   INPUT:
%       x - 2D matrix, x components of initial points
%       y - 2D matrix, y components of initial points
%       u - 2D matrix, x components of arrows
%       v - 2D matrix, y components of arrows
%       maxNumArrows - a positive integer (non-integer should work as well)
%           number limiting the maximum number of plotted arrows. Since vectors
%           of length zero aren't plotted and the matrices are sampled
%           uniformly, it's possible that fewer arrows will be visible in the
%           end. If maxNumArrows is not given its value is set to 1000.
% 
%   OUTPUT:
%       none
% 
%   WARNING!: Using large datasets in combination with choosing maxNumArrows 
%       very large might result in this script running forever.
% 
% --------------------------------------------------------------------------------------
% 
%   EXAMPLE: 
%       [x,y] = meshgrid(linspace(0,10,100),linspace(0,10,100));
%       u = exp(-0.2*(x-5).^2 - 0.2*(y-5).^2);
%       v = -u;
%       quiverC2D(x,y,u,v,1000);
%   
% --------------------------------------------------------------------------------------
% 
%   See also: QUIVER, LINESPEC, COLORMAP.
% 
% 

%% prearrangements

narginchk(4,5);
n_in = nargin;

if ~isequal(size(x),size(y),size(u),size(v))
    error('X,Y,U,V have to be matrices of the same size.');
end


%% assignments

% maximum number of arrows if necessary
if n_in == 4
    maxNumArrows = 1000;
end
% Line width
lw = 2;
% Maximum of arrow head size
hs = 2;


%% initialization
if numel(u) > maxNumArrows
       
    N = ceil(sqrt(numel(u)/maxNumArrows));
    
    x = x(1:N:end,1:N:end);
    y = y(1:N:end,1:N:end);
    u = u(1:N:end,1:N:end);
    v = v(1:N:end,1:N:end);
end

%% taking care of possible issues

x = issues(x);
y = issues(y);
u = issues(u);
v = issues(v);

nu = norm(u);
nv = norm(v);
if nu > 0
    u = u./nu;
end
if nv > 0
    v = v./nv;
end



%% colormap definition
C = colormap;
ncolors = size(C, 1);
I = sqrt(u.^2 + v.^2);
% assume that the colormap matrix contains the colors in its rows
Ic = round(I/max(max(I))*ncolors);
Ic(Ic == 0) = 1;

%% plotting
hold on

if numel(u) > ncolors
    % let's take an intelligent approach: group all values which belong to
    % the same color value
    for k = 1:ncolors
        mask = (Ic == k);    
        quiver(x(mask), y(mask), u(mask), v(mask), 0, 'Color', C(k, :), ...
            'LineWidth', lw, 'maxheadsize', hs);
    end
else
    % if there are so few values, it will be more efficient to plot each
    % arrow individually
    % linear indexing!
    for k = 1:numel(u)
        quiver(x(k), y(k), u(k), v(k), 'Color', C(Ic(k),:), ...
            'LineWidth', lw, 'maxheadsize', hs);
    end
end

hold off

end