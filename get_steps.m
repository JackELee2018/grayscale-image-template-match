function steps = get_steps(P, varargin)
%GET_STEPS Gets a cell array, STEPS, which specified templates rotating steps in each
%pyramid layers.
%
%   STEPS = GET_STEPS(P) creates a array whose elements specify template
%   rotating angle steps of every layer of image pyramid.
%
%   STEPS = GET_STEPS(P, STEP_RECIPE) makes a array of template rotating
%   angle-steps for each pyramid layers according to STEP_RECIPE. Because
%   there may be several rotating step strategy for a specified P. The
%   default value of STEP_RECIPE is 1.
% 
%  Example: steps = get_steps(P, varargin);

% Verify the correct number of input.
error(nargchk(1, 2, nargin));

% Get inputs and set default values.
layer = P;
if nargin > 1
    step_recipe = varargin{1};
else
    step_recipe = 1;
end

% Inputs should be positive integers.
if layer ~= round(layer) || step_recipe ~= round(step_recipe)
    error('Inputs should be integers--get_steps');
end

if layer < 0 || step_recipe < 0
    error('Inputs should be positive integers');
end

% Clear data in variable, STEPS.
steps = [];

switch layer
    case 1
        switch step_recipe
            case 1
                steps(1) = 1;
                steps(2) = 0.1;
            case 2
                steps(1) = 2;
                steps(2) = 0.2;
            otherwise
                error('The step_recipe should be no more than 2 according to the given P.');
        end
        
    case 2
        switch step_recipe
            case 1
                steps(1) = 4;
                steps(2) = 0.5;
                steps(3) = 0.2;
            otherwise
                error('The step_recipe should be no more than 1 according to the given P.');
        end
        
    case 3
        switch step_recipe
            case 1
                steps(1) = 10;
                steps(2) = 2;
                steps(3) = 0.8;
                steps(4) = 0.4;
            otherwise
                error('The step_recipe should be no more than 1 according to the given P.');
        end
        
    case 4
        switch step_recipe
            case 1
                steps(1) = 4;
                steps(2) = 2;
                steps(3) = 1;
                steps(4) = 0;
                steps(5) = 0.5;
            case 2
                steps(1) = 9;
                steps(2) = 1;
                steps(3) = 0.2;
                steps(4) = 0.1;
                steps(5) = 0.05;
            case 3
                steps(1) = 10;
                steps(2) = 2;
                steps(3) = 0.8;
                steps(4) = 0.4;
                steps(5) = 0.2;
            case 4
                steps(1) = 10;
                steps(2) = 2;
                steps(3) = 0.4;
                steps(4) = 0.2;
                steps(5) = 0.1;
            case 5
                steps(1) = 10;
                steps(2) = 2;
                steps(3) = 0.4;
                steps(4) = 0.1;
                steps(5) = 0;
            otherwise
                error('The step_recipe should be no more than 5 according to the given P.');
        end
    case 5
        switch step_recipe
            case 1
                steps(1) = 6;
                steps(2) = 1;
                steps(3) = 0.1;
                steps(4) = 0.04;
                steps(5) = 0;
                steps(6) = 0;
            otherwise
                error('The step_recipe should be no more than 5 according to the given P.');
        end
    otherwise
        error('No P that greater than 5 is available for now.');
end
        
