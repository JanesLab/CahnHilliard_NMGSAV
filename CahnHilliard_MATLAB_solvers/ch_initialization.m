function phi0 = ch_initialization(nx,ny,varargin)
%This function sets initial conditions for Cahn-Hilliard simulations in
%terms of qualitative chemical states in the [-1,+1] interval.
%
%INPUTS
    %nx = Number of grid points in x. No default value but must be even.
    %ny = Number of grid points in y. No default value but must be even.
%
%NAME-VALUE PAIRS
    %inittype = Type of initial conditions:
    %   'spinodal' (default) - assortment of ±1 states with equal probability.
    %   'unirand' - uniform random distribution over [-1,+1].
    %   'droplet' - circle of radius R0 with equilibrated interface profile.
    %   'geometric' - full-height rectangle of width W0 plus a circle of diameter D0.
    %   'file' - delimited file of nx by ny entries.
    %seed = Pseudorandom number generator seed for 'spinodal' and 'unirand'. Default = 'none'.
    %R0 = Initial radius for 'droplet' as a fraction of the domain. Default = 0.1 of domain.
    %epsilon = Interfacial transition distance for the interface profile of 'droplet'. Default = 0.01.
    %W0 = Width of the full-height rectangle for 'geometric' in grid points. Default = 4 grid points.
    %D0 = Diameter of the circle for 'geometric' in grid points. Default = 20 grid points.
    %initfile = Path/file containing values to be mapped to the nx-ny grid. Default = ''.
    %chemstates = Specifies the concentration of soluble and droplet states in 'file'. Default = [-1 1].
%
%OUTPUT
    %phi0 = Initialized chemical states in nx by ny grid points.



%% Set option defaults and parse inputs

%Set parameter defaults
default_inittype = 'spinodal';
default_seed = 'none';
default_R0 = 0.1;
default_epsilon = 0.01;
default_W0 = 4;
default_D0 = 20;
default_initfile = '';
default_chemstates = [-1 1];

ch_initialization_parser = inputParser;

%Set general criteria for inputs and name-value pairs
valid_even = @(x) mod(x,2) == 0;
valid_integer = @(x) x-floor(x) == 0;
valid_inittype = @(x) strcmpi(x,'spinodal') || strcmpi(x,'unirand') || ...
    strcmpi(x,'droplet') || strcmpi(x,'geometric') || strcmpi(x,'file');
valid_pos_num = @(x) isnumeric(x) && (x > 0);
valid_fraction = @(x) isnumeric(x) && (x > 0) && (x < 1);
valid_char = @(x) ischar(x);
valid_chemstate = @(x) length(x) == 2 && (x(1) < x(2));

%Set parser options and valid input criteria
addRequired(ch_initialization_parser,'nx',valid_even);
addRequired(ch_initialization_parser,'ny',valid_even);
   
addParameter(ch_initialization_parser,'inittype',default_inittype,valid_inittype);
addParameter(ch_initialization_parser,'seed',default_seed,valid_pos_num);
addParameter(ch_initialization_parser,'R0',default_R0,valid_fraction);
addParameter(ch_initialization_parser,'epsilon',default_epsilon,valid_pos_num);
addParameter(ch_initialization_parser,'W0',default_W0,valid_integer);
addParameter(ch_initialization_parser,'D0',default_D0,valid_integer);
addParameter(ch_initialization_parser,'initfile',default_initfile,valid_char);
addParameter(ch_initialization_parser,'chemstates',default_chemstates,valid_chemstate);

parse(ch_initialization_parser,nx,ny,varargin{:});

%Extract parsed inputs
nx = ch_initialization_parser.Results.nx;
ny = ch_initialization_parser.Results.ny;
inittype = ch_initialization_parser.Results.inittype;
seed = ch_initialization_parser.Results.seed;
R0 = ch_initialization_parser.Results.R0;
epsilon = ch_initialization_parser.Results.epsilon;
W0 = ch_initialization_parser.Results.W0;
D0 = ch_initialization_parser.Results.D0;
initfile = ch_initialization_parser.Results.initfile;
chemstates = ch_initialization_parser.Results.chemstates;

%% Execute initializations

if ~ischar(seed) %if not 'none' default
    rng(seed);
end

switch inittype
    case 'spinodal'
        phi0 = (-1).^round(rand(nx,ny));        
    case 'unirand'
        phi0 = 1 - 2*rand(nx,ny);
    case 'droplet'
        [xx,yy] = meshgrid(1/nx*(1:nx),1/ny*(1:ny));
        R = sqrt((xx-0.5).^2 + (yy-0.5).^2);
        phi0 = tanh((R0-R)/(sqrt(2)*epsilon)); %Equilibrium interface profile
    case 'geometric'
        phi0 = -1*ones(nx,ny);
        for i = 1:nx
            for j = 1:ny
                if sqrt((i-nx/2)^2+(j-ny/2)^2) <= D0/2 %if within circle diameter
                    phi0(i,j) = 1;
                elseif i > round(nx/2)-W0 && i < round(nx/2)+W0 %if within full-height rectangle
                    phi0(i,j) = 1;
                end
            end
        end
    case 'file'
        c0 = readmatrix(initfile);
        phi0 = (2*c0-chemstates(2)-chemstates(1))/(chemstates(2)-chemstates(1)); %Convert chemical states
end

if ~strcmpi(inittype,'file') %if x-y arrangement is not prespecified
    phi0 = phi0'; %Transpose so that rows are the x-axis and columns the y-axis **CONVENTION OR NOT?**
end

end