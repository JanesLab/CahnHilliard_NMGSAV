function ch_movie_from_file(phi_file,t_out,ny,varargin)
%This function creates a red-white-blue video trajectory of
%chemical states in the current directory
%
%INPUTS
    %phi_file = File name for multidimensional array of chemical states.
    %t_out = Vector of time steps corresponding to the third dimension of phi_t
    %ny = Number of mesh points in the y direction
%
%NAME-VALUE PAIRS
    %dtframes = Number of dt time steps per frame. Default = 1 (all time
    %   steps recorded)
    %filename = Number of movie file to be saved. Default = 'ch_movie'.
    %filetype = Any permissible movie profile in VideoWriter. Most useful:
    %   'MPEG-4' (default) - MPEG-4 file with H.264 encoding; best
    %       compression but may cause artifacts with image fields
    %       of high spatial frequency
    %   'Motion JPEG AVI' - AVI file using Motion JPEG encoding; max
    %       quality JPEG compression
    %   'Uncompressed AVI' - Uncompressed AVI file with RGB24 video;
    %       largest file size
    % colorbar_type = Type of colorbar to be used. Default = "default".
    %   'default': ranges from -1 to +1
    %   'initial_range': ranges is set using initial phi
    %   'variable': can change on each iteration
%
%OUTPUT
    %None except for the saved red-white-blue video

warning off

%% Set option defaults and parse inputs

%Set parameter defaults
default_dtframes = 1;
default_filename = 'ch_movie';
default_filetype = 'MPEG-4';
default_colorbar = 'default';

ch_movie_parser = inputParser;

%Set general criteria for inputs and name-value pairs
valid_vector = @(x) isvector(x);
valid_integer = @(x) mod(x,1) == 0;
valid_filename = @(x) ischar(x) || isstring(x);
valid_filetype = @(x) strcmpi(x,'MPEG-4') || strcmpi(x,'Motion JPEG AVI') ...
    || strcmpi(x,'Uncompressed AVI');

%Set parser options and valid input criteria
addRequired(ch_movie_parser,'phi_file',valid_filename);
addRequired(ch_movie_parser,'t_out',valid_vector);
addRequired(ch_movie_parser,'ny',valid_integer);

addParameter(ch_movie_parser,'dtframes',default_dtframes,valid_integer);
addParameter(ch_movie_parser,'filename',default_filename,valid_filename);
addParameter(ch_movie_parser,'filetype',default_filetype,valid_filetype);
addParameter(ch_movie_parser,'colorbar_type',default_colorbar);

parse(ch_movie_parser,phi_file,t_out,ny,varargin{:});

%Extract parsed inputs
phi_file = ch_movie_parser.Results.phi_file;
t_out = ch_movie_parser.Results.t_out;
ny = ch_movie_parser.Results.ny;
dtframes = ch_movie_parser.Results.dtframes;
filename = ch_movie_parser.Results.filename;
filetype = ch_movie_parser.Results.filetype;
colorbar_type = ch_movie_parser.Results.colorbar_type;

% phi_movie = VideoWriter(strcat(cd,'/',filename),filetype); %Extension will be automatically appended
phi_movie = VideoWriter(filename,filetype);
if strcmpi(filetype,'MPEG-4') || strcmpi(filetype,'Motion JPEG AVI')
    phi_movie.Quality = 100; %Highest quality video compression
end
open(phi_movie);

% loop over t_out but still skip frames according to dtframes
% read in only the necessary frame from the file using
%A = readmatrix('yourtxtFileName', 'FileType', 'text', 'Range', [firstRow, 1, lastRow, 1])
% need to know correct size of phi_t, dtframes, dt_out

for i = 1:dtframes:size(t_out,2)
    phi_temp = readmatrix(phi_file, 'Range', [(i-1)*ny+1, 1, i*ny, 128]);
% for i = 1:dtframes:size(phi_t,3)
    h = figure('visible','off');
    image(phi_temp,'CDataMapping','scaled');
    colorbar; 
    axis square;
    if colorbar_type == "default"
        clim([-1, 1]);
    elseif colorbar_type == "initial_range"
        clim([min(phi_temp,[],'all') max(phi_temp,[],'all')]); %Set color axis to initial conditions
    end
    g = gca;
    %Scale and center display roughly to the size of the mesh,
    %with extra space horizontally for the color bar
    % g.Position = [0.5-size(phi_t,2)/1000 0.5-size(phi_t,1)/1000 ...
    %     2.2*size(phi_t,2)/1000 2*size(phi_t,1)/1000];
    axis([0 ny 0 ny]);
    title(strcat('t = ',num2str(t_out(i))));
    colormap(interp1(1:100:1100,redbluecmap,1:1001)); %Expand redbluecmap to 1000 elements
    colorbar;
    writeVideo(phi_movie,getframe(h));
    if mod(i/(size(t_out,2)-1)*100,5) == 0
        fprintf('%3.0f percent complete\n',i/size(t_out,2)*100)
    end
end
close(phi_movie);