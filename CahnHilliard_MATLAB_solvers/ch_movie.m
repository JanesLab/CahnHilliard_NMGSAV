function ch_movie(phi_t,t_out,varargin)
%This function creates a red-white-blue video trajectory of
%chemical states in the current directory
%
%INPUTS
    %phi_t = Multidimensional array of chemical states.
    %t_out = Vector of time steps corresponding to the third dimension of phi_t
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
valid_3darray = @(x) ismatrix(x(:,:,1));
valid_vector = @(x) isvector(x);
valid_integer = @(x) mod(x,1) == 0;
valid_filename = @(x) ischar(x) || isstring(x);
valid_filetype = @(x) strcmpi(x,'MPEG-4') || strcmpi(x,'Motion JPEG AVI') ...
    || strcmpi(x,'Uncompressed AVI');

%Set parser options and valid input criteria
addRequired(ch_movie_parser,'phi_t',valid_3darray);
addRequired(ch_movie_parser,'t_out',valid_vector);
   
addParameter(ch_movie_parser,'dtframes',default_dtframes,valid_integer);
addParameter(ch_movie_parser,'filename',default_filename,valid_filename);
addParameter(ch_movie_parser,'filetype',default_filetype,valid_filetype);
addParameter(ch_movie_parser,'colorbar_type',default_colorbar);

parse(ch_movie_parser,phi_t,t_out,varargin{:});

%Extract parsed inputs
phi_t = ch_movie_parser.Results.phi_t;
t_out = ch_movie_parser.Results.t_out;
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

for i = 1:dtframes:size(phi_t,3)
    h = figure('visible','off');
    image(phi_t(:,:,i),'CDataMapping','scaled');
    colorbar; 
    axis square;
    if colorbar_type == "default"
        clim([-1, 1]);
    elseif colorbar_type == "initial_range"
        clim([min(phi_t(:,:,1),[],'all') max(phi_t(:,:,1),[],'all')]); %Set color axis to initial conditions
    end
    g = gca;
    %Scale and center display roughly to the size of the mesh,
    %with extra space horizontally for the color bar
    % g.Position = [0.5-size(phi_t,2)/1000 0.5-size(phi_t,1)/1000 ...
    %     2.2*size(phi_t,2)/1000 2*size(phi_t,1)/1000];
    axis([0 size(phi_t,2) 0 size(phi_t,1)]);
    title(strcat('t = ',num2str(t_out(i))));
    colormap(interp1(1:100:1100,redbluecmap,1:1001)); %Expand redbluecmap to 1000 elements
    colorbar;
    writeVideo(phi_movie,getframe(h));
    if mod(i,20)==0
        fprintf("%f%% \n",round(100*i/size(phi_t,3),2));
    end
end
close(phi_movie);