%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define plot variables and functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Functions

% Function to make matrix one long column
makeLong = @(x) x(:);

% Function to calculate color matrices
rgb     = @(x) round(x./255,2);

%% Plot variables

axfntsz         = 14;
titlefntsz      = 16;
lgndfntsz       = 16;
labfntsz        = 20;

lnwid           = 2;
mksz            = 16;

set(0,'DefaultAxesFontName', 'Helvetica');
set(0,'DefaultTextFontname', 'Helvetica');

colz(1,:)   = rgb([0, 160, 0]);%rgb([221,28,26]);
colz(2,:)   = rgb([86, 186, 220]);%rgb([95 100 127]);
colz(3,:)   = rgb([255,188,31]); % yellow
colz(4,:)   = rgb([209,50,33]);
colzedge    = colz*.8;

% Neural network colors
nncolz(1,:)     = rgb([49 135 200]);%rgb([95 100 127]);
nncolz(2,:)     = rgb([255,60,40]);%rgb([221,28,26]);
nncolzedge      = nncolz*.8;

% ERP colors
erpcol      = [rgb([215,48,39]);...
                rgb([252,141,89]);...
                rgb([254,224,144]);...
                rgb([161,170,162]);...
                rgb([145,191,219]);...
                rgb([69,117,180])];