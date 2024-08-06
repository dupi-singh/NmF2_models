function [fof2,nmf2,geoind] = call_nmf2(station,time,fism,ap3)
% SYNTAX
% [fof2,nmf2] = call_nmf2(station,time)
% [fof2,nmf2] = call_nmf2(station,doy,fism,ap3)
% call_nmf2()
% 
% call NmF2 model for a given station and time. Value of FISM2 and ap3 index
% can be kept constant. Call Without Arguments to display list of available
% Ionosonde models.
%
% INPUTS:
%      station: First 4 letter of station name
%         time: UT date and time in MATLAB datetime format 
%               [can be vector or nd array, can be mixed UT]
%               OR
%               Fractional DOY when making experimental call. When providing time
%               as doy, fixed value of fism and ap3 must be passed.
%
% INPUTS FOR EXPERIMENTAL RUNS: LEAVE THESE INPUTS EMPTY TO USE USE REAL
% FISM2 AND AP3 VALUES. PASS A VALUE TO KEEP THE PREDICTORS FIXED.
%         fism: fism index with 48 hour lag [scalar]
%          ap3: Vector of ap3 index with lags = [0,3,6,9,12,24,48,72]. A
%               scalar value can also be passed
% outputs
%      fof2 : Critical frequency of F2 layer
%      nmf2 : NmF2 in 1/m^3. Same size as the input time array.
%    geoind : Structure with following geophysiical indices fields [optional]
%             fism2: fism2 index values
%              f107: F10.7 (sfu) index values
%               ap3: 3 hourly ap index values
%             geoind will have the same size as time
%
%
% EXAMPLES (see examples.m for detailed usage)
% 1. Get Jicamarca NmF2/foF2 values for the year 2010
%    station = 'JICA';
%    sdate = datetime([2010 1 1 0 0 0]);
%    edate = datetime([2010 12 31 23 0 0]);
%    time = (sdate:hours(1):edate)';
%    [fof2,nmf2] = call_nmf2(station,time);
%
% 2. Get seasonal variation of NmF2 over JICAMRACA for fixed solar and
%    geomagnetic activity.
%    station = 'JICA';
%    time = (1/24:1/24:365)';
%    fism = 0.05;
%    ap3 = 6;
%    [fof2,nmf2] = call_nmf2(station,time,fism,ap3)
%
% -------------------------------------------------------------------------
% Dupinder Singh (dupinder@mit.edu)
% MIT Haystack Obserrvatory
% Release Date: 25 Oct 2023 
% Version: 3.0
% -------------------------------------------------------------------------

persistent norms
if isempty(norms)
    load('norms.mat','norms');
end

% validate inputs
if nargin~=0 && isnumeric(time) 
    if nargin<4
        error('Please specify FISM2 and ap3 when calling with DOY');
    else
        if fism2f107(fism)>260
            warning('Solar flux value is too high. Provinding model output equivalent to F10.7 = 260 sfu')
            fism = f1072fism(260);
        elseif fism2f107(fism)<65
            warning('Solar flux value is too low. Provinding model output equivalent to F10.7 = 65 sfu')
            fism = f1072fism(65);
        end
        if ap3>100
            warning('Ap3 value is too high. Provinding model output for ap3 = 100')
        end
    end
end
if nargin==0
    display_stations();
    return
end

% Set internal parameters
mag_proxy = 'ap3';
flux_proxy = 'fism2';
fism_lag = 48;
ver = 'fism2ap3_48';
sw_norm = 1;

% Parse Inputs and predictors
station = upper(station);
siz = size(time);
time = time(:);
nt = length(time);
if isnumeric(time)  % User provided geophysical indices; Experimental call
    geoind.ap3 = ap3;
    geoind.fism2 = fism;
    geoind.f107 = fism2f107(fism);
    % Do normalisation as was done during the model development
    fism = (fism-norms.fism2_mean)/norms.fism2_std;
    ap3 = (ap3-norms.ap3_mean)/norms.ap3_std;
    solp = repmat(fism,nt,1);
    solp2 = solp.^2;
    if length(ap3)==1
        ap3 = repmat(ap3,1,8);
    end
    magp = repmat(ap3(1),nt,1);
    magp_3 = repmat(ap3(2),nt,1);
    magp_6 = repmat(ap3(3),nt,1);
    magp_9 = repmat(ap3(4),nt,1);
    magp_12 = repmat(ap3(5),nt,1);
    magp_24 = repmat(ap3(6),nt,1);
    magp_48 = repmat(ap3(7),nt,1);
    magp_72 = repmat(ap3(8),nt,1);
    DOY = time;
    ut2 = round(24*(DOY-floor(DOY)),2);
else
    solp = get_indices_dist(time,flux_proxy,fism_lag,sw_norm);
    solp2 = solp.^2;
    magp = get_indices_dist(time,mag_proxy,0,sw_norm);
    magp_3 = get_indices_dist(time,mag_proxy,3,sw_norm);
    magp_6 = get_indices_dist(time,mag_proxy,6,sw_norm);
    magp_9 = get_indices_dist(time,mag_proxy,9,sw_norm);
    magp_12 = get_indices_dist(time,mag_proxy,12,sw_norm);
    magp_24 = get_indices_dist(time,mag_proxy,24,sw_norm);
    magp_48 = get_indices_dist(time,mag_proxy,48,sw_norm);
    magp_72 = get_indices_dist(time,mag_proxy,72,sw_norm);
    DOY = day(time,'DayofYear');
    ut2 = time.Hour+time.Minute/60;
end

% Load model coefficients
matfile = ['mat' filesep station '_nmf2_model_v' ver '.mat'];
if ~exist(matfile,'file')
    warning('Station code not found. Type call_nmf2() to display list of avaialble model codes.')
    disp('  ')
    user_input = input('Enter y to display list of available model codes now... ',"s");
    if strcmp(user_input,'y')
        display_stations();
    end
    [fof2,nmf2]=deal(nan);
    return
end
load(matfile,'coeffs_arr','ut','meta');
% make sure the model file used the same predictors
if strcmp(mag_proxy,meta.magnetic_proxy) || strcmp(flux_proxy,meta.flux_proxy)
end
ind1 = isnan(coeffs_arr(1,:));
coeffs_arr(:,ind1)=[];

% Calculate matrix
doysin = sin(2*pi*DOY/365);   % Annual
doycos = cos(2*pi*DOY/365);
doysin2 = sin(4*pi*DOY/365);  % Semi-annual
doycos2 = cos(4*pi*DOY/365);
doysin3 = sin(6*pi*DOY/365);  % Ter-annual 
doycos3 = cos(6*pi*DOY/365);
crossterm1 = solp.*doysin;
crossterm2 = solp.*doycos;
crossterm3 = solp.*doysin2;
crossterm4 = solp.*doycos2;
crossterm5 = solp.*doysin3;
crossterm6 = solp.*doycos3;
crossterm7 = solp2.*doysin;
crossterm8 = solp2.*doycos;
crossterm9 = solp2.*doysin2;
crossterm10 = solp2.*doycos2;
datamat = [ones(nt,1) solp,solp2,magp,magp_3,magp_6,magp_9,magp_12,magp_24,magp_48,magp_72,...
        doysin,doycos,doysin2,doycos2,doysin3,doycos3,...
        crossterm1,crossterm2,crossterm3,crossterm4,crossterm5,crossterm6,...
        crossterm7,crossterm8,crossterm9,crossterm10];

ut2 = floor(ut2 * 2)/2; % round to nearest .5 for matching with ut 
uut2 = unique(ut2);
par_mod = nan(nt,1);
for ii=1:length(uut2)
    [~,ind1] = min(abs(uut2(ii)-ut)); % find nearest UT model
    ind2 = uut2(ii)==ut2;
    par_mod(ind2) = datamat(ind2,:)*coeffs_arr(ind1,:)';
end
nmf2 = 10.^par_mod;
fof2 = nmf2_fof2(nmf2);
nmf2 = reshape(nmf2,siz);
fof2 = reshape(fof2,siz);
% Get un normalized indices with zero lag
if nargout == 3 && nargin==2
    fism2 = get_indices_dist(time,'fism2',0,0);
    ap3 = get_indices_dist(time,'ap3',0,0);
    f107 = get_indices_dist(time,'f107',0,0);
    geoind.fism2 = reshape(fism2,siz);
    geoind.ap3 = reshape(ap3,siz);
    geoind.f107 = reshape(f107,siz);
end
end % EOF call_model_nmf2

% LOCAL FUNCTIONS

function display_stations
% Display 4 letter codes and coordinates of available Ionosonde models
if ~exist('station_list.mat','file')
    station_list = generate_station_list;
    disp(station_list)
else
    load('station_list.mat','station_list')
    disp(station_list);
end
end

function station_list = generate_station_list
% Generate station list file and save it to work directory
files = dir(['mat' fsep '*.mat']);
n = length(files);
station_code = cell(n,1);
[geo_lat,geo_lon] = deal(nan(n,1));
for ii=1:n
    f = fullfile(files(ii).folder,files(ii).name);
    load(f,'meta');
    station_code{ii} = meta.stname;
    geo_lon(ii) = meta.lon;
    geo_lat(ii) = meta.lat;
end
station_list = table(station_code,geo_lat,geo_lon);
save('station_list.mat','station_list')
end


