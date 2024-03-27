function var = get_indices_dist(time,type,lag,stdz)
% Syntax: var = get_indices(time,type,lag,stdz)
% get_mad_indices returns the geomagnetic/solar indices. For time lagged 
% indices, specify lag in hours.
%
% INPUT
%      time: date and time in matlab time format when indices is wanted
%      type: Indices catergory. Following string are accepted
%            kp, ap3, ap, f107, fbar, fism2, fism22 (fism square), fism23
%       lag: Time lag in hours [default=0]
%      stdz: Standardize data if set to 1. [default=0]
% OUTPUT
%      var: output index at hourly resolution with given lag
%
%
% NOTE
% - data gaps are filled using linear interpolation
%
% -------------------------------------------------------------------------
% Dupinder Singh (dupinder@mit.edu)
% MIT Haystack Obserrvatory
% Release Date: 25 Oct 2023 
% Version: --
% -------------------------------------------------------------------------


% Input parsing
if ~isdatetime(time)
    error('First input must be of datetime type')
end
time.TimeZone = '';
type = lower(type);
if ~ismember(type,{'kp','ap3','ap','f107','f1072','f1073','fbar','fism2',...
        'fismr1','fismr2','fismr3','fismr12','fismr22','fismr32',...
        'fismr13','fismr23','fismr33',...
        'fism22','fism23','fism2fl','hp60','ap60','dst',...
        'hpa70','hpa50','hpa40','hpa30','hpa20','hpa15','hpa10'})
    error('Geophysical index not found');
end
if nargin==2
    lag=0;
    stdz=0;
elseif nargin==3
    stdz=0;
end

persistent geoindices norms
if isempty(geoindices)
    geo_file = 'geoindices.mat';
    norm_file = 'norms.mat';
    load(geo_file,'geoindices');
    load(norm_file,'norms');
end
time2 = geoindices.time;
switch type
    case 'f1072'
        var = geoindices.f107.^2;
    case 'f1073'
        var = geoindices.f107.^3;
    case 'fism22'
        var = geoindices.fism2.^2;
    case 'fism23'
        var = geoindices.fism2.^3;
    case 'fismr12'
        var = geoindices.fismr1.^2;
    case 'fismr22'
        var = geoindices.fismr2.^2;  
    case 'fismr32'
        var = geoindices.fismr3.^2;
    case 'fismr13'
        var = geoindices.fismr1.^3;
    case 'fismr23'
        var = geoindices.fismr2.^3;  
    case 'fismr33'
        var = geoindices.fismr3.^3;       
    otherwise
        var = geoindices.(type);
end

% Introduce lag
if lag
    var = circshift(var,lag);
end
var=interp1(time2,var,time,'linear','extrap');

% Standardize data. This is important if these indexes are to be used in
% empirical modelling. This normalisation is based on mean and std of these
% parameter calculated from all the available indices data. 
if stdz
    switch type
        case 'kp'
            % do nothing
        case 'ap3'
            var = (var-norms.ap3_mean)/norms.ap3_std;
        case 'ap'
            % do nothing
        case 'f107'
            var = (var-norms.f107_mean)/norms.f107_std;
        case 'fbar'
            % do nothing
        case 'fism2'
            var = (var-norms.fism2_mean)/norms.fism2_std;
        case 'fismr1'
            var = (var-norms.fismr1_mean)/norms.fismr1_std;
        case 'fismr2'
            var = (var-norms.fismr2_mean)/norms.fismr2_std;
        case 'fismr3'
            var = (var-norms.fismr3_mean)/norms.fismr3_std;
        case 'fismr12'
            var = (var-norms.fismr12_mean)/norms.fismr12_std;
        case 'fismr22'
            var = (var-norms.fismr22_mean)/norms.fismr22_std;
        case 'fismr32'
            var = (var-norms.fismr32_mean)/norms.fismr32_std;    
        case 'fism22'
            var = (var-norms.fism22_mean)/norms.fism22_std;
        case 'fism23'
            var = (var-norms.fism23_mean)/norms.fism23_std;
        case 'fism2fl'
            var = (var-norms.fism2fl_mean)/norms.fism2fl_std;
        case 'hp60'
            var = (var-norms.hp60_mean)/norms.hp60_std;
        case 'ap60'
            var = (var-norms.ap60_mean)/norms.ap60_std;
        case 'dst'
            var = (var-norms.dst_mean)/norms.dst_std;
        case 'hpa30'
            var = (var-norms.hpa30_mean)/norms.hpa30_std;
    end
end

end