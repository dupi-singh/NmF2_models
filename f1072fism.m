function fism = f1072fism(f107)
% Returns fism2 value corresponding to given f107 based in quadratic
% relation between two.
% -------------------------------------------------------------------------
% Dupinder Singh (dupinder@mit.edu)
% MIT Haystack Obserrvatory
% Release Date: 25 Oct 2023 
% Version: --
% -------------------------------------------------------------------------
ifile = 'f107_2_fism.mat';
load(ifile,'pp')
fism = polyval(pp,f107);
end
