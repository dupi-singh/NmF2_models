function f107 = fism2f107(fism)
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
% Use inverse quadratic of f107 to fism coeff
a = pp(1); b=pp(2); c=pp(3);
d = b/(2*a); e = c-b^2/(4*a);
f107 = -sqrt((fism-e)/a)-d;
end
