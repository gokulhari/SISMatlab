function out = integ(in)
% INTEG
% Part of SISMatlab
% Takes input of function in phys-space, does indefinite integration and 
% gives the resultant function in phys-space

% SISMatlab is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% SISMatlab is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with SISMatlab.  If not, see <https://www.gnu.org/licenses/>.
%
% Written by Gokul Hariharan
in = phys2cheb(in);
n = length(in)-1;
neve = (n+1)/2;
nodd = neve;
out = zeros(n,1);
ueve = zeros(neve,1);
uodd = zeros(nodd,1);
% split into eve and odd parts
for i = 0:(neve-1)
    k = i+1;
    ueve(k) = in(2*i+1);
end
for i = 0 : (nodd-1)
    k = i+1;
    uodd(k) = in(2*i+2);
end
out(1) = uodd(1)/4.0;

for i= 1:(neve-1)
    out(2*i+1) =(uodd(i)- uodd(i+1))/(4.0*i);
end
for i = 0:(nodd-2)
    out(2*i+2) = (ueve(i+1)- ueve(i+2))/(4*i+2.0);
end
i = nodd-1;
out(2*i+2) = ueve(i+1)/(4.0*i+2.0);
out = cheb2phys(out);

