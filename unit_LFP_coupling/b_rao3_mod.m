function [Z,p1,U,p2,mrl,ftm] = b_rao3_mod(x)
%RAO    Rao's Spacing Test.
%   [Z,P1,U,P2] = RAO(X) calculates Rayleigh's Z-Statistic (Z) and Rao's
%   U-Statistic (U) with the corresponding p-values (P1 and P2). Input
%   argument X is the vector of phase angles in radians.
%
%   See also WATSON and WATSONTWO.

%   Balazs Hangya
%   Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   New York 11724, USA
%   balazs.cshl@gmail.com
%
%   Institute of Experimental Medicine
%   Szigony street 43, Budapest 1083, Hungary
%   hangyab@koki.hu

%   Updated by Johanna Petra Szab�, 2024

% Input argument check
error(nargchk(1,1,nargin))

if isempty(x)
    [Z,p1,U,p2] = deal([]);
    return
end
% Angle
bang = x;

% bang = mod(x,2*pi);     % from 0 to 2*pi (from -pi to pi: angle(cos(x)+i*sin(x)))

% Mean resultant length
n = length(bang);
ftm = sum(exp(1).^(1i*bang)) / n;    % first trigonometric moment
mrl = abs(ftm);     % mean resultant length

% Rayleigh's Z-Statistic
Z = n * (mrl ^ 2);
p1 = exp(1) ^ (-1 * Z) * (1 + (2 * Z - Z ^ 2) / ...
    (4 * n) - (24 * Z - 132 * Z ^ 2 + 76 * Z ^ 3 - 9 * Z ^ 4) / (288 * n ^ 2));

% Rao's Spacing Test
xdeg = x * 180 / pi;
xdegs = sort(mod(xdeg,360));
n = length(xdegs);
spacings = [diff(xdegs), xdegs(1)-xdegs(n)+360];
U = 0.5 * sum(abs(spacings-360/n));     % Test Statistic
if n < 4
    %     error('Too small sample size.')
    p1 = []; p2 = [];
    return
elseif n <= 30
    trow = n - 3;
elseif n <= 32
    trow = 27;
elseif n <= 37
	trow = 28;
elseif n <= 42
	trow = 29;
elseif n <= 47
	trow = 30;
elseif n <= 62
	trow = 31;
elseif n <= 87
	trow = 32;
elseif n <= 125
	trow = 33;
elseif n <= 175
	trow = 34;
elseif n <= 250
	trow = 35;
elseif n <= 350
	trow = 36;
elseif n <= 450
	trow = 37;
elseif n <= 550
	trow = 38;
elseif n <= 650
	trow = 39;
elseif n <= 750
	trow = 40;
elseif n <= 850
	trow = 41;
elseif n <= 950
	trow = 42;
else
    trow = 43;
end
rao_table = raotable;
if U > rao_table(trow,1)
	p2(1) = 0;
    p2(2) = 0.001;
elseif U > rao_table(trow,2)
	p2(1) = 0.001;
    p2(2) = 0.01;
elseif U > rao_table(trow,3)
	p2(1) = 0.01;
    p2(2) = 0.05;
elseif U > rao_table(trow,4)
	p2(1) = 0.05;
    p2(2) = 0.10;
else
    p2(1) = 0.10;
    p2(2) = 1;
end

% -------------------------------------------------------------------------
function RT = raotable
RT = [247.32, 221.14, 186.45, 168.02;...
    245.19, 211.93, 183.44, 168.66;...
    236.81, 206.79, 180.65, 166.30;...
    229.46, 202.55, 177.83, 165.05;...
    224.41, 198.46, 175.68, 163.56;...
    219.52, 195.27, 173.68, 162.36;...
    215.44, 192.37, 171.98, 161.23;...
    211.87, 189.88, 170.45, 160.24;...
    208.69, 187.66, 169.09, 159.33;...
    205.87, 185.68, 167.87, 158.50;...
    203.33, 183.90, 166.76, 157.75;...
    201.04, 182.28, 165.75, 157.06;...
    198.96, 180.81, 164.83, 156.43;...
    197.05, 179.46, 163.98, 155.84;...
    195.29, 178.22, 163.20, 155.29;...
    193.67, 177.08, 162.47, 154.78;...
    192.17, 176.01, 161.79, 154.31;...
    190.78, 175.02, 161.16, 153.86;...
    189.47, 174.10, 160.56, 153.44;...
    188.25, 173.23, 160.01, 153.05;...
    187.11, 172.41, 159.48, 152.68;...
    186.03, 171.64, 158.99, 152.32;...
    185.01, 170.92, 158.52, 151.99;...
    184.05, 170.23, 158.07, 151.67;...
    183.14, 169.58, 157.65, 151.37;...
    182.28, 168.96, 157.25, 151.08;...
    181.45, 168.38, 156.87, 150.80;...
    177.88, 165.81, 155.19, 149.59;...
    174.99, 163.73, 153.82, 148.60;...
    172.58, 162.00, 152.68, 147.76;...
    170.54, 160.53, 151.70, 147.05;...
    163.60, 155.49, 148.34, 144.56;...
    159.45, 152.46, 146.29, 143.03;...
    154.51, 148.84, 143.83, 141.18;...
    151.56, 146.67, 142.35, 140.06;...
    148.06, 144.09, 140.57, 138.71;...
    145.96, 142.54, 139.50, 137.89;...
    144.54, 141.48, 138.77, 137.33;...
    143.48, 140.70, 138.23, 136.91;...
    142.66, 140.09, 137.80, 136.59;...
    142.00, 139.60, 137.46, 136.33;...
    141.45, 139.19, 137.18, 136.11;...
    140.99, 138.84, 136.94, 135.92];