function [p, R] = polypredcicall_mod(x,y,lev,str,res)
%POLYPREDCICALL   Plot regression line with confidence interval.
%   POLYPREDCICALL(X,Y,LEV) calls POLYPREDCI with input coordinates X and
%   Y. Regression line is overlayed on scatter plot. Confidence interval is
%   plotted at a confidence level defined by LEV.
%
%   POLYPREDCICALL(X,Y,LEV,'ROBUST') uses MATLAB's robust regression tool
%   for generating the regression line and associated t-test.
%
%   Dependence:
%   Star Strider (2020). polypredci (https://www.mathworks.com/matlabcentral
%   /fileexchange/57630-polypredci), MATLAB Central File Exchange. Retrieved 
%   November 20, 2020.
%
%   See also POLYPREDCI and ROBUSTFIT.

%   Balazs Hangya, 20-Nov-2020
%   Institute of Experimental Medicine
%   hangya.balazs@koki.hu 

% UPDATED by Johanna Petra Szab�, 10.2024

% Fit regression

if nargin<5
    res = 1;
end;
x = x(:);  % force column vectors
y = y(:);
n = 1; % polynomial order
alfa = lev; % desired significance
xv = min(x):res:max(x); % optional high-resolution vector
xv = xv(:);

[p, yhat, ci] = polypredci(x, y, n, alfa,xv); % define both ‘alfa’ & ‘xv’

% % Plot regression with confidence interval
% % H1 = figure;
% plot(x, y, 'ko')
% hold on
% plot(xv, yhat, '--r')
% plot(xv, yhat+ci, '--g')
% plot(xv, yhat-ci, '--g')
% hold off
% grid

% Plot with regress
X = [ones(length(x),1) x];
[b,bint,r,rint,stats] = regress(y(:),X);
p = stats(3);
pR = corrcoef(y,X(:,2));
R = pR(3);

if strcmp(str,'robust')  % robust regression (optional)
    [b,stats] = robustfit(x,y);
    p = stats.p(2);   % test for the slope being non-zero
    outliers_ind = abs(stats.resid)>stats.mad_s; %#ok<NASGU>
    nonoutliers_ind = abs(stats.resid)<=stats.mad_s;
end

if p<0.05
    col = 'r';
else
    col = 'k';
end

% H2 = figure; % Regression plot
scatter(x,y,[],'k');
axis square
icp = b(1);   % intercept
gr = b(2);   % gradient
yv = xv .* gr + icp;
hold on
errorshade(xv',yv',ci',[0.7 0.7 0.7])
plot(xv,yv,'Color','k','LineWidth',2)   % overlay regression line
text('Units','normalized','Position',[0.7 0.7],...
    'String',{['p = ' num2str(p)] ['R = ' num2str(R)]},'Color',col)
if strcmp(str,'robust')
    title('Robust regression')
else
    title('Regression (original, not robust)')
end

% Plot regression with confidence interval of non-outliers (identified by
% robustfit)
% [p, yhat, ci] = polypredci(x(nonoutliers_ind), y(nonoutliers_ind), n, alfa, xv); % fit line w confidence interval
% H3 = figure;
% plot(x, y, 'ko')
% hold on
% plot(x(nonoutliers_ind), y(nonoutliers_ind), 'k.', 'MarkerSize',18)
% plot(xv, yhat, '--r')
% plot(xv, yhat+ci, '--g')
% plot(xv, yhat-ci, '--g')
% hold off
% grid
% X = [ones(length(x(nonoutliers_ind)),1) x(nonoutliers_ind)];
% [b,bint,r,rint,stats] = regress(y(nonoutliers_ind),X);  % regression on non-outliers
% p = stats(3);
% pR = corrcoef(y(nonoutliers_ind),X(:,2));
% R = pR(3);
% text('Units','normalized','Position',[0.7 0.7],...
%     'String',{['p = ' num2str(p)] ['R = ' num2str(R)]})
% title('Regression on non-outliers')