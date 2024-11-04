function [polim sign p_vals mrl Z ftm] = PC_polar(allsp_ph,binnr,stat,isfig,varargin)
%PC_POLAR   Distribution of phases associated with unit spikes
%   [polim sign p_vals mrl Z ftm] = PC_POLAR(allsp_ph,binnr,stat,isfig,...)
%   draws a polar histogram of phase values in ALLSP_PH. The histogram is
%   divided to BINNR nr of bins. STAT statistical test is performed. 
% 
% Required inputs:
%   ALLSP_PH    vector of phase values associated with unit spikes
%   BINNR       numeric value, nr of bins in histogram
%   STAT        type of stat test to perform 
%     'ray' | 'rao'
%   ISFIG       if true, figure is generated
%     true | false
%
% Optional inputs:
%   PTEXT      true | false, if true the result of stat testing is written
%               on the plot (def. value: true)
%   COLOR      color of polar plot (by default: 'purple')
%

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


narginchk(4,6);
if nargin==4
    ptext = true;
    color = [];
elseif nargin==5
    ptext = varargin{1};
    color = [];
elseif nargin==6
    ptext = varargin{1};
    color = varargin{2};
end

[Z, ray_p, U, rao_p,mrl,ftm] = b_rao3_mod(allsp_ph);

switch stat
    case 'ray'
        if ray_p<0.05
            sign = 1;
        else
            sign = 0;
        end
        if ~isempty(ray_p)
            p_vals = ray_p;
        else
            p_vals =NaN;
        end
        
    case 'rao'
        if rao_p(2)<0.05
            sign = 1;
        else
            sign = 0;
        end
        p_vals = rao_p;
end


edges=[(0:binnr:360)/180*pi];


if isfig
    if isempty(color)
        color = rgb('purple');
    end
    
    
    polim =polarhistogram(allsp_ph,edges);
    polim.FaceColor = color;
    
    set(gca,'ThetaAxisUnits','radians');
    rL = rlim;
    
    if ptext
        switch stat
            case 'ray'
                if ray_p<0.05
                    c = 'r';
                    text(pi*1.2,rL(2)*1.5,['p = ' num2str(ray_p)],'Color',c);
                else
                    c = 'k';
                    text(pi*1.2,rL(2)*1.5,['p = ' num2str(ray_p)],'Color',c);
                end
                
            case 'rao'
                if rao_p(2)<0.05
                    c = 'r';
                    text(pi*1.2,rL(2)*1.5,['p<' num2str(rao_p(2))],'Color',c);
                else
                    c = 'k';
                    text(pi*1.2,rL(2)*1.5,['p> ' num2str(rao_p(1))],'Color',c);
                end
        end
    end
else
    polim = [];
end