function [r_p,wat_pup,ftm,mu,bestmod,mufit,kappafit,pfit] = ...
    PC_figstat3(sp_hilb,smwin,act_cellid,event,figdir,isplot,fignm,plottype)
%PC_FIGSTAT3    Population phase histogram plot with statistics
% 
% Input parameters:
%     SP_HILB     vector of complex or numeric values; 
%                 if input data is complex, include phase + mrl values in plot
%                 if input data is numeric (phase values), only phase
%                 distr. is plotted
%     SMWIN       1x2 vector, limits of time period (s) rel. to EVENT
%                 from where phase values derive
%     ACT_CELLID  char. array, unit ID
%     EVENT       char. array, event label
%     FIGDIR     path to save figure
%     ISPLOT     true|false, if true figure is generated, if false only
%               statistics are performed
%     FIGNM     file name to save figure
%     PLOTTYPE  type of polar plot
%       'scatter' | 'hist'

% See also: B_RAO3_MOD, B_WATSON, PC_POLAR 


% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

PC_win = length(smwin)-1;
[sign, r_p,ftm,mu,wat_pup,kappa,bestmod] = deal(nan(1,PC_win));
[AIC,BIC,mufit,kappafit,pfit] = deal(cell(1,PC_win));
wat_p = cell(1,PC_win);


if isplot
    HCfig = figure(1);
    set(HCfig,'Visible','off')
    binnr = 20; % bin nr of polar histogram
    rownr = 2;
    colnr = PC_win;
end

% Hilbert phases (Polar histogram) + statistics

for ww = 1:PC_win
    
    if ~isreal(sp_hilb{ww}) % if input data is complex -> include phase + mrl values in plot
        phasvals = angle(sp_hilb{ww}); %spike-coupled phase values of 1 small window of every trial
        mrlvals = abs(sp_hilb{ww});
    else
        phasvals = sp_hilb{ww};
        mrlvals = [];
    end
    
    if isempty(phasvals)
        continue
    end
    
    stat = 'ray';
    
    % Polarhistogram + ray. stat
    phasvals(isnan(phasvals)) = [];
    mrlvals(isnan(mrlvals)) = [];
    
    
    
    if isplot
        spnr = ww;
        subplot(rownr,colnr,spnr)
        
        if strcmp(plottype,'scatter') % only if mrl values are available
            color = rgb('purple');
            polim = polarscatter(phasvals,mrlvals,50,color,'filled');
            maxrL = 0.3;
            [~, r_p(ww), ~, ~,~,ftm(ww)] = b_rao3_mod(phasvals);
            
            
        elseif strcmp(plottype,'hist')
            [polim, sign(ww), r_p(ww), ~, ~,ftm(ww)] = PC_polar(phasvals,binnr,stat,isplot,false);
            %             maxrL = 15;
            LL = rlim;
            maxrL = LL(2);
        end
        
        
        
        rlim([0 maxrL])
        if ww==1
            annotation('textbox',[0.02 0.8 0.2 0.2],'String',['first win: [' num2str(smwin(ww)) ' ' num2str(smwin(ww+1)) '] s' ],...
                'FitBoxToText','on','LineStyle','none')
        elseif ww==2
            annotation('textbox',[0.02 0.75 0.2 0.2],'String',['second win:[' num2str(smwin(ww)) ' ' num2str(smwin(ww+1)) '] s' ],...
                'FitBoxToText','on','LineStyle','none')
        end
    else
        
        [~, ray_p, ~, rao_p,~,ftm(ww)] = b_rao3_mod(phasvals);
        
        switch stat; case 'ray';if ray_p<0.05; sign(ww) = 1; else; sign(ww) = 0; end; r_p(ww) = ray_p;
            case 'rao'; if rao_p(2)<0.05; sign(ww) = 1; else; sign(ww) = 0; end; r_p(ww) = rao_p;
        end
        
    end
    
    
    % Compass + watson
    
    hold on;
    
    [mu(ww),kappa(ww),~,wat_p{ww},~] = b_watson(phasvals);
    rL = rlim;
    polarplot([mu(ww) mu(ww)],[0 rL(2)],'k','LineWidth',3)

    
    % Write stat results to plot
    sp = subplot(rownr,colnr,spnr);
    if r_p(ww)<0.05;c = 'r'; else; c = 'k'; end
    posi = sp.Position;
    txtposi = posi; txtposi(2) = posi(2)*0.4;
    annotation('textbox',txtposi , 'string', ...
        ['Rayleigh stat: p = ' num2str(r_p(ww))],'Color',c,'LineStyle','none')
    
    
    if wat_p{ww}(2)<=0.05;c = 'r'; else; c = 'k'; end
    
    txtposi2 = posi; txtposi2(2) = posi(2)*0.3;
    annotation('textbox', txtposi2, 'string',...
        ['Watson1 stat: p = ' num2str(wat_p{ww}(1)) '-' num2str(wat_p{ww}(2))],'Color',c,'LineStyle','none')
    
    
    wat_pup(ww) = wat_p{ww}(2);
    
    
    spnr = colnr+ww;
%     % Fitting mixture of von Mises distributuions
%     if ~isempty(components)
%         [AIC{ww},BIC{ww},mufit{ww},kappafit{ww},pfit{ww}] = phaselock_modes_simulation_param_mod(phasvals',components,0);
%         
%         figure(1);
%         subplot(rownr,colnr,spnr)
%         plot(AIC{ww})
%         hold on
%         plot(BIC{ww})
%         legend({'AIC','BIC'});
%         setmyplot_balazs(gca)
%         
%         [~, bestmodb] = min(BIC{ww});
%         [~, bestmoda] = min(AIC{ww});
%         if bestmodb~=bestmoda
%             %         keyboard;
%             bestmod(ww) = max([bestmoda bestmodb]);
%         else
%             bestmod(ww) = bestmodb;
%         end
%         %         keyboard;
%     else
        [AIC{ww},BIC{ww},mufit{ww},kappafit{ww},pfit{ww}] = deal([]);
%     end
    
    
    
end

cellidtit = act_cellid; cellidtit(ismember(cellidtit,'_')) = '-';
% suptitle({event,cellidtit})
annotation('textbox',[0.02 0.6 0.2 0.2],'String',{event,cellidtit},...
    'FitBoxToText','on','LineStyle','none')
% set(HCfig,'Position',get(0,'Screensize'));

saveas(HCfig,fullfile(figdir,[fignm '.jpg']));
saveas(HCfig,fullfile(figdir,[fignm '.fig']));
close(HCfig)

close all;
end