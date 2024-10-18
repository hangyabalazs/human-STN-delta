function allsessions = getdata2analyse(rootdir, varargin)

%{
  allsessions = getdata2analyse(rootdir,...) makes a struct containing all
  necessary information for sessions matching the experimental setting defined
  by input parameters.
  
  Required input:
  rootdir     root directory with patient folders

  Optional inputs (name-value pairs):
  'rectype'     specifies type of recording
              'LFP' | 'EEG' | 'BEHAV' (default: 'EEG')

  'rectime'     time of recording relative to surgery (DBS implantation)
              'preop' | 'intraop' | 'postop' (default: 'postop')

  'patients'   patients to analyse
               'allpatients' | cell array of character vectors, ex.: {'pd01','pd02'}
               (default: 'allpatients')

  'side'       side contralateral to the hand used in the task
              'right' | 'left' | 'bothside' (default: 'bothside')

  'condition'  condition of DBS stimulation 
              'stimoff' | 'stimon' | 'bothcond' | 'nocond' (default: 'bothcond')
              ('stimon' and 'bothcond' are used only in postoperative settings)

%}


prs = inputParser;
addRequired(prs,'rootdir',@isdir);
addParameter(prs,'rectype','EEG',@ischar)
addParameter(prs,'rectime','postop',@ischar)
addParameter(prs,'patients','allpatients',@(x) iscell(x)|ischar(x))
addParameter(prs,'side','bothside',@(x) iscell(x)|ischar(x))
addParameter(prs,'condition','bothcond',@(x) iscell(x)|ischar(x));
parse(prs,rootdir,varargin{:})
g = prs.Results;




if strcmp(g.rectype, 'EEG') && strcmp(g.rectime, 'postop')
    
    if strcmp(g.patients,'allpatients')
        g.patients = {'pd01','pd02','pd03','pd04','pd05','pd08','pd10','pd12','pd13','pd14','pd15','pd16'};
    end
    
    patinx = find(cell2mat(cellfun(@(x) ismember(x,{'pd01','pd02','pd03','pd04','pd05','pd08','pd10','pd12','pd13','pd14','pd15','pd16'}),g.patients,'UniformOutput',0)));
    if ~isempty(patinx)
        g.patients = g.patients(patinx);
    end
    
    allsessions = get_sessions(g.rootdir,g.patients,g.side,g.condition); 

    
elseif strcmp(g.rectype, 'BEHAV') && strcmp(g.rectime, 'postop')   
    
    allpats = {'pd01','pd02','pd03','pd04','pd05','pd07','pd08','pd10','pd12','pd13','pd14','pd15','pd16'};
     if strcmp(g.patients,'allpatients')
        g.patients = allpats;
    end
    
    patinx = find(cell2mat(cellfun(@(x) ismember(x,allpats),g.patients,'UniformOutput',0)));
    if ~isempty(patinx)
        g.patients = g.patients(patinx);
    end
    
    allsessions = get_sessions(g.rootdir,g.patients,g.side,g.condition); 

    
    
elseif strcmp(g.rectype, 'EEG') && strcmp(g.rectime, 'intraop')
    
    allpats = {'pd02','pd03','pd04','pd05','pd07','pd08','pd10','pd12','pd13','pd14','pd15','pd16'};
    if strcmp(g.patients,'allpatients')
        g.patients = allpats;
    end

    
      patinx = find(cell2mat(cellfun(@(x) ismember(x,allpats),g.patients,'UniformOutput',0)));
  
    if ~isempty(patinx)
        g.patients = g.patients(patinx);
    end
    
%     sessioninfos: patient name, session, trigger channel, EEG channels
    sessioninfos = {'pd02','181123','CH2','F3','F4';
        'pd03','181130','CH2','F3','F4';
        'pd04','190118','CH2','F3','F4';
        'pd05','190201','CH1','F3','F4';
        'pd07','190405','CH1','F3','F4';
        'pd08','190412','CH1','F3','F4';
        'pd10','191009','CH2','F3','F4';
        'pd12','191025','CH2','F3','F4';
        'pd13','200121','CH1','F3','F4';
        'pd14','200129','CH2','F3','F4';
        'pd15','200205','CH2','F3','F4';
        'pd16','200221','CH2','F3','F4'};
    
    
    
    allsessions = get_sessions(g.rootdir, g.patients, g.side, g.condition,sessioninfos);
    
    
    
    
    
elseif strcmp(g.rectype, 'LFP') && strcmp(g.rectime, 'intraop')
    
    allpats = {'pd01','pd02','pd03','pd04','pd05','pd07','pd08','pd10','pd12','pd13','pd14','pd15','pd16'};
    if strcmp(g.patients,'allpatients')
        g.patients = allpats;
    end
    
    
       patinx = find(cell2mat(cellfun(@(x) ismember(x,allpats),g.patients,'UniformOutput',0)));
  
    if ~isempty(patinx)
        g.patients = g.patients(patinx);
    end
    
    % sessioninfos: patient name, session
    sessioninfos = {'pd01','181116';
        'pd02','181123';
        'pd03','181130';
        'pd04','190118';
        'pd05','190201';
        'pd07','190405';
        'pd08','190412';
        'pd10','191009';
        'pd12','191025';
        'pd13','200121';
        'pd14','200129';
        'pd15','200205';
        'pd16','200221';};
    
    
    allsessions = get_sessions(g.rootdir, g.patients, g.side, g.condition,sessioninfos);
    
    
       
elseif strcmp(g.rectype, 'BEHAV') && strcmp(g.rectime, 'intraop')
    
    allpats = {'pd01','pd02','pd03','pd04','pd05','pd07','pd08','pd10','pd12','pd13','pd14','pd15','pd16'};
    if strcmp(g.patients,'allpatients')
        g.patients = allpats
    end
    

     patinx = find(cell2mat(cellfun(@(x) ismember(x,allpats),g.patients,'UniformOutput',0)));
    if ~isempty(patinx)
        g.patients = g.patients(patinx);
    end
    
    % sessioninfos: patient name, session
    sessioninfos = {'pd01','181116';
        'pd02','181123';
        'pd03','181130';
        'pd04','190118';
        'pd05','190201';
        'pd07','190405';
        'pd08','190412';
        'pd10','191009';
        'pd12','191025';
        'pd13','200121';
        'pd14','200129';
        'pd15','200205';
        'pd16','200221';};
   
    
    allsessions = get_sessions(g.rootdir, g.patients, g.side, g.condition,sessioninfos);
    
    
    
    
    elseif strcmp(g.rectime, 'preop')
    
        allpats = {'pd01','pd02','pd03','pd04','pd05','pd07','pd08','pd10','pd12','pd13','pd14','pd15','pd16'};
    if strcmp(g.patients,'allpatients')
       g.patients = allpats;
    end
    
     patinx = find(cell2mat(cellfun(@(x) ismember(x,allpats),g.patients,'UniformOutput',0)));
    if ~isempty(patinx)
        g.patients = g.patients(patinx);
    end
    
    
    % sessioninfos: patient name, session
    sessioninfos = {'pd01','181115';
        'pd02','181122';
        'pd03','181129';
        'pd04','190117';
        'pd05','190131';
        'pd07','190404';
        'pd08','190411';
        'pd10','191008';
        'pd12','191024';
        'pd13','200107';
        'pd14','200128';
        'pd15','200204';
        'pd16','200220'};
    
    
    allsessions = get_sessions(g.rootdir, g.patients, g.side, g.condition,sessioninfos);
    
    
end

for si = 1:length(allsessions)
    allsessions(si).rectype = g.rectype;
    allsessions(si).rectime= g.rectime;
    
    if strcmp(g.rectype,'LFP')
        allsessions(si).folder = [allsessions(si).folder '_LFP'];
        if ~isfolder(allsessions(si).folder); mkdir(allsessions(si).folder); end;
    end
end

if exist(fullfile(g.rootdir,['sessioninfos_' g.rectype '_' g.rectime '.mat']))~=2 && exist('sessioninfos')==1
    save(fullfile(g.rootdir,['sessioninfos_' g.rectype '_' g.rectime '.mat']),'sessioninfos');
end

%--------------------------------------------------------------------------
function allsessions = get_sessions(rootdir,sespatient, sesside, sescond,sessioninfos)


narginchk(4,5)

if nargin<5
    sessioninfos= [];
end
patfolds = dir([rootdir filesep 'PD_CellBase\patients']); patfolds = patfolds(~ismember({patfolds.name},{'.','..'}));
patfolds = patfolds(find([patfolds.isdir]==1));
if ~strcmp(sespatient,'allpatients')
    patfolds = patfolds(ismember({patfolds.name},sespatient));
    
end


switch sesside
    case {'left','right'}
        jmax =1;
        sidenm= {sesside};
    case 'bothside'
        jmax = 2;
        sidenm = {'left','right'};
end

switch sescond
    case {'stimoff','stimon','nocond'}
        comax =1;
        tagnm = {sescond};
    case 'bothcond'
        comax = 2;
        tagnm = {'stimoff','stimon'};
end



snr =0;
for i = 1:size(patfolds,1)
    patnm = patfolds(i).name;
    currpat = [rootdir filesep 'PD_CellBase\patients' filesep patnm];
    if  ~isempty(sessioninfos)
        patrow = find(strcmp(sessioninfos,patnm));
        posessRL = [str2num(sessioninfos{patrow,2}),str2num(sessioninfos{patrow,2})];
    else
        posessRL = find_postopsess(currpat);
    end
    for j = 1:jmax % 2 sessions from each patient (left & right)
        side = sidenm{j};

            currsess = [currpat filesep strcat(num2str(posessRL(j)),side(1))];
            sessnm = strcat(num2str(posessRL(j)),side(1));
            
            if ~isdir(fullfile(currpat,sessnm))
                fprintf('%s- no files\n',currpat(end-4:end));
                continue
            end
    cd(currsess);
        for condi = 1:comax
            
            snr = snr+1;
            tag = tagnm{condi};
            curr_resdir = [currsess filesep tag '_' side];
            if ~isdir(curr_resdir); mkdir(curr_resdir); end;
            allsessions(snr).sessfolder = currsess;
            allsessions(snr).folder = curr_resdir;
            allsessions(snr).side = side;
            allsessions(snr).tag = tag;
            allsessions(snr).patient = patnm;
            
        end
    end
end

T = struct2table(allsessions);
sortedT = sortrows(T, 'patient');
allsessions = table2struct(sortedT);


function posessRL = find_postopsess(currpat)

patsess = dir(currpat);  patsess = patsess(~ismember({patsess.name},{'.','..'}));
sessnrs = [];
for jp = 1:size(patsess,1)
    sessnrs =[sessnrs str2num(patsess(jp).name(1:end-1))];
end
postopsess = max(sessnrs);
if sum(ismember(sessnrs, postopsess)) == 2
    posessRL = [postopsess, postopsess];
else
    sessnrs(ismember(sessnrs,postopsess)) = []; % if days of right and left sesssions don't match
    postopsess2 = max(sessnrs);
    posessRL = [postopsess postopsess2];
end




