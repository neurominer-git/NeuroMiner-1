%==========================================================================
% FORMAT function [M, cases_m] = nk_MatchID(Sid, S, Did, D, act, infind)
%==========================================================================
% 
% Inputs:
% -------
% Match two m x n arrays (D, S) according to their m x 1 column vector IDs 
% (Did, Sid) and perform some action (act) using the matching info:
% act can be one of the following strings:
% 'intersect' / 'intersectnan'   : in source (S) AND destination (D) array,
% nan option fills entries not found in source with NaN
% 'src_not_dst' : in source array (s) but NOT in destination (D)
% 'dst_not_src' : in destination array (D) but NOT in source (S)
%
% Outputs:
% --------
% M             : matched / extracted arrays
% cases_m       : matched IDs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 06/2012

function [M, cases_m, Sind, Dind] = nk_MatchID(Sid, S, Did, D, act, infind, verbose)

if ~exist('infind','var') || isempty(infind), infind = false; end
if infind, funcomp = 'strfind'; else funcomp = 'strcmp'; end
if ~exist('verbose','var') || isempty(verbose), verbose = 2; end    
mSid=[]; mDid=[];
if ~isempty(S)
   Sfl = true;
   if size(Sid,1) ~= size(S,1)
       error('Source ID and Source Matrix have to have same number of rows')
   end
else
    Sfl = false;
end

if exist('D','var') && ~isempty(D)
   Dfl = true;
   if size(Did,1) ~= size(D,1)
       error('Destination ID and Destination Matrix have to have same number of rows')
   end
else
    Dfl = false;
end

if ~exist('act','var') || isempty(act)
    act = 'intersect';
end
cntnf = 0; 
switch act
            
    case {'intersect', 'intersectnan'}
        
        nanDid = [];
        for i = 1 : numel(Did)

            dstr = strtrim(deblank(Did{i}));
            if verbose == 2, fprintf('\nSearching match for %s ... ',dstr); end
            fndfl = false;

            for j = 1 : numel(Sid)

                sstr = strtrim(deblank(Sid{j}));
                a = feval(funcomp,sstr,dstr); 
                b = feval(funcomp,dstr,sstr);
                if (~isempty(a) && a) ||  (~isempty(b) && b)
                    if verbose == 2, fprintf('match found (%s).',sstr); end
                    mSid = [mSid; j]; 
                    mDid = [mDid; i];
                    fndfl = true;
                    break
                end

            end
            if ~fndfl, 
                if verbose == 2, cprintf('red*','NO MATCH FOUND !!!'); end
                cntnf = cntnf+1; 
                nfid{cntnf} = dstr;    
                if verbose == 2 && strcmp(act,'intersectnan'), cprintf('green', ' filled with NaN!'), end
            end
            nanDid = [nanDid; fndfl];
        end
        nanDid = logical(nanDid);
        
        if Sfl && Dfl
            
            switch act
                case 'intersect'
                    M = [ S(mSid,:) D(mDid,:) ] ;
                case 'intersectnan'
                    M = nan(numel(nanDid),size(S,2)+size(D,2));
                    M(nanDid,1:size(S,2)) = S(mSid,:);
                    M(nanDid,size(S,2)+1:end) = D(mDid,:);
            end
            
        elseif Sfl

            switch act
                case 'intersect'
                     M = S(mSid,:);
                case 'intersectnan'
                     M = nan(numel(nanDid),size(S,2));
                     M(nanDid,:) = S(mSid,:);
            end
            
        elseif Dfl

            
            switch act
                case 'intersect'
                    M = D(mDid,:);
                case 'intersectnan'
                    M = nan(numel(nanDid),size(D,2));
                    M(nanDid,:) = D(mDid,:);
            end

        else

            M = [ mSid mDid ] ;

        end
        
        cases_m = Did(mDid);        
         
    case 'src_not_dst'
        
        for j = 1 : numel(Sid)
            
            sstr = strtrim(Sid{j});
            if verbose == 2, fprintf('\nSearching match for %s ... ',sstr); end
            fndfl = false;
            
            for i = 1 : numel(Did)
                dstr = strtrim(Did{i});
                a = feval(funcomp,sstr,dstr); b = feval(funcomp,dstr,sstr);
                if (~isempty(a) && a) ||  (~isempty(b) && b)
                    if verbose == 2, fprintf('match found (%s).',dstr); end
                    fndfl = true;
                    break
                end
            end
            
            if ~fndfl, if verbose == 2, fprintf('store source ID (%s).',sstr); end; mSid = [ mSid j ]; end
            
        end
        if Sfl
            M = S(mSid,:) ;    
        else
            M = mSid;
        end
        
        cases_m = Sid(mSid);
        
    case 'dst_not_src'
        
        for i = 1 : numel(Did)
            
            dstr = strtrim(Did{i});
            fprintf('\nSearching match for %s ... ',dstr);
            fndfl = false;
            
            for j = 1 : numel(Sid)
                sstr = strtrim(Sid{i});
                a = feval(funcomp,sstr,dstr); b = feval(funcomp,dstr,sstr);
                if (~isempty(a) && a) ||  (~isempty(b) && b)
                    fprintf('match found (%s).',sstr);
                    fndfl = true;
                    break
                end
            end
            
            if ~fndfl, fprintf('store destination ID (%s).',dstr); mDid = [ mDid i ]; end
            
        end
        if Dfl
            M = D(mDid,:) ;    
        else
            M = mDid;
        end
        cases_m = Did(mDid);    
        
end
if cntnf>0,
    cprintf('red','\n==========================================================')
    fprintf('\n'); cprintf('red*','NO MATCHES FOUND FOR %g SUBJECTS!!! ', cntnf);
    for i = 1:cntnf
       cprintf('red','\n\t%s',nfid{i}); 
    end
    cprintf('red','\n==========================================================\n')
end
Sind = mSid; Dind = mDid;