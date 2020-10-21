function [out, fnd, out2, out3, out4] = nk_CheckLoadFile(pth, filetyp, f, d, ovrwrt, nclass)
global FUSION

% Check file type
a = regexp(filetyp,'CVdatamat','ONCE'); if ~isempty(a), a=1; else a=0; end

nvar = size(pth,1); % get number of variates

out = []; fnd = false; out2=[]; out3 = [];
if ~exist('ovrwrt','var') || isempty(ovrwrt), ovrwrt  = false; end

for i=1:nvar

    px = deblank(pth(i,:));
    
    if exist(px,'file')
        fnd = true;
        if ovrwrt
            fprintf('\nFound %s file, CV2 partition [%g,%g].',filetyp, f, d)
        else
            [p,n] = fileparts(px);
            if ~a,
                fprintf('\nFound %s file for modality #%g, CV2 partition [%g,%g].',filetyp, i, f, d)
            else
                fprintf('\nFound %s file, CV2 partition [%g,%g].',filetyp, f, d)
            end
            loadstr = sprintf('\nLoading %s file:\n%s',filetyp, n);
            fprintf(loadstr);
            try
                load(px);
            catch
                cprintf('red','\nCould not open file. May be corrupt. Recompute CV2 partition [%g,%g].',f,d);
                fnd = false;
                return
            end
            if exist('GD','var')
                out = GD;
            elseif exist('pGD','var')
                out = pGD;
            elseif exist('mapY','var')
                
                [iy,jy] = size(mapY.Tr);
                
                switch FUSION.flag % Concatenate modality data into single structure
                    case 2
                        if i == 1, 
                            out = mapY; 
                            out.Tr = cell(iy,jy);
                            out.CV = cell(iy,jy);
                            out.Ts = cell(iy,jy);
                            out.VI = cell(iy,jy);
                            
                            if iscell(mapY.Tr{1,1})
                                for k=1:iy
                                    for l=1:jy
                                        out.Tr{k,l} = cell(nclass,1);
                                        out.CV{k,l} = cell(nclass,1);
                                        out.Ts{k,l} = cell(nclass,1);
                                        out.VI{k,l} = cell(nclass,1);
                                    end
                                end
                            end
                            
                        end
                        
                        if i>1, fprintf('\nAdding data of modality #%g to single data matrix.',i); end
                        for k=1:iy
                            for l=1:jy
                                for j=1:nclass
                                    if i>1
                                        % Create mixtures of data shelves,
                                        % if modality concatenation is
                                        % activated
                                        cnt = 1;
                                        nZo = size(out.Tr{k,l}{j},1);
                                        nZp = size(mapY.Tr{k,l}{j},1);
                                        MixCount = nZo * nZp;
                                        tOut.Tr = cell(MixCount,1);
                                        tOut.CV = cell(MixCount,1);
                                        tOut.Ts = cell(MixCount,1);
                                        tOut.VI = cell(MixCount,1);
                                        for zp = 1:nZp
                                            for zo = 1:nZo
                                                tOut.Tr{cnt} = [ out.Tr{k,l}{j}{zo}, mapY.Tr{k,l}{j}{zp} ];
                                                tOut.CV{cnt} = [ out.CV{k,l}{j}{zo}, mapY.CV{k,l}{j}{zp} ];
                                                tOut.Ts{cnt} = [ out.Ts{k,l}{j}{zo}, mapY.Ts{k,l}{j}{zp} ];
                                                tOut.VI{cnt} = [ out.VI{k,l}{j}{zo}; i*ones(size(mapY.Tr{k,l}{j}{zp},2),1) ]; 
                                                cnt = cnt+1;
                                            end
                                        end
                                        out.Tr{k,l}{j} = tOut.Tr;
                                        out.CV{k,l}{j} = tOut.CV;
                                        out.Ts{k,l}{j} = tOut.Ts;
                                        out.VI{k,l}{j} = tOut.VI;
                                        clear tOut;
                                    else
                                        out.Tr{k,l}{j} = mapY.Tr{k,l}{j};
                                        out.CV{k,l}{j} = mapY.CV{k,l}{j};
                                        out.Ts{k,l}{j} = mapY.Ts{k,l}{j};
                                        for m=1:numel(mapY.Tr{k,l}{j})
                                            out.VI{k,l}{j}{m} = ones(size(mapY.Tr{k,l}{j}{m},2),1);
                                        end
                                    end
                                end   
                            end
                        end
                       
                    case {0,1,3,4} % Create structure storing modality information separately for MKL
                        
                        if i == 1, 
                            out = mapY; 
                            out.Tr = cell(iy,jy,nvar);
                            out.CV = cell(iy,jy,nvar);
                            out.Ts = cell(iy,jy,nvar);
                            if size(mapY.Tr{1,1},2) == nclass
                                for g=1:nvar
                                    for k=1:iy
                                        for l=1:jy
                                            out.Tr{k,l,g} = cell(1,nclass);
                                            out.CV{k,l,g} = cell(1,nclass);
                                            out.Ts{k,l,g} = cell(1,nclass);
                                        end
                                    end
                                end
                            end
                        end
                        for k=1:iy
                            for l=1:jy
                                if size(mapY.Tr{k,l},2) == nclass
                                    for j=1:nclass % loop through binary
                                        % Concatenate CV1 training data
                                        out.Tr{k,l,i}(:,j) = mapY.Tr{k,l}(:,j);
                                        % Concatenate CV1 test data
                                        out.CV{k,l,i}(:,j) = mapY.CV{k,l}(:,j);
                                        % Concatenate CV2 test data
                                        out.Ts{k,l,i}(:,j) = mapY.Ts{k,l}(:,j);
                                    end
                                else
                                   out.Tr{k,l,i} = mapY.Tr{k,l};
                                   % Concatenate CV1 test data
                                   out.CV{k,l,i} = mapY.CV{k,l};
                                   % Concatenate CV2 test data
                                   out.Ts{k,l,i} = mapY.Ts{k,l};
                                end
                                
                            end
                        end
                        
                end
            elseif exist('mapYocv','var') 
                [iy,jy] = size(mapYocv.Ts);

                switch FUSION.flag % Concatenate modality data into single structure
                    case 2
                        if i == 1, 
                            out = mapYocv; out.Ts = cell(iy,jy);
                            if iscell(mapYocv.Ts{1,1})
                                for k=1:iy
                                    for l=1:jy
                                        out.Ts{k,l} = cell(nclass,1);
                                    end
                                end
                            end
                        end
                        if i >1, fprintf('\nAdding data of modality #%g to single data matrix.',i); end
                        for k=1:iy
                            for l=1:jy
                                if iscell(mapYocv.Ts{k,l})
                                    for j=1:nclass % loop through binary   
                                        % Concatenate CV2 test data
                                        out.Ts{k,l}{j} = [out.Ts{k,l}{j} mapYocv.Ts{k,l}{j}];
                                    end
                                else
                                   % Concatenate CV2 test data
                                   out.Ts{k,l} = [out.Ts{k,l} mapYocv.Ts{k,l}];
                                end
                            end
                        end
                        
                    case {0,1,3,4} % Create structure storing modality information separately for MKL
                        if i == 1, 
                            out = mapYocv; 
                            out.Ts = cell(iy,jy,nvar);
                        end
                        for k=1:iy
                            for l=1:jy
                                if iscell(mapYocv.Ts{k,l})
                                    for j=1:nclass % loop through binary 
                                        % Concatenate CV2 test data
                                        out.Ts{k,l,i}{j} = mapYocv.Ts{k,l}{j};
                                    end
                                else
                                   % Concatenate CV2 test data
                                   out.Ts{k,l,i} = mapYocv.Ts{k,l};
                                end
                            end
                        end
                        
                end
            end
        end
        if exist('MD','var')
            out2=MD;
        end
        if exist('xpos','var')
            out3=xpos;out4=ypos;
        end
    end
end