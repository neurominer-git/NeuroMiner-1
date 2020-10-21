function out = nk_mapY2Struct(mapY, permutefl)

if ~exist('permutefl','var') || isempty(permutefl), permutefl = true; end

global FUSION

switch FUSION.flag 
    
    case {0,1,3}
        
        out = mapY; 
        
    case 2 % Concatenate modality structure into single multi-modal structure
        
        nvar    = numel(mapY);
        out     = mapY{1};
        
        if isfield(mapY{1},'Tr')
            [iy,jy] = size(mapY{1}.Tr);
            nclass  = numel(mapY{1}.Tr{1,1});
            out.Tr = cell(iy,jy);
            out.CV = cell(iy,jy);
            out.Ts = cell(iy,jy);
            out.VI = cell(iy,jy);

            for k=1:iy
                for l=1:jy
                    out.Tr{k,l} = cell(nclass,1);
                    out.CV{k,l} = cell(nclass,1);
                    out.Ts{k,l} = cell(nclass,1);
                    out.VI{k,l} = cell(nclass,1);
                end
            end
            
            for i=1:nvar
                 fprintf('\nAdding data of modality #%g to single data matrix.',i);
                 for k=1:iy
                    for l=1:jy
                        for j=1:nclass
                            if i>1
                                % Create mixtures of data shelves,
                                % if modality concatenation is
                                % activated
                                cnt = 1;
                                nZo = size(out.Tr{k,l}{j},1);
                                nZp = size(mapY{i}.Tr{k,l}{j},1);
                                if permutefl
                                    MixCount = nZo * nZp;
                                    tOut.Tr = cell(MixCount,1);
                                    tOut.CV = cell(MixCount,1);
                                    tOut.Ts = cell(MixCount,1);
                                    tOut.VI = cell(MixCount,1);
                                    for zp = 1:nZp
                                        for zo = 1:nZo
                                            tOut.Tr{cnt} = [ out.Tr{k,l}{j}{zo}, mapY{i}.Tr{k,l}{j}{zp}];
                                            tOut.CV{cnt} = [ out.CV{k,l}{j}{zo}, mapY{i}.CV{k,l}{j}{zp}];
                                            tOut.Ts{cnt} = [ out.Ts{k,l}{j}{zo}, mapY{i}.Ts{k,l}{j}{zp}];
                                            tOut.VI{cnt} = [ out.VI{k,l}{j}{zo}; i*ones(size(mapY{i}.Tr{k,l}{j}{zp},2),1) ]; 
                                            cnt = cnt+1;
                                        end
                                    end
                                    out.Tr{k,l}{j} = tOut.Tr;
                                    out.CV{k,l}{j} = tOut.CV;
                                    out.Ts{k,l}{j} = tOut.Ts;
                                    clear tOut;
                                else
                                    if nZp ~= nZo, 
                                        error('Modalities do not have the equal number of data shelves'); 
                                    end
                                    for zp=1:nZp
                                        out.Tr{k,l}{j}{zp} = [ out.Tr{k,l}{j}{zp} mapY{i}.Tr{k,l}{j}{zp} ];
                                        out.CV{k,l}{j}{zp} = [ out.CV{k,l}{j}{zp} mapY{i}.CV{k,l}{j}{zp} ];
                                        out.Ts{k,l}{j}{zp} = [ out.Ts{k,l}{j}{zp} mapY{i}.Ts{k,l}{j}{zp} ];
                                        out.VI{k,l}{j}{zp} = [ out.VI{k,l}{j}{zp}; i*ones(size(mapY{i}.Tr{k,l}{j}{zp},2),1) ]; 
                                    end
                                end
                            else
                                out.Tr{k,l}{j} = mapY{i}.Tr{k,l}{j};
                                out.CV{k,l}{j} = mapY{i}.CV{k,l}{j};
                                out.Ts{k,l}{j} = mapY{i}.Ts{k,l}{j};
                                for m=1:numel(mapY{i}.Tr{k,l}{j})
                                    out.VI{k,l}{j}{m} = ones(size(mapY{i}.Tr{k,l}{j}{m},2),1);
                                end
                            end
                        end   
                    end
                 end
            end
        else
            [iy,jy] = size(mapY{1}.Ts);
            nclass  = numel(mapY{1}.Ts{1,1});
            out.Ts = cell(iy,jy); for k=1:iy, for l=1:jy, out.Ts{k,l} = cell(nclass,1); end; end
            for i=1:nvar
                 fprintf('\nAdding data of modality #%g to single data matrix.',i);
                 for k=1:iy
                    for l=1:jy
                        for j=1:nclass
                            if i>1
                                 % Create mixtures of data shelves,
                                % if modality concatenation is
                                % activated
                                cnt = 1;
                                nZo = size(out.Ts{k,l}{j},1);
                                nZp = size(mapY{i}.Ts{k,l}{j},1);
                                if permutefl
                                    MixCount = nZo * nZp;
                                    tOut.Ts = cell(MixCount,1);
                                    for zp = 1:nZp, for zo = 1:nZo, tOut.Ts{cnt} = [ out.Ts{k,l}{j}{zo}, mapY{i}.Ts{k,l}{j}{zp}]; cnt = cnt+1; end; end
                                    out.Ts{k,l}{j} = tOut.Ts;
                                    clear tOut;
                                else
                                    if nZp ~= nZo, 
                                        error('Modalities do not have the equal number of data shelves'); 
                                    end
                                    for zp=1:nZp
                                        out.Ts{k,l}{j}{zp} = [ out.Ts{k,l}{j}{zp} mapY{i}.Ts{k,l}{j}{zp} ];
                                    end
                                end
                            else
                                out.Ts{k,l}{j} = mapY{i}.Ts{k,l}{j};
                            end
                        end   
                    end
                 end
             end
        end
end