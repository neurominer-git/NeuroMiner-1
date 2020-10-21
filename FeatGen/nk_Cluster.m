function out = nk_Cluster(W, ACT)
global VERBOSE

% Find clusters (connected areas) in an image and return new image containing the
% number for each cluster sorted by size. The largest cluster has number
% one and so on.

%spm_progress_bar('Init',n,'Cluster',' ');

if VERBOSE,
    fprintf('\n\n===============================\n')
    fprintf('\n==     CLUSTER EXTRACTION    ==')
    fprintf('\n===============================\n')
end

% First bring W into 3D format
if ~exist(regexprep(ACT.brainmask,',1',''),'file')
    error(['Space-defining image ' ACT.brainmask ' not found!']);
end

% Read-in space-defining image
Vm = spm_vol(ACT.brainmask);
[Vmx.dims, Vmx.indvol, ~, Vmx.vox] = nk_ReadMaskIndVol(Vm, []);

% Create an empty 3D image
V = zeros(Vmx.dims); out = V;   
% Map label vector into this 3D image
V(Vmx.indvol) = W;              
XYZ = [];

% Compute coordinate map of non-zero label image entries
for j = 1:Vm.dim(3)
    % Get slice
    d = V(:,:,j);  d = d > 0;
    [Qc, Qr] = find(d ~= 0);
    if size(Qc,1)
        XYZ = [XYZ; [Qc Qr j*ones(size(Qc))]];
    end
end

XYZ = XYZ';

% Get non-connected components
A = spm_clusters(XYZ);

% Determine number non-connected components
j = max(A);

sz = zeros(j,1);
if VERBOSE, fprintf('\n%g clusters found.',j); end

% Get size of each component
for i = 1:j
    M = find(A == i);
    sz(i) = size(M,2);	
end

% Remove 1-voxel clusters
%     indg = (sz > 1);
%     A = A(indg);
%     sz = sz(indg);
%     j = max(A);
%     fprintf('\nAfter pruning one-voxel clusters: %g remaining clusters.',sum(indg));
%     
[sz2,i2] = sort(sz);


switch ACT.clustflag

    case 1 % All clusters

        for i = 1:length(i2)
            M    = A == i2(end-i+1);
            ML   = XYZ(:,M);
            for l=1:size(ML,2)
                out(ML(1,l),ML(2,l),ML(3,l)) = i;
            end   
        end

    case 2 % N clusters of sorted list according to cluster extent
        %sf='';
        mn = min(ACT.nclust,j);
        for i = 1:mn
            M    = A == i2(end-i+1);
            ML   = XYZ(:,M);
            numML = sum(M);
            for l=1:size(ML,2)
                out(ML(1,l),ML(2,l),ML(3,l)) = i;
            end
            if VERBOSE, fprintf('\nCluster %g: %g features',i,numML); end
        end
        if VERBOSE,
            fprintf('\n======================================');
            fprintf('\n%g clusters extracted.',mn);
        end
    case 2 % N clusters of sorted list according to cluster mass (not implemented yet)

        % ... to be done

    case 3 % Minimum cluster extent [voxels]
        mn = find( sz2 >= ACT.minclustextent,1,'first');
        cnt=1;
        for i = j:-1:mn
            M    = find(A == i2(i));
            ML    = XYZ(:,M);
            for l=1:size(ML,2)
                out{k}(ML(1,l),ML(2,l),ML(3,l)) = cnt;
            end
            %fprintf('\nextracted cluster #%g: %g voxels',cnt,sz2(i))
            cnt=cnt+1;
        end
        if VERBOSE, fprintf('\n%g clusters extracted using extent threshold of %g voxels.\n', ...
            cnt-1, ACT.minclustextent); end
        %sf = num2str(cnt-1);

    case 4 % Percentile cluster extent

        ... to be done

    case 5 % Percentile cluster mass (not implemented yet)

        ... to be done
end

out = out(Vmx.indvol)';
	
