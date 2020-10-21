function mess = CheckImageInfoConsistency(Vinfo, Vvox, mess)

if ~exist('mess','var'), mess=[]; end

% Check translational 
H = size(unique(Vinfo(:,1:3),'rows'),1);
if H>1, mess = GenerateMessageEntry(mess, 'WARNING: Images do not have the same translation parameters.', 2, 'SystemCommands'); end

% Check rotation
R = size(unique(Vinfo(:,4:6),'rows'),1);
if R>1, mess = GenerateMessageEntry(mess,'WARNING: Images do not have the same rotation parameters.', 2, 'SystemCommands'); end

% Check scaling
S = size(unique(Vinfo(:,7:9),'rows'),1);
if S>1, mess = GenerateMessageEntry(mess,'WARNING: Images do not have the same scaling (zooming) parameters.', 2, 'SystemCommands'); end
%Check affine
A = size(unique(Vinfo(:,10:12),'rows'),1);
if A>1, mess = GenerateMessageEntry(mess,'WARNING: Images do not have the same affine (shearing) parameters.', 2, 'SystemCommands'); end

%Check vox dimensions
V = size(unique(Vvox,'rows'),1);
if V>1, mess = GenerateMessageEntry(mess,'WARNING: Images do not have the same voxel resolution.', 2, 'SystemCommands'); end

