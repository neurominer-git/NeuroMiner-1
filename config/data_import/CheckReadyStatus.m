function disallow = CheckReadyStatus(mess)

disallow = false;
if isempty(mess) || ~isstruct(mess), return; end
nM = numel(mess);

for i=1:nM
    if isempty(mess(i).flag), continue; end
    if mess(i).flag == 1; disallow = true; break; end
end
