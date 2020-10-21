function [act, GridAct] = nk_CVGridSelector(ix, jx, GridAct, proceedfl)

if ~exist('GridAct','var') || isempty(GridAct)
	GridAct = false(ix,jx);
end
if ~exist('proceedfl','var') || isempty(proceedfl)
    proceedfl = true; end

print_grid(ix, jx, GridAct)

menustr = ['Select all CV2 positions|' ...
	'Select single single perm / fold position|' ...
	'Select single perm / fold range|' ...
	'Deselect all CV2 positions|' ...
	'Deselect single single perm / fold position|' ...
	'Deselect single perm / fold range|' ...
	'Toggle all CV2 positions|' ...
	'Toggle single perm / fold position|' ...
	'Toggle perm / fold range|' ...
    'Enter index matrix of CV2 partitions'];

menuact = 1:10;

if sum(GridAct(:)) && proceedfl, 
    menustr = [menustr '|PROCEED >>>']; menuact = [ menuact 11];
end

act = nk_input('Select CV2 partitions',0,'mq',menustr, menuact);

switch act
    
	case 1
		GridAct = true(ix,jx);
    case 2
        selpos = nk_input('Select position at [ perm, fold ]',0,'i',[],2);
        GridAct(selpos(1),selpos(2)) = 1;
	case 3
        selrange = nk_input('Select range [ perm start, perm end, fold start, fold end ]',0,'i',[],4);
        GridAct(selrange(1):selrange(2),selrange(3):selrange(4)) = 1;
	case 4
        GridAct = false(ix,jx);
	case 5
        delpos = nk_input('Deselect position at [ perm, fold ]',0,'i',[],2);
        GridAct(delpos(1),delpos(2)) = 0;
	case 6
        delrange = nk_input('Deselect range [ perm start, perm end, fold start, fold end ]',0,'i',[],4);
        GridAct(delrange(1):delrange(2),delrange(3):delrange(4)) = 0;
	case 7
        GridAct = ~GridAct; 
    case 8
        togpos = nk_input('Toggle position at [ perm, fold ]',0,'i',[],2);
        GridAct(togpos(1),togpos(2)) = ~GridAct(togpos(1),togpos(2));
	case 9
        togrange = nk_input('Toggle range [ perm start, perm end, fold start, fold end ]',0,'i',[],4);
        GridAct(togrange(1):togrange(2),togrange(3):togrange(4)) = ~GridAct(togrange(1):togrange(2),togrange(3):togrange(4));
    case 10
        GridAct = logical(nk_input('Enter index matrix',0,'e',[],[ix jx]));
        
end

% __________________________________
function print_grid(ix, jx, GridAct)

clc
fprintf('******************************\n') 
fprintf('*** CV2 Partition Selector ***\n')
fprintf('******************************\n\n')
fprintf('Select CV2 partitions to be analyzed:\n')
rowindent = 10;
rowindentspc = repmat(' ',1,rowindent);
for i=1:ix
    rowhdr = ['perm ' num2str(i) ':'];
    rowhdrspc = rowindent - length(rowhdr);
    rowhdrspc = repmat(' ',1,rowhdrspc);
	fprintf('\n%s%s',rowhdr,rowhdrspc)
	for j=1:jx
		if j > 1 && j <= jx
			fprintf('-')
		end
		if GridAct(i,j)
			fprintf('*')
		else
			fprintf('o')
		end
    end
%     if i<ix
%         fprintf('\n%s',rowindentspc)
%         for j=1:jx
%             fprintf('|  ')
%         end
%     end
end

fprintf('\n')