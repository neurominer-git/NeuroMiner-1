function set_SubInd(handles, act)

SubInd = true(size(handles.label));

switch act
    
    case 'set'
        SubInd = get(handles.txtSubgroupSet,'Value');
        if ~isempty(SubInd) && ~strcmp(SubInd,'')
            try
                SubInd = evalin('base',SubInd);
            catch
                
            end
            handles.SubInd = SubInd;
        end
    case 'reset'
        
end