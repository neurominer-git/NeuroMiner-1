function NM = nk_InitNMwindowColors(NM, cl)

if ~isfield(NM.defs,'JTextArea') || isempty(NM.defs.JTextArea)    
    cmdWinDoc = com.mathworks.mde.cmdwin.CmdWinDocument.getInstance;
    listeners = cmdWinDoc.getDocumentListeners;
    for i=1:numel(listeners)
        listener_str = listeners(i).toString;
        if ~isempty(strfind(listener_str,'JTextArea'))
            NM.defs.JTextArea = listeners(i); 
            break
        end
    end
end

NM.defs.JTextArea.setBackground(java.awt.Color(cl(1),cl(2),cl(3)))
NM.defs.JTextArea.setForeground(java.awt.Color(0,0,0));

   