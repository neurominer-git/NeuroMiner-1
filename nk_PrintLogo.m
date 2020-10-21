function nk_PrintLogo(maindlg)
global NMinfo EXPERT SPMAVAIL

clc
if ~SPMAVAIL
    cl = [0.1,0.5,0]; mode = 'non-imaging mode';
else
    cl = NMinfo.cllogo; mode = [];
end
if exist('maindlg','var') && maindlg
    cprintf(cl,'\n\t*******************************************')
    cprintf(cl,'\n\t****\\                                 /****')
    cprintf(cl,'\n\t*****\\     ');cprintf([1 .5 0],'~~~~~~~~~~~~~~~~~~~~~ ');cprintf(cl,'    /*****');
    cprintf(cl,'\n\t******\\     ');cprintf('*red','N E U R O M I N E R ');cprintf(cl,'    /******');
    cprintf(cl,'\n\t*******\\   ');cprintf([1 .5 0],'~~~~~~~~~~~~~~~~~~~~~ ');cprintf(cl,'  /*******');
    cprintf(cl,'\n\t*******/                           \\*******')
    cprintf(cl,'\n\t******/     ');cprintf('*black','pattern recognition ');cprintf(cl,'    \\******');
    cprintf(cl,'\n\t*****/      ');cprintf('*black','for neurodiagnostic ');cprintf(cl,'     \\*****');
    cprintf(cl,'\n\t****/           ');cprintf('*black','applications ');cprintf(cl,'         \\****');
    cprintf(cl,'\n\t***/                                   \\***')
    cprintf(cl,'\n\t*******************************************')
    if ~isempty(mode),cprintf(cl,'\n\t%s',mode);end
else
    cprintf('*red','\t~~~~~~~~~~~~~~~~~~~~~~~ \n');
    cprintf('*red','\t  N E U R O M I N E R \n');
    cprintf('*red','\t~~~~~~~~~~~~~~~~~~~~~~~ ');
end
if EXPERT
    fprintf('\n')
    cprintf('*black','\t>>> EXPERT MODE <<< ')
end
cprintf(cl,'\n\t%s \n', NMinfo.info.ver); fprintf('\n')
cprintf('*black','(c) %s, %s ', NMinfo.info.author, NMinfo.info.datever)
%disp(['<a href="matlab: sendmail(''' NMinfo.info.email ''',''NeuroMiner'')">' NMinfo.info.email '</a>'])
fprintf('\n    nm@pronia.eu \n')