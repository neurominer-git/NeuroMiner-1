function NM = DeleteModalityInNM(NM, varind)

NM.datadescriptor(varind)   = [];
NM.Y(varind)                = [];
NM.files(varind)            = [];
NM.brainmask(varind)        = [];
NM.badcoords(varind)        = [];
NM.featnames(varind)        = [];
