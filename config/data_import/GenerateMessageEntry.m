function mess = GenerateMessageEntry(mess, text, flag, format)

if ~exist('format','var') || isempty(format), format = 'red*'; end
if ~exist('flag','var') || isempty(flag), flag = 1 ; end
e = numel(mess);
if e==1 && isempty(mess), e=0; end 
mess(e+1).format = format;
mess(e+1).flag = flag;
mess(e+1).text = text;

