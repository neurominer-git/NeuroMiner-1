function Data=dload(filename)
% DLOAD: loads a column-oriented ascii-data file into memory
% synopsis
%		Data=dload(filename)
% Description
%		the file filename has to be tab or space delimited all-numeric or character datafile
%		first row has to be a header file with valid variable names.
%       if there is an output argument, dload gives it to a struture with
%       field names as the columns
%       Otherwise it assigns them to the workspace.
%       If transforming to numericals is unsuccessful, it leaves them as
%       strings
% v 1.1: checks now where error in file reading occurs
% (o.white@bangor.ac.uk)

fid=fopen(filename,'r');
if (fid==-1)
    fprintf(['Error: Could not find ' filename '\n']);
    Header=[];
    Data=[];
    return;
end;

Header=fgetl(fid);
H={};
Data=[];
Head=Header;
while length(Head)>0
    [r,Head]=strtok(Head);
    if (~isempty(r))
        H={H{1:end} r};
    end;
end;

try
    A=textread(filename,'%f','headerlines',1);
catch

    A=textread(filename,'%s','headerlines',1);
end;

Indx=[1:length(H):length(A)]';

% assign into variable names
for col=1:length(H)
	if (max(Indx)>length(A))
		warning(['File format does not have rows*columns entries (e.g. empty cells):' filename]);
		disp('Error is detected at row...');
		error_line=getIndexError(filename,length(H));
		disp(error_line)
		Data=[];
        return;
    end;
    if (iscell(A))
	V={A{Indx}}';
        d=str2num(char(V));
        if (~isempty(d))
            V=d;
        end;
    else
        V=A(Indx);
    end;
    name=strrep(H{col},'-','_');
    if (nargout==0)        % no output variable: just assign to global scope as single variables
        assignin('caller',name,V);
    else
        try 
            Data=setfield(Data,name,V);
        catch
            error('problem with field names/n make sure variable names are valid and text file is saved as tab-delimited text');
        end;
    end;
    Indx=Indx+1;
end;
fclose(fid);

function error_line=getIndexError(filename,ncol)
% we know the file can be read
	fid=fopen(filename,'r');

    % skip a line 
    fgetl(fid); 
    % now read each line at a time
	tline=0;
	curr_line=0;
	while tline~=-1
		curr_line=curr_line+1;
		tline=fgetl(fid);
		v=sscanf(tline, '%f', inf);

		if (length(v) ~= ncol)
			error_line=curr_line;
			fclose(fid);
			return;
		end
	end;
