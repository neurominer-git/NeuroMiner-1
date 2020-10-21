function D=addstruct(D,A,type,force)
% function D=addstruct(D,A,type,force)
% Adds the fields of structure A to the fields of structure D
% adds the field as a row or column depending on type
% last concatinates 2-d along columns, N-d along the last dimension
% type = 'row' / 'column' / 'last'
%   row: (DEFAULT) add as rows
%   column: add as columns
%   last: add as the last dimension
% 'force': Forces all fields to be added in same row/column
%          This means, if a field is not exisitent or shorter, it will
%          be padded with NaNs
%          If a field is no exisitent in the added structure, it will also
%          be padded with NaNs;
% Joern Diedrichsen
% v.1.1 09/18/05: added support for cell arrays
% v.1.2 12/14/05: taken out reference to header
% v.1.3 06/07/06: force option that includes NaN for
%               so far non-existing fields
% v.1.4 06/07/08: Force option corrected to that a missing last field is
%               done
% -------------------------------------------------------
names=fieldnames(A);
Dnames=fieldnames(struct(D));
if (nargin <3 | strcmp(type,'row'))
    dim=1;
elseif (strcmp(type,'column'))
    dim=2;
elseif (strcmp(type,'last'))
    dim=length(size(eval(['A.' names{1} ';'])));
else
    error('unknown option (row/column/last)');
end;

if (nargin >3 & strcmp(force,'force'))
    force=1;
    if (~isempty(Dnames))
        Dlength=size(getfield(D,Dnames{end}),dim);
    else
        Dlength=0;
    end;
else
    force=0;
end;

for (i=1:length(names))
    if (~isfield(D,names{i}) & ~force)
        eval(['D.' names{i} ' = A.' names{i} ';']);
    elseif (~isfield(D,names{i}) & force)
        if (strcmp(type,'row'))
            eval(['Alength= size(A.' names{i} ',2);']);
            eval(['D.' names{i} '= ones(Dlength,Alength).*NaN;']);
        elseif (strcmp(type,'column'))
            eval(['Alength= size(A.' names{i} ',1);']);
            eval(['D.' names{i} '= ones(Alength,Dlength).*NaN;']);
        else
            error ('force does only work with column / row');
        end;
        eval(['D.' names{i} '= cat(dim,D.' names{i} ' ,A.' names{i} ');']);
    else
        eval(['D.' names{i} '= cat(dim,D.' names{i} ' ,A.' names{i} ');']);
    end;
end;

if (force)
    for i=1:length(Dnames)
        if (isempty(strmatch(Dnames{i},names,'exact')))
            if (strcmp(type,'row'))
                eval(['cases= size(A.' names{1} ',1);']);
                eval(['width= size(D.' Dnames{i} ',2);']);
                eval(['D.' Dnames{i} '= cat(dim,D.' Dnames{i} ' ,ones(cases,width).*NaN);']);
            elseif (strcmp(type,'column'))
                eval(['cases= size(A.' names{1} ',2);']);
                eval(['leng= size(D.' Dnames{i} ',1);']);
                eval(['D.' Dnames{i} '= cat(dim,D.' Dnames{i} ' ,ones(leng,cases).*NaN);']);
            else
                error ('force does only work with column / row');
            end;
        end;
    end;
end;