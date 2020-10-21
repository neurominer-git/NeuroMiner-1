function idx=findIdx(Data,Con,flag);
% Function findIdx:
% returns all the index entrees whic match a certian condition
% function idx=findIndx(Data,Con);
% Requires Two cell arrays, one with the Data, the other with the condition
% Condition can be vectors, then it's an OR withing this vector
% flag: 'and': vectors are anded 
%       'or' : vectors are ored 
% EXAMPLE:
% findIdx({SN,BN,HCond},{[1 2 3],2,[1 3 2],'or'}): returns data index from subject 1-3, BN=2, all Handconditions
% findIdx({SN,BN,HCond},{[1 2 3],2,[3 2 1],'and'}); : returns data index subject 1, HCond 3, subject 2 HCond 2 etc.

l=length(Data);
if(l~=length(Con));
    fprintf('findIndx:Data and Conditional dont have same dimension!\n');
    return;
end;
andcheck=[];
for i=1:l
    if(length(Con{i})==1)
        T(:,i)=(Data{i}==Con{i});
    else 
        if (nargin<3 | strcmp(flag,'or'))
            idx=(Data{i}==Con{i}(1));
            for j=2:length(Con{i})
                idx=idx| Data{i}==Con{i}(j);
            end;
            T(:,i)=idx;
        elseif (strcmp(flag,'and'))
            T(:,i)=ones(length(Data{1}),1);
            andcheck(end+1)=i;
        end;
    end;
end;
if(~isempty(andcheck))
    lcheck=length(Con{andcheck(1)});
    idx=zeros(length(Data{andcheck(1)}),1);
    for j=1:lcheck
        T2(:,j)=ones(length(Data{andcheck(1)}),1);
        for a=1:length(andcheck)
            T2(:,j)=T2(:,j) & (Data{andcheck(a)}==Con{andcheck(a)}(j));
        end;
        idx=idx | T2(:,j);
    end;
    T(:,andcheck(1))=idx;
end;
if (l>1)
    idx=all(T');
else 
    idx=T;
end;
idx=find(idx==1);