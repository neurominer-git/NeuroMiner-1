function [H, P, PT, PerfT, Chi2, Power, Rx] = nk_CompareClass(Cl, L, Clnames, Fun)

alpha = 0.05;
teststr = 'exact';
altstr = 'less';
[~,n] = size(Cl); pfl = false; unbfl = false;
if ~exist('Clnames','var') || isempty(Clnames) || numel(Clnames) ~= n
    Clnames = cellstr([repmat('Clnames',n,1) num2str((1:n)')]);
end
R = sum(L==1)/sum(L==-1);
if R>1.5 || R<0.5, unbfl= true; altstr = 'unequal'; teststr = 'asymptotic'; end

if ~exist('Fun','var') || isempty(Fun), Fun='mcnemar'; end
H = nan(n,n); P = nan(n,n); Chi2 = nan(n,n); Power = nan(n,n);
if nargout>3, Perf = nan(1,n) ; pfl = true; end
for i=1:n-1
    for j=i+1:n
        C1 = Cl(:,i); C2 = Cl(:,j);
        indn = sum(isnan([C1 C2]),2)>0;
        switch Fun
            case 'mcnemar'
                C1 = C1 ~= L; C2 = C2 ~= L;
                C1(indn) = []; C2(indn) = [];
                m00 = sum(C1 == 0 & C2 == 0);
                m01 = sum(C1 == 0 & C2 == 1);
                m10 = sum(C1 == 1 & C2 == 0);
                m11 = sum(C1 == 1 & C2 == 1);
                X = [m00 m01 ; m10 m11];
                [~, Chi2(i,j), P(i,j), Power(i,j)] = mcnemar(X, alpha);
                Rx(i,j).X = X;
                Rx(i,j).desc = sprintf('%s vs. %s',Clnames{i},Clnames{j});
            case 'testcholdout'
                C1(indn) = []; C2(indn) = []; iL = L; iL(indn) = [];
                switch unbfl
                    case true
                        if R > 1.5,
                            Cost = [ 0 R; 1 0 ];
                        else
                            Cost = [ 0 1; 1/R 0];
                        end
                        ClassNames = categorical([1 -1]);
                        [H(i,j), P(i,j)] = testcholdout(C1,C2, iL, 'Cost', Cost, 'ClassNames', ClassNames, 'CostTest', 'likelihood');
                    case false
                        [H(i,j), P(i,j)] = testcholdout(C1,C2, iL, 'Test', teststr, 'Alternative', altstr);
                end
        end
    end
end

PT = mat2cell(P, ones(1,n), ones(1,n));
PT = cell2table(PT,'Rownames',Clnames,'Variablenames',Clnames);

if pfl
    for i=1:n
        Perf(i) = BAC(L,Cl(:,i));
    end
    PerfT = array2table(Perf, 'Variablenames', Clnames);
end