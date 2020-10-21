function S = nk_DefineMLParamStr(P, Pdesc, curclass)

S = [];
if ~isempty(Pdesc{curclass})
    for n = 1:numel(P)
        if iscell(P), Pn = P{n}; else, Pn=P(n); end
        S = sprintf('%s, %s: %1.6f',S, Pdesc{curclass}{n}, Pn);
    end
    S = S(3:end);
end

end