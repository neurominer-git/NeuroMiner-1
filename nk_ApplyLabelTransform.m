function [ inp ] = nk_ApplyLabelTransform( PREPROC, MODEFL, inp )

[ inp.label, inp.targscale, inp.minLbCV, inp.maxLbCV, ~, inp.PolyFact ] = nk_LabelTransform(PREPROC, MODEFL, inp.labels(:,inp.curlabel));
if isfield(inp,'labelOOCV')
    if inp.targscale, 
        IN.minY = inp.minLbCV; IN.maxY = inp.maxLbCV; 
        inp.labelOOCV = nk_PerfScaleObj(inp.labelOOCV(:, inp.curlabel), IN); 
    end
    if ~isempty(inp.PolyFact), inp.labelOOCV = inp.labelOOCV .^ (1/inp.PolyFact); end 
end



