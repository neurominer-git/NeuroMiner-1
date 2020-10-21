function fixedOrder = nk_CheckFixedFeatOrderStackInput(stacking, analyses)

fixedOrder = true;
if stacking == 1 
    for a = 1:numel(analyses)
        if analyses{a}.params.TrainParam.GRD.NodeSelect.mode >1, 
            fixedOrder = false; break
        end
    end
end

