function O = nk_OcclusionMapping(Y, model, IN)

if ~exist('IN','var') || isempty(IN)
   IN.iter = 100;
   IN.featnum = 100;
end

% Apply model to Y
[rs, ds] = nk_GetTestPerf(

%Initialize occlusion map