function Res = nk_AnalyzeMatrix(Y, cases, features)
%========================================================================== 
%FORMAT function Res = nk_AnalyzeMatrix(Y)
%==========================================================================
%Checks Y for NaNs or Infs
[m,n] = size(Y);

%% Check for NaN
Res.NaN.Mat          = isnan(Y);
Res.NaN.inY          = any(Res.NaN.Mat(:));
if Res.NaN.inY
    Res.NaN.Rows.sum     = sum(Res.NaN.Mat,2);
    Res.NaN.Rows.ind     = Res.NaN.Rows.sum > 0;
    Res.NaN.Rows.perc    = Res.NaN.Rows.sum * 100 / n;
    if exist('cases','var') && ~isempty(cases)
        Res.NaN.cases = cases(Res.NaN.Rows.ind);
    end
    Res.NaN.Cols.sum     = sum(Res.NaN.Mat);
    Res.NaN.Cols.ind     = Res.NaN.Cols.sum > 0;
    Res.NaN.Cols.perc    = Res.NaN.Cols.sum * 100 / m;
    if exist('features','var') && ~isempty(features)
        Res.NaN.features = features(Res.NaN.Cols.ind);
    end;
end

%% Check for Inf
Res.Inf.Mat          = isinf(Y);
Res.Inf.InY          = any(Res.Inf.Mat(:));
if Res.Inf.InY
    Res.Inf.Rows.sum     = sum(Res.Inf.Mat,2);
    Res.Inf.Rows.ind     = Res.Inf.Rows.sum > 0;
    Res.Inf.Rows.perc    = Res.Inf.Rows.sum * 100 / n;
    if exist('cases','var') && ~isempty(cases)
        Res.Inf.cases = cases(Res.Inf.Rows.ind);
    end
    Res.Inf.Cols.sum     = sum(Res.Inf.Mat);
    Res.Inf.Cols.ind     = Res.Inf.Cols.sum > 0;
    Res.Inf.Cols.perc    = Res.Inf.Cols.sum * 100 / m;
    if exist('features','var') && ~isempty(features)
        Res.Inf.features = features(Res.Inf.Cols.ind);
    end;
end

%% Check for zero-variance rows and columns
Res.Zeros.Rows.ind = ~any(Y,2);
if exist('cases','var') && ~isempty(cases) && any(Res.Zeros.Rows.ind)
    Res.Zeros.cases = cases(Res.Zeros.Rows.ind);
end
Res.Zeros.Cols.ind = ~any(Y);
if exist('features','var') && ~isempty(features) && any(Res.Zeros.Cols.ind)
    Res.Zeros.features = features(Res.Zeros.Cols.ind);
end

%% Analyze variance
Y = nk_PerfScaleObj(Y);
Res.Variance.Rows = nm_nanvar(Y,1,2);
Res.Variance.Cols = nm_nanvar(Y,1);