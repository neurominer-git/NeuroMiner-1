function [h,p,stats] = nk_CompareClassResample(Cl1, Cl2, L, scoretype, alternative, Test, alpha)

modes = {'5x2t', '5x2F', 'rkCV', 'rank', 'rank_cardillo'};

if ~isequal(size(Cl1), size(Cl2))
    error('Score grid of the two models must have same dimensions')
end

if ~exist('alternative','var') || isempty(alternative),
    error('The type of the classifier score has to be defined: ''loss'' or ''skill''');
end

if ~exist('alternative','var') || isempty(alternative),
    alternative = 'unequal';
end

if ~exist('alpha','var') || isempty(alpha) || ~isscalar(alpha) || ~isfinite(alpha),
    alpha = 0.05;
end

[R,K] = size(Cl1);

if ~exist('Test','var') || isempty(Test)
    if R==5 && K==2, 
        mode = modes{1};
    else
        mode = modes{2};
    end
else
    ind = strcmp(modes, Test);
    if any(ind)
        mode = modes{ind};
    else
        error('Requested statistical test not found');
    end
end

if (~exist('L','var') || isempty(L)) && strcmp(mode,'rkCV')
    error('The label vector has to be provided')
end

switch scoretype
    case 'loss'
        % For loss-type scores
        delta = Cl1(:) - Cl2(:);
    case 'skill'
        % For accuracy-type scores
        delta = Cl2(:) - Cl1(:);
end

switch mode
    case '5x2t'        
        mdelta_r = mean(delta,2);
        s2_r = sum(bsxfun(@minus,delta,mdelta_r).^2,2);
        s2 = sum(s2_r);
        t = delta(1,1)/sqrt(s2/5);
        
        switch alternative
            case 'unequal'
                p = 2*tcdf(-abs(t),5);
            case 'less'
                % delta has a large positive value under H1
                p = tcdf(t,5,'upper');
            case 'greater'
                % delta has a large negative value under H1
                p = tcdf(t,5);
        end
        if nargout == 3
           stats.p = p;
           stats.t = t;
           stats.s2 = s2;
        end
    
    case '5x2F'        
        mdelta_r = mean(delta,2);
        s2_r = sum(bsxfun(@minus,delta,mdelta_r).^2,2);
        s2 = sum(s2_r);
        F = sum(delta(:).^2)/(2*s2);
        p = fcdf(F,10,5,'upper'); % computed only for 'unequal' H1
        if nargout == 3
           stats.p = p;
           stats.F = F;
           stats.s2 = s2;
        end
        
    case 'rkCV' 
         % This is the corrected repeated CV paired t-test 
        % which has an adjusted variance component to account for the 
        % dependence of training and test samples in repeated k-fold CV     
        % Define training and test sample sizes
        n = size(L,1); n2 = n/K; n1 = n-n2;         
        m = mean(delta(:));
        s2 = var(delta(:));
        t = m / sqrt( (1/(R*K) + n2/n1 ) * s2);
        switch alternative
            case 'unequal'
                p = 2*tcdf(-abs(t),K);
            case 'less'
                % delta has a large positive value under H1
                p = tcdf(t,K,'upper');
            case 'greater'
                % delta has a large negative value under H1
                p = tcdf(t,K);
        end
        if nargout == 3
            stats.p = p;
            stats.t = t;
            stats.s2 = s2;
        end

    case {'rank','rank_cardillo'}
        Cl1 = Cl1(:); Cl2 = Cl2(:);
        switch alternative
            case 'unequal'
                tail = 'both';
            case 'greater'
                tail = 'right';
            case 'less'
                tail = 'left';
        end
        if license('test','statistics_toolbox') && strcmp(mode,'rank')
            if size(delta,1)>15
                method = 'exact';
            else
                method = 'approximate';
            end
            [p,~, stats] = signrank(Cl1, Cl2, 'alpha', alpha, 'tail', tail, 'method', method);
        else
            stats = wilcoxon(Cl1, Cl2, alpha);
            p = stats.p;
        end
end

h = p<alpha;