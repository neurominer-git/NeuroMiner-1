function F = nk_ComputeRegulFunction(C, meanflag, RegulType)

if ~exist('meanflag','var') || isempty(meanflag), meanflag = false; end

[m,n,nclass] = size(C);

if nclass > 1  
   if meanflag
       C = mean(C,3); nclass = 1;
   else
       F = zeros(size(C)); 
   end
end

for curclass = 1 : nclass
    switch RegulType
        case {'cos_invert','sin_invert','sin*cos_invert','sin+cos/2_invert','cos','sin','sin*cos','sin+cos/2'}
           sC = scaledata( C(:,:,curclass), [], -pi, pi);
           msC = mean(sC(:)); 
           cmsC = sC - msC;
           switch RegulType
               case 'cos'
                   Fc = scaledata( cos(cmsC),[], 0, 100);
               case 'sin'
                   Fc = scaledata( sin(cmsC),[], 0, 100);
               case 'sin*cos'
                   Fc = scaledata( sin(cmsC) .* cos(cmsC), [], 0, 100 );
               case 'sin+cos/2'
                   Fc = scaledata( (sin(cmsC) + cos(cmsC)) / 2, [], 0, 100 );
               case 'cos_invert'
                   Fc = scaledata( 1 - cos(cmsC),[], 0, 100);
               case 'sin_invert'
                   Fc = scaledata( 1 - sin(cmsC),[], 0, 100);
               case 'sin*cos_invert'
                   Fc = scaledata( (1 - sin(cmsC)) .*  (1 - cos(cmsC)), [], 0, 100);
               case 'sin+cos/2_invert'
                   Fc = scaledata( (1 - sin(cmsC)) + (1 - cos(cmsC)) ./ 2, [] , 0, 100);
           end
        case 'none'
            Fc = scaledata(C(:,:,curclass), [], 0, 100);
            Fc(~isfinite(Fc)) = 100;
        case 'invert'
            Fc = 100 - scaledata(C(:,:,curclass), [], 0, 100);
            Fc(~isfinite(Fc)) = 0;
    end
    if nclass > 1 
       F(:,:,curclass) = Fc;
    else
       F = Fc;
    end
end
if sum(isnan(F))>0
    fprintf('Prob')
end
end