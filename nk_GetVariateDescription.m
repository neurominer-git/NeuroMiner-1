function vargout = nk_GetVariateDescription(res, descrip)

if exist('descrip','var') && ~isempty(descrip)
    vargout = descrip;
end

if isfield(res,'Y') && iscell(res.Y)
    vargout.num_variates = length(res.Y);
    vargout.num_variates_str = ['Number of Variates: ' num2str(length(res.Y))];
    if isfield(res,'datadescriptor')
        for i=1:length(res.Y)
            vargout.variate_description{i} = ...
                ['Variate ' num2str(i) ': '];
            switch res.datadescriptor{i}.type 
                case 1
                    vargout.variate_description{i} = ...
                        [ vargout.variate_description{i} 'Volumetric Imaging Data'];
                case 2
                    vargout.variate_description{i} = ...
                        [ vargout.variate_description{i} 'Diffusion Tensor Imaging'];
                case 3
                    vargout.variate_description{i} = ...
                        [ vargout.variate_description{i} 'fMRI'];
                case 4
                    vargout.variate_description{i} = ...
                        [ vargout.variate_description{i} 'Magnetic Resonance Spectroscopy'];
                case 5
                    vargout.variate_description{i} = ...
                        [ vargout.variate_description{i} 'Non-imaging data'];
            end
            vargout.variate_description{i} = ...
                        [ vargout.variate_description{i} ' (' res.datadescriptor{i}.desc ')'];
        end
    end
else
    vargout.num_variates = 1;
    vargout.num_variates_str = 'One variate (no description available)';
end

return