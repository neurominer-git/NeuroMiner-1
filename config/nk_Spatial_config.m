function [ SPATIAL, PX ] = nk_Spatial_config(SPATIAL, PX, defaultsfl, parentstr)
% Compute spatial consistency of discriminative effects

if ~exist('defaultsfl','var') || isempty(defaultsfl),  defaultsfl = 0; end;
if ~exist('PX','var'), PX =  []; end;
cubetype    = 1;
cubefwhm    = 8;
cubevoxres  = 3;

if ~defaultsfl
   
    if ~isempty(SPATIAL) && isfield(SPATIAL,'cubetype'), cubetype = SPATIAL.cubetype; else SPATIAL.cubetype = cubetype; end
    if ~isempty(SPATIAL) && isfield(SPATIAL,'cubefwhm'), cubefwhm = SPATIAL.cubefwhm; else SPATIAL.cubefwhm = cubefwhm; end
    if ~isempty(SPATIAL) && isfield(SPATIAL,'cubevoxres'), cubevoxres = SPATIAL.cubevoxres; else SPATIAL.cubevoxres = cubevoxres; end
    
    switch cubetype
        case 1
            cubetypestr = 'No filtering';
        case 2
            cubetypestr = 'Absolute difference filtering (6 neighbors)';
        case 3
            cubetypestr = 'Cube variance filtering (27 neighbors)';
        case 4
            cubetypestr = 'Gaussian smoothing';
            cubefwhmstr = nk_ConcatParamstr(cubefwhm);
        case 5
            cubetypestr = 'Resampling';
            cubevoxresstr = nk_ConcatParamstr(cubevoxres);
    end
   
    menustr = sprintf('Select spatial operation [ %s ]', cubetypestr ); menuact = 1;
    switch SPATIAL.cubetype 
        case 4
            menustr = sprintf('%s|Specify Gaussian filter width [ %s ]', menustr, cubefwhmstr ); menuact = [ menuact 2 ];
        case 5
            menustr = sprintf('%s|Specify voxel resolution [ %s ]', menustr, cubevoxresstr ); menuact = [ menuact 3 ];
    end
    
    nk_PrintLogo

    mestr = 'Spatial operations'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\n\nYou are here: %s >>> ',parentstr); 
    act = nk_input(mestr,0,'mq', menustr, menuact);
    
    switch act
        case 1
            SPATIAL.cubetype = uint8(nk_input('Select spatial operation',0,'m',...
                ['No filtering|' ...
                'Absolute difference filtering (6 neighbors)|' ...
                'Cube variance filtering (27 neighbors)|' ...
                'Gaussian smoothing|' ... %'Resampling (=>Voxel size)']
                ],1:4,cubetype));
            switch SPATIAL.cubetype 
                case 4
                    SPATIAL.cubefwhm = 8;
                    PX = nk_AddParam(SPATIAL.cubefwhm, 'FWHM', 0, PX, 'replace');
                case 5
                    SPATIAL.cubevoxres = 3;
                    PX = nk_AddParam(SPATIAL.cubefwhm, 'VOX', 0, PX, 'replace');
                otherwise
                    PX = nk_AddParam(SPATIAL.cubefwhm, 'FWHM', 0, PX, 'reset');
            end
        case 2
            SPATIAL.cubefwhm = nk_input('Specify Gaussian filter width range [mm]',0,'e', cubefwhm);
            PX = nk_AddParam(SPATIAL.cubefwhm,'FWHM', 0, PX, 'replace'); 
        case 3
            SPATIAL.cubevoxres =  nk_input('Specify voxel resolution range [mm]',0,'e', cubevoxres);
            PX = nk_AddParam(SPATIAL.cubevoxres,'VOX', 0, PX, 'replace'); 
    end
else
    act = 0;
end    

if act, [ SPATIAL, PX ] = nk_Spatial_config(SPATIAL, PX, [], parentstr); end
