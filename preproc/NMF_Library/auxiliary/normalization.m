function [ Vo ] = normalization( Vo_org, max_gray_level )
%NOMALIZATION ���̊֐��̊T�v�������ɋL�q
%   �ڍא����������ɋL�q

    Vo_max = max(max(Vo_org));
    
    Vo = Vo_org * (max_gray_level / Vo_max);
end

