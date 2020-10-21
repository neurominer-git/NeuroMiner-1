% Usage:  ret = GetEnvironment
%         Get all 'Environment Variables' from Operating System.
%
%         see also: getenv
%
% Release:  V1.03 07-04-2014
%
% Author:   Jan G. de Wilde   email: jan.dewilde@nl.thalesgroup.com
%           Thales Nederland B.V.
% V1.00 19-02-2014  First release.
% V1.01 20-02-2014  Bug fixes for fieldnames containing '(x64)' and '__'
% V1.02 22-02-2014  Support for Windows and Linux
% V1.03 07-04-2014  Handle lines without assignments (no '=')
function ret = getenvironment2
[ err, sys ] = system('set');
if ~err,
    EOL = [ 0 find( sys == 10 ) ];
    for ind = 2 : length( EOL ),
        lijn = sys(EOL(ind-1)+1 : EOL(ind)-1 );
        pos = find( lijn == '=', 1, 'first' );
        
        if isempty( pos ),  continue;   end     % V1.03
        
        fldname = lijn( 1 : pos-1 );
        fldname = strrep( fldname, '(x86)', '_x86' );   % V1.01
        fldname = strrep( fldname, '__', '' );          % V1.01
        if fldname == '_',  continue;   end             % V1.02
        
        try
            ret.( fldname ) = lijn( pos+1 : end );
        catch ME
            fprintf('%s :: %s   (%s)\n', mfilename, ME.message, lijn );
        end
    end
end
address = java.net.InetAddress.getLocalHost;
ret.HOSTADDRESS	= char( address.getHostAddress );
if nargout == 0,
    disp( ret );   clear ret;
end
return