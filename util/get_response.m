function R=get_response()
% waits for user input 
  while(1==1)
       INP=input('RE>:');
        if (isempty(INP))
            R=2;
            break;
        end;
        if (INP==0 | isnan(INP))
            R=0;
            break;
        end;
        if (INP==1)
            R=1;
            break;
        end;      
    end; 