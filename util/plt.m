% plt.m:   An alternative to plot and plotyy (version 11May10)
% Author:  Paul Mennen (paul@mennen.org)
%          Copyright (c) 2010, Paul Mennen

function [Volv1,Volv2] = plt(varargin);

global copia;

Volv1 = 0;

if ~nargin
  eval('evalin(''base'',''pltdef'','''')');
  w = plt('select','base','who');
  lix = length(w);  b = zeros(1,lix);
  opcion = [];
  wstruct = {}; ws=0;
  rnsion = 6;
  for k=1:lix
    a = plt('select','base',w{k});
    if length(w{k})>6 & findstr(w{k},'pltopt') opcion = [opcion a]; end;
    if strcmp(w{k},'TraceIDlen') rnsion = a;  end;
    if isnumeric(a) & length(a)>3 b(k)=1;
    elseif isstruct(a) & numel(a)==1
      f = fieldnames(a);
      for fn=1:length(f)
        g = getfield(a,f{fn});
        if isnumeric(g) & length(g)>3 ws=ws+1; wstruct{ws} = [w{k} '.' f{fn}]; end;
      end;
    end;
  end;
  w = [w(logical(b)); wstruct'];  nov = length(w);
  if nov<2 disp('Not enough vectors to plot'); return; end;
  spix = get(0,'screenpix');
  Sfuente = (196-spix)/10;
  spix = spix/80 - .2;
  espaciV = round(20 * spix);
  espaciH = round(200 * spix);
  set(0,'units','pix');  ssz = get(0,'screensize');
  fct = fix((ssz(4)-110)/espaciV);
  nov1 = nov+3;  fct = ceil(nov1/fct);
  nr = ceil(nov1/fct);
  wid = fct*espaciH+20; hei = nr*espaciV+40;
  fig = figure('menu','none','number','off','name','PLT select','color','black',...
               'Invert','on','pos',[ssz(3)-wid-5 ssz(4)-hei-65 wid hei]);
  axes('xli',[0 fct],'yli',[0 nr]-.5,'pos',[.02/fct 0 1 1],...
       'color','black','xcol','black','ycol','black');
  hb = hei-25;
  set([uicontrol('str','Plot',    'pos',[  5 hb 40 21],'call','plt select plot;');
       uicontrol('str','CloseAll','pos',[ 60 hb 60 21],'call','plt close;');
       uicontrol('str','Help',    'pos',[135 hb 40 21],'call','plt(''helpv'',0);');
      ],'units','nor','fontsi',Sfuente);
  absc=nr-3; x=0;
  set([text(0,absc+1.1,'Right click on desired x vectors');
       text(0,absc+.45, 'Left click on desired y vectors');
       text(0,absc-.2, 'Double click for right hand axis');
      ],'fonta','ita','fontsi',9,'color',[.5 .5 1]);
  txt = zeros(nov,1); sz1=txt; sz2=txt; len=txt; nvec=txt; sxy=txt; yderecho=txt;
  for k=1:nov
    absc = absc-1;  if absc<0 absc=nr-1; x=x+1;  end;
    sz = size(plt('select','base',w{k}));
    sz1(k)  = sz(1);
    sz2(k)  = sz(2);
    len(k)  = max(sz);
    nvec(k) = min(sz);
    txt(k) = text(x,absc,' ');
    set(txt(k),'ButtonD',['plt(''select'',' int2str(k) ');']);
  end;
  set(fig,'use',[txt sz1 sz2 len nvec sxy yderecho]);
  setappdata(fig,'w',w);
  setappdata(fig,'opt',opcion);
  setappdata(fig,'IDlen',rnsion);
  plt('select','text');
  return;
end;

y1 = varargin{1};
if ~ischar(y1) y1='none'; end;
if y1(1)=='g' y1(1)='G';  end;
switch sum(y1)
case 640
  v2 = varargin{2};
  if isnumeric(v2) | v2(1) ~= 'b'
      w = getappdata(gcf,'w');  e = get(gcf,'use');
      len = e(:,4); sxy = e(:,6); nov = length(len);
      xk = []; for k=1:nov if sxy(k)<0 xk = [xk k]; end; end;
  end;
  switch sum(v2)
    case 411, Volv1 = eval(['evalin(''base'',''' varargin{3} ''');']);
    case 453,
      lnxk = length(xk);
      for k=1:nov
        s = w{k};  if length(s)>20 s=s(1:20); end;
        s = strrep(s,'_','\_');
        s = [s '  (' int2str(e(k,2)) ',' int2str(e(k,3)) ')'];
        if sxy(k) s = [s ' \leftarrow ']; s2 = int2str(abs(sxy(k)));  end;
        c = [1 1 1];
        if lnxk
          if lnxk==1 s2 = ''; end;
          if ~sxy(k)
            if ~length(find(len(xk)==len(k))) c = .5*c; end;
          elseif sxy(k)<0                 c = [.9 .3 .3]; s = [s 'x' s2];
          else
            if e(k,7)               c = [ 1 .6 .2]; s = [s 'y' s2 'R'];
                      else                c = [ 1  1  0]; s = [s 'y' s2];
                      end;
          end;
        else if e(k,5) ~= 1 c = .5*c; end;
        end;
        if c(3)==.5 | c(3)==1 bold = 'norm'; else bold = 'bold';  end;
        set(e(k,1),'color',c,'fontw',bold,'str',s);
      end;
    case 447
      if ~max(sxy) return; end;
      opcion       = getappdata(gcf,'opt');
      rnsion = getappdata(gcf,'IDlen');
      opt = [];
      for k=1:length(opcion)
        t = opcion{k};
        if   isnumeric(t) v = num2str(t);
             if length(t)>1 v = ['[' v ']']; end;
        else v = ['''''' t ''''''];
        end;
        opt = [opt v ','];
      end;
      if length(xk)==1 x = w{find(sxy<0)}; else x = 'X axis'; end;
      derecho = [];
      ny = 0;
      nvL = 0; nvR = 0;
      func = ['plt(''''LabelX'''',''''' x ''''','];
      id = '''''TraceID'''',[''''';
      for k=1:nov
        if sxy(k)>0
          nvec = e(k,5);
          x = w{find(sxy==-sxy(k))};  absc = w{k};
          s2 = strrep(absc,'_','\_');
          ss = strrep(absc,'_','');
          nn = rnsion - (nvec>1);  s = [];
          if length(ss) > nn
            sp = findstr(ss,'.');
            if length(sp)
              nn1 = min(sp-1,floor((nn-1)/2));
              nn2 = min(nn-nn1-1,length(ss)-sp);
              nn1 = nn - nn2 - 1;
              ss = ss([1:nn1-1 sp-1:sp+nn2-1 end]);
            else
              ss = [ss(1:nn-1) ss(end)];
            end;
          end;
          while length(ss) < nn ss = [ss ' ']; end;
          st = ss;
          for v = 1:nvec
            if nvec>1 ss = [st '0'+mod(v,10)]; end;
            s = [s ss ''''';'''''];
          end;
          if e(k,7) derecho = [derecho ny+(1:nvec)];  s2R=s2; nvR=nvR+1;  else  s2L=s2; nvL=nvL+1;  end;
          ny = ny + nvec;
          func = [func x ',' absc ','];
          id = [id s];
        end;
      end;
      if nvL==1 func = [func '''''LabelY'''','''''  s2L ''''',']; end;
      if nvR==1 func = [func '''''LabelYr'''',''''' s2R ''''',']; end;
      if nvR>0  func = [func '''''Right'''',' plt('vtoa','%d',derecho) ',']; end;
      if ny>1 func = [func id(1:end-3) '],']; end;
      func = [func opt];
      func = [func(1:end-1) ');'];
      plt('select','base',func);
    otherwise,
      clic = get(gcf,'SelectionT'); clic = clic(1);  sk = sxy(v2);
      if clic=='a' | ~length(xk)
        if sk>=0
             if e(v2,5)==1 e(v2,6) = min(sxy)-1; end;
        else e(v2,6) = 0;
             e(find(sxy==-sk),6)=0;
        end;
      elseif sk>=0
        if clic=='o' & sk sk = 0; end;
        xki = sort(-sxy(xk(find(len(xk)==len(v2)))));
        if length(xki)
          if sk liy = find(xki==sk)+1;
                if liy>length(xki) sk=0; e(v2,7)=0; else sk=xki(liy);  end;
          else  sk = xki(1);
          end;
          e(v2,6) = sk;
          if clic=='o' & sk  e(v2,7)=1; end;
        end;
      end;
      set(gcf,'use',e);
      plt('select','text');
  end;

case 756
    x = varargin{2};
    liy=0;  nv=x;
    if ~x Volv1=''; Volv2=1; return; end;
    while nv>900  liy=liy+1;  nv=.001*nv; end;
    while nv<.9   liy=liy-1;  nv=1000*nv; end;
    if (liy==1 & x<10000) | (liy==-1 & x>=.1 ) liy=0; nv=x; end;
    Volv2 = nv/x;
    if     liy>5  Volv1 = plt('ftoa','%4w',1/Volv2);
    elseif liy<-6 Volv1 = plt('ftoa','%5w',1/Volv2);
    else        qs=reshape('Atto FemtoPico Nano MicroMilli     Kilo Mega Giga Tera Peta ',5,12)';
                Volv1 = deblank(qs(liy+7,:));
    end;
    if length(Volv1) Volv1 = [Volv1 '-']; end;

case 426

  Csena = varargin{2};  val = varargin{3};  val = val(:);
  if length(val) > 1
    if nargin < 4 arado = ' '; else arado = varargin{4}; end;
    Volv1 = '';
    for k=1:length(val) Volv1 = [Volv1 plt('ftoa',Csena,val(k)) arado]; end;
    Volv1 = Volv1(1:end-length(arado));
    return;
  end;
  atir = length(Csena);  vrto = Csena(atir);
  while vrto==' ' atir=atir-1; vrto = Csena(atir); end;
  if vrto=='w' | vrto=='v'
     v = s2i(Csena(2:atir-1));
     absc = abs(val);
     if isnan(val) Volv1 = 'NaN'; return; end;
     if ~val Volv1='0'; if vrto=='v' Volv1=[blanks(v) Volv1]; end; return;
     elseif val<0 vv = v-1;
     else vv = v;
     end;
     if     absc < 1.0e-99  pf = [7 1 2 2];   elseif absc < 1.0e-9   pf = [6 2 2 3];
     elseif absc < 0.01     pf = [5 1 4 4];   elseif absc < 10.0^vv  pf = [1 1 1 1];
     elseif absc < 1.0e10   pf = [4 1 5 5];   elseif absc < 1.0e100  pf = [5 2 3 4];
     else                pf = [6 3 1 3];
     end;
     pf = pf - 1;
     fp = v-pf(1);  if val<0 fp=fp-1; end;
     if fp<0  Volv1='0'; return; end;
     if pf(1) tipo = 'e'; else tipo = 'g'; end;
     if ~fp  pf = pf + [0,1,-1,-1]; end;
     Volv1 = sprintf(sprintf('%%1.%d%c',fp,tipo),val);
     lix = length(Volv1);
     if pf(1) pf = [1:v-pf(2) v+pf(3):v+pf(4)];
     else     pf = 1:min(v+1,lix);  end;
     if max(pf)<=lix Volv1 = Volv1(pf); else Volv1 = '0'; end;
     if vrto=='v'
       zp = v-length(Volv1);
       if isempty(findstr('.',Volv1))
          if zp>0  Volv1 = [Volv1 '.']; end;
       else zp = zp + 1;
       end;
       for k=1:zp  Volv1 = [Volv1 '0']; end;

     end;
  else
     Volv1 = sprintf(Csena,val);
  end;

case 442
  tipo = varargin{2};  v = varargin{3};
  if nargin < 4 arado = ' '; else arado = varargin{4}; end;
  Volv1 = plt('ftoa',tipo,v,arado);
  if length(v)>1 Volv1 = ['[' Volv1 ']']; end;
case 759
  c = datevec(varargin{2});
  if isempty(c) | isnan(c(1)) Volv1 = []; return; end;
  if nargin<3 tipo=0; else tipo=varargin{3}; end;
  frac = sprintf('%4.3f',mod(c(6),1));
  c(6) = floor(c(6));
  if strcmp(frac,'1.000')
    if c(6)==59 frac = '0.999'; else c(6) = c(6)+1;  end;
  end;
  Volv1 = [datestr(datenum(c),tipo) frac(2:end)];
  if strcmp(Volv1(7:9),'-20') Volv1(8:9) = []; end;
case 543
  if varargin{2} f = get(gcbo,'use'); else  f = ''; end;
  if isempty(f)
     if ispc f = 'plt.chm'; else f = 'plt.htm'; end;
  end;
  if isempty(findstr(f,filesep)) 
    foobar = 0;
    if exist('foobar')
         f = feval('which',f);
    else f = fullfile(fileparts(GetExe),f);
    end;
  end;
  if exist(f)
    if strcmpi(f(end-2:end),'chm') dos(['hh ' f ' &']);
    elseif exist('foobar')         feval('web',['file:///' f ],'-browser');
    else                           ibrowse(f);
    end;
  end;

case 335

  k = 2;  a1 = varargin{2};
  if ischar(a1)
       g = gca;
       a = [0 .3 .4];
       a = axes('vis','of','XtickL',' ','YtickL',' ','TickLen',[0 0]','color',a,'xcol',a,'ycol',a);
       pd = {a   [.1 1 .9]   [0 0]   'none'   []     0         ''      1       1};
       Volv1 = text(0,0,'','use',pd,'interp','none','color',[1 1 .4]);
       set(Volv1,'ButtonD',['plt(''pop'',' int2str(Volv1*8192) ',''open'',0);']);
       axes(g);
  else lix = length(a1);  Volv1 = a1;
       if lix==1 if Volv1>8192 Volv1=Volv1/8192; end;
               k=k+1; pd=get(Volv1,'use'); a=pd{1};
       else    paraN = varargin{3};  param = varargin{4};  s = size(param,1);
               if s for k=1:lix plt('pop',Volv1(k),paraN,param(min(k,s),:)); end;
               else for k=1:lix plt('pop',Volv1(k),paraN,[]); end;
               end;
               return;
       end;
  end;
  while k <= length(varargin)
    paraN  = lower(varargin{k});  param = varargin{k+1}; k=k+2;
    switch sum(paraN)
    case 857, if max(param)>1 s='pixels'; else s='normal'; end;
                    if length(param)<2 fp=get(a,'pos');  param=[fp(1:3) param]; end;
                    set(a,'pos',param,'units',s);
    case 734,  lix = length(param);  pd{6} = param;
                    set(Volv1,'str',param{1},'pos',[.08 .5-lix],'use',pd);
                    set(a,'yli',[-lix 0]);
    case 434,     if ~pd{8} return; end;
                    c = get(gcf,'SelectionT'); c = c(1);
                    if c=='n'
                      uh = findobj(pd{5},'vis','on');
                      setappdata(Volv1,'Uhide',uh);
                      set([Volv1; uh],'vis','of');  ch = pd{6};  lix = length(ch);
                      of = pd{3};  p=get(a,'pos'); pt = p(4)+p(2)+of(2);
                      s = get(a,'Units');
                      if s(1)=='n' & pt>1
                        of(2)=of(2)-pt+1; pd{3}=of; set(Volv1,'use',pd);
                      end;
                      set(a,'vis','on','pos',p+[of 0 0]);  axes(a);
                      s = ['plt(''pop'',' int2str(8192*Volv1) ',''index'','];
                      for liy=1:lix 
                          ht = text(.08,.5-liy,ch{liy});
                          set(ht,'interp',pd{4},'color',pd{2},'ButtonD',[s int2str(liy) ');']);
                          if liy==pd{9} set(ht,'fontw','bol'); end;
                      end;
                      cr = get(a,'color');  z = zeros(1,lix-1);
                      x = [z;z+1;z];  z = [z;z;z+NaN];   absc = 1:lix-1;  absc = -[absc;absc;absc];
                      line('z',z(:),'y',absc(:),'x',x(:),'color',cr-.2+.4*(cr<.5));
                    else
                      lix = length(pd{6});  v = pd{9};
                      if v==lix v=1; else v=v+1; end;
                      plt('pop',Volv1,'index',v);
                    end;
    case 410,     pd{5} = param;   set(Volv1,'use',pd);
    case 647,   if length(param)<2 param=[0 param]; end;  pd{3}=param; set(Volv1,'use',pd);
    case 759,  pd{2} = param;   set(Volv1,'use',pd);
    case 748,  set(a,'color',param,'xcol',param,'ycol',param);
    case 658,   pd{4} = param;  set(Volv1,'use',pd,'interp',param);
    case 536,    pd{9} = abs(param);  set(Volv1,'use',pd);  v = get(a,'vis');
                    if v(2)=='n'
                      of = pd{3};
                      if any(of) set(a,'pos',get(a,'pos')-[of 0 0]); end;
                      set(a,'vis','of');  ch = get(a,'child');  ch(find(ch==Volv1)) = [];
                      delete(ch);
                      uh = getappdata(Volv1,'Uhide');
                      set([Volv1; uh],'vis','on');
                      ch = gcf;
                      matVer = version;  mver = (matVer([1 3])-'0') * [10; 1];
                      if mver<75 & matVer(2)=='.'
                        set(ch,'Share','off');
                        set(ch,'Share','on');
                      end;
                    end;
                    kk = pd{6}{abs(param)};  set(Volv1,'str',kk);
                    if param>0 evalRep(pd{7},{'@STR',kk,'@IDX',int2str(param)}); end;
    case 617,   pd{7} = param;  set(Volv1,'use',pd);
    case 320,      Volv1 = pd{6}; xk = pd{9};
                    switch sum(param)
                      case 536,   Volv1 = xk;
                      case 663,  Volv1 = Volv1{xk};
                      case 734,
                      otherwise, disp(['Invalid argv=' param ' in plt(''pop'',H,''get'',argv)']);;
                    end;
    case 615,   pd{8} = param;  set(Volv1,'use',pd);
    otherwise,      set(Volv1,paraN,param);
    end;
  end;

case 422

  k = 2;  a1 = varargin{2};
  if ischar(a1)
       pd = {''     1     1   -1e9  1e9    1    1    ''    ''  '%7w'};
       Volv1 = text(0,0,'1','horiz','cent','color',[1 1 .4],'use',pd);
       set(Volv1,'ButtonD',['plt(''edit'',' int2str(Volv1*8192) ',''click'',0);']);
  else lix = length(a1);  Volv1 = a1;
       if lix==1 if Volv1>8192 Volv1=Volv1/8192; end;
               k=k+1; pd=get(Volv1,'use'); a=pd{1};
       else    paraN = varargin{3};  param = varargin{4};  s = size(param,1);
               if s for k=1:lix plt('edit',Volv1(k),paraN,param(min(k,s),:)); end;
               else for k=1:lix plt('edit',Volv1(k),paraN,[]); end;
               end;
               return;
       end;
  end;
  while k <= length(varargin)
    paraN  = lower(varargin{k});  param = varargin{k+1}; k=k+2;
    switch sum(paraN)
    case 518,
      if param
        c = get(gcf,'CurrentChar');
        if ~length(c) c=666;  end;
        if param==2 c=27; end;
        t = get(Volv1,'str');   u = get(Volv1,'color');
        p = findstr(t,'_');  lix = length(t);   termi = 0;  valuar = 0;
        switch abs(c)
        case 666,
        case 8,    if p>1 t(p-1)=[]; end;
        case 127,  if p==lix t='_';
                   else    t(p+1)=[];
                   end;
        case 127,  if p==lix t='_'; else t(p+1)=[];  end;
        case 28,   if p>1 t(p)=t(p-1); t(p-1)='_'; end;
        case 29,   if p<lix t(p)=t(p+1); t(p+1)='_'; end;
        case 27,   termi = 1;  t = pd{8};
        case 13,   termi = 1;
                   t(p) = [];
                   if pd{7}
                     tt=t;  t = 'error';
                     try, s = eval(tt);
                        if pd{7}==-1 | pd{7}==length(s)
                          if length(s)==1
                            if isreal(s)
                                 s = max(min(s,pd{5}),pd{4});
                                 t = ['  ' plt('ftoa',pd{10},s) '  '];
                            else pd{6} = imag(s);
                                 t = pd{8};
                            end;
                          else  t = ['[' num2str(s) ']'];
                          end;
                          if isreal(s)  pd{3} = s;  valuar = 1; end;
                        end;
                     end;
                   else pd{3}=0; valuar=1;
                   end;
        otherwise, t = strrep(t,'_',[c '_']);
        end;
        set(Volv1,'str',t);
        if termi
          set(Volv1,'use',pd,'color',u([3 1 2]),'erase','norm','interp',pd{9});
          set(gcf,'keypress','');
          if valuar evalRep(pd{1},{'@VAL',t,'@OBJ',int2str(8192*Volv1)}); end;
        end;
      else
        kp = get(gcf,'keypress');  f = findstr(kp,'click');
        if length(f)==1  kp(f+7) = '2';  eval(kp); return; end;
        if ~pd{2} return; end;
        c = get(gcf,'SelectionT'); c = c(1);
        if c ~= 'a'
          v = pd{3}; w = pd{6};
          if pd{7}==1 & abs(round(v/w)*w-v) < 1e-9
            cp = get(gca,'curre'); cp = cp(1,1);
            imX = get(gca,'xli');  rp = get(Volv1,'pos');
            cp = (cp-imX(1))/diff(imX) > rp(1);
            cp = w * (2*cp-1);
            s = v + cp;
            if s>pd{5} s=pd{4}; elseif s<pd{4} s=pd{5}; end;
            pd{3} = s;  s = ['  ' plt('ftoa',pd{10},s) '  '];
            set(Volv1,'str',s,'use',pd);
            evalRep(pd{1},{'@VAL',s,'@OBJ',int2str(8192*Volv1)});
            return;
          else c='a';
          end;
        end;
        if c=='a'
          t = get(Volv1,'str');
          if ~strcmp(t,'error') pd{8} = t; end;
          t = deblank(t);  while(t(1)==' ') t(1)=[]; end;
          viehoc = get(Volv1,'color');
          pd{9} = get(Volv1,'interp');
          set(Volv1,'interp','none','str',[t '_'],'color',viehoc([2 3 1]),'erase','xor','use',pd);
          set(gcf,'keypress',['plt(''edit'',' int2str(Volv1*8192) ',''click'',1);']);
        end;
      end;
    case 320,     Volv1 = pd{3};
    case 541,   pd{3}  = param;  set(Volv1,'str',['  ' plt('ftoa',pd{10},param) '  '],'use',pd);
    case 617,  pd{1} = param;  set(Volv1,'use',pd);
    case 324,     pd{4}    = param;  set(Volv1,'use',pd);
    case 326,     pd{5}    = param;  set(Volv1,'use',pd);
    case 428,    pd{6}   = param;  set(Volv1,'use',pd);
    case 642,  pd{7}    = param;  set(Volv1,'use',pd);
    case 615,  pd{2}    = param;  set(Volv1,'use',pd);
    case 649,  pd{10}    = param;  set(Volv1,'use',pd);
    otherwise,     set(Volv1,paraN,param);
    end;
  end;

case 643
  Hin = varargin{2};
  if length(Hin)>1
     if nargin>2  en1 = varargin{3}; else en1 = [50 0 100]; end;
     if nargin>3  en2 = varargin{4}; else en2 = '';    end;
     if nargin>4  en3 = varargin{5}; else en3 = '';    end;
     if nargin>5  en4 = varargin{6}; else en4 = 1;     end;
     if nargin>6  en5 = varargin{7}; else en5 = '%6w'; end;
     if nargin>7  en6 = varargin{8}; else en6 = [.75 .75 .75; 0 1 1]; end;

     if length(en1)<5 en1 = [en1 -1e99 1e99]; end;
     if ischar(en5) & length(en5(:,1))==1 en5 = {'%2w'; en5; '%3w'}; end;
     if iscell(en5) en5 = char(en5); end;
     if length(en6(:,1))<3 en6 = [en6; 0 0 0; 0 0 0]; end; 

     t = get(findobj('tag','PltSlider'),'Val');
     if iscell(t) Volv1 = 1+max([t{:}]);
     else Volv1 = 1+length(t);
     end;
     hs = zeros(1,8);
     hs(2)  = uicontrol('sty','text','str',en2,'tag','PltSlider','Val',Volv1,'horiz','cent');
     if isempty(en2) set(hs(2),'vis','off'); end;
     hs(4) = uicontrol('sty','text','str','a','use',en5,'horiz','left');
     hs(5) = uicontrol('sty','text','str','a','use',en3,'horiz','right');
     hs(3)  = uicontrol('sty','slid','use',en1(1),'call',sprintf('plt(''slider'',%d,''CBslide'');',Volv1));
     hs(6)  = uicontrol('sty','edit','backg',en6(2,:),'foreg',en6(4,:),'use',en1(4:5),'horiz','cent',...
                      'call',sprintf('plt(''slider'',%d,''CBedit'');',Volv1));
     set(hs(2:5),'backg',en6(1,:),'foreg',en6(3,:));
     en4 = [en4 10];  hs([7 8]) = en4(1:2);
     set(hs(2),'use',hs);
     plt('slider',Volv1,'set','minmax',en1(2:5),en1(1));
     plt('slider',Volv1,'set','position',Hin);
     return;
  end;

  hs = findobj('tag','PltSlider','Val',Hin);
  if length(hs) ~= 1 disp('Bad PltSlider handle'); return; end;
  hs   = get(hs,'use');
  fmat = get(hs(4),'use');

  accion = varargin{3};
  if nargin>3  en0 = varargin{4};  else en0 = ''; end;
  if nargin>4  en1 = varargin{5};  else en1 = []; end;
  if nargin>5  en2 = varargin{6};  else en2 = []; end;

  switch sum(accion)
  case 662
     ov = get(hs(3),'use');    nv = get(hs(3),'Val');
     smax   = get(hs(3),'max');   smin   = get(hs(3),'min');
     dv = nv-ov;  sdv = sign(dv);
     switch hs(7)
     case 4,
        v = hs(8)*round(nv/hs(8));
        if v==ov v = v + hs(8)*sdv; end;
        nv = v;
     case 2,
        v = round(nv);
        if v==ov v = v + sdv; end;
        nv = v;
     case 5,
        paso = exp(log(smax/smin)*abs(dv)/(smax-smin));
        if dv>0 nv=ov*paso; elseif dv<0 nv=ov/paso; end;
     case 3,
        if nv >= get(hs(3),'use')-eps  nv=2^nextpow2(nv);
        else                                  nv=2^nextpow2(nv/2);
        end;
     case 6,
        pdv = abs(dv)/(smax-smin);
        if pdv > .11     nv = hs(8)*round(nv/hs(8));
        elseif pdv > .02 nv = ov + hs(8)*sdv;
        else             nv = ov + sdv;
        end;
     end;
     if nv<smin nv = smin; end;
     if nv>smax nv = smax; end;
     set(hs(6),'str',plt('ftoa',fmat(2,:),nv));
     set(hs(3),'Val',nv,'use',nv);
     evalRep(get(hs(5),'use'),{'@VAL',sprintf('%g',nv)});

  case 555
     nv = s2d(get(hs(6),'str'));
     if isempty(nv)
         set(hs(6),'str',plt('ftoa',fmat(2,:),get(hs(3),'Val')));
     else
        switch hs(7)
        case 2,  nv = round(nv);
        case 3,  nv = 2 ^ nextpow2(nv/1.414);
        case 4, nv=hs(8)*round(nv/hs(8));
        end;
        minmax = get(hs(6),'use');
        nv = min(max(nv,minmax(1)),minmax(2));
        Volv1 = nv;
        set(hs(6),'str',plt('ftoa',fmat(2,:),nv));
        nv = max(min(get(hs(3),'max'),nv),get(hs(3),'min'));
        set(hs(3),'Val',nv,'use',nv);
        evalRep(get(hs(5),'use'),{'@VAL',sprintf('%g',nv)});
     end;

  case 320
     switch sum(en0)
     case 338, Volv1 = sum(get(hs(2),'vis')) == 221;
     case 308, Volv1 = sum(get(hs(2),'ena')) == 221;
     case 315, Volv1 = hs(2:6);
     otherwise Volv1 = s2d(get(hs(6),'str'));
     end;
  case 332
     if isnumeric(en0) en1 = en0;  en0 = 'value'; end;
     switch sum(en0)
     case 541
        set(hs(6),'str',num2str(en1));
        Volv1 = plt('slider',get(hs(2),'Val'),'CBedit');
     case 323
        set(hs(6),'str',num2str(en1));
        sv = get(hs(5),'use');  set(hs(5),'use','');
        Volv1 = plt('slider',get(hs(2),'Val'),'CBedit');
        set(hs(5),'use',sv);
     case 495, set(hs(2:6),'vis','on')
                   if isempty(get(hs(2),'str')) set(hs(2),'vis','of'); end;
     case 557, set(hs(2:6),'vis','of');
     case 527, set(hs([3 6]),'ena','of');
     case 465,  set(hs([3 6]),'ena','on');
     case 885,
        a = .1*get(0,'screenp')-9.6;
        p = get(gcf,'pos');
        if length(en1<3) en1 = [en1 .13]; end;
        if en1(1)<1 en1(1) = en1(1)*p(3); end;
        if en1(2)<1 en1(2) = en1(2)*p(4); end;
        if en1(3)<1 en1(3) = en1(3)*p(3); end;
        lmax = 7 + 7*fix(s2d(fmat(3,2:(length(deblank(fmat(3,:)))-1))));
        lmin = 7 + 7*fix(s2d(fmat(1,2:(length(deblank(fmat(1,:)))-1))));
        if isnan(lmax) lmax = 28; end;
        if isnan(lmin) lmin = 14; end;
        ha = a + 17;       ya = en1(2) - a - 17;
        set(hs(2:6),'units','pix',{'pos'},...
               {[en1(1:3) 16+a/2];
                [en1(1),en1(2)-33,en1(3),17-a-1];
                [en1(1),ya,lmin,ha];
                [en1(3)+en1(1)-lmax,ya,lmax,ha];
                [en1(1)+lmin,ya-1,en1(3)-(lmin+lmax),ha+2]},...
            'units','nor');
     case 421
        hs(7) = en1(1);
        if length(en1)>1 hs(8) = en1(2); end;
        set(hs(2),'use',hs);
        Volv1 = plt('slider',get(hs(2),'Val'),'CBedit');
     case 650
        if length(en1)<4 en1 = [en1 -1e99 1e99]; end;
        set(hs(5),'str',[plt('ftoa',fmat(3,:),en1(2)) ' ']);
        set(hs(4),'str',[' ' plt('ftoa',fmat(1,:),en1(1))]);
        set(hs(6),'use',en1(3:4));
        v = s2d(get(hs(6),'str'));
        if ~isempty(en2) v = en2; end;
        v = min(max(v,en1(3)),en1(4));
        set(hs(6),'str',plt('ftoa',fmat(2,:),v));
        set(hs(3),'min',en1(1),'max',en1(2),'Val',min(max(en1(1),v),en1(2)));
        Volv1 = s2d(get(hs(6),'str'));
     case 512,  set(hs(2),'str',en1);
                   if isempty(en1) set(hs(2),'vis','of'); end;
     end;
  otherwise disp(sprintf('Invalid Action in plt(''slider''), Action=%s',accion));
  end;

case 390
  v2 = varargin{2};
  if nargin==2 accion='update'; else  accion = varargin{3};  end;
  if sum(accion) == 436
    if nargin<4 grdclr = [.3,.3,.3]; else grdclr = varargin{4}; end;
    if nargin<5 erMode = 'xor'; else erMode = varargin{5}; end;
    axes(v2);  f = findobj(v2,'type','line');
    if length(f) f=get(f(1),'ButtonD'); end;
    rejaH = line('color',grdclr,'Tag','SkipCur','ButtonD',f,'LineStyle',':','use','grid','erase',erMode);
    accion = 'on';
    if nargin >= 6 set(rejaH,varargin{6},varargin{7}); end;
    if nargin >= 8 set(rejaH,varargin{8},varargin{9}); end;
    set(v2,'child',flipud(get(v2,'child')));
  else  rejaH = findobj(v2,'use','grid');  z = get(rejaH,'z');
        Volv1 = length(z);
        if Volv1 Volv1 = z(1) | Volv1>=2; end;
  end;

  switch sum(accion)
  case 643, if Volv1 accion='on'; end;
                fixMark;
  case 642, if Volv1 accion='off'; else accion='on'; end;
  end;

  x=0; absc=0;
  switch sum(accion)
  case 221,  z=1; Volv1=1;
             l = [get(v2,'xli') get(v2,'yli')];
             ckY = [l(3) get(v2,'YTICK') l(4)];  liy = length(ckY);
             ckX = [l(1) get(v2,'XTICK') l(2)];  lix = length(ckX);
             if ckY(2) <= ckY(1)   ckY = ckY(2:liy);   liy=liy-1; end;
             if ckY(liy) <= ckY(liy-1) liy=liy-1;  ckY = ckY(1:liy);  end;
             if ckX(2) <= ckX(1)   ckX = ckX(2:lix);   lix=lix-1; end;
             if ckX(lix) <= ckX(lix-1) lix=lix-1;  ckX = ckX(1:lix);  end;
             t = liy + lix - 4;
             if t s = ones(3,1);  absc = [1 liy liy];  x = [1 lix lix];  lix=lix-1;  liy=liy-1;
                  z = [0;0;NaN] * ones(1,t);           z=z(:);
                  absc = [s*ckY(2:liy)  ckY(absc)'*ones(1,lix-1)]; absc=absc(:);
                  x = [ckX(x)'*ones(1,liy-1) s*ckX(2:lix)];  x=x(:);
             end;
  case 315, z=0; Volv1=0;
  end;
  set(rejaH,'z',z,'y',absc,'x',x);

case 436
  accion = varargin{2};
  if nargin > 2  en0 = varargin{3};
                 if nargin > 3  en1 = varargin{4}; end;
  end;
  switch sum(accion)
  case 436
     copia(12)= en0;
     set(0,'units','pix');  sz = get(0,'screens');
     szw = sz(3) - 350 - 4;
     ppos  = get(gcf,'pos');
     if ppos(4)<1 sp = sz(3:4); ppos = ppos .* [sp sp]; end;
     xp = min(ppos(1)+ppos(3)+6,szw);
     yp = ppos(2)+ppos(4)-105-30*length(findobj('type','fig','use',gcf));
     if nargin ==4  renderer = en1; else renderer = '-painters';  end;
     copia(1)=figure('menu','none','number','off','Back','off','resize','off','pos',[xp yp 350 105],'color',[0,.4,.4],'name','Hardcopy','use',gcf);
     foobar = 0;
     if exist('foobar')
          cpath = feval('which','plt.m');
     else cpath = GetExe;
     end;
     filen = [fileparts(cpath) filesep 'pltHcpy.mat'];
     f = fopen(filen);
     if f>=0 fclose(f); load(filen);
     else
       z = [0 0];  t = z+1;  w = [0 1];
       pS1=[z;z;[2 1];w;t;w;t;z;t;w;z];
       pS2 = 'nofile.000';
     end;
     copia(3)   = uicontrol('sty','pop','str',['Win Meta|Bit Map|HPGL|LaserJet IIp|Post Script|Encaps PS|Windows'],...
                             'call','plt(''hcpy'',''ModePUcb'');');
     copia(5)  = uicontrol('sty','radio','str','color','call','plt(''hcpy'',''rbg1'',''C'');');
     copia(4)     = uicontrol('sty','radio','str','B&W','call','plt(''hcpy'',''rbg1'',''BW'');');
     copia(6) = uicontrol('sty','radio','str','Clip Board','call','plt(''hcpy'',''rbg2'',''cb'');');
     copia(8)    = uicontrol('sty','radio','str','Device','call','plt(''hcpy'',''rbg2'',''dev'');');
     copia(7)   = uicontrol('sty','radio','str','File','call','plt(''hcpy'',''rbg2'',''file'');');
     copia(11) = uicontrol('sty','text','str',pS2,'horiz','left');
     copia(10)  = uicontrol('sty','text','str', 'Path/File:','horiz','left');
     copia(9)  = uicontrol('str','Print','call','plt(''hcpy'',''print'');','use',renderer);
     copia(2)   = uicontrol('str','Select File','call',['plt(''hcpy'',''OpenPBcb'');']);
     cncl            = uicontrol('str','Cancel','call','close(gcf)');
     h1 = copia([3 5 4 6 8 7]);
     h2 = copia([11 10]);
     h3 = [copia([9 2]) cncl];
     set(h1,{'Val'},{pS1(3,1); pS1(5,1); pS1(4,1); pS1(6,1); pS1(8,1); pS1(7,1)});
     set([h1 h2],'backg',[.8,.8,.9],'foreg','black');
     set([h1 h2 h3],{'pos'},{[115 80 100 15]; [ 10 85 65 15]; [10 70  65 15]; [245 85 95 15];
                           [245 55  95 15]; [245 70 95 15]; [10 35 330 15]; [ 10 50 65 15];
                           [110  5  55 20]; [ 10  5 85 20]; [180 5  55 20]});

     for kk=1+1:12-1
       if pS1(kk,2)==1  set(copia(kk),'ena','on'); else set(copia(kk),'ena','of'); end;
     end;

     if get(copia(8),'Val') set(copia(2),'str','Select Dev');
     else                       set(copia(2),'str','Select File');
     end;
     set(copia(1),'vis','on');
  case 364
     switch sum(en0)
     case 67,  set(copia(5),'Val',1); set(copia(4),'Val',0);
     case 153, set(copia(4),'Val',1);    set(copia(5),'Val',0);
     end;
  case 365
     switch sum(en0)
     case 319,  set(copia(6),'Val',0);  set(copia(7),'Val',0);
                 set(copia(9),'ena','on'); set(copia(11),'ena','of');
                 set(copia(8),'Val',1);      set(copia(2),'str','Select Dev','ena','on');
     case 197,   set(copia(6),'Val',1);   set(copia(7),'Val',0);
                 set(copia(9),'ena','on'); set(copia(11),'ena','of');
                 set(copia(8),'Val',0);     set(copia(2),'ena','of');
     case 416, set(copia(6),'Val',0);  set(copia(7),'Val',1);
                 set(copia(8),'Val',0);     set(copia(2),'str','Select File','ena','on');
                 cher = get(copia(11),'str');
                 if length(cher)<5 cher='none'; set(copia(11),'str',cher); end;
                 if strcmp(cher(1:4),'none') set(copia(9),'ena','of'); end;
                 set(copia(11),'ena','on');
     end;
  case 751
     en = [1 1 0 1 1; 1 1 0 1 0; 1 1 0 1 1; 1 1 0 0 1; 0 1 0 1 1; 0 1 0 1 1; 0 0 1 1 1];
     va = [0 1 0 1 0; 0 1 0 1 0; 0 1 0 0 1; 0 1 0 0 1; 0 1 0 1 0; 0 1 0 1 0; 0 0 1 1 0];
     s = get(copia(3),'Val');
     ipo=reshape('WMFBMPHGLJETPS EPSxxxxxx',3,8)';
     if s==7  set(copia(2),'str','Select Dev', 'ena','on'); set(copia(11),'ena','of');
     else     set(copia(2),'str','Select File','ena','on'); set(copia(11),'ena','on');
              cher = get(copia(11),'str');
              if length(cher)<5 cher='none'; set(copia(11),'str',cher); end;
              nfn=length(cher);
              if cher(nfn-3) == '.'
                   set(copia(11),'str',[cher(1:nfn-3) ipo(s,:)]);
              else set(copia(11),'str',[cher ipo(s,:)]);
              end;
     end;
     ena = {'off' 'on'};
     set(copia(6),'ena',ena{1+en(s,1)},'Val',va(s,1));
     set(copia(7),  'ena',ena{1+en(s,2)},'Val',va(s,2));
     set(copia(8),   'ena',ena{1+en(s,3)},'Val',va(s,3));
     set(copia(5), 'ena',ena{1+en(s,4)},'Val',va(s,4));
     set(copia(4),    'ena',ena{1+en(s,5)},'Val',va(s,5));
  case 745
     if strcmp('Select File',get(copia(2),'str'))
        ipo=reshape('WMFBMPHGLJETPS EPSxxxxxx',3,8)';
        [cher,pathN]= uiputfile(['*.',ipo(get(copia(3),'Val'),:)],'Open (new) File');
        if cher
           fi = cher;
           if isempty(findstr(',',fi)) fi = [fi '.' ipo(get(copia(3),'Val'),:)]; end;
           set(copia(11),'str',[pathN fi]);  set(copia(9),'ena','on');
        end;
     elseif strcmp('Select Dev',get(copia(2),'str')) print -dsetup
     else   disp('Invalid OpenPBcb string');
     end;
  case 557
     trvS= get(copia(12),'invert');
     set(copia(12),'invert','off');
     figure(copia(12)); drawnow;
     PrintMode = get(copia(3),'Val');
     colorFlg  = get(copia(5),'Val');
     nirS = sprintf('print -f%d',copia(12));
     renderer = get(copia(9),'use');
     opcion={[' ' renderer ' -dmeta '];' -dbitmap ';' -dhpgl ';' -dljet2p '};
     if colorFlg opcion = [opcion; {' -dpsc ';' -depsc ';' -dwinc '}];
     else        opcion = [opcion; {' -dps ' ;' -deps  ';' -dwin  '}];
     end;
     nirS=[nirS opcion{PrintMode,:} ' '];
     if get(copia(7),'Val')==1
        PathFilen = get(copia(11),'str');  nirS=[nirS '''' PathFilen ''''];
     end;
     nirS = [nirS ' -noui'];  axh = [];  axh = findobj(copia(12),'type','axes');  nvx = length(axh);
     if ~colorFlg
        figC=get(copia(12),'color');
        nino = [];   axCol = zeros(nvx,3);  imar = zeros(nvx,1);
        for kk=1:nvx
           x = get(axh(kk),'color');
           if strcmp('none',x) imar(kk) = 1; else axCol(kk,:)  = x; end;
           nino = [nino; get(axh(kk),'child')];
           xCol(kk,:) = get(axh(kk),'xcol');  yCol(kk,:) = get(axh(kk),'ycol');
           tCol(kk,:) = get(get(axh(kk),'title'),'color');
           if xCol(kk,:) == figC  imar(kk) = imar(kk) + 2; end;
        end;
        nnino = length(nino);
        for kk=1:nnino
           if sum(get(nino(kk),'type'))==528   kidCol(kk,:) = get(nino(kk),'facecolor');
                                                    set(nino(kk),'facecolor',[.25 .25 .25]);
           else kidCol(kk,:) = get(nino(kk),'color'); set(nino(kk),'color','black');
           end;
        end;
        set(copia(12),'color','white');
        for kk= 1:nvx
           if imar(kk)==1 | imar(kk)==3
           else set(axh(kk),'color','white');
           end;
           if imar(kk)>=2 set(axh(kk),'xcol','white','ycol','white' );
           else            set(axh(kk),'xcol','black','ycol','black' );
           end;
        end;
        for kk=1:nvx
           if get(get(axh(kk),'title'),'color')==figC set(get(axh(kk),'title'),'color','white');
           else                                    set(get(axh(kk),'title'),'color','black');
           end;
           if imar(kk) >=2
                set(get(axh(kk),'xlabel'),'color','white'); set(get(axh(kk),'ylabel'),'color','white');
           else set(get(axh(kk),'xlabel'),'color','black'); set(get(axh(kk),'ylabel'),'color','black');
           end;
        end;
     end;
     set(copia(1),'pointer','watch'); drawnow; pause(1);
     if PrintMode == 2
        set(copia(1),'vis','of');
        refresh
        eval(nirS);
        set(copia(1),'vis','on');
     else  drawnow discard;  eval(nirS);
     end;
     if ~colorFlg & (PrintMode ~=4 | PrintMode ~=7)
        set(copia(12),'color',figC);
        for kk=1:nvx
           if ~rem(imar(kk),2) set(axh(kk),'color',axCol(kk,:)); end;
           set(axh(kk),'xcol',xCol(kk,:)); set(axh(kk),'ycol',yCol(kk,:));
           set(get(axh(kk),'title'),'color',tCol(kk,:));
        end;
        for kk=1:nnino
           if sum(get(nino(kk),'type'))==528 set(nino(kk),'facecolor',kidCol(kk,:));
           else                                   set(nino(kk),'color',kidCol(kk,:));
           end;
        end;
        if PrintMode==2 drawnow; else drawnow discard; end;
     end;
     pS1=zeros(12-1,2);
     for kk=1+1:12-1
        pS1(kk,1) = get(copia(kk),'Val');
        if sum(get(copia(kk),'ena'))==221 pS1(kk,2)=1; end;
     end;
     pS2=get(copia(11),'str');
     foobar = 0;
     if exist('foobar')
          cpath = feval('which','plt.m');
     else cpath = GetExe;
     end;
     save([fileparts(cpath) filesep 'pltHcpy.mat'],'pS1','pS2');
     close(copia(1));
     set(copia(12),'invert',trvS);
     clear copia
  otherwise disp([accion ' invalid Action in plt(hcpy)']);
  end;

case 670

    if nargin<3 disp('Not enough arguments plt(''cursor'',...)'); return;  end;
    cador = varargin{2};
    accion = varargin{3};
    if nargin > 3
      en0 = varargin{4};
      if nargin > 4
        en1 = varargin{5};
        if nargin > 5
          en2 = varargin{6};  if iscell(en2) en2 = char(en2); end;
          if nargin > 6
            en3 = varargin{7};
            if nargin > 7
              en4 = varargin{8};
              if nargin > 8
                en5 = varargin{9};  if iscell(en5) en5 = char(en5); end;
                if nargin > 9
                  en6 = varargin{10};
                  if nargin > 10
                    en7 = varargin{11};
                    if nargin > 11 en8 = varargin{12}; end;
                  end;
                end;
              end;
            end;
          end;
        end;
      end;
    end;

    CurPri = getappdata(0,'CurMain');
    sact = sum(accion);
    if sact ~= 436
      if isempty(cador) return; end;
      Hc = get(CurPri(cador),'use');
      v2 = Hc(13);  ax2 = Hc(14);
      aconex = 0;
      if ax2
           if sum(get(ax2,'vis')) == 315  ax2=0;
           else aconex = [get(get(ax2,'ylabel'),'str') '  '];
                aconex = aconex(1)~=92 | aconex(2) ~= 'd';
           end;
      end;
      nw1 = Hc(4);  nw2 = Hc(6);
      tm = get(nw1,'use');
      actv = tm(4);
      iact = 14 + actv;
      hact = Hc(iact);
      LMda = [get(v2,'xli') get(v2,'yli')];
      xyli = LMda;
      if get(hact,'par')==ax2 xyli(3:4) = get(ax2,'yli'); end;
      DuaC = getappdata(v2,'DualCur');
    end;
    nargin2=nargin;
    switch sact
    case 436
      v2 = cador(1);
      DuaC = getappdata(v2,'DualCur');
      if isempty(DuaC) DuaC = 0;  setappdata(v2,'DualCur',0); end;
      Hc(13) = v2;
      if length(cador)==2  ax2 = cador(2); else ax2 = 0; end;
      Hc(14) = ax2;
      hf = get(v2,'par');
      if nargin2>=10 vis=en6; else vis = 'on'; end;
      Cn = length(CurPri);
      f = find(~CurPri);
      if length(f) cador = f(1);
      else CurPri = [CurPri 0]; cador = Cn+1;
      end;
      if cador>200 disp('Warning: Clear actions missing for 200 cursor inits'); end;
      Volv1 = cador;
      rasHT = sprintf('plt(''cursor'',%d,',cador);
      hRastro = getappdata(v2,'Lhandles');
      if isempty(hRastro)
        hRastro = flipud(findobj(v2,'type','line'));
      end;
      fct = findobj(hRastro,'Tag','SkipCur');
      if length(fct) [c kk] = setdiff(hRastro,fct);
                    c = sortrows([c kk],2);
                    hRastro = c(:,1);
      end;
      axes(v2);
      if isempty(en1) en1 = [.7 .7 .7; 0 0 0; 1 1 .5; 1 0 0; 0 0 0]; end;
      trk = ~sum(en1(2,:));
      if length(hRastro)
         en1n = length(en1(:,1));
         cHilo = length(hRastro);
         cli4 = cHilo + 4;
         if en1n < cli4
           en1 = en1(min(1:cli4,en1n*ones(1,cli4)),:);
         end;
         en3n = length(en3);
         if en3n<cHilo
           en3 = en3(min(1:cHilo,en3n*ones(1,cHilo)));
         end;
         en1n     = length(en1(:,1));
         objectoC = en1(5:en1n,:);
         tact = 0;
         if ax2  hLines2 = findobj(ax2,'type','line'); else hLines2 = 0; end;
         for kk=1:cHilo
            hi = hRastro(kk);
            argo = [get(hi,'x'); get(hi,'y')];
            if sum(objectoC(kk,:)) curColor = objectoC(kk,:);
            else                  curColor = get(hi,'color');
            end;
            Hc(14+kk)=line('x',argo(1,1),'y',argo(2,1),'erase','xor','color',curColor,...
                            'clipping','on','vis',vis,'use',[hi trk]);
            set(Hc(14+kk),'Marker',en3(kk),'MarkerSize',en4);
            set(hi,'ButtonD',sprintf('%s''lineCB'',%d);',rasHT,kk));
            if ~tact & (sum(get(hi,'vis'))==221) tact = kk; end;
            if length(find(hi==hLines2)) set(Hc(14+kk),'par',ax2); end;
         end;
      else disp('no lines to attach cursors to')
      end;
      if max(max(en0))>2 unitt = 'Pixels'; else unitt = 'Normal'; end;
      Sfuente = (196-get(0,'screenpix'))/10;
      for kk=2:3
        if kk==2  cordo  = '''x'');'; else cordo  = '''y'');'; end;
        cordo= [rasHT '''scaleAxis'',' cordo];
        if isempty(en2) en2 = ['x';'y']; end;
        Hc(kk) = uicontrol(hf,'Style','text','fontsi',Sfuente,'vis',vis,'str',deblank(en2(kk-2+1,:)),...
                'Units',unitt,'pos',en0(kk-1,:),'backg',en1(1,:),'horiz','cent','ButtonD',cordo,'call',cordo);
      end;
      for kk=4:7
        if trk==1 objectoC = [1 1 1]; else objectoC = en1(2,:); end;
        cordo  = sprintf('%s''editCB'',%d);',rasHT,kk-4);
        if kk==5 | kk==7   objectoC=get(gcf,'color'); end;
        Hc(kk)  = uicontrol(hf,'Style','edit','fontsi',Sfuente,'vis',vis,'str',' ','Units',unitt,'pos',en0(kk-1,:),...
                           'foreg',get(v2,'color'),'backg',objectoC,'horiz','left','call',cordo);
      end;
      if ~DuaC set(Hc(7),'vis','of'); end;
      set(Hc(5),'vis','of');
      nw1 = Hc(4);  nw2 = Hc(6);
      ih = [8 9 10];  c = [173 175 68];  u = {-inf inf ''};
      cbs = {'''peakval'',0);', '''peakval'',1);', '''markCB'','''');'};
      for k = 1:3
        liy = ih(k);
        Hc(liy) = uicontrol(hf,'Units',unitt,'pos',en0(liy-1,:),'vis',vis,'horiz','cent',...
                  'str',char(c(k)),'fontname','symbol','fontw','bol',...
                  'call',[rasHT cbs{k}],'use',u{k});
      end;
      set(Hc(8), 'ui',uicontextmenu('call',[rasHT '''peakval'',2);']));
      set(Hc(9),'ui',uicontextmenu('call',[rasHT '''peakval'',3);']));
      set(Hc(10),'fontsi',12,'ena','of');
      Hc(12) = line('x',[],'y',[],'erase','xor','color',en1(4,:),'vis','of');
      set(Hc(12),'Marker','+','MarkerSize',5*en4,'tag','DeltaC');
      if nargin2<9 en5 = ''; end;
      if isempty(en5) en5 = ['%7w';'%7w']; end;
      set(Hc(3),'use',en5);
      if length(en0(:,1)) > 9
        uicontrol(hf,'style','slider','Units',unitt,'pos',en0(10,:),...
             'backg',[.3 .3 .3],'Min',0,'Max',1000,'Val',500,'use',500,...
             'vis',vis,'call',[rasHT '''xslider'');'],'tag','xslider');
      end;
      set(get(v2,'xlabel'),'ButtonD',[rasHT '''xincr'');']);
      kk = [rasHT '''axisCB'');'];
      Hc(11) = line('x',[],'y',[],'erase','xor','vis','of','color',en1(3,:),'ButtonD',kk);
      set([hf v2],'ButtonD',kk);
      if nargin2>=11 & length(en7) jst1 = en7; else jst1 = -1; end;
      if nargin2>=12 set(Hc(5),'use',en8); else set(Hc(5),'use',''); end;
      tm = zeros(1,8);
      if ax2 set(ax2,'ButtonD',[rasHT '''AxisCBaux'');']); end;
      tm(5) = jst1;
      tm(4) = max(1,tact);
      set(nw1,'use',tm);
      set(Hc(7),'use','','Tag','');
      set(Hc(2),'use',Hc);
      CurPri(cador) = Hc(2);
      setappdata(0,'CurMain',CurPri);

    case 719
      if nargin<8
        en4 = [];
        if nargin<7
          en3 = [];
          if nargin<5 en1 = []; end;
        end;
      end;
      if isempty(en1)
        RastH = get(hact,'use');
        x = get(RastH(1),'x');  absc = get(RastH(1),'y');
        if en0<1 en0 = tm(1); end;
        lix = length(x);
        if en0>lix  lix = round(lix/2);  en1 = x(lix);    en2 = absc(lix);
        else                       en1 = x(en0);  en2 = absc(en0);
        end;
      end;
      if isempty(en3) en3 = sum(get(Hc(12),'vis'))==221; end;

      lmx = xyli(1:2);   ylim = xyli(3:4);
      xd = 1.01 * (en1-lmx);
      yd = 1.01 * (en2-ylim);
      if xd(1)<=0 lmx = lmx + xd(1); elseif xd(2)>=0 lmx = lmx+xd(2); end;
      if yd(1)<=0 ylim = ylim + yd(1); elseif yd(2)>=0 ylim = ylim+yd(2); end;
      if ~isequal(xyli,[lmx ylim])
         set(v2,'xli',lmx);
         if ax2 set(ax2,'xli',lmx); end;
         set(get(hact,'par'),'yli',ylim);
         if length(findobj(v2,'use','grid')) plt('grid',v2,'update'); end;
         evalQ(get(Hc(5),'use'));
      end;

      set(hact,'y',en2,'x',en1);
      tipo = get(Hc(3),'use');
      if length(en4)
         fr = repmat(max(en4.*[.9 1 .5])<=.5,1,3);
         set([nw1 nw2],'backg',en4,'foreg',fr);
         allC = Hc(14+1:end);
         if DuaC allC(DuaC) = []; end;
         set(allC,'vis','of');
         set(hact,'vis','on');
         if ax2 | en3  set(Hc(7),'backg',en4,'foreg',fr); end;
         set(Hc([2:3]),'ena','inact'); set(Hc(10),'ena','on');
      end;
      set(nw2,'str',plt('ftoa',tipo(2,:),en2));
      set(nw1,'str',plt('ftoa',tipo(1,:),en1));
      if en3
         bk = get(nw1,'backg');
         fr = get(nw1,'foreg');
         set(Hc(7),'str',plt('ftoa',tipo(2,:),en2-get(Hc(12),'y')),'backg',bk,'foreg',fr,'vis','on');
         set(Hc(5),'str',plt('ftoa',tipo(1,:),en1-get(Hc(12),'x')),'backg',bk,'foreg',fr,'vis','on');
      end;
      tm([1 2 3]) = [en0 en1 en2];
      set(nw1,'use',tm);
      s = get(Hc(7),'Tag');
      if length(s) eval(s); end;
      sx = getappdata(v2,'xstr');  sy = getappdata(v2,'ystr');  s = [sx sy];
      rep1 = {'@XVAL','real(@XY)',...
              '@YVAL','imag(@XY)',...
              '@XY',  'plt("cursor",@CID,"get","position")',...
              '@IDX', 'plt2nd({"cursor",@CID,"get","position"})',...
              '@LNUM','plt("cursor",@CID,"get","activeLine")',...
              '@HAND','plt2nd({"cursor",@CID,"get","activeLine"})',...
              '@XU'  ,['get8192(' int2str(8192*sx) ',"user")'],...
              '@YU'  ,['get8192(' int2str(8192*sy) ',"user")'],...
              '@CID', sprintf('%d',cador)};
      for k=1:length(s)
        set(s(k),'str',evalRep2(getappdata(s(k),'evl'),rep1));
      end;
      if length(getappdata(v2,'moveCBext')) | length(find(v2==findobj(gcf,'type','axes')))
        s = get(Hc(7),'use');
        if length(s) evalRep(s,rep1); end;
      end;

    case 632
      tipo  = get(Hc(3),'use');
      kvs = tm(1);
      RastH = get(Hc(14+DuaC),'use'); lh1 = RastH(1);
      argo = [get(lh1,'x'); get(lh1,'y')];
      tm(6) = argo(2,kvs);  set(nw1,'use',tm);
      largo = length(argo(1,:)); if kvs>largo kvs=largo; end;
      set(Hc(14+DuaC),'x', max(xyli(1),min(xyli(2),argo(1,kvs))),...
                      'y', max(xyli(3),min(xyli(4),argo(2,kvs))));
      if sum(get(Hc(12),'vis'))==221 return; end;
      if RastH(2)==1
         set(Hc(7),'backg',get(lh1,'color'));
         set(Hc(5),'vis','of');
      end;
      set(Hc(7),'str',plt('ftoa',tipo(2,:),argo(2,kvs)));
    case 472
      p0 = getappdata(get(v2,'xlab'),'OldCur');
      sc = get(v2,'yscale');  liny = sc(2)=='i';
      if ischar(en0)
          Cpunta = get(ax2,'curre');  cy = Cpunta(1,2);
          absc = get(ax2,'yli');
          if liny  absc = absc + p0(2) - cy;
          else     absc = absc * p0(2) / cy;
          end;
          set(ax2,'yli',absc);
      else
        Cpunta = get(v2,'curre');
        cx = Cpunta(1,1);  cy = Cpunta(1,2);
        if en0<1 x = LMda(1:2);
                 sc = get(v2,'xscale');
                 if sc(2)=='i' x = x + p0(1) - cx;
                 else          x = x * p0(1) / cx;
                 end;
                 if ax2 axb = [v2 ax2]; else axb = v2; end;
                 set(axb,'xli',x);
        end;
        if en0   absc = LMda(3:4);
                 if liny  yi = p0(2) - cy;  absc = absc + yi;
                 else     yi = p0(2) / cy;  absc = absc * yi;
                 end;
                 set(v2,'yli',absc);
                 if aconex
                   yr = get(ax2,'yli');
                   if liny  yr = yr + yi * diff(yr)/diff(absc);
                   else     yr = yr * (yr(2)/yr(1))^(log(yi)/log(absc(2)/absc(1)));
                   end;
                   set(ax2,'yli',yr);
                 end;
        end;
      end;
      fixMark;
    case 606
      p0 = getappdata(get(v2,'xlab'),'OldCur');
      if ischar(en0)
          Cpunta = get(ax2,'curre');
          absc = get(ax2,'yli');  y0 = absc(1);  y1 = absc(2);
          cy = Cpunta(1,2)-y0;  if ~cy cy=1e-06; end;
          set(ax2,'yli',[y0 y0 + abs((p0(2)-y0)*(y1-y0)/cy)]);
      else
        Cpunta = get(v2,'curre');
        if en0<1 x0 = LMda(1);  x1 = LMda(2);
                 cx = Cpunta(1,1)-x0;  if ~cx cx=1e-06; end;
                 if ax2 axb = [v2 ax2]; else axb = v2; end;
                 set(axb,'xli',[x0 x0 + abs((p0(1)-x0)*(x1-x0)/cx)]);
        end;
        if en0   y0 = LMda(3);  y1 = LMda(4);  dy = y1-y0;
                 cy = Cpunta(1,2)-y0;  if ~cy cy=1e-06; end;
                 cy = dy/cy;  p2 = p0(2)-y0;
                 set(v2,'yli',[y0 y0 + abs(p2*cy)]);
                 if aconex
                    sc = get(v2,'yscale');
                    if sc(2)=='i'
                       absc = get(ax2,'yli');  y0 = absc(1);  y1 = absc(2);
                       p2 = (y1-y0) * p2 / dy;
                       set(ax2,'yli',[y0 y0 + abs(p2*cy)]);
                    else
                       yr = get(ax2,'yli');  yr0 = yr(1);  yr1 = yr(2);
                       Cpunta = get(ax2,'curre');
                       cy = Cpunta(1,2)-yr0;  if ~cy cy=1e-06; end;
                       p0 = yr0 * exp(log(p0(2)/y0) * log(yr1/yr0) / log(y1/y0));
                       set(ax2,'yli',[yr0 yr0 + abs((p0-yr0)*(yr1-yr0)/cy)]);
                     end;
                 end;
        end;
      end;
      fixMark;
    case 662
      Cpunta  = get(v2,'curre');
      urX = max(LMda(1),min(LMda(2),Cpunta(1,1)));  yCur = max(LMda(3),min(LMda(4),Cpunta(1,2)));
      Cpunta   = get(Hc(11),'use');
      xAntiq  = Cpunta(1);  yAntiq  = Cpunta(2);
      if sum(get(v2,'yscale')) == 322
            vz = abs(log(LMda(3)/LMda(4))) < 50 * abs(log(yAntiq/yCur));
      else  vz = diff(LMda(3:4)) < 50 * abs(yAntiq-yCur);
      end;
      if vz
         if sum(get(v2,'xscale')) == 322
               vz = abs(log((.0001*LMda(2)+LMda(1))/LMda(2))) < 50 * abs(log(xAntiq/urX));
         else  vz = diff(LMda(1:2)) < 50 * abs(xAntiq-urX);
         end;
      end;
      if vz
          set(Hc(11),'x',[xAntiq urX urX xAntiq xAntiq],'y',[yAntiq yAntiq yCur yCur yAntiq]);
          tipo=get(Hc(3),'use');
          if sum(get(Hc(11),'vis')) == 315
             set(Hc(11),'vis','on');
             set(nw2,'str',plt('ftoa',tipo(2,:),yAntiq));
             set(nw1,'str',plt('ftoa',tipo(1,:),xAntiq));
             set(Hc([5 7 4 6]),'ena','on');
             set(Hc(10),'ena','of');
          end;
          bk = get(nw1,'backg');
          fr = get(nw1,'foreg');
          set(Hc(7),'str',plt('ftoa',tipo(2,:),yCur),'backg',bk,'foreg',fr,'vis','on');
          set(Hc(5),'str',plt('ftoa',tipo(1,:),urX),'backg',bk,'foreg',fr,'vis','on');
      end;
    case {570, 557},
      Mtipo = sum(get(gcf,'SelectionT'));
      if (sum(get(Hc(11),'vis'))==315) & ~tm(7) & ~tm(8)
         rasHT = sprintf('plt(''cursor'',%d,',cador);
         smv = '';
         PareMovi = 'set(gcf,{''WindowButtonMotionFcn'',''WindowButtonUpFcn''},{'''',''''});';
         Cpunta  = get(v2,'curre');  cx = Cpunta(1,1);  cy  = Cpunta(1,2);
         if     cx<LMda(1) | cx>LMda(2)   dragy =  '1);';
         elseif cy<LMda(3) | cy>LMda(4)   dragy =  '0);';
         else                               dragy = '-1);';
         end;
         urX = max(LMda(1),min(LMda(2),cx));
         dxy = diff(LMda);
         if isempty(getappdata(gcf,'NoEdit')) &  urX-cx > .08*dxy(1)
           switch Mtipo
             case 649, plt('click','EDIT','1');
             case 321,    plt('click','EDIT','2');
           end;
           return;
         end;
         switch Mtipo
         case 649
            yCur = max(LMda(3),min(LMda(4),cy));
            if length(get(hact,'ButtonD'))
              setappdata(hact,'ilast',getappdata(hact,'i'));
              plt('click','EDIT','0');
            end;
            ovt = 0;
            if sact==570
               nargin2 = 4;
               RastH = get(hact,'use');
               if sum(get(RastH(1),'vis')) ~= 221
                  for kk=14+1:length(Hc)
                      RastH = get(Hc(kk),'use');
                      if sum(get(RastH(1),'vis')) == 221 break; end;
                  end;
                  actv = kk-14;
               end;
               en0 = actv;
               if ax2
                 ylm = get(ax2,'yli');  dylm = diff(ylm);
                 cvert = (cy - LMda(3)) / dxy(3);
                 mdst = 1e+99;
                 idst = 0;
                 for kk=14+1:length(Hc)
                   if get(Hc(kk),'par') ~= ax2 continue; end;
                   RastH = get(Hc(kk),'use'); RastH = RastH(1);
                   if sum(get(RastH,'vis')) == 315 continue; end;
                   x = get(RastH,'x'); absc = get(RastH,'y');
                   if length(x)>999 | all(diff(x)>0)
                     [toss,j] = min(abs(urX-x));
                     dvert = (absc(j) - ylm(1)) / dylm;
                     acv = abs(cvert-dvert);
                   else
                     acv = min(abs((x-urX)/dxy(1)) + abs((absc-ylm(1))/dylm - cvert));
                   end;
                   if acv < mdst  idst = kk;  mdst = acv; end;
                 end;
                 if mdst < .02
                   ovt = 1;
                   ii = idst-14;
                   if ii ~= DuaC en0 = ii; end;
                 end;
               end;
            end;
            RastH = get(Hc(14+en0),'use');
            x = get(RastH(1),'x'); absc = get(RastH(1),'y');
            if nargin2==4
               set(Hc(8), 'use',-inf);
               set(Hc(9),'use',inf);
               if tm(5)==-1
                  df = diff(x);
                  jst1 = all(df>0) | all(df<0);
               else jst1 = tm(5);
               end;
               if RastH(2)==1 clRastro = get(RastH(1),'color');
               else        clRastro = get(nw2,'backg');
               end;
               indic = sum(get(Hc(12),'vis'))==221;
            else clRastro=[];  jst1=en1;  indic=en2;
            end;
            if jst1 [junk,kvs] = min(abs(urX-x));
            elseif get(RastH(1),'par')==ax2
                  ylm = get(ax2,'yli');  dylm = diff(ylm);
                  cvert = (cy - LMda(3)) / dxy(3);
                  [toss,kvs] = min(abs((x-urX)/dxy(1)) + abs((absc-ylm(1))/dylm - cvert));
            else  [toss,kvs] = min(abs((x-urX)/dxy(1)) + abs(yCur-absc)/dxy(3));
            end;
            actv = en0;  tm(4) = actv;
            iact = 14 + actv;
            hact = Hc(iact);
            set(nw1,'use',tm);
            plt('cursor',cador,'mainCur',kvs,x(kvs),absc(kvs),indic,clRastro);
            if DuaC plt('cursor',cador,'auxCur'); end;
            if (sact==557) | ovt
              if nargin2==4
                smv = sprintf('%s''lineCB'',%d,%d,%d);',rasHT,en0,jst1,indic);
              end;
            else setappdata(get(v2,'xlab'),'OldCur',[cx cy LMda]);
                 PareMovi = [rasHT '''svHist'');' PareMovi];
                 smv = [rasHT '''panAX'',' dragy];
            end;
         case 321
            setappdata(get(v2,'xlab'),'OldCur',[cx cy LMda]);
            PareMovi = [rasHT '''svHist'');' PareMovi];
            smv = [rasHT '''zoomAX'',' dragy];
         otherwise,
            set(Hc(11),'use',[max(LMda(1),min(LMda(2),cx)) max(LMda(3),min(LMda(4),cy))]);
            smv = [rasHT '''expbox'');'];
         end;
         if length(smv) set(gcf,'WindowButtonMotionFcn',smv,'WindowButtonUpFcn',PareMovi); end;
      elseif (sum(get(Hc(11),'vis'))==221) | tm(7) | tm(8)
        switch Mtipo
        case 649, plt('cursor',cador,'scale','new');
        case 321,    plt('cursor',cador,'restore');
        end;
        tm([7 8]) = 0;  set(nw1,'use',tm);
      end;
    case 872
      rasHT = sprintf('plt(''cursor'',%d,',cador);
      Cpunta  = get(ax2,'curre');  setappdata(get(v2,'xlab'),'OldCur',[Cpunta(1,1) Cpunta(1,2)]);
      switch sum(get(gcf,'SelectionT'));
      case 321,   smv = [rasHT '''zoomAX'',''R'');'];
      otherwise,   smv = [rasHT '''panAX'',''R'');'];
      end;
      set(gcf,'WindowButtonMotionFcn',smv,'WindowButtonUpFcn','set(gcf,{''WindowButtonMotionFcn'',''WindowButtonUpFcn''},{'''',''''});');
    case 555
      x   = str2double(get(Hc(4+en0),'str'));
      tipo = get(Hc(3),'use'); tipo = tipo((1+fix(en0/2)),:);
      xpx = (sum(get(Hc(11),'vis'))==221) | tm(7) | tm(8);
      if xpx  argo = Hc(11); argo = [get(argo,'x') get(argo,'y')];
              if en0 xold=argo(3*en0); else xold=argo(1); end;
      else    argo = [get(hact,'x') get(hact,'y')];
              if ~en0 xold=argo(1); elseif en0==2 xold=argo(2); else xold=0; end;
      end;
      if length(x)
         if xpx  ixy = [1 4 5; 2 3 3; 6 7 10; 8 9 9];  argo(ixy(en0+1,:)) = x;
                 set(Hc(11),'x',argo(1:5),'y',argo(6:10));
         else
            editd = length(get(hact,'ButtonD'));
            if imag(x) x = imag(x); setappdata(hact,'alt',0); end;
            if ~en0
               RastH = get(hact,'use'); RastH = RastH(1);
               xd = get(RastH,'x');  yd = get(RastH,'y');
               if editd set(hact,'x',x); plt('click','EDIT','/');
               else [v xk] = min(abs(xd-x));
                    xv = xd(xk);
                    if xv<xyli(1) | xv>xyli(2) plt('cursor',cador,'set','xlim',xv+[-.5 .5].*diff(xyli(1:2))); end;
                    plt('cursor',cador,'mainCur',xk);
               end;
            elseif en0==2  set(hact,'y',x);
                           if x<xyli(3) | x>xyli(4)
                              plt('cursor',cador,'set','ylim',x+[-.5 .5].*diff(xyli(3:4)));
                           end;
                           if editd plt('click','EDIT','/'); end;
            else set(Hc(4+en0),'str',plt('ftoa',tipo,tm(6)));
            end;
         end;
      else set(Hc(4+en0),'str',plt('ftoa',tipo,xold));
      end;
    case 641
      set(Hc(11),'x',[LMda(1:2) 0 0 0],...
                      'y',[LMda(3) 0 LMda(4) 0 0]);
      p0 = getappdata(get(v2,'xlab'),'OldCur');
      set(v2,'xli',p0(3:4),'yli',p0(5:6));
      plt('cursor',cador,'scale','new',0,0);
    case 520
      toria = get(nw2,'use');
      if isempty(toria) toria = [LMda 1]; end;
      lExp   = length(toria(:,1));
      gasto = find(toria(:,5)==1);
      skip   = [];
      switch sum(en0)
      case 330
         argo = Hc(11); argo = [get(argo,'x') get(argo,'y')];
         apXY = [sort(argo(1:2)) sort(argo([6 8])) 1];
         if isequal(LMda,apXY(1:4))
            if length(gasto) toria(gasto,5)=1; end;
         else
            if length(gasto)
               toria(gasto,5)=0;
               if max(gasto,lExp) < 4
                  toria = [toria; apXY];
               else
                  if gasto==4  toria(2:4,:)=[toria(3:4,:); apXY];
                  elseif gasto >= 1  toria(gasto+1,:)=apXY;
                  end;
               end;
            else
               if lExp < 4  toria = [toria; apXY];
               else               toria = [toria(lExp-1,:); apXY];
               end;
            end;
         end;
         escars=0;
      case 319
          if sum(get(gcf,'SelectionT'))==321 escars=1;
          else
            if isempty(gasto) gasto = lExp;
            else               toria(gasto,5)=0;  gasto = gasto-1;
            end;
            if gasto  escars=0;  apXY=toria(gasto,:); toria(gasto,5)=1;
            else       escars=1;
            end;
          end;
      case 441
          switch sum(en1)
          case 120,    escars=2;
          case 121,    escars=3;
          case 429, escars=1;
          otherwise, disp([en1 ' is not a valid action in plt(cursor,CurID,scale,auto,In2)']);
          end;
      otherwise, disp([en0 ' is not a valid action in plt(cursor,CurID,scale,In1)']);
      end;
      if escars == 0
          set(v2,'xli',sort(apXY(1:2)),'yli',sort(apXY(3:4)));
          wii = 0;
      else
         Nrastro=length(Hc)-14;  ximo=+inf;  xima=-inf;
         mc = zeros(Nrastro,2);  Rlista = [];   wii = 1;
         for kk=1:Nrastro
            rastro=get(Hc(14+kk),'use');
            if sum(get(rastro(1),'vis')) == 221
               xdatos = get(rastro(1),'x');
               ximo  = min(ximo,min(xdatos));
               xima  = max(xima,max(xdatos));
               if wii ydatos = get(rastro(1),'y'); end;
               Rlista = [kk Rlista];  wii = 0;
            end;
         end;

         if wii msgbox('Possible autoscale without any visible lines','Warning','warn','modal');
         else
            for kk=1:Nrastro
              k = Hc(14+kk);
              mc(kk,:) = [get(k,'x') get(k,'y')];
              set(k,'x',xdatos(1),'y',ydatos(1));
            end;
            skip = findobj(v2,'Tag','SkipCur','vis','on'); set(skip,'vis','of');
            switch escars
            case 1,  set(v2,'YlimM','auto');
                         apXY = [ximo xima get(v2,'yli') 1];
                         set(v2,'xli',apXY(1:2));
            case 2, apXY = [ximo xima get(v2,'yli') 1];
                         set(v2,'xli',apXY(1:2));
            case 3,
               tm=get(nw1,'use');
               if ismember(actv,Rlista)
                  rastro = get(hact,'use');  ydatos = get(rastro(1),'y');
                  ymax = max(ydatos);  ymin = min(ydatos);  dy = 0.25*(ymax-ymin);
                  ymin = ymin - dy;   ymax = ymax + dy;
                  if ymin ~= ymax  set(v2,'YlimM','man'); set(v2,'yli',[ymin ymax]); drawnow;
                  else             set(v2,'YlimM','auto');  drawnow;
                  end;
               else set(v2,'YlimM','man'); drawnow; set(v2,'YlimM','auto');
               end;
               apXY = [LMda 1];
            end;
            for kk=1:Nrastro set(Hc(14+kk),'x',mc(kk,1),'y',mc(kk,2)); end;
         end;
      end;
      if ~wii
         if ax2 set(ax2,'xli',apXY(1:2)); end;
         set(nw2,'use',toria);
         if aconex & nargin<6
           yr = get(ax2,'yli');
           set(ax2,'yli',yr(1) + (get(v2,'yli') - LMda(3)) * diff(yr) / (LMda(4)-LMda(3)));
         end;
         plt('cursor',cador,'restore');
         axes(v2); evalQ(get(Hc(5),'use'));
      end;
      set(v2,'YlimM','man');
      set(skip,'vis','on');
      plt('grid',v2,'update');
    case 548,
      RastH = get(hact,'use');
      lx = length(get(RastH(1),'x'));   xk = tm(1);
      ty = get(gcf,'SelectionT');  ty = ty(1);
      if ty=='o' & getappdata(gcbo,'ty') ty = 'a'; end;
      ty = ty=='a';
      if ty xk = max(xk-1,1);
      else  xk = min(xk+1,lx);
      end;
      setappdata(gcbo,'ty',ty);
      plt('cursor',cador,'mainCur',xk);
      if DuaC plt('cursor',cador,'auxCur'); end;
    case 763,
      RastH = get(hact,'use'); h1 = RastH(1);
      x = get(h1,'x');   absc = get(h1,'y'); lx = length(x);
      lmx = LMda(1:2);   ylim = LMda(3:4);  dxlim = diff(lmx);
      xk = tm(1);    v = get(gcbo,'Val');
      x1 = min(x);  x2 = max(x);   xvalo = x2 - x1;
      if xvalo/dxlim < 2
        dx = round(v - 500);
        if abs(dx)==10 pmove = .01; else pmove = .05; end;
        xk = max(min(lx,xk + sign(dx)*round(lx*pmove)),1);
        x = x(xk);  xd = 1.01 * (x-lmx);
        if xd(1)<0 lmx = lmx + xd(1); elseif xd(2)>0 lmx = lmx+xd(2); end;
        v = 500;
      else
        dx = v - get(gcbo,'use');
        switch round(abs(dx))
          case 10,   lmx = lmx + sign(dx)*dxlim/10;
          case 100,  lmx = lmx + sign(dx)*dxlim;
          otherwise, lmx = x1 + xvalo*v/1000 + dxlim*[-.5 .5];
        end;
        if     lmx(1)<x1  lmx = [x1 x1+dxlim] - dxlim/50;
        elseif lmx(2)>x2  lmx = [x2-dxlim x2] + dxlim/50;
        end;
        xlc = mean(lmx);
        v = 1000*(xlc-x1)/xvalo;
        v = max(min(v,1000),0);
        [dmy xk] = min(abs(xlc-x));
      end;
      set(gcbo,'Val',v,'use',v);
      absc = absc(xk);  yd = 1.01 * (absc-ylim);
      if yd(1)<0 ylim = ylim + yd(1); elseif yd(2)>0 ylim = ylim+yd(2); end;
      if sum(LMda - [lmx ylim])
         set(v2,'xli',lmx,'yli',ylim);
         if ax2 set(ax2,'xli',lmx); end;
         if length(findobj(v2,'use','grid')) plt('grid',v2,'update'); end;
         evalQ(get(Hc(5),'use'));
      end;
      plt('cursor',cador,'mainCur',xk);
      if DuaC plt('cursor',cador,'auxCur'); end;
    case 740
      b = gcbo;  set(b,'vis','of'); drawnow; set(b,'vis','on');
      switch en0
        case 2, en0=0; set(Hc(8), 'use',-inf);
        case 3, en0=1; set(Hc(9),'use',inf);
      end;
      RastH = get(hact,'use'); h1 = RastH(1);  absc = get(h1,'y'); x = get(h1,'x');
      x12 = find(x <= LMda(2) & x >= LMda(1));
      if isempty(x12) x = LMda(1); x12 = 1:length(x); disp('You must select a line for the min/max finder'); end;
      absc = absc(x12);  nabsc = length(absc);
      if en0
         xk=find(absc < [absc(2:nabsc) inf] & ...
                 absc < [inf  absc(1:nabsc-1)] & ...
                 absc > get(Hc(9),'use'));
         if isempty(xk)
              [absc xk] = min(absc);
         else [absc kk]  = min(absc(xk));
              xk = xk(kk);
         end;
         set(Hc(9),'use',absc);
      else
         xk=find(absc > [absc(2:nabsc) -inf] & ...
                 absc > [-inf absc(1:nabsc-1)] & ...
                 absc < get(Hc(8),'use'));
         if isempty(xk)
              [absc xk] = max(absc);
         else [absc kk] = max(absc(xk));
              xk = xk(kk);
         end;
         set(Hc(8),'use',absc);
      end;
      xk = xk + x12(1) - 1;
      if sum(get(Hc(12),'vis')) == 221
         indic = 1;
         if RastH(2)==1 clRastro = get(RastH(1),'color');
         else        clRastro = get(nw2,'backg');
         end;
      else indic = 0;
           clRastro = [];
      end;
      plt('cursor',cador,'mainCur',xk,x(xk),absc,indic,clRastro);
      if DuaC plt('cursor',cador,'auxCur'); end;
    case 772
      tipo = get(Hc(3),'use');
      set(nw2,'str',plt('ftoa',tipo(2,:),tm(3)));
      set(nw1,'str',plt('ftoa',tipo(1,:),tm(2)));
      set(Hc(11),'vis','of');
      if DuaC
         RastH = get(Hc(14+DuaC),'use');
         if RastH(2)==1
              clRastro = get(RastH(1),'color');
              set(Hc(7),'str',plt('ftoa',tipo(2,:),tm(6)),'ena','on','backg',clRastro);
         else set(Hc(7),'str',plt('ftoa',tipo(2,:),tm(6)),'ena','on')
         end;
         set(Hc(5),'str','','ena','of','vis','of');
         plt('cursor',cador,'auxCur');
      else
         if sum(get(Hc(12),'vis')) == 315
            set(Hc([5 7]),'vis','of','str','');
         end;
      end;
      set(hact,'y',max(xyli(3),min(xyli(4),tm(3))),...
               'x',max(xyli(1),min(xyli(2),tm(2))));
      set(Hc([2:3]),'ena','inact');  set(Hc(10),'ena','on');
    case 560
      set(Hc(10),'vis','of','backg',1-get(Hc(10),'backg'));
      drawnow; set(Hc(10),'vis','on');
      if sum(get(Hc(12),'vis')) == 315
            set(Hc(12),'vis','on','par',get(hact,'par'),'use',actv,...
                'x',get(hact,'x'),'y',get(hact,'y'));
      else  set(Hc(12),'vis','of'); set(Hc(8:9),'ena','on');
            plt('cursor',cador,'restore');
      end;
    case 493
      fg = get(nw1,'par');
      a = flipud(findobj(get(findobj(fg,'use','TraceID'),'child'),'type','text'));
      if length(a)<actv  s = 'Yval';
      else s = deblank(get(a(actv),'str'));
      end;
      set(findobj(fg,'use','idcur'),{'str','color'},{s,get(hact,'color')});
    case 320
      switch sum(en0)
      case 315
         Volv1 = Hc(2:12);
      case 885
         if nargin2==5 rasNo=en1; else rasNo=actv; end;
         k = Hc(14+rasNo);
         Volv1 = complex(get(k,'x'),get(k,'y'));
         Volv2 = tm(1);
      case 1028
         RastH = get(hact,'use');
         Volv1 = actv;
         Volv2 = RastH(1);
      case 625
         toria = get(nw2,'use');
         if isempty(toria)
            Volv1  = [[LMda 1]; [zeros(3,4), -ones(3,1)]];
         else aa = 4 - length(toria(:,1));
              Volv1 = [ toria; [ zeros(aa,4), -ones(aa,1) ] ];
         end;
         if ax2 Volv1 = [Volv1; [get(ax2,'xli') get(ax2,'yli') 2]];
         else   Volv1 = [Volv1; [zeros(1,4) -1]];
         end;
      otherwise, disp('error in plt(cursor,get)');
      end;
    case 332
      switch sum(en0)
      case 557
         h = [Hc([2:10 14+1:end]) ...
              findobj(gcf,'tag','xstr') ...
              findobj(gcf,'tag','ystr') ...
              findobj(gcf,'use','idcur') ...
              findobj(gcf,'tag','xslider')];
         h = findobj(h,'vis','on');
         set(h,'vis','of');
         setappdata(Hc(2),'hid',h);
      case 495, set(getappdata(Hc(2),'hid'),'vis','on');
      case 885
         for kk=2:10 set(Hc(kk),'pos',en1(kk-1,:)); end;
      case 334,      set(Hc(7),'vis',en1);
                      if DuaC set(Hc(DuaC),'vis',en1); end;
      case 423,     set(Hc(2),'str',en1);
      case 424,     set(Hc(3),'str',en1);
      case 570,   set(Hc(5),'use',en1);
      case 572,   if length(en1) & ischar(en1) & en1(1)==';'
                         setappdata(v2,'moveCBext',1);  en1(1) = [];
                      end;
                      set(Hc(7),'use',en1);
      case 622,  set(Hc(7),'Tag',en1);
      case {563,443,442},
         oldEh = get(nw2,'use');
         toria = oldEh;
         if length(toria)
             gasto = find(toria(:,5)==1);
             if length(gasto) toria(gasto,5)=0; end;
         end;
         switch sum(en0)
         case 563, xyl = en1;  toria(1,:) = [xyl 1];
         case 443
            if length(oldEh)  toria(1,:) = [LMda(1:2) en1,1];  xyl = [LMda(1:2) en1];
            else              toria=[];                         xyl = [LMda(1:2) en1];
            end;
         case 442
            if length(oldEh)  toria(1,:) = [en1 LMda(3:4) 1]; xyl  = [en1 LMda(3:4)];
            else              toria=[];                        xyl = [en1 LMda(3:4)];
            end
         end;
         set(v2,'xli',sort(xyl(1:2)),'yli',sort(xyl(3:4)));
         set(nw2,'use',toria);
         if ax2
            if nargin2==6 argo = en2; else  argo = get(ax2,'yli'); end;
            set(ax2,'xli',xyl(1:2),'yli',argo);
            set(nw1,'use',tm);
         end;
         if length(findobj(v2,'use','grid')) plt('grid',v2,'update'); end;
         evalQ(get(Hc(5),'use'));
      case 540, set(Hc(8), 'use',-inf); set(Hc(9),'use',inf);
      case 625,
         if length(en1(1,:)) < 5
            plt('cursor',cador,'set','xylim',en1);
         elseif isequal(size(en1),[5 5])
           toria = [];
           for kk=1:5
             if kk<=4
               if en1(kk,5) >= 0
                 toria=[toria; en1(kk,:)];
                 if en1(kk,5) == 1
                    LMda = en1(kk,1:4);
                    set(v2,'xli',en1(kk,1:2),'yli',en1(kk,3:4));
                    if length(findobj(v2,'use','grid'));
                       plt('grid',v2,'update');
                    end;
                 end;
               end
             else
               if en1(kk,5)==2 set(ax2,'xli',en1(kk,1:2),'yli',en1(kk,3:4)); end;
             end;
           end;
           set(nw2,'use',toria);
         else disp('error in plt(cursor,CurId,set,expHis,xxx), xxx is wrong shape');
         end;
      case 1028
        if en1
          set(hact,'vis','of');
          actv = en1;
          tm(4) = actv;
          iact = 14 + actv;
          hact = Hc(iact);
          set(hact,'vis','on');
          set(nw1,'use',tm);
        end;
        if nargin<6 en2 = 1; end;
        plt('cursor',cador,'mainCur',en2,[],[],[],get(hact,'color'));
      end;
    case 519
      delete(Hc([2:12 14+1:end]));
      CurPri(cador) = 0;
      if ~sum(CurPri) CurPri = []; end;
      setappdata(0,'CurMain',CurPri);
    case 925
      IsY = strcmp(en0,'y');
      Mtipo = get(gcf,'SelectionT');
      if strcmp(Mtipo,'normal')
         if tm(7+IsY)
            plt('cursor',cador,'scale','new');  tm(7)=0;  tm(8)=0;  set(nw1,'use',tm);
         else
            tipo = get(Hc(3),'use'); tipo = tipo(IsY+1,:);
            of2 = IsY*2;
            set(Hc(4 +of2),'str',plt('ftoa',tipo,LMda(of2+1)),'ena','on' );
            set(Hc(5+of2),'str',plt('ftoa',tipo,LMda(of2+2)),'ena','on' ,...
                                 'backg',get(nw1,'backg'),'foreg',get(nw1,'foreg'),'vis','on');
            tm(7+IsY)=1;
            if ~tm(7) | ~tm(8)
               set(Hc(11),'x',LMda([1 2 2 1 1]),'y',LMda([3 3 4 4 3]),'vis','on');
            end;
            set(nw1,'use',tm);
         end;
      elseif  strcmp(Mtipo,'alt')
          plt('cursor',cador,'scale','auto',en0);
      end;
    otherwise disp([accion ' is not a valid action in plt(cursor)']);
    end;

case 785
  switch varargin{2}
  case 0,
     p = gcbo;  e = getappdata(p,'edt');  f = get(e,'use');
     obj = getappdata(p,'obj');  obj = obj(1);
     c = get(p,'str');   v = get(p,'Val');  prop = deblank(c(v,:));
     if strcmp(prop,'Delete') & ishandle(obj) delete(obj); end;
     if ishandle(obj)
       s = get(obj,prop);
       if isnumeric(s)  set(e,'str',plt('vtoa','%6w',s,'  '));
                        if f & length(s)==3 plt('ColorEdit',2,f,2); end;
       else             set(e,'str',s);
       end;
     else set(e,'str','Deleted');
     end;
  case 1,
     e = gcbo;  p = getappdata(e,'pop'); obj = getappdata(p,'obj');
     c = get(p,'str'); v = get(p,'Val');   prop = deblank(c(v,:));
     s = get(e,'str'); f = get(e,'use');
     if strcmp(prop,'Delete')
       if strcmpi(s,'all')
         if length(c(:,1))==9 c='line'; else c='text'; end;
         delete(findobj(get(get(e,'par'),'use'),'type',c,'tag','mark'));
       end;
     else if isnumeric(get(obj(1),prop))
             s = str2num(s);
             tid = getappdata(p,'tid');
             if length(s)==3 set(tid,prop,s); end;
          end;
          set(obj,prop,s);
          if f plt('ColorEdit',2,f,2); end;
     end;
   case 2,
     p = gcbo;  e = getappdata(p,'edt');  v = get(p,'Val');  f = get(e,'use');
     h = get(p,'use');  h = h{v}(1);
     if v==3 c = get(h,'xcol');
     else    c = get(h,'color');
     end;
     set(e,'str',plt('vtoa','%6w',c,'  '));
     if f plt('ColorEdit',2,f,2); end;
   case 3,
     e = gcbo;  p = getappdata(e,'pop');  v = get(p,'Val');
     c = str2num(get(e,'str'));
     h = get(p,'use');  h = h{v};
     if v==3 set(h,'xcol',c,'ycol',c);
     else    set(h,'color',c);
     end;
  end;

case 642
  t = gco;
  Cpunta  = get(get(t,'par'),'curre');
  if nargin==1
    switch sum(get(gcf,'SelectionT'));
    case 649
      setappdata(t,'Dxy',diff([Cpunta(1,:); get(t,'pos')]));
      set(gcf,'WindowButtonMotionFcn','plt(''marker'',0);','WindowButtonUpFcn','set(gcf,{''WindowButtonMotionFcn'',''WindowButtonUpFcn''},{'''',''''});');
    case 321
      callc = 'plt(''ColorEdit'',0);';
      callp = 'plt(''MarkEdit'',0);';
      calle = 'plt(''MarkEdit'',1);';
      set(0,'units','pix');  sz = get(0,'screens');
      szw = sz(3) - 302 - 4;
      ppos  = get(gcf,'pos');
      if ppos(4)<1 sp = sz(3:4); ppos = ppos .* [sp sp]; end;
      xp = min(ppos(1)+ppos(3)+6,szw);
      yp = ppos(2)+ppos(4)-85-30*length(findobj('type','fig','use',gcf));
      figure('menu','none','number','off','Back','off','resize','off','pos',[xp yp 302 85],'color',[0,.4,.4],'name','Edit Marker','use',gcf,...
                     'closereq','plt(''click'',''mark'',4);');
      ed1 = uicontrol('sty','edit','pos',[8 5 130 22],'call',calle,'ButtonD',callc,'use',0);
      pu1 = 'Delete|Color|LineStyle|LineWidth|Marker|MarkerSize|Xdata|Ydata|Zdata';
      pu1 = uicontrol('sty','pop','str',pu1,'pos',[8 35 130 20],'call',callp,'Val',5);
      uicontrol('sty','text','str', 'Marker properties:','pos',[8 64 130 17]);
      ed2 = uicontrol('sty','edit','pos',[145 5 150 22],'call',calle,'ButtonD',callc,'use',0);
      pu2 = uicontrol('sty','pop','str','Delete|Color|FontAngle|FontName|FontSize|FontWeight|HorizontalAlign|Position|Rotation|String|VerticalAlign',...
                   'pos',[145 35 150 20],'call',callp,'Val',10);
      uicontrol('sty','text','str', 'String properties:','pos',[145 64 150 17]);
      set([ed1 ed2 pu1 pu2],'backg',[.8,.8,.9],'foreg','black');
      l = get(t,'use');
      setappdata(ed1,'pop',pu1); setappdata(pu1,'edt',ed1); setappdata(pu1,'obj',l);
      setappdata(ed2,'pop',pu2); setappdata(pu2,'edt',ed2); setappdata(pu2,'obj',t);
      if ishandle(l) l=get(l,'Marker'); else l='Deleted'; end;
      if ishandle(t) t=get(t,'str');      else t='Deleted'; end;
      set(ed1,'str',l);  set(ed2,'str',t);
    end;
  else  dxy = Cpunta(1,:) + getappdata(t,'Dxy'); set(t,'pos',dxy(1:2));
  end;

case 901
  switch varargin{2}
  case 0,
    e = gcbo;
    c = str2num(get(e,'str'));
    if length(c) ~= 3 return; end;
    f = get(e,'use');
    if ~f
      par = get(e,'par');
      f = figure('menu','none','number','off','Back','off','pos',get(par,'pos')+[0 -150 0 120],'color',[0,.4,.4],...
                 'name','Color select','use',get(par,'use'),...
                 'closereq','plt(''ColorEdit'',4);');
      s1 = [0 100];
      s4 = 'plt(''ColorEdit'',1,%d,%d);';  d = c*100;
      s = [plt('slider',[.025 .870 110],[d(1) s1 s1],'Red (%)'  ,sprintf(s4,1,f),2);
           plt('slider',[.025 .545 110],[d(2) s1 s1],'Green (%)',sprintf(s4,2,f),2);
           plt('slider',[.025 .220 110],[d(3) s1 s1],'Blue (%)' ,sprintf(s4,3,f),2);
          ]';
      v2 = axes('xli',[.8 12.15],'yli',[.8 15.4],'color',[0 0 0],'xcol',[0,.4,.4],'ycol',[0,.4,.4],'XtickL',' ','YtickL',' ','TickLen',[0 0]',...
                'units','nor','pos',[.408 .02 .57 .96]);
      ph = zeros(11,11);
      cb = sprintf('plt(''ColorEdit'',2,%d);',f);
      for fila = 1:11
        for col = 1:11
          ph(fila,col) = patch(col+[0 1 1 0],fila+[0 0 1 1],[0 0 0],'use',[fila col],'ButtonD',cb);
        end;
      end;
      pat = patch([1 12 12 1],12+[.5 .5 3 3],c,'use',ph,'ButtonD',sprintf('plt(''ColorEdit'',2,%d,1);',f));
      set(v2,'use',[s e pat c]);  set(e,'use',f);
    end;
    plt('ColorEdit',2,f,2);
  case 1,
    s = varargin{3};
    f = varargin{4};
    v2 = findobj(f,'type','axes');
    u = get(v2,'use');
    c = [0 0 0];
    for k=1:3
      h = u(k);
      c(k) = plt('slider',h,'get')/100;
      if k==s  bk = [1 1 0]; else  bk = [0 1 1]; end;
      obj = plt('slider',h,'get','obj');
      set(obj(5),'backg',bk);
    end;
    ph = get(u(5),'use');
    phc = [0 0 0];
    phc(s) = c(s);
    v = mod(s,3)+1;  w = mod(v,3)+1;
    for fila = 1:11
      phc(v) = (fila-1)/10;
      for col = 1:11
        phc(w) = (col-1)/10;
        set(ph(fila,col),'FaceColor',phc);
      end;
    end;
    set(u(5),'FaceColor',c);
    es = get(u(4),'str');
    if es(1)~='[' return; end;
    set(u(4),'str',plt('vtoa','%6w',c,'  '));
    p = getappdata(u(4),'pop');
    s = get(p,'str');  v = get(p,'Val');  obj = getappdata(p,'obj');
    if length(obj)
      prop = deblank(s(v,:));
      set(obj,prop,c);
      tid = getappdata(p,'tid');
      if tid set(tid,prop,c); end;
    else  h = get(p,'use');  h = h{v};
          if v==3 set(h,'xcol',c,'ycol',c);
          else    set(h,'color',c);
          end;
    end;
  case 2,
    f = varargin{3};
    v2 = findobj(f,'type','axes');
    u = get(v2,'use');
    if nargin<4 in4 = 0; else in4 = varargin{4}; end;
    switch in4
      case 0, c = get(gcbo,'FaceColor');
      case 1, c = u(6:8);
      case 2, c = str2num(get(u(4),'str'));
              if length(c)~=3 return; end;
    end;
    for k=1:3
      obj = plt('slider',u(k),'get','obj');   bk =get(obj(5),'backg');
      if bk(1) break;  end;
    end;
    v = mod(k,3)+1;  w = mod(v,3)+1;
    plt('slider',u(v),'set',c(v)*100);
    plt('slider',u(w),'set',c(w)*100);
    plt('slider',u(k),'set',c(k)*100);
  case 4,
    u = get(findobj(gcf,'type','axes'),'use');
    if ishandle(u(4)) set(u(4),'use',0); end;
    closereq;
  end;

case 434
  if nargin<2 [fi pth] = uigetfile('plt.plt','Select plt figure to open');
              if isnumeric(fi) return; end;
              fi = [pth fi];
  else        fi = varargin{2};
  end;
  ydat = [];
  feval('load',fi,'-mat');
  if isempty(ydat)
    p = find(fi=='.');
    if length(p) f = [fi(1:p(end)) 'mat']; end;
    dos(['copy ' fi ' ' f ' > NUL:']);
    load(f); delete(f);
  end;
  xdat = [xdat'; ydat'];
  plt(xdat{:},params{:});

case 431
  if nargin<2 [fi pth] = uiputfile('plt.plt','Save plt figure as');
              if isnumeric(fi) return; end;
              fi = [pth fi];
  else        fi = varargin{2};
  end;
  AX    = findobj(gcf,'Tag','click');
  cador = get(AX,'use');
  AX2   = findobj(gcf,'YAxisLoc','right','use',cador);
    RastH = getappdata(gcf,'Lhandles');
    xdat = get(RastH,'x');  ydat = get(RastH,'y');
    xlm =  get(AX,'xli');   ylm =  get(AX,'yli');
    xymult = getappdata(gcf,'xymult');
    if xymult(1) ~= 1
       mult = 1/xymult(1);
       for k = 1:length(xdat) xdat{k} = mult * xdat{k}; end;
       xlm = xlm * mult;
    end;
    ym = 1;
    for k = 1:length(ydat)
      mult = xymult(k+1);
      if mult ~= 1  
         mult = 1/mult;
         if ym ylm=ylm*mult; ym=0;  end;
         ydat{k} = mult * ydat{k};
      end;
    end;
    v = get(RastH,'vis');  v = [v{:}];  
    params = [getappdata(gcf,'params') { ...
             'position' get(gcf,'pos') ...
             'DIStrace' v(find(v=='o')+1)=='f' ...
             'Xlim'     xlm ...
             'Ylim'     ylm     }];
    if AX2 params = [params {'YlimR' get(AX2,'yli')}]; end;
    ver = '11May10';
    save(fi,'xdat','ydat','params','ver');

case 518
  y2    = varargin{2};
  AX    = findobj(gcf,'Tag','click');
  cador = get(AX,'use');
  AX2   = findobj(gcf,'YAxisLoc','right','use',cador);
  AXrl  = [AX AX2];
  if isempty(AX2) AX2=0; end;
  CurPri = getappdata(0,'CurMain');
  Hc = get(CurPri(cador),'use');
  if ischar(y2)
    switch sum(y2)
    case 430
      s = get(gcbo,'str');
      if s(1)==92 & s(2)=='d'  s = s(6:end-5);
      else                             s = ['\div ' s ' \div'];
      end;
      set(gcbo,'str',s);
    case 733
      imX = get(AX,'xli'); imY = get(AX,'yli');
      d = [-.2 .2]; dr = 'o';
      clic = get(gcf,'SelectionT'); clic = clic(1);
      if clic=='a' | (clic=='o' & getappdata(AX,'dir')=='i')
         d = d/-1.4; dr = 'i';
      end;
      set(AXrl,'xli',imX + diff(imX)*d);
      set(AX  ,'yli',imY + diff(imY)*d);
      if AX2  axl = get(get(AX2,'ylabel'),'str');
              if axl(1)~=92 | axl(2) ~= 'd';
              imY = get(AX2,'yli'); set(AX2,'yli',imY + diff(imY)*d);
              end;
      end;
      plt('grid',AX,'update');
      setappdata(AX,'dir',dr);
      axes(AX); evalQ(get(Hc(5),'use'));
    case 427
      if nargin>2 in3 = varargin{3}; else in3 = -1; end;
      switch in3
        case 3,
          cFIGbk = get(gcf,'color');
          cPLTbk = get(AX,'color');
          if isstr(cPLTbk) cPLTbk = get(AX2,'color'); end;
          cXYax  = get(AX,'xcol');
          cXYlbl = get(get(AX,'xlab'),'color');
          cGRID  = get(findobj(gcf,'type','line','use','grid'),'color');
          cDELTA = get(findobj(gcf,'type','line','tag','DeltaC'),'color');
          cTRACE = get(getappdata(gcf,'Lhandles'),'color');
          cTRACE = reshape([cTRACE{:}],3,length(cTRACE))';
          cFile  = get(findobj(gcf,'style','push','str','D'),'tag');
          if isempty(cFile) [cFile pth] = uiputfile('*.mat','Select file for saving colors');
                            cFile = [pth cFile];
          end;
          if sum(cFile)
             save(cFile,'cFIGbk','cPLTbk','cXYax','cXYlbl','cGRID','cDELTA','cTRACE');
          else disp('No file was selected'); end;
          return;
        case 4,
           h = findobj(gcf,'style','edit');
           for k=1:length(h)
             f = get(h(k),'use');
             if f & ishandle(f) close(f); end;
           end;
           closereq;
           return;
      end;
      set(0,'units','pix');  sz = get(0,'screens');
      szw = sz(3) - 302 - 4;
      ppos  = get(gcf,'pos');
      if ppos(4)<1 sp = sz(3:4); ppos = ppos .* [sp sp]; end;
      xp = min(ppos(1)+ppos(3)+6,szw);
      yp = ppos(2)+ppos(4)-85-30*length(findobj('type','fig','use',gcf));
      callc = 'plt(''ColorEdit'',0);';
      if in3==2
        g = gcf;
        figure('menu','none','number','off','Back','off','resize','off','pos',[xp yp+25 302 60],'color',[0,.4,.4],'name','Edit figure colors','use',g,...
                      'closereq','plt(''click'',''mark'',4);');
        ps = 'Figure background|Plot background|Axis color|Axis labels|Grid color|Delta cursor';
        pu = uicontrol('sty','pop','str',ps,'pos',[80 35 140 20],'call','plt(''MarkEdit'',2);');
        ed = uicontrol('sty','edit','pos', [80  5 140 22],'call','plt(''MarkEdit'',3);','ButtonD',callc,'use',0);
        set([ed pu],'backg',[.8,.8,.9],'foreg','black');
        setappdata(ed,'pop',pu); setappdata(pu,'edt',ed); setappdata(pu,'obj',[]);
        set(ed,'str',plt('vtoa','%6w',get(g,'color'),'  '));
        v2 = findobj(g,'type','axes');  a = v2(1); ap = a;
        if isstr(get(ap,'color')) ap = v2(end); end;
        set(pu,'use',{[g v2(2)];
                     ap;
                     a;
                     [get(a,'xlab'); get(a,'ylab'); get(v2(2),'child')];
                     findobj(g,'type','line','use','grid');
                     findobj(g,'type','line','tag','DeltaC');
                    });
        return;
      end;
      nw1 = Hc(4);
      tm = get(nw1,'use');
      actv = tm(4);
      iact = 14 + actv;
      hact = Hc(iact);
      clic = get(gcf,'SelectionT'); clic = clic(1);
      if clic=='a' | nargin>2
        h = findobj(gcf,'use','TraceID');  tx = [];  hb = [];  tid = [];
        if length(h) h = flipud(get(h,'child'));
                     tx = findobj(h,'type','text');  tid = tx(actv);
                     h = findobj(h,'type','line');
                     bt = get(h,'button');
                     if iscell(bt)
                       hb = h(find(cellfun('length',bt)));
                       if length(hb)>=actv tid = [tid hb(actv)]; end;
                     end;
        end;
        c = hact;
        h = get(c,'use'); h = h(1);
        if in3>=0  alll = in3;
        else       alll = get(Hc(11),'vis');
                   alll = alll(2)=='n';
        end;
        if alll c = [];  h = [];
                for kk=14+1:length(Hc)
                  RastH = get(Hc(kk),'use');
                  c = [c Hc(kk)];  h = [h RastH(1)];
                end;
                tid = [tx; hb];
                fname = 'Edit all lines';
        else    fname = sprintf('Edit Line %d',actv);
        end;
        callp = 'plt(''MarkEdit'',0);';
        calle = 'plt(''MarkEdit'',1);';
        figure('menu','none','number','off','Back','off','resize','off','pos',[xp yp 302 85],'color',[0,.4,.4],'name',fname,'use',gcf,...
               'closereq','plt(''click'',''mark'',4);');
        props = 'Color|LineStyle|LineWidth|Marker|MarkerSize|Xdata|Ydata|Zdata';
        ed1 = uicontrol('sty','edit','pos',[8 5 130 22],'call',calle,'ButtonD',callc,'use',0);
        pu1 = uicontrol('sty','pop','str',props,'pos',[8 35 130 20],'call',callp);
        uicontrol('sty','text','str', 'Line properties:','pos',[8 64 130 17]);
        ed2 = uicontrol('sty','edit','pos',[145 5 150 22],'call',calle,'ButtonD',callc,'use',0);
        pu2 = uicontrol('sty','pop','str',props,'pos',[145 35 150 20],'call',callp,'Val',5);
        uicontrol('sty','text','str', 'Cursor properties:','pos',[145 64 150 17]);
        set([ed1 ed2 pu1 pu2],'backg',[.8,.8,.9],'foreg','black');
        setappdata(ed1,'pop',pu1); setappdata(pu1,'edt',ed1); setappdata(pu1,'obj',h); setappdata(pu1,'tid',tid);
        setappdata(ed2,'pop',pu2); setappdata(pu2,'edt',ed2); setappdata(pu2,'obj',c); setappdata(pu2,'tid',[]);
        set(ed1,'str',plt('vtoa','%6w',get(h(1),'color'),'  '));
        set(ed2,'str',get(c(1),'MarkerSize'));
        axes(AX);
      else
        nw2 = Hc(6);
        axes(AX);
        x = get(hact,'x'); absc = get(hact,'y');
        p = get(hact,'par');
        rlim = [];
        if p ~= AX
          rlim = get(p,'yli');
          ylim = get(AX,'yli');
          absc = ylim(1) + diff(ylim) * (absc - rlim(1)) / diff(rlim);
        end;
        l = line(x,absc,'marker','s');
        t = text(x,absc,['   (' get(nw1,'str') ', ' get(nw2,'str') ')'],...
                 'fontsi',get(p,'fontsi'),'use',l,'ButtonD','plt(''marker'');');
        set([t l],'color',get(hact,'color'),'tag','mark');
        if length(rlim)
          set(l,'tag','markR','use',[t AX p ylim rlim]);
        end;
      end;

    case 674
      hl = findobj(gcf,'ButtonD','plt(''click'',''TGLlogy'');');
      if strcmp(get(AX,'Yscale'),'log')
           sc='linear'; st='LinY';
      else sc='log';    st='LogY';
           absc = get(AX,'yli');    if absc(1)<=0 set(AX,'yli',absc(2)*[.001 1]); end;
           if AX2
             absc = get(AX2,'yli'); if absc(1)<=0 set(AX2,'yli',absc(2)*[.001 1]); end;
           end;
      end;
      set(AXrl,'Yscale',sc);  set(hl,'str',st);
    case 673
      hl = findobj(gcf,'ButtonD','plt(''click'',''TGLlogx'');');
      if strcmp(get(AX,'Xscale'),'log')
           sc='linear'; st='LinX';
      else sc='log';    st='LogX';
           x = get(AX,'xli'); if x(1)<=0 set(AXrl,'xli',x(2)*[.001 1]); end;
      end;
      set(AXrl,'Xscale',sc);  set(hl,'str',st);
    case 653, set(AX,'TickLen',(1-plt('grid',AX,'toggle'))*[.01 .025]);
                   a = getappdata(gcf,'axis');
                   if AX2 a(end) = []; end;
                   for k=2:length(a)
                     set(a(k),'TickLen',(1-plt('grid',a(k),'toggle'))*[.01 .025]);
                   end;

    case 668
      if strcmp(get(gcf,'SelectionT'),'normal')
        f = get(gcf,'menu');
        if   f(1)=='f' delete(findobj(gcf,'type','uimenu')); set(gcf,'menu','none');
        else set(gcf,'menu','fig');
             v = get(0,'ShowHidden');
             set(0,'ShowHidden','on');
             a = findobj(gcf,'label','&File');
             uimenu(a,'Label','p&lt  save','separator','on','call','plt(''save'');');
             uimenu(a,'Label','pl&t  open','call','plt(''open'');');
             uimenu(a,'Label','plt  &hardCopy','call','plt(''hcpy'',''init'',gcf);');
             c = get(a,'child');
             set(a,'child',c([4:end-3 1:3 end-2:end]));
             set(0,'ShowHidden',v);
             liy = uimenu('label','&Color');
             uimenu(liy,'label','&Edit line','accel','e','call','plt(''click'',''mark'',0);');
             uimenu(liy,'label','&Edit all lines',       'call','plt(''click'',''mark'',1);');
             uimenu(liy,'label','&Edit figure colors',   'call','plt(''click'',''mark'',2);');
             uimenu(liy,'label','&Save figure colors',   'call','plt(''click'',''mark'',3);','separator','on');
        end;
      else
        k = get(gcbo,'use');  if isempty(k) k=1; end;
        mrk = 'none';  sty = '-';
        switch k
          case 1, mrk = 'o';  sty = 'none';
          case 2, mrk = 'o';
          case 3, k = 0;
        end;
        set(gcbo,'use',k+1);
        set(getappdata(gcf,'Lhandles'),'marker',mrk,'linestyle',sty);
      end;
    case 242
      Hc = get(CurPri(cador),'use');
      tm = get(Hc(4),'use');
      lH = get(Hc(14+tm(4)),'use');
      nw1 = Hc(4);  nw2 = Hc(6);  hix2 = Hc(5);  hiy2 = Hc(7);
      x = get(lH(1),'x');  absc = get(lH(1),'y');
      lmx = get(Hc(13),'xli');
      x12 = find(x >= lmx(1)  &  x <= lmx(2));
      absc = absc(x12); absc = absc(~isnan(absc)); nabsc = length(absc);
      p = findobj(gcf,'use','idcur');  ps = get(p,'str');
      switch sum(ps)
        case 286,  s = 'RMS'; r = plt('ftoa','%7w',sqrt(sum(absc.^2)/nabsc));
        case 242,  s = 'y/x';
                    xr = s2d(get(nw1,'str'));
                    r = getappdata(p,'idcur');
                    r = plt('ftoa','%7w',s2d(r{2})/xr);
                    if sum(get(hiy2,'vis'))==221 & sum(get(hix2,'vis'))==221
                       xr = str2num(get(hix2,'str'));
                       set(hiy2,'str',s2d(get(hiy2,'str'))/xr);
                    end;
        case 288,  s = '\surdx^2+y^2';
                    xr = s2d(get(nw1,'str'));
                    q = getappdata(p,'idcur');
                    r = plt('ftoa','%7w',abs(xr + s2d(q{2})*1j));
                    if sum(get(hiy2,'vis'))==221 & sum(get(hix2,'vis'))==221
                       xr = str2num(get(hix2,'str'));
                       set(hiy2,'str',abs(xr + s2d(q{3})*1j));
                    end;
        case 1110, r = getappdata(p,'idcur');  s = r{1}; set(hiy2,'str',r{3});  r = r{2};
        otherwise,  setappdata(p,'idcur',{ps get(nw2,'str') get(hiy2,'str')});
                    s = 'Avg';
                    r = plt('ftoa','%7w',sum(absc)/nabsc);
      end;
      set(p,'str',s); set(nw2,'str',r);
    case 294
      v3 = varargin{3};  en2 = v3 - '0';
      if en2>5
         Cpunta = get(get(gco,'par'),'curre');  en2 = en2-5;
         if bitget(en2,1) set(gco,'x',Cpunta(1,1)); end;
         if en2>1         set(gco,'y',Cpunta(1,2)); end;
         return;
      end;
      if en2>2
         set(gcf,'WindowButtonMotionFcn',['plt click EDIT ' char(v3+3) ';'],...
                 'WindowButtonUpFcn',['set(gcf,{''WindowButtonMotionFcn'',''WindowButtonUpFcn''},{'''',''''});' ' plt click EDIT /;']);
         return;
      end;
      Hc = get(CurPri(cador),'use');
      tm = get(Hc(4),'use');
      hact = Hc(14 + tm(4));
      mksz = get(hact,'markersize');
      sz = (196-get(0,'screenpix'))/10;
      if en2==-1
        RastH = get(hact,'use');  RastH = RastH(1);
        x = get(RastH,'x');  absc = get(RastH,'y');   kk = tm(1);
        x1 = get(hact,'x');  y1 = get(hact,'y');
        if strcmp(get(gcf,'SelectionT'),'alt') | length(getappdata(hact,'alt'))
             setappdata(hact,'alt','');
             iUlti = getappdata(hact,'ilast');  x0 = x(iUlti);  y0 = absc(iUlti);
             inc = 1 - 2*(iUlti>kk);  passo = max(1,abs(iUlti-kk));
             dx = (x1-x0)/passo; dy = (y1-y0)/passo;
             for k = iUlti:inc:kk x(k)=x0; absc(k)=y0; x0=x0+dx; y0=y0+dy; end;
        else
             if mksz==sz x(kk)=x1; absc(kk)=y1;
             else
                ylim = get(Hc(13),'yli');
                if y1<ylim(1) x(kk) = [];  absc(kk) = [];
                else          x = [x(1:kk) x1 x(kk+1:end)];
                              absc = [absc(1:kk) y1 absc(kk+1:end)];
                end;
             end;
        end;
        set(RastH,'x',x,'y',absc);
        setappdata(hact,'i',kk);
        setappdata(gcf,'NewData',1);
        return;
      end;
      liy = get(hact,'marker');
      p = get(hact,'par');
      v2 = findobj(gcf,'type','axes','tag','click');
      if en2==2 sz=sz+4; end;
      switch liy(1)
        case '^',  k = 1;
        case '>',  k = 2;
        case 'd',  k = 3;
        otherwise, k = 0;
      end;
      switch k + 4*isempty(getappdata(gcf,'DataEdit')) + 8*~en2
        case {1,5},   c = '>';  b = 'plt click EDIT 3;';
        case {3,4,7}, c = '^';  b = 'plt click EDIT 4;';
        case {0,2,6}, c = 'd';  b = 'plt click EDIT 5;';
        otherwise,              b = '';
      end;
      if p ~= v2  set(v2,'vis','of'); end;
      if length(b) ~= length(get(hact,'ButtonD'))
        ch = get(p,'child');
        if length(b)  a = find(ch==hact);
                      setappdata(hact,'swapp',{a liy mksz});
        else          a = getappdata(hact,'swapp');
                      c = a{2};  sz = a{3};  a = a{1};
                      if p ~= v2 set(v2,'vis','on'); end;
        end;
        sv = ch(a); ch(a) = ch(1); ch(1) = sv;  set(p,'child',ch);
      end;
      set(hact,'marker',c,'markersize',sz,'ButtonD',b);
      axes(AX);
    end;
    plt('grid',AX,'update');
  else y8 = y2/8192;  j = y8(2);
       Mtipo = get(gcf,'SelectionT');
       p = {'color'; 'marker'; 'linestyle'; 'linewidth'}; q = {[0 .3 .3] 'none' '-' 9};
       if strcmp(Mtipo,'normal')
          k = y8(1);  s = get(k,'vis');
          mk = getappdata(j,'mk');
          if s(2)=='f' set(k,'vis','on'); set(j,'fonta','nor','fontw','bol'); set(mk,p,get(k,p));
          else         set(k,'vis','of'); set(j,'fonta','ita','fontw','nor'); set(mk,p,q);
          end;
       else
          Mtipo = strcmp(Mtipo,'open');
          for k = get(get(j,'par'),'child')'
             t = get(k,'ButtonD');  bk = findstr(t,'[');  w = get(k,'type');
             if w(1)=='t' & length(bk)
               t = s2i(t(bk+1:findstr(t,' ')))/8192;
               mk = getappdata(k,'mk');
               if Mtipo | k==j  set(t,'vis','on');  set(k,'fonta','nor','fontw','bol');  set(mk,p,get(t,p));
               else               set(t,'vis','of');  set(k,'fonta','ita','fontw','nor');  set(mk,p,q);
               end;
             end;
          end;
       end;
       DuaC = getappdata(AX,'DualCur');
       Lh = getappdata(AX,'Lhandles');
       if DuaC v = get(Lh(DuaC),'vis'); else v = 'off'; end;
       plt('cursor',cador,'set','aux',v);
       ls = findobj(Lh,'par',AX);
       v = 'off';
       for k=1:length(ls)
         if strcmp(get(ls(k),'vis'),'on') v = 'on';  break; end;
       end;
       set(get(AX,'Ylabel'),'vis',v);
       if AX2
         ls = findobj(Lh,'par',AX2);
         v = 'off';  izquC = get(AX2,'color');
         for k=1:length(ls)
           if strcmp(get(ls(k),'vis'),'on') v = 'on';  izquC = 'none';
                                              set(AX2,'xli',get(AX,'xli'));
                                              break;
           end;
         end;
         set(AX2,'vis',v);  set(AX,'color',izquC);
       end;
       TIDback = getappdata(get(j,'par'),'TIDcback');
       if isempty(TIDback) return; end;
       evalRep(TIDback,{'@TID',int2str(y2(2)),'@LINE',int2str(y2(1))});
       axes(AX);
  end;

case 708,
  v2 = getappdata(gcf,'axis');  v2 = v2(1);
  cid = getappdata(gcf,'cid');  fct = length(cid);
  ylbl = get(v2,'Ylabel');
  obj = getappdata(ylbl,'obj');
  if length(obj)
    setappdata(ylbl,'obj',[]);
    for k=1:fct plt('cursor',cid(k),'set','visON'); end;
    if obj(1) set(obj,'vis','on'); end;
  else
    for k=1:fct plt('cursor',cid(k),'set','visOFF'); end;
    a = findobj(gcf,'str','Zout');
    if length(a) a = get(a(1),'par');
                 obj = [a; get(a,'child')];
                 set(obj,'vis','of');
    else         obj = 0;
    end;
    setappdata(ylbl,'obj',obj);
  end;

case 944,
  h = findobj(gcf,'use','TraceID');
  if isempty(h) return; end;
  h = flipud(findobj(get(h,'child'),'type','text'));
  lix = length(h);
  e = varargin{2};
  e(find(e>lix)) = [];
  v = zeros(lix,1);
  v(e) = 1;
  v = find(xor(v,3-cellfun('length',get(getappdata(gcf,'Lhandles'),'vis'))));
  if isempty(v) return; end;
  b = get(h,'button');
  eval([b{v}]);

case {959, 534},   for f = findobj('type','fig')'
                               b = get(f,'ButtonD');
                               if length(b)>3 & strcmp(b(1:4),'plt(')
                                 set(f,'closereq','closereq');
                                 close(f);
                               end;
                             end;
case 632,  set(flipud(findobj(get(findobj(gcf,'use','TraceID'),'child'),'type','text')), ...
                      {'str'},varargin{2}(:));
case 774, Volv1 = '11May10'; Volv2 = 0;
case 425,    eval('plt(''helpv'',0);');

otherwise,
cTRACE  = [0  1  0;   1  0  1;   0  1  1;   1  0  0;  .2 .6  1;
           1  1  1;   1 .6 .2;   0  0  1;   1 .2 .6;  .2  1 .6;
          .6  1 .2;  .6 .2  1;   1  1  0;   0 .6  0;  .6  0 .6;
           0 .6 .6;  .6 .6  0;  .7 .7 .7;  .6  0  0;  .2 .2 .7;
          .5 .5 .5;  .7 .2 .2;  .2 .7 .2;   0  0 .6;  .3 .3 .3;
           0 .9 .4;   0 .4 .9;  .9 .4  0;  .4 .9  0;  .9  0 .4;
          .4  0 .9;  .8 .5 .5;  .5 .8 .5;  .5 .5 .8;];
posFIG  = [9 45 700 525];
posFIGd = 1;
o4      = ones(1,4);
axisPOS = o4;
idPOS   = o4;
POScid  = [3,-.075];
cFIGbk  = [.25 .15 .15];
cPLTbk  = [0   0   0  ];
cXYax   = [1   1   1  ];
cXYlbl  = [.7 .8 .95];
CURcDEF = [1   1   .50];
cCRSR = CURcDEF;
cDELTA  = [1   0   0  ];
cGRID   = [.3  .3 .3  ];
for k=35:99 cTRACE = [cTRACE; .75 * cTRACE(k-34,:)]; end;
cDEFAULT = cTRACE;
LabelX = 'X axis';
LabelY = 'Y axis (Left)';
LabelYr = 'Y axis (Right)';
Titulo = '';
LiNueva = '';
NombreF = 'plt';
Xlim   = 'default';
Ylim   = 'default';
YlimR  = 'default';
Xscale = 'linear'; Xsc = 'LinX';
Yscale = 'linear'; Ysc = 'LinY';
Quiv = 0;
cabeQ = [.3 .3];
Reja = 'on';
Mbar = 0;
Xslide = 0;
ayudaF = '';
cFile = '#';
ENApre = ones(1,2);
ENAcur = ones(1,99);
DISras = 0;
estilo  = 0;
marca = 0;
AXr = 0;
Derecho = [];
DuaC = 0;
TIDcback = '';
Ridcol = 0;
Mcaja = [1 1 1 1 0 1 1 1 1];
RASid = reshape(sprintf('Line%2d',1:99),6,99)';
TRACEmk = 0;
Xstring = '';
Ystring = '';
moveCB = '';
axisCB = '';
aConex = 1;
SubPlot = 0;
NoCursor = 0;
GridEr = 'norm';
aLp = {}; aLv = {};
aRp = {}; aRv = {};
lLp = {}; lLv = {};
lRp = {}; lRv = {};
lXp = {}; lXv = {};
sz = posFIG([3 4 3 4]);
posAX   = [.1429  .0933  .8329  .8819];
posPICO = [.0071  .0441  .0286  .0343];
posVALY = [.0071  .0060  .0286  .0343];
posDEL  = [.0386  .0194  .0300  .0457];
posSLDR = [.0071  .0080  .1250  .0200];
posCXL  = [.1386  .0095  .0200  .0410];
posC1X  = [.1643  .0076  .1000  .0450];
posC2X  = [.2686  .0076  .1000  .0450];
posCYL  = [.7609  .0095  .0200  .0410];
posC1Y  = [.7865  .0076  .1000  .0450];
posC2Y  = [.8908  .0076  .1000  .0450];
posAX2  = [.0070  .0870  .0580  .0000];

Sfuente = (196-get(0,'screenpix'))/10;
if sum(lower(varargin{1})) == 310
     FIG = varargin{2};  figure(FIG);  set(FIG,'vis','of');
else FIG = figure('menu','none','number','off','Back','off','vis','of');
end;
set(FIG,'PaperPositionMode','auto','invert','off','PaperOrient','land',...
        'PaperUnits','norm','DoubleBuf','on','Invert','on');
AX = axes('units','nor','fontsi',Sfuente,'Tag','click');
Volv1 = [];  nt = 0;
kparam = [];
k = 1;
pp = 1;
while k<=nargin
  absc  = varargin{k};  k=k+1;
  if ischar(absc)
    if k>nargin
       disp('Error using ==> plt.  Not enough input arguments.');
       disp('For help on using plt, type "help plt"');
       eval('plt(''helpv'',0);');
       return;
    end;
    kparam = [kparam k-1 k];
    absc = lower(absc);  yy = varargin{k};   k=k+1;
    pfx = zeros(1,5);
    while 1
      b = findstr(absc(1),'+-<>.');
      if isempty(b) break; end;
      pfx(b) = 1;  absc(1) = [];
    end;
    switch sum(absc)
      case 546,     Titulo    = yy;
      case 442,      Xlim     = yy;
      case 443,      Ylim     = yy;
      case 557,     if absc(1)=='y'  YlimR = yy;  AXr = 1; else cXYax = yy; end;
      case 632,    LabelX   = yy;
      case 633,    LabelY   = yy;
      case 747,   LabelYr  = yy;  AXr = 1;
      case 542,     Derecho    = yy;  AXr = 1;
      case 752,   DuaC  = yy;
      case 727,   NombreF  = yy;
      case 640,    cPLTbk   = yy;
      case 614,    cFIGbk   = yy;
      case 626,    cTRACE   = yy;
      case 621,    cDELTA = yy;
      case 654,    cXYlbl   = yy;
      case 769,   cCRSR  = yy;
      case 521,     if any(yy<0) GridEr = 'xor'; end;
                       cGRID    = abs(yy);
      case 676,    if iscell(yy) estilo = char(yy);   else estilo = yy;   end;
      case 757,   if iscell(yy) marca = char(yy);  else marca = yy;  end;
      case 732,   if iscell(yy) RASid = char(yy); else RASid = yy; end;
      case 743,   if length(yy)==1 & yy  yy = [yy (yy+.9)/2 .9]; end;
                       TRACEmk = yy;
      case 638,    ENAcur   = yy;
      case 847,  DISras = yy;
      case 885,  if ~yy(3) yy(3) = yy(4)/.944; elseif ~yy(4) yy(4) = yy(3)*.944; end;
                       posFIG   = yy;
                       posFIGd  = 0;
      case 775,   axisPOS  = yy(1:4);
                       switch length(yy) case 5, idPOS(3)=yy(5); case 8, idPOS=yy(5:8); end;
      case 821,  TIDcback = yy;
      case 975, Ridcol = yy;
      case 783,   Xstring  = yy;
      case 784,   Ystring  = yy;
      case 873,  LiNueva = yy;
      case 635,    ENApre   = yy;
      case 668,    Quiv     = yy;  if min(yy)<2 disp('No quiver tail position '); return; end;
      case 515,     cabeQ    = yy;
      case 841,  ayudaF = yy;
      case 959, cFile    = yy;
      case 846,  cPLTbk = get(0,'defaultaxescolor'); cGRID = .4*cPLTbk + .3;
                       cFIGbk = get(0,'defaultfigurecolor');
                       cXYax  = get(0,'defaultaxesxcolor');  cXYlbl = cXYax;
                       cTRACE = get(0,'defaultaxescolororder');
                       if length(yy(1,:))==3 cTRACE = [cTRACE; yy]; end;
      case 636,    moveCB = yy;
      case 634,    axisCB = yy;
      case 867,  aConex = yy;
      case 777,   SubPlot = yy;
      case 310,
      case 780,   kq = 0;
        while kq < length(yy)
          kq = kq + 1;
          switch yy(kq)
            case 'T', Reja = 'off';
            case 'M', Mbar = 1;
            case 'X', Xsc = 'LogX';  Xscale = 'Log';
            case 'Y', Ysc = 'LogY';  Yscale = 'Log';
            case 'S', Xslide = 1;
            case 'N', NoCursor = 1;
            case '-', kq = kq + 1;  km =findstr(yy(kq),'HXYGPFMZRA');
                      if length(km) if km==10 Mcaja=0; else Mcaja(km)=0; end; end;
            case '+', kq = kq + 1;  km =findstr(yy(kq),'HXYGPFMZRA');
                      if length(km) if km==10 Mcaja=ones(1,9);
                                    else if pp & km~=5 Mcaja=0; pp=0; end;
                                         Mcaja(km)=1;
                                    end;
                      end;
          end;
        end;
      otherwise,
         if sum(pfx)
           if pfx(1) aLp = [aLp {absc}]; aLv = [aLv {yy}]; end;
           if pfx(2) aRp = [aRp {absc}]; aRv = [aRv {yy}]; end;
           if pfx(3) lLp = [lLp {absc}]; lLv = [lLv {yy}]; end;
           if pfx(4) lRp = [lRp {absc}]; lRv = [lRv {yy}]; end;
           if pfx(5) lXp = [lXp {absc}]; lXv = [lXv {yy}]; end;
         else
           if iscell(yy)
             if length(yy)==length(Volv1)
               absc = {absc};
               yy = yy(:);
             else fprintf('Warning: For parameter %s, found %d elements but expected %d\n',absc,length(yy),length(Volv1));
             end;
           end;
           set(Volv1,absc,yy);
         end;
    end;
  else
     if k<=nargin         yy = varargin{k}; else yy = 'a'; end;
     if ~isreal(absc)        yy = imag(absc);  absc = real(absc);
     elseif isnumeric(yy) k = k+1;
     else                 yy=absc; absc=1:length(yy);
     end;
     nxt = length(Volv1) + 1;
     if length(find(Quiv==nxt))
        a = Volv1(colaQ);
        x0 = get(a,'x');  x0 = transpose(x0(:));
        y0 = get(a,'y');  y0 = transpose(y0(:));
        x1 = transpose(absc(:));
        y1 = transpose(yy(:));
        rNaN = repmat(NaN,size(y1));
        x01 = x0+x1;  y01 = y0+y1;
        a = x01-cabeQ(1)*(x1+cabeQ(2)*y1);
        b = x01-cabeQ(1)*(x1-cabeQ(2)*y1);
        absc  = [x0; x01; rNaN; a; x01; b; rNaN];
        a = y01-cabeQ(1)*(y1-cabeQ(2)*x1);
        b = y01-cabeQ(1)*(y1+cabeQ(2)*x1);
        yy = [y0; y01; rNaN; a; y01; b; rNaN];
        absc = absc(:);  yy = yy(:);
        ENAcur(nxt) = 0;
        nxt = 0;
     end;
     H = line(absc,yy);  Volv1 = [Volv1; H];   nt = nt + length(H);
     if nxt colaQ = nt; end;
  end;
end;

if iscell(Ylim) YlimR = Ylim{2}; Ylim = Ylim{1}; AXr = 1; end;
nSP = length(SubPlot) - 1;
if (iscell(LabelX) | iscell(Xlim)) & ~nSP & nt>1 SubPlot = [100 -50 100]; nSP = 2; end;
if iscell(LabelX) LabelXr = LabelX{2};  LabelX = LabelX{1}; else LabelXr = 'X axis'; end;
if iscell(Xlim)   Xlimr   = Xlim{2};    Xlim   = Xlim{1};   else Xlimr   = 'default'; end;

nSPl = nSP;
nSPr = 0;
SPw = 1;
if nSP
  SubPlot = SubPlot/100;
  k = find(SubPlot<0);
  if length(k) SPw = -SubPlot(k);
               SubPlot(k) = [];      nSP = nSP-1;
               nSPr = nSP + 2 - k;   nSPl = nSP - nSPr;
  end;
end;
nID = nt - nSP;
if nID>99 & RASid disp(sprintf('Max # of traceIDs = %d',99)); return; end;
if ~iscell(LabelY) LabelY = {LabelY}; end;
if length(LabelY) >= nSP+2  LabelYr = LabelY{nSP+2}; AXr = 1; end;
setappdata(FIG,'params',varargin(kparam));
k = sum(Ridcol);
if k  Ridcol = [nID-k Ridcol]; else Ridcol = nID; end;
ntid = max(Ridcol);
ncol = length(Ridcol);
if ncol>1 & all(axisPOS == o4) & all(idPOS == o4)
  idPOS(3) = ncol;
  axisPOS = [.4 + ncol/2, 1, (210-11*ncol-ncol^2)/200, 1];
end;
if all(idPOS == o4) & (nt<6 | isempty(LabelY{1})) & ~nSP
   idPOS(3) = 1.2;
end;
fsep = length(findstr(cFile,filesep));
if length(cFile) & ~fsep
  foobar = 0;
  if exist('foobar')
       liy = feval('dbstack');
       if length(liy)
          lix = liy(end).name;
          nq = findstr('(',lix);
          if length(nq) lix = lix(1:nq(1)-2); end;
          np = findstr(filesep,lix);
          if length(np) np=np(end); else np=0; end;
          if length(lix)-np>30 & length(liy)>1  lix = liy(end-1).name; end;
          liy = feval('which',lix);
       end;
  else liy = GetExe;
  end;
  [pth, name] = fileparts(liy);
  if cFile(1) == '#' cFile = [name 'Color']; end;
  cFile = fullfile(pth,[cFile '.mat']);
end;
if fsep cFile = [cFile '.mat']; end;
if exist(cFile) == 2
  load(cFile);
end;
if length(ENAcur)<nt  ENAcur = [ENAcur ones(1,nt-length(ENAcur))]; end;
if length(DISras)<nt  DISras = [DISras zeros(1,nt-length(DISras))]; end;
nC = length(cTRACE(:,1));
if length(Titulo)
   axisPOS = axisPOS .* [1 1 1 .96];
   ntx = findstr('[TexOff]',Titulo);
   if length(ntx)==1
        Titulo(ntx:ntx+7) = [];
        title(Titulo,'color',cXYlbl,'HandleV','on','interp','none');
   else title(Titulo,'color',cXYlbl,'HandleV','on');
   end;
end;
setappdata(AX,'DualCur',DuaC);
if AXr & nID>1
   if isempty(Derecho) Derecho = nID; end;
   axisPOS = axisPOS .* [1 1 .94 1];
   if ischar(YlimR)
     mn = inf;  mx = -inf;
     for k=Derecho
       absc = get(Volv1(k),'y');  mn = min(mn,min(absc));  mx = max(mx,max(absc));
     end;
      df=(mx-mn)/20; YlimR=[mn-df mx+df];
      if ~diff(YlimR) YlimR = [mn mn+max(1e-12,mn*1e-12)];  end;
   end;
   if length(Derecho)==1 yclr = cTRACE(mod(Derecho-1,nC)+1,:); else yclr = cXYax; end;
   AXr = axes('units','nor','fontsi',Sfuente,'YAxisLoc','right','yli',YlimR,...
              'color',cPLTbk,'xcol',cXYax,'ycol',yclr,'xtick',[]);
   if ~aConex LabelYr = ['\div ' LabelYr ' \div']; end;
   ylabel(LabelYr,'color',yclr,'HandleV','on','ButtonD','plt click link;');
   set(Volv1(Derecho),'par',AXr);
   axes(AX);
   AXrl = [AX AXr];
else AXrl = AX; set(AX,'Box','On');
end;
mrk = repmat('+',1,nt);
mrk(Derecho) = 'o';
if DuaC
  if length(find(DuaC==Derecho)) mrk(DuaC)='s'; else mrk(DuaC) = '*'; end;
end;
ceq = isequal(cCRSR,CURcDEF);  cEXcaja = cCRSR;
if ceq & sum(cPLTbk)>2 cEXcaja=1-cEXcaja; cCRSR=1-cCRSR; end;
curclr = [.7 .7 .7; 0 0 0; cEXcaja; cDELTA];
ENAcurS = sum(ENAcur(1:nt));
if ceq & ENAcurS>1 cCRSR = [0 0 0]; end;
for k=1:nt  set(Volv1(k),'color',cTRACE(mod(k-1,nC)+1,:));
            if ENAcur(k) curclr=[curclr; cCRSR]; else  set(Volv1(k),'Tag','SkipCur'); end;
end;
if estilo  if length(estilo(:,1)) < nt estilo=estilo'; end;
          for k=1:nt
            if length(findstr(estilo(k,1),'+o*.xsd^v<>ph'))
                  set(Volv1(k),'LineStyle','none','Marker',estilo(k,1));
            else  set(Volv1(k),'LineStyle',estilo(k,:));  end;
          end;
end;
if marca if length(marca(:,1)) < nt marca=marca'; end;
          for k=1:nt set(Volv1(k),'Marker',marca(k,:)); end;
end;
for k=1:nt if DISras(k) set(Volv1(k),'vis','of'); end; end;

Ret1a = Volv1;
axS = [];
if nSP
  if posFIGd posFIG([3 4]) = [840 570]; end;
  ySP = .11;  ySPr = ySP;
  ySP  = cumsum([ySP   (1-ySP)*SubPlot(1:nSPl+1)]);
  hSP  = diff(ySP)  - .03;  ySP  = ySP  + .012;
  dx = ((AXr>0)*.025 - .07) * (nSPr>0);
  posAX = [posAX(1) ySP(1) (posAX(3)+dx)*SPw hSP(1)];
  dx = posFIG(3);  dx1 = 64/dx;  dx2 = 14/dx;  dx3 = 67/dx;  dx4 = 19/dx;
  dy = posFIG(4);  dy1 = 20/dy;  y1a = .0076;  y2a = y1a + 24/dy;
  x1 = .1386;  x2 = x1+dx4;  x3 = x2+dx3;  x4 = x2+nSPl*dx3;
  posCXL = [x1 y2a dx2 dy1];
  posCYL = [x1 y1a dx2 dy1];
  posC1X = [x2 y2a dx1 dy1];
  posC2X = [x3 y2a dx1 dy1];
  posC1Y = [x4 y1a dx1 dy1];
  posC2Y = posC1Y + [dx3 0 0 0];
  if nSPr
    ySPr = cumsum([ySPr (1-ySPr)*SubPlot(nSPl+2:end)]);
    hSPr = diff(ySPr) - .03;  ySPr = ySPr + .012;
    SPx = sum(posAX([1 3])) + .075;
    x1 = SPx-.03;  x2 = x1+dx4;  x3 = x2+dx3;  x4 = x2+(nSPr-1)*dx3;
    posCXLr = [x1 y2a dx2 dy1];
    posCYLr = [x1 y1a dx2 dy1];
    posC1Xr = [x2 y2a dx1 dy1];
    posC2Xr = [x3 y2a dx1 dy1];
    posC1Yr = [x4 y1a dx1 dy1];
    posC2Yr = posC1Yr + [dx3 0 0 0];
    p2r = [posCXLr;posCYLr;posC1Xr;posC2Xr;posC1Yr;posC2Yr];
  end;
  for k=1:nSP
    l = Volv1(nt+k-nSP);
    Ret1a(nt+k-nSP) = 0;
    c = get(l,'color');
    a = axes('ycol',c,'xcol',c);
    setappdata(a,'Lhandles',l);
    set(l,'par',a);
    if length(LabelY)>k ylabel(LabelY{k+1}); end;
    if k == nSPl+1 xlabel(LabelXr); end;
    axS = [axS a];
  end;
end;
set(AXrl,'pos',posAX.*axisPOS,'Xscale',Xscale,'Yscale',Yscale);
if ischar(Ylim)
   Ylim = get(AX,'yli');
end;
Ret1a(find(~Ret1a)) = [];
setappdata(AX,'Lhandles',Ret1a);
setappdata(FIG,'Lhandles',Volv1);
axData = [AX axS];
if AXr axData = [axData AXr]; end;
setappdata(FIG,'axis',axData);

pb = [posPICO;posVALY;posDEL];
if Xslide
   posAX2(2) = posAX2(2) + .028;  pb(:,2) = pb(:,2) + .028;  pb = [pb; posSLDR];
end;
cador = plt('cursor',AXrl,'init',...
  [posCXL;posCYL;posC1X;posC2X;posC1Y;posC2Y;pb],...
  curclr,'', mrk, 0.8*Sfuente,'','on',[],LiNueva);

set(findobj(gcf,'style','push','str','D'),'use',cDEFAULT,'tag',cFile);
rasHT = sprintf('plt(''cursor'',%d,',cador);
set(AXrl,'use',cador);
Izqu = setdiff(1:nt-nSP,Derecho);
nIzqu = length(Izqu);
if (AXr | nSP) & nIzqu==1
  yclr = cTRACE(mod(Izqu-1,nC)+1,:); leftclr=yclr;
else yclr = cXYax;  leftclr = cXYlbl;
end;
set(AX,'xcol',cXYax,'ycol',yclr);
if ischar(Xlim)
  if nIzqu Xlim = get(AX,'xli'); else Xlim = get(AXr,'xli'); end;
end;
[prefix xmult] = plt('metricp',max(abs((Xlim))));
if ENApre(1) & xmult~=1
  for k=1:nt set(Volv1(k),'x',xmult*get(Volv1(k),'x')); end;
  LabelX = [prefix LabelX];
else xmult = 1;
end;
set(AXrl,'xli',Xlim*xmult);
[prefix mult] = plt('metricp',max(abs((Ylim))));
ymult = ones(1,nt);
if ENApre(2) & mult~=1
  for k=1:nIzqu
     kk = Volv1(Izqu(k));
     set(kk,'y',mult*get(kk,'y'));
     ymult(Izqu(k)) = mult;
  end;
  LabelY{1} = [prefix LabelY{1}];
else mult = 1;
end;
set(AX,'yli',Ylim*mult);
setappdata(FIG,'xymult',[xmult ymult]);
xlabel(LabelX,'color',cXYlbl,'HandleV','on');
hYlab = ylabel(LabelY{1},'color',leftclr,'HandleV','on');
plt('cursor',cador,'set','moveCB2',[rasHT '''MVcur'');']);
plt('grid',AX,'init',cGRID,GridEr);
set(AX,'TickLen',(1-plt('grid',AX,Reja))*[.01 .025]);
axes('pos',posCYL+[0 0 .2 0],'vis','of');
text(-.02,.45,'','fontsi',Sfuente,'horiz','right','ButtonD','plt click RMS;','use','idcur');
if length(Xstring)
  if ischar(Xstring) & Xstring(1) == '?'
        Xstring(1)=[];
        a = uicontrol('sty','edit','units','nor','pos',posC2X.*[1 1 1.7 1],'horiz','cent',...
                 'backg',[.2 .2 .2],'foreg',[1 1 .3]);
  else  a = text(-2.22,.45,'','color',cXYlbl);
  end;
  set(a,'fontsi',Sfuente,'tag','xstr');
  setappdata(a,'evl',Xstring);
  ch = get(FIG,'child');
  ch([1 end-1]) = ch([end-1 1]);  set(FIG,'child',ch);
end;
if length(Ystring)
  if ischar(Ystring) & Ystring(1) == '?'
        Ystring(1)=[];
        a = uicontrol('sty','edit','units','nor','pos',posC2Y,'horiz','cent',...
                 'backg',[.2 .2 .2],'foreg',[1 1 .3]);
  else  a = text(.6,.45,'','color',cXYlbl);
  end;
  set(a,'fontsi',Sfuente,'tag','ystr');
  setappdata(a,'evl',Ystring);
  ch = get(FIG,'child');
  ch([1 end-2]) = ch([end-2 1]);  set(FIG,'child',ch);
end;
nMenu = sum(Mcaja);
aid = 0;
ahi = .035*nMenu;
if nID>1 & RASid
  h = 19*ntid;
  ahip = ahi * sz(2);
  hr = (sz(2)-85) / (h+ahip);
  if hr<1 ahi = ahi*hr;  h = h*hr;  end;
  aidp = idPOS.*[3 sz(2)-4-h 50 h]./sz;
  aid = axes('xli',[0 ncol],'yli',[-ntid 0]-.5,'color',cPLTbk,'xcol',cFIGbk,'ycol',cFIGbk,'XtickL',' ','YtickL',' ','TickLen',[0 0]','use','TraceID',...
        'units','nor','pos',aidp);
  setappdata(aid,'TIDcback',TIDcback);
  cRid = cPLTbk + .16*(2*(cPLTbk<.5)-1);
  fila = 1;  col = 1;
  bln = 0;
  lpl = {'color'; 'marker'; 'linestyle'; 'linewidth'};
  [tn tw] = size(RASid);
  if nID>tn RASid = [RASid; repmat(' ',nID-tn,tw)]; end;
  for k=1:nID
    if length(find(k==Derecho)) cR=cRid; line(col-[1 .05],[0 0]-fila,'color',cR,'LineWidth',9);
    else                      cR=cPLTbk;
    end;
    s = RASid(k,:);
    if all(s==' ') bln=bln+1; continue; end;
    d = text(col-.93,-fila,s);
    ms = sprintf('plt(''click'',[%d %d]);',8192*Volv1(k),8192*d);
    set(d,'fontsi',Sfuente,'fontw','bol','color',cTRACE(mod(k-1,nC)+1,:),'ButtonD',ms);...
    if TRACEmk
      mk = line(TRACEmk,TRACEmk*0-k,lpl,get(Volv1(k),lpl),'ButtonD',ms);
      setappdata(d,'mk',mk);
      if TRACEmk(1)<.25 set(d,'color',cR); end;
    else mk = [];
    end;
    if DISras(k) set(d,'fonta','ita','fontw','nor'); set(mk,lpl,{[0 .3 .3] 'none' '-' 9}); end;
    fila = fila+1;
    if fila>Ridcol(col) col=col+1; fila=1; end;
  end;
  if bln & ncol==1
     dy = aidp(4) * (1 - (nID-bln)/nID);
     set(aid,'yli',[bln-nID 0]-.5,'pos',aidp+[0 dy 0 -dy]);
  end;
end;
if nMenu
  posAX2(4) = ahi;  c = (cXYax + cFIGbk)/2;
  amb = axes('units','nor','pos',posAX2,'yli',[-nMenu 0]-.45,'Box','On','XaxisLoc','top',...
           'color',cFIGbk,'xcol',c,'ycol',c,'XtickL',' ','YtickL',' ','TickLen',[0 0]','Tag','MenuBox');
  b=0; t=[];
  txt = { 'Help',    Xsc,       Ysc,       'Grid',    'Print',   'Menu',    'Mark',    'Zout',   'XY\leftrightarrow'};
  btn = { 'plt(''helpv'',1);'; 'plt(''click'',''TGLlogx'');'; 'plt(''click'',''TGLlogy'');'; 'plt(''click'',''TGLgrid'');'; 'plt(''hcpy'',''init'',gcf);'; 'plt(''click'',''TGLmenu'');'; 'plt(''click'',''mark'');'; 'plt(''click'',''ZoomOut'');'; ''   };
  for k=1:length(Mcaja)
    if Mcaja(k)
       b=b-1;  t = [t text(.5,b,txt{k},'interp','tex')];  bt = btn{k};
       if isempty(bt) bt = [rasHT '''scale'',''old'');']; end;
       set(t(end),'ButtonD',bt);
       if k==1 set(t(1),'use',ayudaF); end;
    end;
  end;
  set(t,'fontsi',Sfuente,'color',cXYlbl,'horiz','cent');
end;
CurPri = getappdata(0,'CurMain');
Hc = get(CurPri(cador),'use');
set(Hc(4),'ButtonD','plt click EDIT 1;');
set(Hc(6),'ButtonD','plt click EDIT 2;');
if posFIG(1)<0 posFIG = abs(posFIG);
else for k=fliplr(findobj('Type','figure')')
      if get(k,'pos') == posFIG  posFIG = posFIG + [30 25 0 0]; end;
     end;
end;
set(FIG,'pos',posFIG,'Name',NombreF,'color',cFIGbk,'CloseReq',...
 [rasHT '''clear'');delete(findobj(''type'',''fig'',''user'',gcf));closereq;']);
if Mbar plt('click','TGLmenu'); end;
ac = plt('cursor',cador,'get','activeLine');
t = Volv1(ac); x=get(t,'x'); absc=get(t,'y'); k=round(length(x)/2);
set(AX,aLp,aLv);
set(get(AX,'Ylabel'),lLp,lLv);
set(get(AX,'Xlabel'),lXp,lXv);
if AXr set(AXr,aRp,aRv);
       set(get(AXr,'Ylabel'),lRp,lRv);
end;
xyll = get(AX,{'xli' 'yli'});
plt('cursor',cador,'mainCur',k,x(k),absc(k),0,get(t,'color'));
setappdata(AX,'xstr',findobj(gcf,'tag','xstr'));
setappdata(AX,'ystr',findobj(gcf,'tag','ystr'));
set(AX,{'xli' 'yli'},xyll);
if DuaC & ~DISras(DuaC) plt('cursor',cador,'auxCur');
else                            plt('cursor',cador,'set','aux','off');
end;
v = 'off';
for k=1:nIzqu
  if strcmp(get(Volv1(Izqu(k)),'vis'),'on') v = 'on'; break; end;
end;
set(hYlab,'vis',v,'ui',uicontextmenu('call','plt hideCur;'));
izquC = cPLTbk;
if AXr
  ls = Volv1(Derecho)';
  set(Hc(14+Derecho),'par',AXr,'erase','xor');
  v = 'off';
  for k=1:length(ls)
    if strcmp(get(ls(k),'vis'),'on')
       v = 'on';
       izquC = 'none';
       break;
    end;
  end;
  set(AXr,'vis',v);
end;
plt('cursor',cador,'MVcur');
if length(moveCB) plt('cursor',cador,'set','moveCB',moveCB); end;
if length(axisCB) plt('cursor',cador,'set','axisCB',axisCB); end;
set(AX,'color',izquC);
cidS = cador;
if nSP
  if aid & all(idPOS([1 2 4]) == [1 1 1])
    set(aid,'units','nor');
    p = get(aid,'pos');
    p(2) = ySP(2)-p(4)-.025;
    if p(2) < sum(posAX2([2 4]))+.015
       p(2) = posAX2(2);  posAX2(2) = p(2) + p(4) + .02;
       set(amb,'pos',posAX2);
    end;
    set(aid,'pos',p);
  end;
  h = plt('cursor',cador,'get','obj');
  q = [-1 .1 .1 .1];  u = [q;q;q;q];
  py1 = get(h(5),'pos');   p1 = [u;py1;u];
  p = {'pos' 'Xlim' 'color' 'Xscale' 'FontSize' 'TickLength' 'Box'};
  v = get(AX,p);
  if ischar(v{3}) v{3} = get(AXr,'color'); end;
  for k=1:nSPl
    v{1}([2 4]) = [ySP(k+1) hSP(k+1)];
    set(axS(k),p,v);
    p1(5,1) = p1(5,1) - dx3;
    cidS = [cidS plt('cursor',axS(k),'init',p1,'','','+',8)];
    set(axS(k),'use',cidS(end));
  end;
  if nSPr
    v{1}([1 3]) = [SPx .98-SPx];  p2r = [p2r;q;q;q];
    if ischar(Xlimr) k = get(Volv1(end),'x');  v{2} = [min(k) max(k)];
    else             v{2} = Xlimr;
    end;
    for k=1:nSPr
      kr = k+nSPl;
      v{1}([2 4]) = [ySPr(k) hSPr(k)];
      set(axS(kr),p,v);
      cidS = [cidS plt('cursor',axS(kr),'init',p2r,'','','+',8)];
      set(axS(kr),'use',cidS(end));
      if k==1 p2r([1 2 3 4 6],:) = [u;q]; end;
      p2r(5,1) = p2r(5,1) - dx3;
    end;
  end;
  creq = [];
  s = 'plt(''cursor'',';
  cidSS = {[s int2str(cidS(1))]};
  for k=1:nSP
    SS = [s int2str(cidS(k+1))];
    cidSS = [cidSS {SS}];
    creq = [creq SS ',''clear'');'];
    plt('grid',axS(k),'init',cGRID,GridEr);
    if length(Reja)==3 plt('grid',axS(k),'off'); end;
  end;
  set(FIG,'closereq',[creq get(FIG,'closereq')]);
  setappdata(FIG,'c',0);
  s1a = 'if getappdata(gcf,"c")==%d setappdata(gcf,"c",0); else setappdata(gcf,"c",getappdata(gcf,"c")+1);';
  s1 = sprintf(s1a,nSPl);
  s2 = ',"set","activeLine",0,@IDX); end';
  s3 = ',"set","xlim",get(gca,"xlim")); end';
  for k=0:nSPl
    if k==nSPl j=1; else j=k+2; end;
    plt('cursor',cidS(k+1),'set','moveCB',[s1 cidSS{j} s2]);
    plt('cursor',cidS(k+1),'set','axisCB',[s1 cidSS{j} s3]);
  end;
  xk = Xlim(1) + diff(Xlim)/9;
  [a xk] = min(abs(x-xk));
  plt('cursor',cador,'mainCur',xk);
  if nSPr
    s1 = sprintf(s1a,nSPr);
    for k=1:nSPr
      kr = k+nSPl+1;
      if k==nSPr j=nSPl+2; else j=kr+1; end;
      plt('cursor',cidS(kr),'set','moveCB',[s1 cidSS{j} s2]);
      plt('cursor',cidS(kr),'set','axisCB',[s1 cidSS{j} s3]);
    end;
    plt('cursor',cidS(nSPl+2),'mainCur',round(length(get(Volv1(end),'x'))/9));
    axS(nSPl+1) = [];
  end;
  set(axS,'xticklabel',[]);
  set(findobj(FIG,'use','idcur'),'vis','of');
end;
setappdata(FIG,'cid',cidS);
if NoCursor plt('cursor',cador,'set','visOFF'); end;
axes(AX);
if length(Reja)==2 plt('grid',AX); end;
set(FIG,'vis','on');
drawnow;
end;

function fixMark()
  markR = findobj(gcf,'tag','markR')';
  if isempty(markR) return; end;
  for l = markR
    u = get(l,'use');
    t = u(1);  v2 = u(2);  axr = u(3);  ylim = u(4:5);  rlim = u(6:7);
    if ishandle(t) p = get(t,'pos');  else p = [0 0]; end;
    absc = [get(l,'y') p(2)];
    absc = rlim(1) + diff(rlim) * (absc - ylim(1)) / diff(ylim);
    rlim = get(axr,'yli');   ylim = get(v2,'yli');
    absc = ylim(1) + diff(ylim) * (absc - rlim(1)) / diff(rlim);
    set(l,'y',absc(1),'use',[u(1:3) ylim rlim]);
    if ishandle(t) p(2) = absc(2);  set(t,'pos',p); end;
  end;

function r2 = plt2nd(v)
  [r1 r2] = plt(v{:});

function r = get8192(h,prop)
  r = get(h/8192,prop);

function evalQ(a)
  if ischar(a)     a = strrep(a,'"','''');
                   eval(a);
  elseif iscell(a) feval(a{:});
  else             feval(a);
  end;

function evalRep(a,rep)
  if ischar(a)     for k=1:2:length(rep) a = strrep(a,rep{k},rep{k+1}); end;
                   a = strrep(a,'"','''');
                   eval(a);
  elseif iscell(a) feval(a{:});
  else             feval(a);
  end;

function r = evalRep2(a,rep)
  if ischar(a)     for k=1:2:length(rep) a = strrep(a,rep{k},rep{k+1}); end;
                   a = strrep(a,'"','''');
                   r = eval(a);
  elseif iscell(a) r = feval(a{:});
  else             r = feval(a);
  end;

function v = s2d(s)
   v = sscanf(s,'%f');

function v = s2i(s)
   v = sscanf(s,'%d');
