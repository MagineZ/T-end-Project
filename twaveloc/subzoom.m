function subzoom(arg)
%SUBZOOM subplot zoom utility: make equal Xlim and Ylim in subplots
%This file should be placed in a directory listed in Matlab search path.
%Click on the figure window containing subplots to make it the current 
%one before calling SUBZOOM, or call SUBZOOM inside a m-file immediately 
%after the creation of the figure.  
%After zooming in a subplot, click in the figure window outside subplot
%frames with the right mouse botton to get a context menu, then click
%with the left mouse botton on 'Equal Xlim' or 'Equal Ylim'.
%
% Author: Qinghua Zhang
% Copyright 2005 INRIA

if nargin<1
  arg = 'init';
end

switch arg
  case 'init'
    ucmh=uicontextmenu;
    set(gcf , 'UIContextMenu', ucmh)
    uimenu(ucmh, 'Label', 'Equal Xlim', 'callback', 'subzoom xlim');
    uimenu(ucmh, 'Label', 'Equal Ylim', 'callback', 'subzoom ylim');
     
  case {'xlim', 'ylim'}
    lasta = gca;
    if isempty(lasta)
      return
    end
    xylim = get(lasta, arg);
    ah = get(gcf, 'child');
    for ka=1:length(ah)
      if strcmp(get(ah(ka),'type'), 'axes')
        set(ah(ka), arg,xylim)
      end
    end 
end