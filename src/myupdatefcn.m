function txt = myupdatefcn(~,event_obj,t)
% Customizes text of data tips

pos = get(event_obj, 'Position');
% disp(pos(1));
txt = {['Coprime Pair: (', num2str(t{pos(1)}(1)), ',', num2str(t{pos(1)}(2)), ')'], ...
       ['Additional Periods: ', num2str(pos(2))], ...
       ['Power of PSL: ', num2str(pos(3))]};
