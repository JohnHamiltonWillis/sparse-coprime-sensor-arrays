function txt = MinProductAnalysis_updatefcn(~,event_obj,t)
% Customizes text of data tips

pos = get(event_obj, 'Position');
% disp(pos(1));
txt = {['Coprime Pair: (', num2str(t{pos(1)}(1)), ',', num2str(t{pos(1)}(2)), ')'], ...
       ['Periods: ', num2str(pos(2))], ...
       ['PSL Power: ', num2str(pos(3))]};