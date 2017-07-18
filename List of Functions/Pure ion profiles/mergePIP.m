function mergePIP(PIPs, stringTitle)
if nargin == 1, stringTitle = ''; end
InputFig = figure(                          ...
    'Visible'          , 'on'             , ...
    'Name'             , stringTitle      , ...
    'Toolbar'          , 'none'           , ...
    'MenuBar'          , 'none'           , ...
    'Units'            , 'normalized'     , ...
    'WindowStyle'      , 'normal'         , ...
    'Resize'           , 'on');

AxisH1 = axes(InputFig ,...
    'Units'            , 'normalized'     , ...
    'OuterPosition'    , [0 0 1/2 1]);

rpts = length(PIPs);
hold on
cFOM = [];
for ii = 1:rpts
    plot(AxisH1, PIPs{ii}.x, PIPs{ii}.y, 'k')
    cFOM = [cFOM; PIPs{ii}.FOM];
end
hold off

xlabel([PIPs{ii}.AxisX.Label, ' / ', PIPs{ii}.AxisX.Unit]);
ylabel([PIPs{ii}.AxisY.Label, ' / ', PIPs{ii}.AxisY.Unit]);

f1 = ['%.', num2str(signFig(std(cFOM(:,1)))), 'f']; 
str1 = ['Max Intensities :'];
f2 = ['%.', num2str(signFig(std(cFOM(:,2)))), 'f'];
str2 = ['Times @ IMax  :'];
f3 = ['%.', num2str(signFig(std(cFOM(:,3)))), 'f'];
str3 = ['Peaks Area      :'];
f4 = ['%.', num2str(signFig(std(cFOM(:,4)))), 'f'];
str4 = ['Peaks centre   :'];
f5 = ['%.', num2str(signFig(std(cFOM(:,5)))), 'f'];
str5 = ['Peaks variance:'];
f9 = ['%.', num2str(signFig(std(cFOM(:,9)))), 'f'];
str9 = ['Accurate mass ='];
for ii = 1:rpts
    str1 = [str1, ' ', f1];
    str2 = [str2, ' ', f2];
    str3 = [str3, ' ', f3];
    str4 = [str4, ' ', f4];
    str5 = [str5, ' ', f5];
    str9 = [str9, ' ', f9];
end

String4Edit{1} = '  FIGURES OF MERITS';
String4Edit{2} = '_____________________';
String4Edit{3} = '';

String4Edit{4} = sprintf([str1, ' %s'], cFOM(:,1), PIPs{1}.AxisZ.Unit);
String4Edit{5} = sprintf([str2, ' %s'], cFOM(:,2), PIPs{1}.AxisX.Unit);
String4Edit{6} = sprintf([str3, ' %s'], cFOM(:,3),  'a.u.');
String4Edit{7} = sprintf([str4, ' %s'], cFOM(:,4), PIPs{1}.AxisX.Unit);
String4Edit{8} = sprintf([str5, ' %s^2'], cFOM(:,5), PIPs{1}.AxisX.Unit);
String4Edit{10} = sprintf([str9 ' a.m.u.'],  cFOM(:,9));

EditH1 = uicontrol(InputFig ,...
    'Style'              , 'Edit'      , ...
    'Visible'            , 'on'        , ...
    'Units'              , 'normalized', ...
    'Max'                , 2           , ...
    'String'             , String4Edit ,...
    'HorizontalAlignment', 'left', ...
    'Position'           , [1/2 0 1/2 1]);
end
