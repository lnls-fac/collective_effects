function plot_impedances(w, budget, mult_beta, escalax, escalay, save, name)
% function plot_impedances(w, budget, mult_beta, escalax, escalay, save, name)
% Funcao que plota impedancias:
%
% Inputs:
%   escalax   : define se escala x sera logaritmica, 'log', ou 'linear'.  Default = 'log';
%   escalay   : define se escala y sera logaritmica, 'log', ou 'linear'.  Default = 'log';
%   mult_beta : plotar Zt ou betat*Zt. true ou false.  Default = true;
%   save      : salvar ou não os gráficos. Default = false;
%   name      : name of the file to save the plots.
%   budget: array de celulas, sendo que cada celula contem uma estrutura
%           com os seguintes campos obrigatorios:
%       name     : nome do elemento ao qual a impedancia pertence
%       Zl       : impedancia longitudinal [Ohm]
%       Zv       : impedancia vertical [Ohm/m]
%       Zh       : impedancia horizontal [Ohm/m]
%       w        : frequencia angular [rad/s]
%       quantity : numero desse elemento no anel
%       betax    : beta horizontal media nas posicoes desses elementos
%       betay    : beta vertical media nas posicoes desses elementos

if ~exist('escalax','var'),   escalax = 'log';  end
if ~exist('escalay','var'),   escalay = 'log';  end
if ~exist('mult_beta','var'), mult_beta = true; end
if ~exist('save','var'),      save = false; name=''; end

Zh = zeros(length(budget),length(w));
Zv = Zh;
Zl = Zh;
if mult_beta
    for i=1:length(budget)
        Zh(i,:) = budget{i}.Zh*budget{i}.quantity*budget{i}.betax;
        Zv(i,:) = budget{i}.Zv*budget{i}.quantity*budget{i}.betay;
        Zl(i,:) = budget{i}.Zl*budget{i}.quantity;
    end
    labelh = '\beta_x Z_{x}';
    labelv = '\beta_y Z_{y}';
else
    for i=1:length(budget)
        Zh(i,:) = budget{i}.Zh*budget{i}.quantity;
        Zv(i,:) = budget{i}.Zv*budget{i}.quantity;
        Zl(i,:) = budget{i}.Zl*budget{i}.quantity;
    end
    labelh = 'Z_{x}';
    labelv = 'Z_{y}';
end

DispNames = getcellstruct(budget,'name',1:length(budget));
labelx = '\omega [rad/s]';

labelL = 'Z_{||}';
create_figures(w, Zl, escalax, escalay, labelx, labelL, DispNames, save,['l_' name]);
create_figures(w, Zh, escalax, escalay, labelx, labelh, DispNames, save,['x_' name]);
create_figures(w, Zv, escalax, escalay, labelx, labelv, DispNames, save,['y_' name]);



function create_figures(w, Z, scalex, scaley, labelx, labely, DispNames, save, name)

units = '\Omega';
if strcmp(scaley,'linear'),
    maxZ = max(max(abs(Z)));
    if maxZ > 1e6,     Z = Z/1e6; units = 'M\Omega';
    elseif maxZ > 1e3, Z = Z/1e3; units = 'k\Omega'; end
end

if strcmp(scalex,'log')
    indx = w > 0 ;
    w = w(:,indx);
    Z = Z(:,indx);
end
w = repmat(w,size(Z,1),1);
indp = true(size(Z));
indn = [];
IZpos = imag(-Z);
IZneg = [];
if strcmp(scaley,'log')
    indn = imag(-Z) < 0;
    indp = imag(-Z) > 0;
    IZpos(indn) = 0;
    IZneg = imag(Z);
    IZneg(indp) = 0;
end

colors = {[0, 0, 1],...
          [0, 0.5, 0],...
          [1, 0, 0],...
          [0, 0.75, 0.75],...
          [0.75, 0, 0.75],...
          [0.75, 0.75, 0],...
          [0.25, 0.25, 0.25],...
          [1, 0.69, 0.39]};
lc = length(colors);

% Create figure
f1 = figure('Position',[1,1,1100,680]);


% Real Axes
axes1 = axes('Parent',f1,'Position',[0.088 0.103 0.885 0.432],...
    'FontSize',16, 'YScale',scaley,'YMinorTick','on',...
                   'XScale',scalex,'XMinorTick','on');
box(axes1,'on'); hold(axes1,'all'); grid(axes1,'on');
plot1 = zeros(1,size(Z,1));
for i=1:size(Z,1)
    plot1(i) = plot(w(i,:),abs(real(Z(i,:))),'Parent',axes1,...
                 'LineWidth',2,'Color',colors{mod(i-1,lc)+1});
end
xlabel(labelx,'FontSize',16);
ylabel(['Re(',labely,') [',units,']'],'FontSize',16);
xlim(axes1,[min(w(:));            max(w(:))]);
ylim(axes1,[min(abs(real(Z(:)))); max(abs(real(Z(:))))]);



% Imaginary Axes
axes2 = axes('Parent',f1,'Position',[0.088 0.536 0.885 0.432],...
    'FontSize',16,'YScale',scaley,'YMinorTick','on',...
                  'XScale',scalex,'XMinorTick','on',...
    'XTickLabel',{''});
box(axes2,'on');hold(axes2,'all');grid(axes2,'on');

for i=1:size(Z,1)
    plot(w(i,:),IZpos(i,:),'Parent',axes2,'LineWidth',2,...
                    'DisplayName',DispNames{i}, 'Color',get(plot1(i),'Color'));
end
legend(axes2,'show','Location','Best');
if strcmp(scaley,'log')
    for i=1:size(Z,1)
        plot(w(i,:),IZneg(i,:),'Parent',axes2,'LineWidth',2,...
            'Color',get(plot1(i),'Color'),'LineStyle','--');
    end
end
ylabel(['-Im(',labely,') [',units,']'],'FontSize',16);
xlim(axes2,[min([w(indp);w(indn)]), max([w(indp);w(indn)])]);
ylim(axes2,[min([IZpos(indp); IZneg(indn)]), ...
            max([IZpos(indp); IZneg(indn)])]);

if save
    if ~exist('name','var')
        name = 'no_name';
    end
    saveas(f1,['Z', name, '.fig']);
end

