function Z = lnls_load_fitted_impedance(w, fold, plano, codigo)
% Carrega impedancias de arquivos texto:
% Inputs:
% fold = caminho completo do diretorio onde encontra-se o arquivo
% kind = 'Longitudinal' para Longitudinal ou 'V' para vertical e
% 'Horizontal'para horizontal
%
% Outputs:
%
% Z    = Impedancia [Ohm ou Ohm/m]
% w    = frequencia angular [rad/s]
%
% Assume-se que os dados estao organizados em duas colunas, sendo a
% primeira de frequencia em GHz e a segunda de Impedancia, no SI.
% Ainda, assume-se um padrao para o nome do arquivo a ser carregado. Para
% mais detalhes, ver o codigo.
%
% Fernando 2013-03-21

if strcmpi(plano,'Longitudinal')
    calc_ressonator = @lnls_calc_impedance_longitudinal_resonator;
else
    calc_ressonator = @lnls_calc_impedance_transverse_resonator;
end

filename = ['Zfit_' plano '_' codigo '.dat'];

try
    fullname = fullfile(fold, filename);
    A = importdata(fullname);
end

if ~exist('A','var')
    disp('Unable to read impedance file');
    return;
end

Z = calc_ressonator(A.data(:,1), A.data(:,3), A.data(:,2)*1e9, w);

