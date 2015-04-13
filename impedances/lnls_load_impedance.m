function [Z w] = lnls_load_impedance(fold, plano, codigo)
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
multi = 1000; % no GdFdl os resultados saem em kiloOhm/m para Zt
if strcmpi(plano,'Longitudinal')
    pl = 'Zlong';
    multi = 1;
elseif strcmpi(plano, 'Vertical')
    pl = 'Zydip';
else
    pl = 'Zxdip';
end

filename = [pl codigo '.txt'];

try
    fullname = fullfile(fold, ['Re' filename]);
    data = load(fullname);
    ReZ = multi*data(:,2);
    Z=ReZ;
    w = 2*pi*1e9*data(:,1);
end

try
    fullname = fullfile(fold, ['Im' filename]);
    data = load(fullname);
    ImZ = multi*data(:,2);
    Z=1i*ImZ;
    w = 2*pi*1e9*data(:,1); % + 0.0/1.7e-6*2*pi*864
end

if (exist('ImZ','var') && exist('ReZ','var')), Z = ReZ + 1i*ImZ;end
if ~exist('Z','var'), disp('Unable to read impedance file');end

w = [-flipud(w); w(2:end)]';
if strcmpi(plano,'Longitudinal')
    Z = [flipud(conj(Z)); Z(2:end)]';
else
    Z = [-flipud(conj(Z)); Z(2:end)]';
end

