% Script to load numeric impedance and fit it with known mathematical
% functions.

% load impedance
[Z w] = lnls_load_impedance('/home/ABTLUS/fernando.sa/MATLAB/Impedancias/Simulacoes/BPM','L');

%% plot impedance
% Create figure
figure1 = figure('Position',[994.0  2.0  927.0  972.0],...
    'XVisual',...
    '0x9b (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)');

% Create subplot
subplot1 = subplot(2,1,1,'Parent',figure1);
% Uncomment the following line to preserve the X-limits of the axes
xlim(subplot1,[0 w(end)]);
box(subplot1,'on');
hold(subplot1,'all');

% Create plot
plotReZ = plot(w(floor(end/2):end),real(Z(floor(end/2):end)),'Parent',subplot1);
% Create xlabel
xlabel('w [rad/s]');
% Create ylabel
ylabel('Re(Z_l) [Ohm]');

% Create subplot
subplot2 = subplot(2,1,2,'Parent',figure1);
% Uncomment the following line to preserve the X-limits of the axes
xlim(subplot2,[0 w(end)]);
box(subplot2,'on');
hold(subplot2,'all');
% Create plot
plotImZ = plot(w(floor(end/2):end),imag(Z(floor(end/2):end)),'Parent',subplot2);
% Create xlabel
xlabel('w [rad/s]');
% Create ylabel
ylabel('Im(Z_l) [Ohm]');


%% escolha dos parametros para o fitting
% Escolha do numero de ressonadores para fittar a impedancia;
prompt='Quantos picos voce quer fitar?';
name='numero de ressonadores';
numlines=1;
defaultanswer={'0'};
options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';
answer = inputdlg(prompt,name,numlines,defaultanswer,options);
num_peaks = str2double(answer);

% Escolha dos parametros iniciais para o fitting
prompt={'wr [Grad/s]', 'Rs [Ohm]','Q'};
name = sprintf('Parametros para o ressonador %d', 0);
numlines=1;
defaultanswer={'5','0','10000'};
options.Resize='on';
options.Position= [500.0  500.0  200  100];
options.WindowStyle='normal';
options.Interpreter='tex';
% Rs = zeros(1,num_peaks);
% wr = zeros(1,num_peaks);
% Q  = zeros(1,num_peaks);
% 
% for m=1:num_peaks
%     name = sprintf('Parametros para o ressonador %d', m);
%     answer = inputdlg(prompt,name,numlines,defaultanswer,options);
%     if isempty(answer)
%         num_peaks = m-1;
%         wr = wr(1:num_peaks);
%         Rs = Rs(1:num_peaks);
%         Q  = Q(1:num_peaks);
%         break;
%     else
%         wr(m)  = str2double(answer(1))*1e9;
%         Rs(m)  = str2double(answer(2));
%         Q(m)   = str2double(answer(3));
%     end
% end


%% faz o fitting;

Zl = lnls_calc_impedance_longitudinal_resonator(Rs, Q, wr, w);
% Create plot
plotReZl = plot(w(floor(end/2):end),real(Zl(floor(end/2):end)),'Parent',subplot1);
plotImZl = plot(w(floor(end/2):end),imag(Zl(floor(end/2):end)),'Parent',subplot2);
drawnow;
sleep(2);

% funcao peso (fp) para as frequencias:
fp = real(lnls_calc_impedance_longitudinal_resonator(1e3./Rs, Q, wr, w))';
% maximo erro relativo permitido:
err = 0.01;
% Define a semente dos numeros aleatorios:
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));

% calcula residuo:
residue = abs(Zl-Z)*fp;

for i=1:num_peaks
    par = [Rs; Q; wr];
    mask = zeros(size(par));
    mask(:,i) = 1;
    improvement = 100;
    fprintf('Ressonator: %d\n',i);
    while improvement > 1
        error = mask.*(err*2*(rand(size(par))-0.5));
        new_par = par.*(1+error);
        new_Zl  = lnls_calc_impedance_longitudinal_resonator(new_par(1,:), new_par(2,:), new_par(3,:), w);
        new_residue = abs(new_Zl - Z)*fp;
        change = residue - new_residue;
        if change > 0
            fprintf('Melhora: %f\n',change);
            par = new_par;
            residue = new_residue;
            improvement = change;
            % Create plot
            delete(plotReZl);
            delete(plotImZl);
            plotReZl = plot(w(floor(end/2):end),real(new_Zl(floor(end/2):end)),'Parent',subplot1);
            plotImZl = plot(w(floor(end/2):end),imag(new_Zl(floor(end/2):end)),'Parent',subplot2);
            drawnow;
        end
    end
    Rs(i) = par(1,i);
    Q(i)  = par(2,i);
    wr(i) = par(3,i);
end
