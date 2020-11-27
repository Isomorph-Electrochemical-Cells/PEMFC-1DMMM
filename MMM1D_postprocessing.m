clear all

%% Compute solution
[I,U,SOL,x, Lsum, Np, Neq, domains] = MMM1D;

%% Postprocess

% PLOT POTENTIALS AND FLUXES
close all
fig_names = {'Potentials', 'Fluxes'};
unit_scale = [1 1 1 1 1 1 1 1;
              1e-4 1e-4 1e-4 1e2 1e2 1e2 1e2 1e2];
quantity = {'$\phi_e [\mathrm{V}]$', '$\phi_p [\mathrm{V}]$', '$T [\mathrm{K}]$', '$\lambda [-]$', '$w_{H_2O} [-]$', '$w_{O_2} [-]$', '$s [-]$', '$P_\mathrm{gas} [\mathrm{Pa}]$';
            '$j_e [\mathrm{A/cm}^2]$', '$j_p [\mathrm{A/cm}^2]$', '$j_T [\mathrm{W/cm}^2]$', '$j_\lambda [\mu\mathrm{mol/cm}^2\mathrm{s}]$', ...
            '$j_{H_2O} [\mu\mathrm{g/cm}^2\mathrm{s}]$', '$j_{O_2} [\mu\mathrm{g/cm}^2\mathrm{s}]$', '$j_s [\mu\mathrm{mol/cm}^2\mathrm{s}]$',...
            '$\rho_\mathrm{gas} \cdot u_\mathrm{gas} [\mathrm{\mu g}/\mathrm{cm}^2\mathrm{s}]$'};
c = jet(Np);
for m = 1:2
    figure('Name', fig_names{m},'units','centimeters','position',[0 20-(m-1)*20 35 12])
    for n = 1:Neq
        subplot(2,4,n)
        box on
        hold on
        us = unit_scale(m,n);
        for k = 1:Np
            plot(SOL{k}.x*1e6, SOL{k}.y(2*(n-1)+m,:)*us, 'Color', c(k,:), 'DisplayName', [num2str(U(k)) ' V'])
        end
        xlim([Lsum(find(domains(n,:),1,'first')) Lsum(find(domains(n,:),1,'last')+1)]*1e6)
        ylim(ylim)
        xlabel('x [um]')
        ylabel(quantity(m,n),'Interpreter','latex')
        for x = Lsum(2:end-1)
            l = line([x x]*1e6, ylim, 'Color', 'k');
            set(get(get(l, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off')
        end
    end
end

% PLOT POLARIZATION CURVE
figure('Name', 'Polarization curve','units','centimeters','position',[35 20 20 15])
hold on
P = I.*U;
yyaxis left
plot(I, U,'b')
ylabel('Cell voltage [V]')
yyaxis right
plot(I, P,'r')
xlabel('Current density [A/cm^2]')
ylabel('Power density [W/cm^2]')

xlim([0 max(I)])
ylim([0 max([U P])])