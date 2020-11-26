clear all
close all

[I,U,SOL,x, Lsum, Np, Neq, domains] = MMM1D;

%%
% PLOT SOLUTION
fig_names = {'Potentials', 'Fluxes'};
unit_scale = [1 1 1 1 1 1 1;
              1e-4 1e-4 1e-4 1e2 1e2 1e2 1e2];
quantity = {'\phi_e [V]', '\phi_p [V]', 'T [K]', '\lambda', 'w_{H_2O}', 'w_{O_2}', 's';
            'j_e [A/cm^2]', 'j_p [A/cm^2]', 'j_T [W/cm^2]', 'j_\lambda [umol/cm^2s]', ...
            'j_{H_2O} [ug/cm^2s]', 'j_{O_2} [ug/cm^2s]', 'j_s [umol/cm^2s]'};
c = jet(Np);
for m = 1:2
    figure('Name', fig_names{m})
    for n = 1:Neq
        subplot(3,3,n)
        box on
        hold on
        us = unit_scale(m,n);
        for k = 1:Np
            plot(SOL{k}.x*1e6, SOL{k}.y(2*(n-1)+m,:)*us, 'Color', c(k,:), 'DisplayName', [num2str(U(k)) ' V'])
        end
        xlim([Lsum(find(domains(n,:),1,'first')) Lsum(find(domains(n,:),1,'last')+1)]*1e6)
        ylim(ylim)
        xlabel('x [um]')
        ylabel(quantity(m,n))
        for x = Lsum(2:end-1)
            l = line([x x]*1e6, ylim, 'Color', 'k');
            set(get(get(l, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off')
        end
    end
end

% PLOT POLARIZATION CURVE
figure('Name', 'Polarization curve')
hold on
P = I.*U;
fnplt(cscvn([I; U]))
fnplt(cscvn([I; P]))
xlabel('Current density [A/cm^2]')
ylabel({'Cell voltage [V]'; 'Power density [W/cm^2]'})
xlim([0 max(I)])
ylim([0 max([U P])])