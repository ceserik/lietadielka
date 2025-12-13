figure(4)
h = bodeplot(podelny_system(3,1));
setoptions(h,'PhaseVisible','on')
h.DisplayName = 'System 1';
hold on

h2 = bodeplot(podelny_system(3,2));
h2.DisplayName = 'System 2';

legend show
