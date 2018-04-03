import matplotlib.pyplot as plt
import numpy as np

colors = ['magenta', 'red','orange','yellow','green','blue', 'cyan','violet','black']

h2o_x = np.asarray([296.0, 350.0, 400.0, 450.0, 500.0, 600.0, 800.0, 1000.0])
h2o_y = np.asarray([1396.8, 1411.6, 1427.4, 1444.9, 1464.1, 1506.8, 1605.8, 1719.8])

fig, ax = plt.subplots()
fit = np.polyfit(h2o_x, h2o_y, deg=2)
ax.plot(h2o_x, fit[0]*h2o_x**2 + fit[1]*h2o_x + fit[2], color=colors[0])
ax.scatter(h2o_x, h2o_y,color=colors[0],label='h2o')



d2o_x = np.asarray([296.0, 350.0, 400.0, 450.0, 500.0, 600.0, 800.0, 1000.0])
d2o_y = np.asarray([940.91, 961.62, 982.93, 1006.1, 1030.9, 1085.1, 1209, 1350])

fit = np.polyfit(d2o_x, d2o_y, deg=2)
ax.plot(d2o_x, fit[0]*d2o_x**2 + fit[1]*d2o_x + fit[2], color=colors[1])
ax.scatter(d2o_x, d2o_y,color=colors[1],label='d2o')


be_x = np.asarray([296, 400, 500, 600, 700, 800, 1000, 1220])
be_y = np.asarray([405.64, 484.22, 568.53, 657.66, 749.69, 843.63, 1035., 1229.3])

fit = np.polyfit(be_x, be_y, deg=2)
ax.plot(be_x, fit[0]*be_x**2 + fit[1]*be_x + fit[2], color=colors[2])
ax.scatter(be_x, be_y,color=colors[2],label='beryllium')


graphite_x = np.asarray([296, 400, 500, 600, 700, 800, 1000, 1200, 1600, 2000])
graphite_y = np.asarray([713.39, 754.68, 806.67, 868.38, 937.64, 1012.7, 1174.9, 1348.2, 1712.9, 2091])

fit = np.polyfit(graphite_x, graphite_y, deg=3)
ax.plot(graphite_x, fit[0]*graphite_x**3 + fit[1]*graphite_x**2 + fit[2]*graphite_x + fit[3], color=colors[3])
ax.scatter(graphite_x, graphite_y,color=colors[3],label='graphite')


benzene_x = np.asarray([296, 350, 400, 450, 500, 600, 800, 1000])
benzene_y = np.asarray([1165.9, 1177.8, 1191.4, 1207.7, 1226, 1268.7, 1373.4, 1497.7])

fit = np.polyfit(benzene_x, benzene_y, deg=2)
ax.plot(benzene_x, fit[0]*benzene_x**2 + fit[1]*benzene_x + fit[2], color=colors[4])
ax.scatter(benzene_x, benzene_y,color=colors[4],label='benzene')


zr_zrh_x = np.asarray([296, 400, 500, 600, 700, 800, 1000, 1200])
zr_zrh_y = np.asarray([317.27, 416.29, 513.22, 611.12, 709.60, 808.43, 1006.8, 1205.7])

fit = np.polyfit(zr_zrh_x, zr_zrh_y, deg=2)
ax.plot(zr_zrh_x, fit[0]*zr_zrh_x**2 + fit[1]*zr_zrh_x + fit[2], color=colors[5])
ax.scatter(zr_zrh_x, zr_zrh_y,color=colors[5],label='zr_zrh')


h_zrh_x = np.asarray([296, 400, 500, 600, 700, 800, 1000, 1200])
h_zrh_y = np.asarray([806.79, 829.98, 868.44, 920.08, 981.82, 1051.1, 1205.4, 1373.4])

fit = np.polyfit(h_zrh_x, h_zrh_y, deg=2)
ax.plot(h_zrh_x, fit[0]*h_zrh_x**2 + fit[1]*h_zrh_x + fit[2], color=colors[6])
ax.scatter(h_zrh_x, h_zrh_y,color=colors[6],label='h_zrh')



beo_x = np.asarray([296, 400, 500, 600, 800, 1000, 1200])
beo_y = np.asarray([596.4, 643.9, 704.6, 775.3, 935.4, 1109.8, 1292.3])

fit = np.polyfit(beo_x, beo_y, deg=2)
ax.plot(beo_x, fit[0]*beo_x**2 + fit[1]*beo_x + fit[2], color=colors[7])
ax.scatter(beo_x, beo_y,color=colors[7],label='beo')


h_ch2_x = np.asarray([296, 350])
h_ch2_y = np.asarray([1222, 1239])     # Values gateff in thermr.f90
h_ch2_y = np.asarray([1204.4, 1215.1]) # Values in the osti.gove report 

fit = np.polyfit(h_ch2_x, h_ch2_y, deg=1)
ax.plot(h_ch2_x, fit[0]*h_ch2_x + fit[1], color=colors[8])
ax.scatter(h_ch2_x, h_ch2_y,color=colors[8],label='h_ch2')






plt.legend()
fig.show()
plt.show()
