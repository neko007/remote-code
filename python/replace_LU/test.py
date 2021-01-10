#%%
import matplotlib.pyplot as plt 

plt.figure()
xx,yy = np.meshgrid(range(100),range(100))
gg = np.sqrt(xx*2+yy*2)
CS = plt.contourf(gg) #, cc, zz_miss)
proxy = [plt.Rectangle((0,0),1,1,fc = pc.get_facecolor()[0]) for pc in CS.collections]

plt.legend(proxy, [str(i) for i in range(10)], frameon=False)
plt.xlabel('gamma')
plt.ylabel('C = 1 / lambda')

plt.show()

