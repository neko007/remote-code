#%%
from Module import *

if __name__ == '__main__':
    ### 站点观测
    df = pd.read_excel('data/纳木错站2017-2018.xlsx', index_col=0)
    obs_NamCo = df.loc[pd.date_range('2017-03-01',
    '2017-12-01')]['降水量']

    # plt.style.use(['science', 'ieee'])
    fig, ax = plt.subplots(figsize=(8,4))
    ax.plot(obs_NamCo.index, obs_NamCo, lw=1.2, c='k', label='OBS')#, marker='o', mfc=None, markersize=1.4)
    ax.set_title('NamCo Station', loc='left', y=0.85, x=0.04, fontsize=14, weight='bold')
    ax.set_ylabel('Precipitation / $\mathrm{mmd^{-1}}$', fontsize=13)
    ax.legend(fontsize=14, loc='upper right', frameon=False)
    ax.axvspan('2017-06-01', '2017-07-01', alpha=0.25, color='grey')

    import matplotlib.dates as mdate    
    ax.xaxis.set_major_formatter(mdate.DateFormatter('%m-%d'))

    fig.savefig('fig-test-9.4/NamCo_station_obs.jpg', dpi=300)