for lt in steps:
    fc_steps          = forecast_anom.sel(step=pd.Timedelta(lt,'D')) #loop through each month
    hc_steps          = hindcast_anom.sel(step=pd.Timedelta(lt,'D'))
    era_steps          = era_anom.sel(step=pd.Timedelta(lt,'D'))         
    
    dim               = 'validation_time.month'    
    
    fc_group          = list(fc_steps.groupby(dim)) 
    hc_group          = list(hc_steps.groupby(dim))
    era_group          = list(era_steps.groupby(dim))
  
    for m,(mf,fcdata) in enumerate(fc_group): #loop through each month
        mh,hcdata          = hc_group[m]
        me,eradata          = era_group[m]
        
        dim                 = 'validation_time.day'
        hc_group_day        = list(hcdata.groupby(dim))
        fc_group_day        = list(fcdata.groupby(dim))
            
        hcmax = []
        hcmin = []
            
        fcmax = []
        fcmin = []
        fc50 = []
        fc75 = []
        fc25 = []
        plotdata_test = []
        plotdata = []
            
        for hcm,(hcmm,hcdata_sel_day) in enumerate(hc_group_day): #loop through each month
            fcm, fcdata_sel_day = fc_group_day[hcm]
            plotdata_test = []
            plotdata = []    
            max_hc = hcdata_sel_day.max('year').max('member').drop('time').drop('step').drop('validation_time')
            min_hc = hcdata_sel_day.min('year').min('member').drop('time').drop('step').drop('validation_time')
            hcmax.append(max_hc)
            hcmin.append(min_hc)
                
            max_fc = fcdata_sel_day.quantile(1,dim='member')
            plotdata_test.append(max_fc)
            min_fc = fcdata_sel_day.quantile(0,dim='member')
            plotdata_test.append(min_fc)
            p50_fc = fcdata_sel_day.quantile(0.5,dim='member')
            plotdata_test.append(p50_fc)
            p25_fc = fcdata_sel_day.quantile(0.25,dim='member')
            plotdata_test.append(p25_fc)
            p75_fc = fcdata_sel_day.quantile(0.75,dim='member')
            plotdata_test.append(p75_fc)
            
            plotdata = xr.concat(plotdata_test,dim='quantile')
            
            
            plt.close()
                 
            varplot = plotdata 
            levels_plot = np.linspace(-10,10,21)
            levels_cbar = np.linspace(-10,10,11)
            plot_title  = 't2m anomaly
            fname       = 't2m_anomaly'
            label_text  = 'K'

            im = varplot.plot(
  x               = 'lon',
  y               = 'lat',
  col              = 'time',
  col_wrap         = 3,
  levels           = levels_plot,
  subplot_kws      = dict(projection=ccrs.PlateCarree()),
  transform        = ccrs.PlateCarree(),
  cbar_kwargs      = {'label': label_text, 'ticks': levels_cbar}
        )

  
for i,ax in enumerate(im.axes.flat):
  ax.coastlines(resolution = '10m', 
                color      = 'black',
                linewidth  = 0.2)
  ax.set_title(dataopen_ymonmean.time[i].dt.month.values)
        
plt.suptitle(plot_title)

plt.savefig(fname+'.png',dpi='figure',bbox_inches='tight')
plt.close()
print('Figure stored at: '+fname+'.png')  
              
              
              
              
            fcmax.append(max_fc)
            fcmin.append(min_fc)
            fc50.append(p50_fc)
            fc75.append(p75_fc)
            fc25.append(p25_fc)
                
                
        hcmax_plot = xr.concat(hcmax,dim='validation_time') 
        hcmin_plot = xr.concat(hcmin,dim='validation_time') 
        fcmax_plot = xr.concat(fcmax,dim='validation_time')
        fcmin_plot = xr.concat(fcmin,dim='validation_time')
        fc50_plot = xr.concat(fc50,dim='validation_time')
        fc75_plot = xr.concat(fc75,dim='validation_time')
        fc25_plot = xr.concat(fc25,dim='validation_time')
        
        
        
        
        #fcdata_sel_df = fcdata_sel_day.drop('step').to_dataframe().reset_index(level = 1,drop=True).reset_index(level=0)
        hcdata_sel_df = hcdata_sel.drop('step').drop('year').to_dataframe().reset_index(level=0,drop=True).reset_index(level=1,drop=True).reset_index(level=0)
        #hcdata_sel_df = hcdata_sel.drop('step').to_dataframe().reset_index(level=0,drop=True).reset_index(level=0).reset_index(level=0)
        eradata_sel_df = eradata.drop('step').to_dataframe().reset_index(level=0,drop=True)
        
        
            plt.close()
            fig,ax=plt.subplots()
            ax.plot(eradata_sel_df.validation_time,np.zeros(eradata_sel_df.validation_time.shape[0]),alpha=0.1)
            hclabel = ax.fill_between(eradata_sel_df.validation_time, hcmax_plot.squeeze(),hcmin_plot.squeeze(),alpha=0.1,zorder=30, facecolor='gray', label='hindcast')
            fclabel = ax.fill_between(eradata_sel_df.validation_time, fcmax_plot.squeeze(),fcmin_plot.squeeze(),alpha=0.1,zorder=30, facecolor='blue', label='forecast')
            fclabel_25_75 = ax.fill_between(eradata_sel_df.validation_time, fc75_plot.squeeze(),fc25_plot.squeeze(),alpha=0.1,zorder=30, facecolor='purple',label='forecast Q1-Q3')
            fclabel50 = ax.plot(eradata_sel_df.validation_time,fc50_plot.squeeze(),color='purple', label='forecast median')
            era = ax.plot(eradata_sel_df.validation_time,eradata_sel_df.t2m, color='red',label='ERA5')
            ax.legend(loc='lower right')
            x_dates = eradata_sel_df['validation_time'].dt.strftime('%m-%d').sort_values().unique()
            ax.set_xticklabels(labels=x_dates, rotation=45, ha='right')
            ax.set_ylim([-8, 8]) 
            figname = 'HC_FC_step_' + str(lt.days) + '_month_' + str(mf) + '_' + reg + '_full_ens.png'
            plt.savefig(figname,dpi='figure',bbox_inches='tight')
            
            
