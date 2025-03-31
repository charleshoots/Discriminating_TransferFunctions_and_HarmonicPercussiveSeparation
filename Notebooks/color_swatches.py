# lines=['#7370cb','royalblue','dodgerblue','black','red']
# lineslabel=['inliers?','inliers?','inliers?','median','outliers']
# fig,axes=plt.subplots(1,3,figsize=(10,10))
# ax = axes[0]
# d=lines
# y=np.linspace(0,1.1,len(d))
# [ax.axhline(yy,c=lc,linewidth=140) for yy,lc in zip(y,d)]
# [ax.text(0.5,yy-0.05,f'{lt}\n{c}',color='white',horizontalalignment='center',fontweight='bold',fontsize=15) for yy,lt,c in zip(y,lineslabel,d)]
# ax.set_title('Scatter',fontweight='bold',fontsize=15)
# ax.tick_params(left=False, right=False, top=False, bottom=False, labelleft=False, labelbottom=False)

# ax = axes[1];x='instrument'
# lineslabel=list(ColorStandard[x].keys())
# d=[ColorStandard[x][i] for i in lineslabel]
# y=np.linspace(0,1.1,len(d))
# [ax.axhline(yy,c=lc,linewidth=100) for yy,lc in zip(y,d)]
# [ax.text(0.5,yy,f'{lt}\n{c}',color='white',horizontalalignment='center',fontweight='bold',fontsize=15,verticalalignment='center') for yy,lt,c in zip(y,lineslabel,d)]
# ax.set_title(x,fontweight='bold',fontsize=15)
# ax.tick_params(left=False, right=False, top=False, bottom=False, labelleft=False, labelbottom=False)

# ax = axes[2];x='network'
# lineslabel=list(ColorStandard[x].keys())
# d=[ColorStandard[x][i] for i in lineslabel]
# y=np.linspace(0,1.1,len(d))
# [ax.axhline(yy,c=lc,linewidth=60) for yy,lc in zip(y,d)]
# [ax.text(0.5,yy,f'{lt}\n{c}',color='white',horizontalalignment='center',fontweight='bold',fontsize=15,verticalalignment='center') for yy,lt,c in zip(y,lineslabel,d)]
# ax.set_title(x,fontweight='bold',fontsize=15)
# ax.tick_params(left=False, right=False, top=False, bottom=False, labelleft=False, labelbottom=False)

# # ax = axes[3];x='components'
# # lineslabel=list(ColorStandard[x].keys())
# # d=[ColorStandard[x][i] for i in lineslabel]
# # y=np.linspace(0,1.1,len(d))
# # [ax.axhline(yy,c=lc,linewidth=90) for yy,lc in zip(y,d)]
# # [ax.text(0.5,yy,lt,color='white',horizontalalignment='center',fontweight='bold',fontsize=15,verticalalignment='center') for yy,lt in zip(y,lineslabel)]
# # ax.set_title(x,fontweight='bold',fontsize=15)
# # ax.tick_params(left=False, right=False, top=False, bottom=False, labelleft=False, labelbottom=False)

# cmap = plt.cm.magma  # Choose the colormap (e.g., 'viridis')
# norm = mpl.colors.Normalize(vmin=0, vmax=1)  # Define the range of the colorbar

# # Create a ScalarMappable and add the colorbar
# sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
# sm.set_array([])  # No data is needed
# # cbar = fig.colorbar(sm, cax=ax)  # Add the colorbar to the axes
# bottom = fig.subplotpars.bottom
# top = fig.subplotpars.top
# # Add a new axis for the colorbar with the same height as the subplots
# cbar_ax = fig.add_axes([0.92, bottom, 0.1, top - bottom]) 
# cbar = fig.colorbar(sm, cax=cbar_ax,label='Spectrograms')  # Add the colorbar to the new axis
# # Customize the colorbar label fontweight and fontsize
# cbar.ax.yaxis.label.set_fontweight('bold')
# cbar.ax.yaxis.label.set_fontsize(15)
# cbar.ax.yaxis.label.set_position((1, 0.5))
# fig.suptitle('Color Swatch',fontsize=20,fontweight='bold')