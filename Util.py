
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs




# IMPORT DATA FROM NETCDF FILE


def open_mfdataset_plasim_monthly(*args, **kwargs):
    data = xr.open_mfdataset(*args, **kwargs)
    vals = data['time'].values
    res = vals - vals.astype(int)
    date_str = np.char.zfill(vals.astype(int).astype(str), 8)
    date = f'{date_str[0][0:4]}-{date_str[0][4:6]}'
    data["time"] = xr.date_range(date, periods=len(data.time), freq="MS")
    return data 





# WEIGHTED AVERAGE/ZONAL MEAN

def area_mean(data):
    return data.weighted(np.cos(data["lat"]/180*np.pi)).mean(["lat","lon"])



# PLOT 2D MAP

def plot_map_2D(data, title, units, nx, ny, nn, fig=None, cmap='viridis'):
    if fig is None:
        fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(ny, nx, nn, projection=ccrs.PlateCarree())
    data.plot(ax=ax, transform=ccrs.PlateCarree(), cmap=cmap, cbar_kwargs={'label': units})
    ax.coastlines()
    ax.set_title(title)
    return fig, ax

#plot.pcolormesh or plot.contour or plot.contourf

def v_plot_map_2D(data, title, units, nx, ny, nn, fig=None, cmap='viridis', **kwargs):
    if fig is None:
        fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(ny, nx, nn, projection=ccrs.PlateCarree())
    
    # Add vmin/vmax or other kwargs to data.plot
    data.plot(ax=ax, transform=ccrs.PlateCarree(), cmap=cmap, 
              cbar_kwargs={'label': units}, **kwargs)
    
    ax.coastlines()
    ax.set_title(title)
    return fig, ax



# Funzione aggiornata per plot_map_2D
def new_plot_map_2D(data, title, units, nx, ny, nn, fig=None, cmap='viridis', show_colorbar=True):
    # Usa la figura principale se fornita, altrimenti creane una nuova
    if fig is None:
        fig = plt.figure(figsize=(15, 20))  # Usa figsize globale per la figura
    ax = fig.add_subplot(ny, nx, nn, projection=ccrs.PlateCarree())
    # Crea il grafico senza la colorbar se 'show_colorbar' è False
    img = data.plot(ax=ax, transform=ccrs.PlateCarree(), cmap=cmap, cbar_kwargs={'label': units} if show_colorbar else {})
    ax.coastlines()
    ax.set_title(title)
    return fig, ax, img


# aggiungere funzione per plot 2D con land sea mask 
# cercare di aggiustare la question kwargs in modo tale da non scrivere dims, cftime, nested etc

def open_mfdataset_plasim_monthly2( concat_dim="time", combine="nested", use_cftime=True, decode_times=False, *args,  **kwargs):
    data = xr.open_mfdataset(*args, concat_dim, combine, use_cftime, decode_times, **kwargs)
    vals = data['time'].values
    res = vals - vals.astype(int)
    date_str = np.char.zfill(vals.astype(int).astype(str), 8)
    date = f'{date_str[0][0:4]}-{date_str[0][4:6]}'
    data["time"] = xr.date_range(date, periods=len(data.time), freq="MS")
    return data 







