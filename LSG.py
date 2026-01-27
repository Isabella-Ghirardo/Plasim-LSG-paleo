
# LSG AREA + LAYER THICKNESS

# Helper function to compute layer depths (works only properly with full z-axis!)


from queue import Full
import numpy as np
import xarray as xr


def layer_weights(dep,d,d_top=0):
    # Format input arrays
    try:
        d=d.values
    except AttributeError:
        pass
    try:
        dep=dep.item()
    except AttributeError:
        pass
    # Return 0 if top is below bottom
    if dep<d_top:
        return xr.DataArray(
        data=np.zeros(len(d)),dims=["depth"],coords={'depth':d},name='layer_depths')
    # Else
    d_sep=(d[1:]+d[:-1])/2
    layer_bnds=np.asarray([d_top]+list(d_sep[np.logical_and(d_sep>d_top,d_sep<dep)])+[dep])
    d_layers=layer_bnds[1:]-layer_bnds[:-1]
    d_layers_all=np.zeros(len(d))
    if dep<=d_sep[-1]:
        d_layers_all[np.argmax(d_sep>d_top):np.argmax(d_sep>=dep)+1]=d_layers
    else:
        if d_top<d_sep[-1]:
            d_layers_all[np.argmax(d_sep>d_top):]=d_layers
        else:
            d_layers_all[-1]=d_layers[-1]
    return xr.DataArray(
        data=d_layers_all,dims=["depth"],coords={'depth':d},name='layer_depths')





# LSG AREA + LAYER THICKNESS (lsg_area_volume.ipynb)

{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "fa38da45-e366-410b-8cef-d8846c3ab10a",
   "metadata": {},
   "outputs": [],
   "source": [
    "import numpy as np\n",
    "import xarray as xr"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "id": "b08af4c0-5540-4206-b74a-2f8ae52dae11",
   "metadata": {},
   "outputs": [],
   "source": [
    "def lsg_area(wet, dlat=2.5, dlon=2.5):\n",
    "    # NB: This function only works correctly on _regular_ lat/lon grid, where the spacing between grid points is equal\n",
    "    deg_length = np.pi*6371000/180 # m/degree (at the equator)\n",
    "    cell_area = wet*np.cos(wet[\"lat\"]/180*np.pi)*dlat*dlon*deg_length**2\n",
    "    return cell_area"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "id": "2263276f-9cd2-48f7-8f6d-ace56fd05fd6",
   "metadata": {},
   "outputs": [],
   "source": [
    "def layer_thickness(depth_axis, upper=0, lower=6500):\n",
    "    thickness_mid = list((depth_axis.values[2:]-depth_axis.values[:-2])/2)\n",
    "    thickness_top = depth_axis.values[0]-upper\n",
    "    thickness_bottom = lower-depth_axis.values[-1]\n",
    "    return xr.DataArray(data=[thickness_top]+thickness_mid+[thickness_bottom], coords={\"depth\": depth_axis})"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "id": "29d86e9a-bca5-4ea0-baf5-9f89bc0ebcac",
   "metadata": {},
   "outputs": [],
   "source": [
    "ds_lsg = xr.open_dataset(\"/home/omehling/work/smalltests/lsg_testout.nc\") # regridded on r144x72 grid"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "id": "bbdb6a28-b901-45ff-b9af-e6bf260f2ec4",
   "metadata": {},
   "outputs": [],
   "source": [
    "area = lsg_area(ds_lsg[\"wet\"].isel(time=0)) # Leave out .isel(time=0) if the mask is time-dependent\n",
    "thickness = layer_thickness(ds_lsg[\"depth\"])\n",
    "volume = area*thickness"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "id": "e21ec35d-f1b3-4058-acfd-03d5c7dc289f",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Surface area: 3.602E+14 m^2\n",
      "Volume: 1.427E+18 m^3\n"
     ]
    }
   ],
   "source": [
    "# Consistency check: ocean surface should be about 3.57e14 m^2, total volume should be about 1.335e18 m^3\n",
    "print(\"Surface area: {:.3E} m^2\".format(area.sum([\"lat\",\"lon\"]).isel(depth=0).item()))\n",
    "print(\"Volume: {:.3E} m^3\".format(volume.sum().item()))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": Full,
   "id": "f01f03f7-0e3d-48ef-bda7-6c1edf0d76be",
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "xrmask",
   "language": "python",
   "name": "xrmask"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.10.10"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
