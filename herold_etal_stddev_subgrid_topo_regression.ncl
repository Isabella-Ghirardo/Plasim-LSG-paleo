; Calculates standard deviation of ETOPO1 topography within each grid cell of a provided topography dataset, then calculates a linear regression to apply to paleotopography.
;
; You need: - NCL (http://www.ncl.ucar.edu/)
;	    - a paleotopography file on a regular grid (see the provided Eocene topography/bathymetry file for a template)
;           - a copy of the ETOPO1 dataset 'ETOPO1_Bed_c_gmt4.nc' renamed from http://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/cell_registered/netcdf/ETOPO1_Ice_c_gmt4.grd.gz
;
; Output:   - netCDF file containing standard deviation of sub grid scale paleotopography.
;
; nherold July 2013 using NCL 6.1.2

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
; ***************************************
; BEGIN USER INPUT
; ***************************************

        input_topography = "eocene_topo_markwick_ANTSCAPE_1x1.nc"
        output_file = "herold_etal_stddev_subgrid_etopo1_to_eocene_1x1_21042014.nc"		; desired output filename

; ***************************************
; END USER INPUT
; ***************************************

        etopo1file = "ETOPO1_Bed_c_gmt4.nc";"ETOPO1_Ice_c_gmt4.nc"             ; DON'T CHANGE. ETOPO1 cell registered (in order to have even number of longitudes for lonFlip), with ice.

; add files and variables
	add0 = addfile(etopo1file,"r")
	add1 = addfile(input_topography,"r")
        system("rm -f "+output_file)
	add2 = addfile(output_file,"c")

        etopo1tmp = add0->z
	etopo1 = int2flt(etopo1tmp)
	etopo1@_FillValue = default_fillvalue(typeof(etopo1))
        etopo1 = lonFlip(etopo1)
	paleotopo = add1->topo

	lat = add1->lat
	lon = add1->lon

; determine latitude and longitude increments
	latinc = (lat(3)-lat(2))/2
	loninc = (lon(3)-lon(2))/2
	print("latinc: "+latinc)
	print("loninc: "+loninc)

; interpolate etopo1 to atm grid resolution and assign attributes
        deg2rad = 0.0174532925
        cosweight = cos(etopo1&y*deg2rad)
	etopo1low = area_hi2lores_Wrap(etopo1&x,etopo1&y,etopo1,True,cosweight,lon,lat,False)
	etopo1low!0 = "lat"
	etopo1low!1 = "lon"
        etopo1low&lat@units = "degrees_north"
        etopo1low&lon@units = "degrees_east"

; Mask out oceans and land ice. Removal of land ice was discarded as recommended by reviewer 2.
        etopo1low = where(etopo1low.lt.0,etopo1low@_FillValue,etopo1low)

        lat2d = conform(etopo1low,lat,0)
        lon2d = conform(etopo1low,lon,1)
;        etopo1low = where(lat2d.lt.-61,etopo1low@_FillValue,etopo1low)          						; Antarctica
;        etopo1low = where(lat2d.lt.90 .and. lat2d.gt.57 .and. lon2d.gt.300 .and. lon2d.lt.350,etopo1low@_FillValue,etopo1low)	; Greenland

	delete([/lat2d,lon2d/])
        lat2d = conform(etopo1,etopo1&y,0)
        lon2d = conform(etopo1,etopo1&x,1)
;        etopo1 = where(lat2d.lt.-61,etopo1@_FillValue,etopo1)          						; Antarctica
;        etopo1 = where(lat2d.lt.90 .and. lat2d.gt.57 .and. lon2d.gt.300 .and. lon2d.lt.350,etopo1@_FillValue,etopo1)	; Greenland

; find the stddev of the high res etopo1 in each low res etopo1 grid cell
        stddev_etopo1 = etopo1low       ; get metadata
        stddev_etopo1 = stddev_etopo1@_FillValue

        do a = 0, dimsizes(lat)-1
                do b = 0, dimsizes(lon)-1
                        if(.not.ismissing(etopo1low(a,b))) then
                        ; create a box from the high res etopo1
                        box = etopo1({lat(a)-latinc:lat(a)+latinc},{lon(b)-loninc:lon(b)+loninc})
                        assignFillValue(box,etopo1)
                        box = where(box.lt.0,box@_FillValue,box)        ; don't include bathymetry in calculation of sub gridcell stddev
                        stddev_etopo1(a,b) = stddev(box)
                        delete(box)
                        end if
                end do
        end do

; Create bins of fixed height intervals (e.g. 100 m) and then find area weighted mean and stddev of each bin.
	numberbins = 61
	bins = fspan(0,6000,numberbins)
	binheights = new((/3,numberbins/),"float")	; mean, stddev and stderror of each bin, respectively
	binheights(:,:) = binheights@_FillValue
	deg2rad = (2*3.14159)/360.
	latwgt = cos(etopo1low&lat*deg2rad)
	do a = 0,dimsizes(bins)-2
		temp = where(etopo1low.ge.bins(a) .and. etopo1low.lt.bins(a+1),etopo1low,etopo1low@_FillValue)	; mask out values outside of bin elevation range
		tempstddev = where(etopo1low.ge.bins(a) .and. etopo1low.lt.bins(a+1),stddev_etopo1,etopo1low@_FillValue)
		binheights(0,a) = wgt_areaave_Wrap(temp,latwgt,1.,0)
		binheights(1,a) = wgt_areaave_Wrap(tempstddev,latwgt,1.,0);stddev(temp)

		if(num(temp).eq.1) then
			binheights(2,a) = stddev(tempstddev)/(num(temp))
		else
			if(num(temp).gt.1) then

		; calculate weighted stddev
				area_ave_stddev = wgt_areaave_Wrap(tempstddev,latwgt,1.,0)
				diffs = new((/num(temp)/),"float")
				weights = tempstddev	; get meta data
				weights = conform(tempstddev,cos(lat*deg2rad),0)
				weights = where(etopo1low.ge.bins(a) .and. etopo1low.lt.bins(a+1),weights,weights@_FillValue)
				onedweights = ndtooned(weights)
				onedstddev = ndtooned(tempstddev)

		                sum_wgts = sum(onedweights)
                		n = num(onedweights)
				y = 0
				do z = 0,dimsizes(onedstddev)-1
					if(.not.ismissing(onedstddev(z))) then
						diffs(y) = onedweights(z)*(onedstddev(z)-area_ave_stddev)^2
						y = y + 1
					end if
				end do
				stddev_wgt = sqrt(sum(diffs)/sum_wgts)
				binheights(2,a) = stddev_wgt/(sqrt(num(onedweights)))
				delete(diffs)
		; non-weighted stddev
		;		binheights(2,a) = stddev(tempstddev)/(sqrt(num(temp)))
			else
				binheights(2,a) = binheights@_FillValue
			end if
		end if
	end do

; remove missing and zero values from arrays, then write out to netCDF and ASCII
	empty = num(ismissing(binheights(2,:)) .or. binheights(2,:).eq.0)
	outvar = new((/3,dimsizes(binheights(2,:))-empty/),"float")
	outvar(0,:) = binheights(0,0:dimsizes(binheights(0,:))-empty-1)
        outvar(1,:) = binheights(1,0:dimsizes(binheights(1,:))-empty-1)
        outvar(2,:) = binheights(2,0:dimsizes(binheights(2,:))-empty-1)

	add2->stderror = outvar(2,:)
	add2->mean_of_stddev = outvar(1,:)
        add2->mean_of_heights = outvar(0,:)

;	asciiwrite("std_error",outvar(2,:))
;        asciiwrite("mean_of_stddev",outvar(1,:))
;        asciiwrite("mean_of_heights",outvar(0,:))

; calculate paleo std dev of sub grid scale topography and write out 
	paleorough = (/paleotopo/)
	copy_VarCoords(paleotopo,paleorough)
	paleorough = where(paleorough.lt.0,paleorough@_FillValue,paleorough)
        rc1 = regline(outvar(0,0:29),outvar(1,0:29))
        rc2 = regline(outvar(0,30:53),outvar(1,30:53))
	paleorough = where(paleotopo.ge.0 .and. paleotopo.lt.3000,paleotopo*rc1+rc1@yintercept,paleorough)
        paleorough = where(paleotopo.ge.3000,paleotopo*rc2+rc2@yintercept,paleorough)

	paleorough@long_name = "Standard deviation of sub grid cell scale paleotopgraphy"
	etopo1low@long_name = "ETOPO1 regridded to 1x1 degrees"

	add2->paleo_stddev_subgrid_topo = paleorough
	add2->paleotopo = paleotopo
	add2->etopo1 = etopo1low
	add2->etopo1_stddev_subgrid_topo = stddev_etopo1

	print("*************************************************")
	print(output_file+" created.")
end
