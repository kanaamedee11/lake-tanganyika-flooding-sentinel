# Data provenance and licences

All files in `data/` were exported by the author from Google Earth Engine with `lake_tanganyika_export.js` and `lake_tanganyika_export2.js`
(the area is a box around the northern end of Lake Tanganyika, Bujumbura and the Rusizi delta). They are derived from:

* Copernicus Sentinel-1 GRD and Sentinel-2 surface reflectance (contains modified Copernicus Sentinel data, 2017 to 2025). Free and open under the Copernicus data terms.
* JRC Global Surface Water (Pekel et al., 2016), used as a reference layer and for the control pixels.
* ESA WorldCover 2020 v100: Zanaga et al., CC BY 4.0, doi:10.5281/zenodo.5571936.
* SRTM elevation (NASA), used for a slope mask.

Check each dataset's Earth Engine catalogue page for current terms before reuse. The exports are small and may be regenerated with the scripts.
