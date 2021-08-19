Explanations to 
cssl_metadata_all/cssl_metadata_all.csv and
cssl_metadata_manuscript_subset/cssl_metadata.csv


- sample_id: unique identifier for each soil samples; project specific
- sample_location: defines sampling location; is the same for all samples for different soil depths but from the same sampling location
- counry_name: country of origin
- country_code: Alpha-3 code, ISO 3166-1 for corresponding country of origin (more details: https://en.wikipedia.org/wiki/List_of_ISO_3166_country_codes)
- province_name: province of origin within the country of origin (source for each raster point: https://gadm.org/)
- subdivision_name: subdivision name (e.g. territory in the Democratic Republic of Congo; source: https://gadm.org/)
- subdivision_type: type of subdivision corresponding to 'subdivision_name' (source: https://gadm.org/)
- sampling_date: date of sampling, format is differing for different sample sets
- sampling_layer: depth of soil sampling on the 'sample_location' in cm
- gps_long: longitude (WGS84)
- gps_lat: latitude (WGS84)
- gps_true: indicating, if the true GPS coordinates (gps_long, gps_lat) habe been taken at sample_location. In case GPS data were not available (gps_true = no), the coordinates of the center of the indicated site were used. 
- altitude: elevation in meters above sea level
- land_use: land use on the 'sample_location'
- crop_type: where available, the crops cultivated on the 'sample_location' are indicated
- forest_type: where available, the forest type on the 'sample_location' is indicated
- savanna_type: where available, the savanna type on the 'sample_location' is indicated

