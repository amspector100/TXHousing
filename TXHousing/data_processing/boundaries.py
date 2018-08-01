"""Paths and methods for boundary files (uas, csa, counties, blocks, etc)"""
from ..utilities import measurements

csa_path = "data/cb_2017_us_csa_500k/cb_2017_us_csa_500k.shp"
cbsa_path =  "data/cb_2017_us_cbsa_500k/cb_2017_us_cbsa_500k.shp"
ua_path = "data/cb_2017_us_ua10_500k/cb_2017_us_ua10_500k.shp"

zip_boundaries_path = "data/cb_2017_us_zcta510_500k/cb_2017_us_zcta510_500k.shp"
county_boundaries_path = "data/cb_2017_us_county_500k/cb_2017_us_county_500k.shp"
texas_blocks_path = "data/ACS_2016_5YR_BG_48_TEXAS.gdb"
texas_places_path = measurements.texas_places_path


# Intersect zoning with zip codes with great precision
def zip_intersect(gdf, zips):
    """
    :param gdf: A geodataframe of some sort (it should probably be in the US, otherwise this is pointless).
    :param zips: The zip codes in the area of interest.
    :return: The gdf but with an extra column for each zip code which specifies the fractional area in the zip code.
    """

    # Make sure zips input and data index is full of strings
    zips = [str(something) for something in zips]
    gdf.index = [str(ind) for ind in gdf.index]
    gdf = gdf.to_crs({'init': 'EPSG:4326'})

    # Get the zip code geodata
    zipdata = helpers.get_zip_boundaries()
    zipdata = zipdata.loc[zipdata['ZCTA5CE10'].isin(zips), ['geometry']]

    # Create spatial index of zoning data
    spatial_index = gdf.sindex

    # Run through zip codes
    print('Starting to calculate which zones are in which zip codes')
    for code, row in tqdm(zipdata.iterrows()):

        # Find intersections
        polygon = row['geometry']
        possible_matches_index = list(spatial_index.intersection(polygon.bounds))
        possible_matches = gdf.iloc[possible_matches_index]
        precise_matches = possible_matches['geometry'].intersection(polygon)

        # Calculate percent of area in the zip code - this is a bit wasteful because most zones are only in one zip code,
        # but it avoids fragmenting the shapes which are quite hard to put back together.
        gdf[code] = precise_matches.area / gdf['geometry'].area
        gdf[code].fillna(value = 0, inplace = True)

    return gdf