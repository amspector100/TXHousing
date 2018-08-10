import TXHousing

def generate_mastermaps(austin_save_path = None, dallas_save_path = None, houston_save_path = None):
    """ Generates the mastermaps for Austin, Dallas, Houston. Takes ~3 minutes."""

    TXHousing.analysis.mastermaps.create_austin_mastermap(save_path = austin_save_path)
    TXHousing.analysis.mastermaps.create_dallas_mastermap(save_path = dallas_save_path)
    TXHousing.analysis.mastermaps.create_houston_mastermap(save_path = houston_save_path)

    print("""Finished generating mastermaps""")

def generate_zoning_graphs():
    """ Generates the graphs based primarily on zoning data. Takes ~8 minutes."""

    TXHousing.analysis.zoning_graphs.plot_minimum_lot_size()
    TXHousing.analysis.zoning_graphs.plot_hd_locations()
    TXHousing.analysis.zoning_graphs.plot_broad_zones_proportion()
    TXHousing.analysis.zoning_graphs.plot_zone_income_histogram(calculate = False)

    # This one takes ~5 minutes
    TXHousing.analysis.zoning_graphs.map_broad_zones_dallas_austin()

    print("""Finished generating zoning_graphs""")

def generate_parcel_graphs():
    """ Generates graphs based on parcel data. """

    TXHousing.analysis.parcel_graphs.plot_singlefamily_lotsizes()
    TXHousing.analysis.parcel_graphs.plot_percent_undeveloped()
    TXHousing.analysis.parcel_graphs.calc_parking_costs()

    print("""Finished generating parcel_graphs""")

def generate_permit_graphs():

    TXHousing.analysis.permit_graphs.plot_permit_scatterplot()
    TXHousing.analysis.permit_graphs.plot_permit_locations()

    print("""Finished generating permit_graphs""")

def generate_suburb_graphs():

    TXHousing.analysis.suburbs.texas_job_centers()
    TXHousing.analysis.suburbs.suburbs_scatterplot()
    TXHousing.analysis.suburbs.analyze_transportation_networks()

    # loop through and plot suburbs choropleths
    for name, zoning_input in zip(['austin', 'dallas', 'houston'], [TXHousing.data_processing.zoning.austin_inputs,
                                                                    TXHousing.data_processing.zoning.dallas_inputs,
                                                                    TXHousing.data_processing.zoning.houston_inputs]):
        TXHousing.analysis.suburbs.suburbs_choropleth(name, zoning_input)

    print("""Finished generating suburb graphs""")

if __name__ == '__main__':

    generate_permit_graphs()
    generate_suburb_graphs()
    generate_parcel_graphs()
    generate_zoning_graphs()
    generate_mastermaps()