import TXHousing

# Generate corrected dallas and houston permit data. It's probably easier just to pull the requesite files off of
# shared_data - this will take a while otherwise.

TXHousing.data_processing.boundaries.calc_percent_residential_in_block_groups()
TXHousing.data_processing.permit.get_corrected_dallas_permit_data()
TXHousing.data_process.permit.scrape_houston_permit_data(kind = 'structural')
TXHousing.data_processing.permit.scrape_houston_permit_data(kind = 'demolition')

# Cache parcel data - this will take several hours
TXHousing.data_processing.parcel.cache_municipal_parcel_data()
TXHousing.data_processing.parcel.cache_all_parcel_data()