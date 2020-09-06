"Precip_and_boot_EOBS_ERA5_5gripdpoints" contains a list of:
-2 matrices, daily_precip_EOBS and daily_precip_ERA5
-4 arrays, boot_sigma_EOBS , boot_sigma_ERA5 , boot_xi_EOBS , boot_xi_ERA5 , which are the seasonal bootstrapped values of the GPD parameters, xi and sigma, obtained from ERA5 and EOBS precipitation EGPD fitting

read the file in R using the command "load(file = "Data/Precip_and_boot_EOBS_ERA5_5gripdpoints")"

Dimension of the matrices: size 14609*5 = number_of_days*number_of_stations
the row name indicates the date, the column name indicates the station name

Dimension of the arrays: size 4*5*200 = number_of_seasons*number_of_stations*bootstrap_size
the row name indicates the season, the column name indicates the station name

Note that there is no bootstrapped values of the parameters in JJA (summer) for Porto, because there
were less that 500 wet days (>1mm) in summer over the whole period at this location, which was a sine qua non condition to fit EGPD.

Location of the stations :
Station Bern Paris Porto Moscow Sofia
Longitude 7.44 2.35 -8.63 37.62 23.32
Latitude 46.95 48.86 41.16 55.76 42.70