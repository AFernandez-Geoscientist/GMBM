# Function to calculate incoming solar radiation
# Calculation of incoming solar radiation at the respective elevation in W/m2*d
# From equations available at:
# Duffie, J. A. and Beckman, W. A.: Solar engineering of thermal processes, John Wiley & Sons, 2013
# Annandale, J., Jovanovic, N., Benade, N., and Allen, R.: Software for missing data error analysis of Penman-Monteith reference evapotranspiration, Irrigation Science, 21, 57–67, 2002.
# Thanks to Alvaro Gonzalez for this: https://cp.copernicus.org/preprints/cp-2019-37/cp-2019-37.pdf

function in_radiation(ct,lat,dates,hyp,tdif)
    # ct: is the solar constant (W/m2). Good to keep it as parameter in order to test the past
    # lat: latitude of the glacier centroid in degrees
    # dates: dates of the climate data (use dates_nc)
    # hyp: hypsometric distribution of the glacier
    # tdif: time series of daily temp difference (Tmax-Tmin)

    # Convert latitude to radians
    lrad=round(lat*(pi/180);digits=3)

    # Calculation of solar declination
    # - Day of the year
    year=Dates.year(dates)
    d_days=dates-Date(year,1,1)
    d=Dates.value(d_days)

    # - Solar declination
    sc=(23.4*pi/180)*sin(2*pi*((284+d)/365))

    # Determine sunrise hour angle (This seems only to work for latitudes 66N to 66S)
    w=round(acos(-tan(lrad)*tan(sc));digits=3)

    # Correction for eccentricity of Earth's orbit
    dr=round(1+0.033*cos((2*pi/365)*d);digits=3)

    # Extraterrestrial radiation
    R=(24*60/pi)*ct*dr*(w*sin(lrad)*sin(sc)+cos(sc)*cos(lrad)*sin(w))

    # Atmospheric transmissivity associated to tmin and tmax (0.16 is a constant)
    # - Tdif is the difference between Tmax and Tmin squared
    Tt=0.16.*(1 .+hyp.*2.7e-5).*(tdif)

    # Calculation of incoming solar radiation
    G=round.(R.*Tt;digits=2)

    #Output variable
    in_radiation=G
    return(in_radiation)
end
