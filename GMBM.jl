# GMBM Model Fernández, Manquehual-Cheuque, Somos-Valenzuela
using ArchGDAL, GeoData, GeoArrays, NCDatasets, CSV, DataFrames, Dates, StatsBase

# Get the time the simulation begins
init_time=now()

# Functions' descriptions
# - Get coordinates to the corresponding glacier from ERA5 elevation field
include("Functions/dem_pixel.jl")
# -  Calculate incoming radiation
include("Functions/in_radiation.jl")
# - Calculate albedo
include("Functions/albedo.jl")

# Defined dates
dates_nc=collect(Dates.Date(2020,1,1):Dates.Day(1):Dates.Date(2095,12,31))

# Extract initial amd end dates of nc corresponding to projection
init_period=Date(2020,1,1)
end_period=Date(2095,12,31)

# Deine years for calculations
years_sol=collect(2020:2095)

#Sequence of elevations from RGI data
seq_elev=collect(25:50:6975)

# Define glacier zones to be processed
zones=[5] # Can add more zones here (e.g., zones = collect(1:8))

# Define climate change experiments
experiment=("ssp245","ssp585","g6Solar")

# Define temperature and precipitation variables (can be adjusted based on variable names)
varT="tas"
varPp="pr"

# Define climate models to be used
models=("IPSL","CNRM-ESM2-1")


# Solar constant in MJ/m2 from 1365.7 W/m2
in_sc=0.0819

# Loop through each climate model
for mods in 1:length(models)
    mO=models[mods]

    # Read annual change in solar constant for the model, using Visioni et al calculations.
    global sc=CSV.read(string("Data_in/Solar_Constant/rsdt_",mO,".csv"),DataFrame)
    global sc_u=sc .* in_sc # Solar constant per year
    global sc_ts=collect(1:length(dates_nc)) .*0 .+in_sc  # Initialize solar constant per day

    # Loop through each experiment
    for ex in collect(1:length(experiment))
        global in_scTs=sc_u[:,ex]

        # Loop to get the solar constant per day
        for y in collect(1:length(years_sol))
            global sc_id=Date(years_sol[y],1,1)
            global sc_ed=Date(years_sol[y],12,31)

            # Find corresponding dates in dates_nc
            global ini_ys=findall(x->x==sc_id,dates_nc)
            global end_ys=findall(x->x==sc_ed,dates_nc)

            # Update solar constant per day in sc_ts
            global sc_ts[ini_ys[1]:end_ys[1]].=in_scTs[y]
        end

        # Loop through each glacier zone
        for gg in zones
            zone=string("G",gg)

            # Define file names for DEM and glacier locations
            dem_in=string("Data_in/Topo_data/",zone,"_DEM.tif")
            elev_in=string("Data_in/Coordinates/",zone,"_coords.csv")
            locs_gl=CSV.read(elev_in,DataFrame)

            # Get glacier elevations from ERA5
            glac_elev=dem_pixel(dem_in,elev_in)

            # Get glacier hypsometry
            in_hyp=CSV.read(string("Data_in/Hypsometry/",zone,"_hypso.csv"),DataFrame)

            # Get parameters for model
            pars=CSV.read(string("Data_in/Optimized_Parameters/",mO,"/",zone,"_optimal_params.csv"),DataFrame)

            # Read nc files (daily data)
            # - Temperature
            T_nc=Dataset(string("Data_in/Climate_data_input/",mO,"/",zone,"_T_",experiment[ex],".nc"))
            v=T_nc[varT]
            T=v[:,:]
            close(T_nc)

            # - Precipitation
            Pp_nc=Dataset(string("Data_in/Climate_data_input/",mO,"/",zone,"_Pp_",experiment[ex],".nc"))
            v=Pp_nc[varPp]
            Pp=v[:,:]./1000 # Pp in mm because gradients are in mm/m
            close(Pp_nc)

            # Remove negative or missing precipitation
            Ppbz=findall(x->ismissing(x),Pp)
            if length(Ppbz) > 0
                Pp[Ppbz].=0
            end
            Ppbz=findall(x->x<0,Pp)
            if length(Ppbz) > 0
                Pp[Ppbz].=0
            end

            # Load Tmin and Tmax to calculate atmospheric transmissivity for solar radiation
            # - Minimum temperature
            Tmin_nc=Dataset(string("Data_in/Climate_data_input/",mO,"/",zone,"_Tmin_",experiment[ex],".nc"))
            v=Tmin_nc[varT]
            Tmin=v[:,:]
            close(Tmin_nc)

            # - Maximum temperature
            Tmax_nc=Dataset(string("Data_in/Climate_data_input/",mO,"/",zone,"_Tmax_",experiment[ex],".nc"))
            #Tmax_nc=Dataset(string("Data_in/T_max/",mO,"/",zone,"Tmax_hist.nc"))
            v=Tmax_nc[varT]
            Tmax=v[:,:]
            close(Tmax_nc)

            # Calculation of temperature difference for atmospheric transmissivity
            To=Tmax-Tmin
            negs=findall(x->x<0,To)
            To[negs].=0
            To_u=sqrt.(To)

            # Remove Tmin and Tmax to free memory
            v=nothing
            Tmin=nothing
            Tmax=nothing

            # Extract initial and end dates of the nc file, corresponding to the climate projection
            ini_nc=findall(x->x==init_period,dates_nc)
            end_nc=findall(x->x==end_period,dates_nc)

            # Loop over each glacier
            for i in 1:length(pars.RGIId)

                # Use only calibrated glaciers
                if pars.selglac[i] == 0
                    continue
                end

                # Get the hypsometry of the modeled glacier
                gl=pars.RGIId[i]
                hyp_in=findall(x->x==gl,in_hyp.RGIId)
                hyp_int=Matrix(in_hyp[hyp_in,4:end])
                hyp_glac=vec(hyp_int)
                hyp_pos=findall(x->x>0,hyp_glac)

                # Skip glacier if hypsometry data is missing
                if length(hyp_pos) < 1
                    continue
                end

                # Select elevations and corresponding hypsometry weights
                seq_elev_u=seq_elev[hyp_pos]
                hyp_weigth=hyp_glac[hyp_pos]./1000

                # Find glacier location in elevation data
                locs_in=findall(x->x==gl,glac_elev.RGIId)
                if length(locs_in) < 1
                    continue
                end

                # Calculate elevation difference relative to ERA5
                glac_diff= seq_elev_u.-glac_elev.Elev[locs_in][1]

                # Extract daily climate data for the glacier location
                T_ts=T[glac_elev.LonIND[locs_in],glac_elev.LatIND[locs_in],ini_nc[1]:end_nc[1]]
                T_ts=vec(convert(Array,T_ts))
                P_ts=Pp[glac_elev.LonIND[locs_in],glac_elev.LatIND[locs_in],ini_nc[1]:end_nc[1]]
                P_ts=vec(convert(Array,P_ts))
                T_os=To_u[glac_elev.LonIND[locs_in],glac_elev.LatIND[locs_in],ini_nc[1]:end_nc[1]]
                T_os=vec(convert(Array,T_os))

                # Apply temperature and precipitation lapse rates
                T_diff=glac_diff.*(pars.LrT[i]/1000) #divide by 1000 to get in m.
                Pp_diff=glac_diff.*(pars.LrP[i]/1000) #divide by 1000 to get in m.

                # Initialize global variables to store mass balance components (to be filled in the loop)
                global mb=Any[]
                global ab=Any[]
                global ac=Any[]
                global IR_out=Any[]

                # Daily Loop for the model
                for var in 1:length(T_ts)

                     # Daily temperature considering lapse rate
                    T_gl=T_ts[var].+T_diff

                    # Ablation due to temperature using a weighted sum by hypsometry
                    abTs=T_gl.*(pars.DDF[i]/1000)
                    abz=findall(x->x<0,abTs)
                    abTs[abz].=0
                    abElev=abTs.*hyp_weigth
                    ab=push!(ab,round(sum(abElev);digits=5))

                    # Daily precipitation considering lapse rate
                    Pp_gl=(P_ts[var].+Pp_diff)./1000
                    Ppbz=findall(x->x<0,Pp_gl) # Find negative precipitation values (if any)
                    if length(Ppbz) > 0
                        Pp_gl[Ppbz].=0
                    end

                    # Exclude precipitation on days exceeding temperature threshold (likely melting)
                    Tt=findall(x->x>2,T_gl)
                    if length(Tt) > 0
                        Pp_gl[Tt].=0 # Set precipitation to zero for days above threshold
                    end

                    # Daily accumulation (consider zero precipitation)
                    if P_ts[var] == 0
                        global ac_elev=0 .*hyp_weigth # No accumulation if no precipitation
                    else
                        global ac_elev=Pp_gl.*hyp_weigth # Accumulation per band weighted by hypsometry
                    end
                    ac_t=round(sum(ac_elev);digits=5) # Accumulate and round daily accumulation (in meters)
                    ac=push!(ac,ac_t)

                    # Calculate Ablation from Incoming Radiation via function
                    # Solar constant: 1365.7 W/m2 = 81942 J/m2/min = 0.0819 MJ/m2
                    IR=in_radiation(sc_ts[var],glac_elev.Lat[i],dates_nc[var],seq_elev_u,T_os[var])

                    # Determine albedo
                    # - Initialize hypsometry for albedo calculation
                    if var == 1
                        # - Reset hyp_ac at start of each year
                        global hyp_ac=(seq_elev_u.*0)
                    else
                        # - After var == 1 determine which bin has accumulation to record that day as snow day
                        local    acbz=findall(x->x>0,ac_elev)
                        if length(acbz) > 0
                            global    hyp_ac[acbz].=var
                        end
                    end

                    # - Calculate albedo of the current time step via function
                    alb_o=albedo(var,hyp_ac,ac_elev)

                    #Calculate absorbed solar radiation
                    IRalb=IR.*((alb_o.-1).*-1)
                    SRF_u=pars.SRF[i]/1000

                    # Calculate ablation from Incoming Radiation
                    IRMelt=((SRF_u.*IRalb).*hyp_weigth)
                    IRMeltsum=round.(sum(IRMelt);digits=5)
                    IR_out=push!(IR_out,IRMeltsum)

                    # If ablation from temperature is zero, set IRMelt to zero
                    if sum(abElev) == 0
                        IRMelt= 0 .*hyp_weigth
                    end

                    # Calculate Mass Balance
                    mbTs=ac_elev.-(abElev.+IRMelt)
                    mbi=round(sum(mbTs).*0.85;digits=5) ./length(hyp_pos)
                    mb=push!(mb,mbi)
                end

                # End of glacier loop

                # Create DataFrame for Mass Balance per glacier (values in mm)
                global gl_ts=DataFrame(MB=mb.*1000.0,ACC=ac.*1000.0,AB=ab.*1000.0,IR=IR_out.*1000.0,
                T_in=round.(T_ts;digits=2),Pp_in=round.(P_ts;digits=5).*1000)

                # Save DataFrame to CSV file
                fnm=string("Data_out/",experiment[ex],"/",mO,"/",zone,"/",gl,".csv")
                CSV.write(fnm,gl_ts)
            end
        end
    end
end

# Print execution time
println("Time taken")
tot_time=now()-init_time
println(Dates.value(tot_time)/60000)
print("minutes")
