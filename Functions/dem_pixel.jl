# Function to get the ERA pixel elevation for each glacier
function dem_pixel(dem,locs)
  # Initialize arrays for data processing and output
  elev=Any[]
  rgi=Any[]
  Lon=Any[]
  Lat=Any[]
  Lon_in=Any[]
  Lat_in=Any[]

  # Open DEM
  dem_use=GeoArrays.read(dem)

  # Get DEM dimensions
  dem_size=size(dem_use)

  # Read glacier location data from CSV file
  locs_in=CSV.read(locs,DataFrame)

  # Loop through each glacier location
  for i in 1:length(locs_in[:,1])
      # Get DEM indices for current glacier location
      ind=indices(dem_use,[locs_in.Lon[i],locs_in.Lat[i]])

       # Check if indices are within DEM bounds
      if ind[1] > dem_size[1]  || ind[2] > dem_size[2] || ind[1] < 1  || ind[2] < 1
          continue  # Skip if out of bounds
      else
          elev=push!(elev,dem_use[ind[1],ind[2]]) # Append elevation value to output array
      end

      # Append RGIId, longitude, latitude, and DEM indices to output arrays
      rgi=push!(rgi,locs_in.RGIId[i])
      Lon=push!(Lon,locs_in.Lon[i])
      Lat=push!(Lat,locs_in.Lat[i])
      Lon_in=push!(Lon_in,ind[1])
      Lat_in=push!(Lat_in,ind[2])
  end

  # Create DataFrame from output arrays
  dem_pixel=DataFrame(RGIId=rgi,Elev=elev,Lon=Lon,Lat=Lat,LonIND=Lon_in,LatIND=Lat_in)

  # Return DataFrame
  return(dem_pixel)
end
