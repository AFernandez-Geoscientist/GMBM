# Function to calculate albedo according to hypsometry
# Following Oerlemans, J. and Knap, W.(10.3189/S0022143000002574)
function albedo(time,lastsnowday,acc)
# time: is the iteration in the main program
# lastsnowday: array indicating the last day of snow in a particular elevation bin
# acc: array of accumulation per elevation bin

    # Albedo parameters
    alb_ice=0.3 # albedo of ice
    alb_fsnow=0.9 # albedo fresh snow
    alb_firn=0.6 # albedo of firn
    d=0.4 # scale coefficient for snowdepth
    tscel=3 # Time between transition of snow albedo to old snow albedo

    # Get albedo of snow
    weight=exp.((lastsnowday.-time)./tscel)
    albsnow=alb_firn.+(alb_fsnow-alb_firn).*weight

    # Get glacier albedo
    weight=exp.(-acc./d)
    alb=albsnow.+(alb_ice.-albsnow).*weight
    snow2much=findall(x->x>alb_fsnow,alb)
    alb[snow2much].=alb_fsnow
    ice2much=findall(x->x<alb_ice,alb)
    alb[ice2much].=alb_ice

    # Output albedo
    albedo=alb
    return(albedo)
end
