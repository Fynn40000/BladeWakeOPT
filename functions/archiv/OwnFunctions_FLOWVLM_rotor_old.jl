#=##############################################################################
# NOTE:
    This code has to be used in the file "src/FLOWVLM_rotor.jl" 
    in the FLOWVLM Repository under the Julia packages
    => it is necessary to overwrite the old within this file by the according functions within this file!
=###############################################################################

#-------------------------------------------------------------------------------------------------------------------
"""
  `Rotor(CW, r, chord, theta, LE_x, LE_z, B, airfoil)`

Object defining the geometry of a rotor/propeller/wind turbine. This class
behaves as an extension of the WingSystem class, hence all functions of
WingSystem can be applied to a Rotor object.

  # Arguments
  * CW::Bool                   : True for clockwise rotation, false for CCW.
  * r::Array{Float64,1}        : Radius position for the following variables.
  * chord::Array{Float64,1}    : Chord length.
  * theta::Array{Float64,1}    : Angle of attack (deg) from the rotor's plane
                                  of rotation.
  * LE_x::Array{Float64,1}     : x-position of leading edge (positive is ahead
                                  of radial axis relative to rotation).
  * LE_z::Array{Float64,1}     : z-position of leading edge (height from plane
                                  of rotation).
  * B::Int64                   : Number of blades.

  # Optional Arguments
  * airfoils::Array{Tuple{Float64, airfoilprep.Polar},1} : 2D airfoil properties
                                  along blade in the form [ (r_i, Polar_i) ]
                                  with Polar_i describes the airfoil at i-th
                                  radial position r_i (both the airfoil geometry
                                  in Polar_i and r_i must be normalized). At
                                  least root (r=0) and tip (r=1) must be given
                                  so all positions in between can be
                                  extrapolated. This properties are only used
                                  when calling CCBlade and for generating good
                                  loking visuals; ignore if only solving the VLM.

NOTE: r here is the radial position after precone is included in the geometry,
hence the need of explicitely declaring LE_z.

  # PROPERTIES
  * `sol` : Contains solution fields specific for Rotor types. They are formated
            as sol[field_name] = Dict(
                            "field_name" => output_field_name,
                            "field_type" => "scalar" or "vector",
                            "field_data" => data
                            )
            where `data` is an array data[i] = [val1, val2, ...] containing
            this field values (scalar or vector) of all control points in the
            i-th blade.

<!-- NOTE TO SELF: r is the y-direction on a wing, hence, remember to build the
               blade from root in the direction of positive y. -->
"""
mutable struct Rotor

  # Initialization variables (USER INPUT)
  CW::Bool                      # True for clockwise rotation
  r::FArrWrap                   # Radius position for the following variables
  chord::FArrWrap               # Chord length
  theta::FArrWrap               # Angle of attack (deg) from the rotor's axis
  LE_x::FArrWrap                # x-position of leading edge
  LE_z::FArrWrap                # z-position of leading edge (Height from plane
                                #                                  of rotation)
  B::IWrap                      # Number of blades
  # Optional inputs
  airfoils::Array{Tuple{FWrap,ap.Polar},1} # 2D airfoil properties along blade
  turbine_flag::Bool            # Whether this is a wind turbine or a propeller
  aoa_bound_min::FArrWrap       # Input Boundaries for angle of attack on each airfoil position 
  aoa_bound_max::FArrWrap                 # (define the boundaries on each pos. an airfoil is defined along the radius)

  # Properties
  RPM::Any                      # Current revs per minute
  hubR::FWrap                   # Hub radius
  rotorR::FWrap                 # Rotor radius
  m::IWrap                      # Number of control points (per blade)
  sol::Dict{String,Any}         # Solution fields for CCBlade (not FLOWVLM)

  # Data storage
  _wingsystem::WingSystem       # Rotor assembly
  _r::FArrWrap                  # Radius of each control point (on one blade)
  _chord::FArrWrap              # Chord length at each control point
  _theta::FArrWrap              # Angle of attack (deg) at each control point
  _LE_x::FArrWrap
  _LE_z::FArrWrap
  _polars::Array{ap.Polar, 1}   # Polar object at each control point (with x,y
                                #  containing the exact geometric airfoil)
  _polarroot::ap.Polar          # Polar at the root
  _polartip::ap.Polar           # Polar at the tip
  _aoa_bound_min::FArrWrap      # minimum angle of attack boundaries on each rediscretized airfoil
  _aoa_bound_max::FArrWrap      # maximum angle of attack boundaries on each rediscretized airfoil
  
  Rotor(
          CW, r, chord, theta, LE_x, LE_z, B,
          airfoils=Tuple{FWrap, ap.Polar}[],
          turbine_flag=false,
          aoa_bound_min=FWrap[],
          aoa_bound_max=FWrap[], 
          RPM=nothing,
            hubR=r[1], rotorR=r[end],
            m=0, sol=Dict(),
          _wingsystem=WingSystem(),
            _r=FWrap[], _chord=FWrap[], _theta=FWrap[],
            _LE_x=FWrap[], _LE_z=FWrap[],
            _polars=ap.Polar[],
              _polarroot=ap.dummy_polar(), _polartip=ap.dummy_polar(), 
              _aoa_bound_min=FWrap[], 
              _aoa_bound_max=FWrap[]
        ) = new(
          CW, r, chord, theta, LE_x, LE_z, B,
          airfoils,
          turbine_flag,
          aoa_bound_min,
          aoa_bound_max,
          RPM,
            hubR, rotorR,
            m, sol,
          _wingsystem,
            _r, _chord, _theta,
            _LE_x, _LE_z,
            _polars, _polarroot, _polartip,
            _aoa_bound_min, _aoa_bound_max
        )
end





#-------------------------------------------------------------------------------------------------------------------
"Calculates the distributed loads from CCBlade. It also stores normal (Np) and
tangential (Tp) components relative to the plane of rotation as given by
CCBlade, if `include_comps==true`.

If `return_performance==true`, it returns propulsive efficiency `eta`,
thrust coefficient `CT`, and torque coefficient `CQ` of each blade.

NOTE: These loads are per unit length of span"
function calc_distributedloads(self::Rotor, Vinf, RPM, rho::FWrap;
                                t::FWrap=0.0, include_comps=true,
                                return_performance=false, Vref=nothing,
                                Uinds=nothing,
                                sound_spd=nothing,
                                _lookuptable::Bool=false, _Vinds=nothing,
                                hubtiploss_correction=hubtiploss_nocorrection,
                                AR_to_360extrap = true,
                                debug=false)
  data = Array{FArrWrap}[]

  if debug
      if _lookuptable
          data_thetaeffdeg = FArrWrap[]
      end
  end
  if include_comps
    data_Np     = FArrWrap[]
    data_Tp     = FArrWrap[]
    data_cd     = FArrWrap[]
    if debug
        data_u      = FArrWrap[]
        data_v      = FArrWrap[]
        data_cl     = FArrWrap[]
        data_cn     = FArrWrap[]
        data_ct     = FArrWrap[]
        data_twistdeg = FArrWrap[]
        data_flowangledeg = FArrWrap[]
        data_Vx_out = FArrWrap[]
        data_Vy_out = FArrWrap[]
        data_F_out  = FArrWrap[]
        data_loc_solidity = FArrWrap[]
        if !_lookuptable
            data_a      = FArrWrap[]
            data_ap     = FArrWrap[]
            data_phi    = FArrWrap[]
            data_alpha  = FArrWrap[]
            data_W      = FArrWrap[]
            data_F      = FArrWrap[]
            data_G      = FArrWrap[]
        end
    end
    data_roR     = FArrWrap[]
  end
  if debug
      ccbrotors, ccbsectionss, ccbopss = [], [], []
  end

  gammas, mus_drag = _lookuptable ? ([], []) : (nothing, nothing)

  turbine_flag = self.turbine_flag

#   turbine_flag = true  # This is a flag for ccblade to swap signs

  # Calculates inflows
  # calc_inflow(self, Vinf, RPM; t=t, Vinds=(_lookuptable ? _Vinds : nothing) )
  calc_inflow(self, Vinf, RPM; t=t, Vinds=_Vinds)

  if return_performance; coeffs = []; end;

  # Distributed load on each blade
  for blade_i in 1:self.B

    # Formats the inflow as CCBlade's Inflow
    inflow = self.sol["CCBInflow"]["field_data"][blade_i]
    inflow_x = [V[1] for V in inflow]
    inflow_y = [V[2] for V in inflow]
    inflow_x = (-1)^(!turbine_flag)*inflow_x # The negative is needed to counteract the
    occbinflow = OCCBInflow(inflow_x, inflow_y, rho) # propeller swapping sign in CCBlade

    # Generates old-CCBlade Rotor object
    occbrotor = FLOWVLM2OCCBlade(self, RPM, blade_i, turbine_flag;
                                                            sound_spd=sound_spd, AR_to_360extrap=AR_to_360extrap)
    # Convert old-CCBlade rotor to current CCBlade rotor type
    ccbrotor, ccbsections, ccbops = OCCB2CCB(occbrotor, turbine_flag,
                                                        occbinflow; pitch=0.0)
    ccboutputs = nothing
    Np, Tp, uvec, vvec, gamma, thetaeffdeg = (nothing for i in 1:6)

    if _lookuptable
      (Np, Tp, uvec, vvec, gamma,
      cn, ct, cl, cd, thetaeffdeg, mu_drag,
      twistdeg, flowangledeg, Vx_out, Vy_out, 
      F_out, loc_solidity) = _calc_distributedloads_lookuptable(occbrotor,
                                                        occbinflow, turbine_flag;
                                                        hubtiploss_correction=hubtiploss_correction)
      push!(gammas, gamma)
      push!(mus_drag, mu_drag)
      zeroarr = zeros(size(Np))
      ccboutputs = ccb.Outputs(Np, Tp, zeroarr, zeroarr, uvec, vvec,
                                zeroarr, zeroarr, zeroarr, cl,
                                cd, cn, ct, zeroarr, zeroarr)
    else
      # Calls CCBlade
      # NOTE TO SELF: Forces normal and tangential to the plane of rotation
      # Np, Tp, uvec, vvec = ccb.distributedloads(ccbrotor, ccbinflow, turbine_flag)
      ccboutputs = ccb.solve.(Ref(ccbrotor), ccbsections, ccbops)
      Np, Tp, a, ap, uvec, vvec, phi, alpha, W, cl, cd, cn, ct, loss, effloss = ccboutputs.Np, ccboutputs.Tp, ccboutputs.a, ccboutputs.ap, ccboutputs.u, ccboutputs.v, ccboutputs.phi, ccboutputs.alpha, ccboutputs.W, ccboutputs.cl, ccboutputs.cd, ccboutputs.cn, ccboutputs.ct, ccboutputs.F, ccboutputs.G
    end

    # Convert forces from CCBlade's c.s. to global c.s.
    ccb_Fs = [ [Np[i], Tp[i], 0.0] for i in 1:get_mBlade(self) ]
    Fs = [ _ccblade2global( get_blade(self, blade_i), ccb_F, self.CW;
                          translate=false) for ccb_F in ccb_Fs ]                     
    #println("ccb_Fs")
    #println(ccb_Fs)
    #println("Fs")
    #println(Fs)


    # Stores the field
    push!(data, Fs)
    if debug
        if _lookuptable
            push!(data_thetaeffdeg, thetaeffdeg)
        end
    end
    if include_comps
      push!(data_Np     , Np)
      push!(data_Tp     , Tp)
      push!(data_cd     , cd)
      if debug
          push!(data_u      , uvec)
          push!(data_v      , vvec)
          push!(data_cl     , cl)
          push!(data_cn     , cn)
          push!(data_ct     , ct)
          push!(data_twistdeg, twistdeg)
          push!(data_flowangledeg, flowangledeg)
          push!(data_Vx_out, Vx_out)
          push!(data_Vy_out, Vy_out)
          push!(data_F_out, F_out)
          push!(data_loc_solidity, loc_solidity)
            if !_lookuptable
                push!(data_a      , a)
                push!(data_ap     , ap)
                push!(data_phi    , phi)
                push!(data_alpha  , alpha)
                push!(data_W      , W)
                push!(data_F      , loss)
                push!(data_G      , effloss)
            end
          push!(ccbrotors, ccbrotor)
          push!(ccbsectionss, ccbsections)
          push!(ccbopss, ccbops)
      end
      push!(data_roR     , self._r/self.rotorR)
    end

    if return_performance
      if Vref==nothing
        error("Performance coefficients requested, but no reference freestream"*
              " has been provided.")
      end
      # T, Q = ccb.thrusttorque(ccbrotor, [ccbinflow], turbine_flag)
      # eta, CT, CQ = ccb.nondim(T, Q, Vref, 2*pi*RPM/60, rho,
      #                           ccbrotor.Rtip, ccbrotor.precone, turbine_flag)
      #NOTE: Define ccboutputs
      T, Q = ccb.thrusttorque(ccbrotor, ccbsections, ccboutputs)
      eta, CT, CQ = ccb.nondim(T, Q, Vref, 2*pi*RPM/60, rho, ccbrotor)
      push!(coeffs, [eta,CT,CQ])
    end

    if Uinds!=nothing; push!(Uinds, [uvec, vvec]); end;
  end

  # Adds solution fields
  field = Dict(
              "field_name" => "DistributedLoad",
              "field_type" => "vector",
              "field_data" => data
              )
  self.sol[field["field_name"]] = field
  if debug
      if _lookuptable
          field = Dict(
                      "field_name" => "ThetaEffDeg",
                      "field_type" => "scalar",
                      "field_data" => data_thetaeffdeg
                      )
          self.sol[field["field_name"]] = field

      end
  end
  if include_comps
    field = Dict(
                "field_name" => "Np",
                "field_type" => "scalar",
                "field_data" => data_Np
                )
    self.sol[field["field_name"]] = field
    field = Dict(
                "field_name" => "Tp",
                "field_type" => "scalar",
                "field_data" => data_Tp
                )
    self.sol[field["field_name"]] = field
    field = Dict(
            "field_name" => "cd",
            "field_type" => "scalar",
            "field_data" => data_cd
            )
    self.sol[field["field_name"]] = field
    if debug
        field = Dict(
            "field_name" => "cl",
            "field_type" => "scalar",
            "field_data" => data_cl
            )
        self.sol[field["field_name"]] = field
        field = Dict(
                "field_name" => "cn",
                "field_type" => "scalar",
                "field_data" => data_cn
                )
        self.sol[field["field_name"]] = field
        field = Dict(
                "field_name" => "ct",
                "field_type" => "scalar",
                "field_data" => data_ct
                )
        self.sol[field["field_name"]] = field
        field = Dict(
                "field_name" => "twistdeg",
                "field_type" => "scalar",
                "field_data" => data_twistdeg
                )
        self.sol[field["field_name"]] = field
        field = Dict(
                "field_name" => "flowangledeg",
                "field_type" => "scalar",
                "field_data" => data_flowangledeg
                )
        self.sol[field["field_name"]] = field
        field = Dict(
                "field_name" => "Vx_out",
                "field_type" => "scalar",
                "field_data" => data_Vx_out
                )
        self.sol[field["field_name"]] = field
        field = Dict(
                "field_name" => "Vy_out",
                "field_type" => "scalar",
                "field_data" => data_Vy_out
                )
        self.sol[field["field_name"]] = field
        field = Dict(
                "field_name" => "F_out",
                "field_type" => "scalar",
                "field_data" => data_F_out
                )
        self.sol[field["field_name"]] = field
        field = Dict(
                "field_name" => "loc_solidity",
                "field_type" => "scalar",
                "field_data" => data_loc_solidity
                )
        self.sol[field["field_name"]] = field
        if !_lookuptable
            field = Dict(
                    "field_name" => "u",
                    "field_type" => "scalar",
                    "field_data" => data_u
                    )
            self.sol[field["field_name"]] = field
            field = Dict(
                        "field_name" => "v",
                        "field_type" => "scalar",
                        "field_data" => data_v
                        )
            self.sol[field["field_name"]] = field
            field = Dict(
                        "field_name" => "a",
                        "field_type" => "scalar",
                        "field_data" => data_a
                        )
            self.sol[field["field_name"]] = field
            field = Dict(
                        "field_name" => "ap",
                        "field_type" => "scalar",
                        "field_data" => data_ap
                        )
            self.sol[field["field_name"]] = field
            field = Dict(
                        "field_name" => "phi",
                        "field_type" => "scalar",
                        "field_data" => data_phi
                        )
            self.sol[field["field_name"]] = field
            field = Dict(
                        "field_name" => "alpha",
                        "field_type" => "scalar",
                        "field_data" => data_alpha
                        )
            self.sol[field["field_name"]] = field
            field = Dict(
                        "field_name" => "W",
                        "field_type" => "scalar",
                        "field_data" => data_W
                        )
            self.sol[field["field_name"]] = field
            field = Dict(
                        "field_name" => "F",
                        "field_type" => "scalar",
                        "field_data" => data_F
                        )
            self.sol[field["field_name"]] = field
            field = Dict(
                        "field_name" => "G",
                        "field_type" => "scalar",
                        "field_data" => data_G
                        )
            self.sol[field["field_name"]] = field
        end
        field = Dict(
                    "field_name" => "ccbrotor",
                    "field_type" => "not-vtk",
                    "field_data" => ccbrotors
                    )
        self.sol[field["field_name"]] = field
        field = Dict(
                    "field_name" => "ccbsections",
                    "field_type" => "not-vtk",
                    "field_data" => ccbsectionss
                    )
        self.sol[field["field_name"]] = field
        field = Dict(
                    "field_name" => "ccbops",
                    "field_type" => "not-vtk",
                    "field_data" => ccbopss
                    )
        self.sol[field["field_name"]] = field
    end
    field = Dict(
            "field_name" => "roR",
            "field_type" => "scalar",
            "field_data" => data_roR
            )
    self.sol[field["field_name"]] = field
  end

  if return_performance
    return coeffs, gammas, mus_drag
  else
    return nothing, gammas, mus_drag
  end
end








#-------------------------------------------------------------------------------------------------------------------
"Calculates the airfoils at each control point"
function _calc_airfoils(self::Rotor, n::IWrap, r::FWrap,
  central, refinement; rediscretize::Bool=true,
  rfl_n_lower::IWrap=15, rfl_n_upper::IWrap=15,
  rfl_r::FWrap=14.0, rfl_central::Bool=true)
#println("_calc_airfoils is executed to calculate the polars and contours at the new control points")
# Erases previous polars
self._polars = ap.Polar[]

# Normalized position of each interval boundary
aux_f(x) = x
if size(refinement)[1]!=0
# Precalulations of complex refinement
nsecs = size(refinement)[1]
ntot = sum([refinement[i][2] for i in 1:nsecs])
ctot = sum([refinement[i][1] for i in 1:nsecs])
ns = [] # Number of lattices in each section
for i in 1:nsecs
if i==nsecs
push!(ns, n-sum(ns))
else
push!(ns, floor(n*refinement[i][2]/ntot))
end
end
points = gt.multidiscretize(aux_f, 0, 1,
[(sec[1], ns[i], sec[3], false) for (i,sec) in enumerate(refinement)])
else
points = gt.discretize(aux_f, 0, 1, n, r; central=central)
end

# Normalized position of control points
norm_CPs = [ (points[i]+points[i-1])/2 for i in 2:size(points)[1] ]

# ------ Iterates over CPs calculating the polar of the airfoil at that CP ---

prev_contour, next_contour = nothing, nothing # Airfoil contours
prev_i, next_i = nothing, 1 # Indices of airfoils previous and next to cur CP

# WARNING: The blending function only blends airfoil properties and contours,
#       but it ignores all other parameters (Re, Ma, etc)

# Iterates over each control point
for this_CP in norm_CPs

# Finds the next airfoil closest to it
while this_CP > self.airfoils[next_i][1]
prev_i, prev_contour = next_i, next_contour
next_i += 1
next_contour = nothing
end

# Rediscretizes any contours as needed
if prev_contour==nothing
if rediscretize
prev_contour = _rediscretize_airfoil(self.airfoils[prev_i][2].x,
                        self.airfoils[prev_i][2].y,
                        rfl_n_lower, rfl_n_upper, rfl_r,
                        rfl_central)
else
prev_contour = (self.airfoils[prev_i][2].x, self.airfoils[prev_i][2].y)
end
end
if next_contour==nothing
if rediscretize
next_contour = _rediscretize_airfoil(self.airfoils[next_i][2].x,
                        self.airfoils[next_i][2].y,
                        rfl_n_lower, rfl_n_upper, rfl_r,
                        rfl_central)
else
next_contour = (self.airfoils[next_i][2].x, self.airfoils[next_i][2].y)
end
end

# Weight of the next airfoil vs the previous
weight = (this_CP-self.airfoils[prev_i][1]
) / (self.airfoils[next_i][1]-self.airfoils[prev_i][1])

# Blends Cd and Cl curves
prev_pypolar = self.airfoils[prev_i][2].pyPolar
next_pypolar = self.airfoils[next_i][2].pyPolar
blended_pypolar = prev_pypolar.blend(next_pypolar, weight)

# Blends airfoil geometry
blended_x = weight*next_contour[1] + (1-weight)*prev_contour[1]
blended_y = weight*next_contour[2] + (1-weight)*prev_contour[2]

# Calculate the maximum and minimum angle of attack boundaries if aoa_bound_min and aoa_bound_max are given
if length(self.aoa_bound_min) != 0
prev_aoa_min = self.aoa_bound_min[prev_i]
next_aoa_min = self.aoa_bound_min[next_i]
aoa_min = prev_aoa_min + weight * (next_aoa_min - prev_aoa_min) 
push!(self._aoa_bound_min, aoa_min)
end
if length(self.aoa_bound_max) != 0
prev_aoa_max = self.aoa_bound_max[prev_i]
next_aoa_max = self.aoa_bound_max[next_i]
aoa_max = prev_aoa_max + weight * (next_aoa_max - prev_aoa_max) 
push!(self._aoa_bound_max, aoa_max)
end

# Blended Polar object
blended_polar = ap._pyPolar2Polar(blended_pypolar; x=blended_x, y=blended_y)

push!(self._polars, blended_polar) # Stores CP airfoils in self._polars
end

# Stores root and tip airfoils
if rediscretize
root_x, root_y = _rediscretize_airfoil(self.airfoils[1][2].x,
                      self.airfoils[1][2].y,
                      rfl_n_lower, rfl_n_upper, rfl_r,
                      rfl_central)
tip_x, tip_y = _rediscretize_airfoil(self.airfoils[end][2].x,
                      self.airfoils[end][2].y,
                      rfl_n_lower, rfl_n_upper, rfl_r,
                      rfl_central)
self._polarroot = ap._pyPolar2Polar(self.airfoils[1][2].pyPolar;
                      x=root_x, y=root_y)
self._polartip = ap._pyPolar2Polar(self.airfoils[end][2].pyPolar;
                      x=tip_x, y=tip_y)
else
self._polarroot = self.airfoils[1][2]
self._polartip = self.airfoils[end][2]
end

end









#-------------------------------------------------------------------------------------------------------------------
"""Calculates the load distribution by using the airfoil lookup table on the
given inflow (this assumes that the inflow already includes all induced velocity
and it is the effective inflow).
"""
function _calc_distributedloads_lookuptable(ccbrotor::OCCBRotor,
  ccbinflow::OCCBInflow,
  turbine_flag::Bool;
  hubtiploss_correction=hubtiploss_nocorrection)

# check if propeller
swapsign = turbine_flag ? 1 : -1

# initialize arrays
n = length(ccbrotor.r)
Np = fill(0.0, n)
Tp = fill(0.0, n)
cn = fill(0.0, n)
ct = fill(0.0, n)
cl = fill(0.0, n)
cd = fill(0.0, n)
uvec = fill(0.0, n)
vvec = fill(0.0, n)
gamma = fill(0.0, n)
mu_drag = fill(0.0, n)
thetaeffdeg = fill(0.0, n)
twistdeg = fill(0.0, n)
flowangledeg = fill(0.0, n)
Vx_out = swapsign*ccbinflow.Vx
Vy_out = ccbinflow.Vy
F_out = fill(0.0, n)
loc_solidity = fill(0.0, n)

# Radial length of every horseshoe
lengths = Float64[2*(ccbrotor.r[1]-ccbrotor.Rhub)]
for i in 2:n
push!(lengths, 2*(ccbrotor.r[i]-ccbrotor.r[i-1])-lengths[end] )
end


for i in 1:n
twist = swapsign*ccbrotor.theta[i]
Vx = swapsign*ccbinflow.Vx[i]
Vy = ccbinflow.Vy[i]

aux1 = 0.5*ccbinflow.rho*(Vx*Vx+Vy*Vy)*ccbrotor.chord[i]

# Effective angle of attack (rad)
thetaV = atan(Vx, Vy)
thetaeff = thetaV - twist

thetaeffdeg[i] = thetaeff*180/pi
twistdeg[i] = twist*180/pi
flowangledeg[i] = thetaV*180/pi
#println("angles = $([thetaeff, thetaV, twist]*180/pi)\tVx,Vy=$([Vx, Vy])")

# airfoil cl/cd
#cl[i], cd[i] = occb_airfoil(ccbrotor.af[i], thetaeff) # THIS IS THE ORIGINAL CODELINE
cl[i], cd[i] = occb_airfoil(ccbrotor.af[i], -1*swapsign*thetaeff)
cl[i] = -1*swapsign*cl[i] 


# Tip and hub correction factor
B, Rtip, Rhub, r = ccbrotor.B, ccbrotor.Rtip, ccbrotor.Rhub, ccbrotor.r[i]
(eh1, eh2, eh3, maxah), (et1, et2, et3, maxat) = hubtiploss_correction

factorhub = B/2.0*(   (r/Rhub)^eh1 - 1   )^eh2/abs(sin(max(maxah*pi/180, abs(thetaV))))^eh3
Fhub = 2.0/pi*acos(exp(-factorhub))

factortip = B/2.0*(  (Rtip/r)^et1 - 1  )^et2/abs(sin(max(maxat*pi/180, abs(thetaV))))^et3
Ftip = 2.0/pi*acos(exp(-factortip))

F = Ftip * Fhub

cl[i] *= F
F_out[i] = F
loc_solidity[i] = (ccbrotor.B*ccbrotor.chord[i])/ (2*pi*ccbrotor.r[i]) # calculate the local solidity
# NOTE: Here we leave cd uncorrected assuming that it is all friction and form drag
# cd[i] *= F

# normal and tangential coefficients
sthtV = sin(thetaV)
cthtV = cos(thetaV)
cn[i] = cl[i]*cthtV + cd[i]*sthtV # ?????
ct[i] = swapsign*(cl[i]*sthtV - cd[i]*cthtV)


# Normal and tangential forces per unit length
Np[i] = cn[i]*aux1
Tp[i] = ct[i]*aux1

# If a turbine is calculated then return the cn and ct coefficients normed by the freestream velocity in x direction
if turbine_flag
magVinfx = 11.4                                             # ENTER THE FREESTREAM VELOCITY HERE
aux2 = 0.5*ccbinflow.rho*(magVinfx*magVinfx)*pi*(ccbrotor.Rtip * ccbrotor.Rtip)  # Wind turbine norm factor

cn[i] = (Np[i]*lengths[i]) / aux2
ct[i] = (Tp[i]*lengths[i]) / aux2
end


# Circulation
gamma[i] = cl[i]*sqrt(Vx*Vx+Vy*Vy)*ccbrotor.chord[i]/2

# Dragging line dipole strength (see Caprace's thesis, 2020)
# mu_drag[i] = cd[i]*sqrt(Vx*Vx+Vy*Vy)*ccbrotor.chord[i]/2
# NOTE: Here I'm correcting DG's model dividing by chord length and later
#       multiplying by the traveled distance to account for the time step
mu_drag[i] = cd[i]*sqrt(Vx*Vx+Vy*Vy)/2
end

# reverse sign of outputs for propellers
# Tp *= swapsign # I already reversed this
vvec *= swapsign

return Np, Tp, uvec, vvec, gamma, cn, ct, cl, cd, thetaeffdeg, mu_drag, twistdeg, flowangledeg, Vx_out, Vy_out, F_out, loc_solidity
end












#-------------------------------------------------------------------------------------------------------------------
function Base.deepcopy_internal(x::Rotor, stackdict::IdDict)
  if haskey(stackdict, x)
      return stackdict[x]
  end

  y = Rotor(  Base.deepcopy_internal(x.CW, stackdict),
              Base.deepcopy_internal(x.r, stackdict),
              Base.deepcopy_internal(x.chord, stackdict),
              Base.deepcopy_internal(x.theta, stackdict),
              Base.deepcopy_internal(x.LE_x, stackdict),
              Base.deepcopy_internal(x.LE_z, stackdict),
              Base.deepcopy_internal(x.B, stackdict),
              Base.deepcopy_internal(x.airfoils, stackdict),
              Base.deepcopy_internal(x.turbine_flag, stackdict),
              Base.deepcopy_internal(x.aoa_bound_min, stackdict),
              Base.deepcopy_internal(x.aoa_bound_max, stackdict),
              Base.deepcopy_internal(x.RPM, stackdict),
              Base.deepcopy_internal(x.hubR, stackdict),
              Base.deepcopy_internal(x.rotorR, stackdict),
              Base.deepcopy_internal(x.m, stackdict),
              Base.deepcopy_internal(x.sol, stackdict),
              Base.deepcopy_internal(x._wingsystem, stackdict),
              Base.deepcopy_internal(x._r, stackdict),
              Base.deepcopy_internal(x._chord, stackdict),
              Base.deepcopy_internal(x._theta, stackdict),
              Base.deepcopy_internal(x._LE_x, stackdict),
              Base.deepcopy_internal(x._LE_z, stackdict),
              Base.deepcopy_internal(x._polars, stackdict),
              Base.deepcopy_internal(x._polarroot, stackdict),
              Base.deepcopy_internal(x._polartip, stackdict),
              Base.deepcopy_internal(x._aoa_bound_min, stackdict),
              Base.deepcopy_internal(x._aoa_bound_max, stackdict))

  stackdict[x] = y
  return y
end