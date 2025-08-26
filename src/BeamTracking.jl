module BeamTracking
using GTPSA,
      ReferenceFrameRotations,
      StaticArrays, 
      SIMD,
      VectorizationBase,
      EnumX,
      Unrolled,
      MacroTools,
      Adapt,
      Accessors,
      SpaceCharge

using KernelAbstractions
      
import GTPSA: sincu, sinhcu
import Base: setproperty!

# Put AtomicAndPhysicalConstants in a box for now for safety
include("Constants.jl")
using .Constants: Constants, Species, massof, chargeof, nameof, C_LIGHT, isnullspecies
export Species

export Bunch, State, ParticleView, sincu, sinhcu, sincuc, expq, quat_mul
export LinearTracking, Linear
export ExactTracking, Exact
export IntegrationTracking, SplitIntegration, DriftKick, BendKick, SolenoidKick, MatrixKick
export SpaceChargeIntegrationTracking, SpaceChargeIntegration
export track!

include("utils.jl")
include("types.jl")
include("kernel.jl")



include("modules/ExactTracking.jl") #; TRACKING_METHOD(::ExactTracking) = Exact
include("modules/LinearTracking.jl") #; TRACKING_METHOD(::LinearTracking) = Linear
include("modules/IntegrationTracking.jl") #; TRACKING_METHOD(::LinearTracking) = SplitIntegration, DriftKick, BendKick, SolenoidKick, MatrixKick
include("modules/SpaceChargeIntegrationTracking.jl") #; TRACKING_METHOD(::SpaceChargeIntegration) = SpaceCharge


# Empty tracking method to be imported+implemented by package extensions
function track! end

# --------------------------------------------------


# Modules separated:
#include("MatrixKick/MatrixKick.jl")
#include("Linear/Linear.jl")


end
