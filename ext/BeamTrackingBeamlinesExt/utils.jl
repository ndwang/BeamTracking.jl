function check_Brho(Brho_ref, bunch::Bunch)
  if isnan(bunch.Brho_ref)
    if isnan(Brho_ref)
      @warn "Both the bunch and beamline do not have any set Brho_ref. If any LineElements have unnormalized fields stored as independent variables, tracking results will be NaNs"
    else
      @info "Setting bunch.Brho_ref = $Brho_ref (from the Beamline)"
      setfield!(bunch, :Brho_ref, typeof(bunch.Brho_ref)(Brho_ref)) #Brho_ref = Brho_ref
    end
  elseif !isnan(Brho_ref)  && !(Brho_ref ≈ bunch.Brho_ref)
    @warn "The reference energy of the bunch does NOT equal the reference energy of the Beamline. 
    Normalized field strengths in tracking ALWAYS use the reference energy of the bunch."
  end
end
