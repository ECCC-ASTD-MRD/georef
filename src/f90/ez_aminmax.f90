SUBROUTINE ez_AMINMAX(FLDMIN, FLDMAX, FLD, NI, NJ)
      REAL, intent(out) :: FLDMIN, FLDMAX
      INTEGER, intent(in) :: NI, NJ
      REAL, dimension(NI,NJ), intent(in) :: FLD
  
      FLDMIN = MINVAL(FLD)
      FLDMAX = MAXVAL(FLD)
  END
  