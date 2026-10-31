MODULE refine_random_mod
!
! Handles the random state 
!
! Needed to ensure that all callc to the macro for the derivatives use the 
! same set of random seeds at startup.
!
INTEGER                              :: refine_nseed
INTEGER, DIMENSION(0:64)             :: refine_seeds  ! Current seeds
INTEGER                              :: refine_nseed_post
INTEGER, DIMENSION(0:64)             :: refine_seeds_post ! Seeds after the macro ran at the centre point
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE refine_save_seeds
!-
!   Get and save the current seeds
!+
!
USE random_state_mod
IMPLICIT NONE
!
refine_nseed = random_nseeds()
!
CALL random_current(refine_nseed, refine_seeds)
!CALL random_current(refine_seeds)
CALL put_seeds(refine_nseed, refine_seeds)   ! Same state, but also drops a cached
                                             ! GASDEV value, so the centre run starts
                                             ! exactly like the restored derivative runs
!
END SUBROUTINE refine_save_seeds
!
!*******************************************************************************
!
SUBROUTINE refine_restore_seeds
!-
!   Get and save the current seeds
!+
!
USE random_state_mod
IMPLICIT NONE
!
CALL put_seeds(refine_nseed, refine_seeds)
!
END SUBROUTINE refine_restore_seeds
!
!*******************************************************************************
!
SUBROUTINE refine_save_seeds_post
!-
!   Save the random state reached after the macro ran at the centre point.
!   Restoring it after the derivative runs leaves the generator exactly where
!   a repeated run at the centre point would have left it.
!+
!
USE random_state_mod
IMPLICIT NONE
!
refine_nseed_post = random_nseeds()
!
CALL random_current(refine_nseed_post, refine_seeds_post)
!
END SUBROUTINE refine_save_seeds_post
!
!*******************************************************************************
!
SUBROUTINE refine_restore_seeds_post
!-
!   Restore the random state saved by refine_save_seeds_post
!+
!
USE random_state_mod
IMPLICIT NONE
!
CALL put_seeds(refine_nseed_post, refine_seeds_post)
!
END SUBROUTINE refine_restore_seeds_post
!
!*******************************************************************************
!
END MODULE refine_random_mod
