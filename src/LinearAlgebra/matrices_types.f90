MODULE matrices_types

  TYPE MAT_CSR_TYP ! A type to store matrices in CSR format
     LOGICAL                        :: start ! keeps track if it is the first solve
     INTEGER                        :: n ! Number of (locally owned) columns in the matrix
     INTEGER                        :: nnz ! Number of (locally owned) non-zeros
     INTEGER, DIMENSION(:), POINTER :: rowptr => NULL() ! Index of first element of each row in cols and vals
     INTEGER, DIMENSION(:), POINTER :: cols => NULL() ! Column of each element
     REAL(8), DIMENSION(:), POINTER :: vals => NULL() ! Value of each element
     INTEGER, DIMENSION(:), POINTER :: loc2glob => NULL() ! Global column number of the local columns
  END TYPE MAT_CSR_TYP

  TYPE RHS_TYP
     INTEGER                        :: n ! Number of (locally owned) columns in the matrix
     REAL(8), DIMENSION(:), POINTER :: vals => NULL() ! Value of each element
     INTEGER, DIMENSION(:), POINTER :: loc2glob => NULL() ! Global column number of the local columns
  END TYPE RHS_TYP

  TYPE MAT_IJV_TYP ! A type to store matrices in IJV format
     INTEGER                        :: n ! Number of (locally owned) columns in the matrix
     INTEGER                        :: nnz ! Number of (locally owned) non-zeros
     INTEGER, DIMENSION(:), POINTER :: rows => NULL() ! Row of each element
     INTEGER, DIMENSION(:), POINTER :: cols => NULL() ! Column of each element
     REAL(8), DIMENSION(:), POINTER :: vals => NULL() ! Value of each element
  END TYPE MAT_IJV_TYP

END MODULE matrices_types
