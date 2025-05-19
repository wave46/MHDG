MODULE interpolation
  USE printutils
CONTAINS

  FUNCTION binarysearch(length, array, VALUE, delta)
    ! Given an array and a value, returns the index of the element that
    ! is closest to, but less than, the given value.
    ! Uses a binary search algorithm.
    ! "delta" is the tolerance used to determine if two values are equal
    ! if ( abs(x1 - x2) <= delta) then
    !    assume x1 = x2
    ! endif

    IMPLICIT NONE
    INTEGER, INTENT(IN)                 :: length
    REAL, DIMENSION(length), INTENT(IN) :: array
    !f2py depend(length) array
    REAL, INTENT(IN)                    :: VALUE
    REAL, INTENT(IN), OPTIONAL          :: delta

    INTEGER                             :: binarysearch

    INTEGER                             :: left, middle, right
    REAL                                :: d

    IF (PRESENT(delta) .EQV. .TRUE.) THEN
       d = delta
    ELSE
       d = 1e-9
    ENDIF

    left = 1
    right = length
    DO
       IF (left > right) THEN
          EXIT
       ENDIF
       middle = NINT((left + right)/2.0)
       IF (ABS(array(middle) - VALUE) <= d) THEN
          binarySearch = middle
          RETURN
       ELSE IF (array(middle) > VALUE) THEN
          right = middle - 1
       ELSE
          left = middle + 1
       END IF
    END DO
    binarysearch = right

  END FUNCTION binarysearch

  REAL FUNCTION interpolate(x_len, x_array, y_len, y_array, f, x, y, delta)
    ! This function uses bilinear interpolation to estimate the value
    ! of a function f at point (x,y)
    ! f is assumed to be sampled on a regular grid, with the grid x values specified
    ! by x_array and the grid y values specified by y_array
    ! Reference: http://en.wikipedia.org/wiki/Bilinear_interpolation
    IMPLICIT NONE
    INTEGER, INTENT(IN)                       :: x_len, y_len
    REAL, DIMENSION(x_len), INTENT(IN)        :: x_array
    REAL, DIMENSION(y_len), INTENT(IN)        :: y_array
    REAL, DIMENSION(x_len, y_len), INTENT(IN) :: f
    REAL, INTENT(IN)                          :: x, y
    REAL, INTENT(IN), OPTIONAL                :: delta
    !f2py depend(x_len) x_array, f
    !f2py depend(y_len) y_array, f

    REAL                                      :: denom, x1, x2, y1, y2
    INTEGER                                   :: i, j

    i = binarysearch(x_len, x_array, x, delta)
    j = binarysearch(y_len, y_array, y, delta)

    IF (i == x_len) THEN
       WRITE (6, *) "Problem IN the binary search"
       WRITE (6, *) "x", x
       WRITE (6, *) "max(x_array)", MAXVAL(x_array)
       WRITE (6, *) "min(x_array)", MINVAL(x_array)
    END IF

    IF (j == y_len) THEN
       WRITE (6, *) "Problem IN the binary search"
       WRITE (6, *) "y", y
       WRITE (6, *) "max(y_array)", MAXVAL(y_array)
       WRITE (6, *) "min(y_array)", MINVAL(y_array)
    END IF

    x1 = x_array(i)
    x2 = x_array(i + 1)

    y1 = y_array(j)
    y2 = y_array(j + 1)

    denom = (x2 - x1)*(y2 - y1)

    interpolate = (f(i, j)*(x2 - x)*(y2 - y) + f(i + 1, j)*(x - x1)*(y2 - y) + &
         f(i, j + 1)*(x2 - x)*(y - y1) + f(i + 1, j + 1)*(x - x1)*(y - y1))/denom

  END FUNCTION interpolate

  FUNCTION nodesearch(x, y, xy_len, x_array, y_array)
    ! Given a  point (x,y), returns the index of the closest node IN 2D
    IMPLICIT NONE
    INTEGER,INTENT(IN)               :: xy_len
    REAL,INTENT(IN)                  :: x, y
    REAL,DIMENSION(xy_len),INTENT(IN):: x_array, y_array
    REAL*8                           :: d(xy_len)

    INTEGER :: nodesearch

    d = SQRT( (x_array - x)**2 + (y_array - y)**2 )
    nodesearch = MINLOC(d, DIM = 1)

  END FUNCTION nodesearch
#ifndef PARALL
  SUBROUTINE lineintegration(qp_len, x_vec, y_vec, f, Xc, T, nodes2D, nli)
#else
   SUBROUTINE lineintegration(qp_len, x_vec, y_vec, f, Xc, T, nodes2D, ghost_elems, nli)
#endif
   ! Given a set of points (x, y) along a line of sight, returns the line integration of
   ! f. The function f is evaluated IN the closest nodes to points (x, y).
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: qp_len, nodes2D, T(:,:)
   REAL*8, INTENT(IN)  :: x_vec(:), y_vec(:), f(:), Xc(:,:)
   REAL*8, INTENT(OUT) :: nli
   INTEGER             :: i, inp, n_ind, f_ind(2), np_len
   REAL*8              :: x, y, x_prev, y_prev, f_prev, dl, f_cur, d(nodes2D)
   LOGICAL             :: is_inside
   INTEGER             :: element_index
#ifdef PARALL
   INTEGER, INTENT(IN) :: ghost_elems(:)
#endif

   ! Initialize result
   nli = 0.0D0
   x_prev = 0.0D0
   y_prev = 0.0D0
   f_prev = 0.0D0
   dl = 0.0D0
   f_cur = 0.0D0
   d = 0.0D0

   DO i = 1, qp_len
      x = x_vec(i)
      y = y_vec(i)
      np_len = SIZE(Xc, 1)

      ! Check if the point (x, y) is inside the mesh
      CALL point_in_mesh(x, y, Xc, T, is_inside, element_index)

      IF (.NOT. is_inside) THEN
         x_prev = x
         y_prev = y
         f_prev = 0.
         CYCLE  ! Skip this point if it is outside the mesh
      ENDIF
#ifdef PARALL
      ! Check if the element is a ghost element
      IF ((ghost_elems(element_index) .EQ. 1)) THEN
         x_prev = x
         y_prev = y
         f_prev = 0.
         CYCLE
      END IF
#endif

      ! Compute distances to the point (x, y)
      d = SQRT((Xc(T(element_index, :), 1) - x)**2 + (Xc(T(element_index, :), 2) - y)**2)

      ! Find the index of the closest node in the element
      inp = MINLOC(d, DIM = 1)
      f_cur = f((element_index-1)*nodes2D + inp)

      ! Compute the line segment length and update the summation
      IF (i > 1) THEN
         dl = SQRT((x - x_prev)**2 + (y - y_prev)**2) * 1.901e-3
         nli = nli + 0.5D0 * (f_cur + f_prev) * dl
      END IF

      ! Update previous values
      x_prev = x
      y_prev = y
      f_prev = f_cur
   END DO

  END SUBROUTINE lineintegration

  SUBROUTINE point_in_mesh(x, y, Xc, T, is_inside, element_index)
     IMPLICIT NONE
     REAL*8, INTENT(IN) :: x, y       ! Coordinates of the point
     REAL*8, INTENT(IN) :: Xc(:,:)    ! Node coordinates
     INTEGER, INTENT(IN) :: T(:,:)    ! Connectivity matrix
     LOGICAL, INTENT(OUT) :: is_inside
     INTEGER, INTENT(OUT) :: element_index
     INTEGER :: i, j
     REAL*8 :: x_nodes(3), y_nodes(3)
  
     is_inside = .FALSE.
     element_index = -1
  
     ! Loop over all elements in the mesh
     DO i = 1, SIZE(T, 1)
         ! Get the coordinates of the nodes of the current element
         DO j = 1, 3
             x_nodes(j) = Xc(T(i, j), 1)
             y_nodes(j) = Xc(T(i, j), 2)
         END DO
        
         ! Check if the point (x, y) is inside the current element
         IF (point_in_triangle(x, y, x_nodes, y_nodes)) THEN
             is_inside = .TRUE.
             element_index = i
             RETURN
         END IF
     END DO
  END SUBROUTINE point_in_mesh

  LOGICAL FUNCTION point_in_triangle(x, y, x_nodes, y_nodes)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: x, y          ! Coordinates of the point
    REAL*8, INTENT(IN) :: x_nodes(3), y_nodes(3)  ! Triangle vertices
    REAL*8 :: area, area1, area2, area3

    ! Compute the area of the triangle
    area = 0.5D0 * ABS(x_nodes(1)*(y_nodes(2)-y_nodes(3)) + &
                       x_nodes(2)*(y_nodes(3)-y_nodes(1)) + &
                       x_nodes(3)*(y_nodes(1)-y_nodes(2)))

    ! Compute the areas of the sub-triangles formed with the point
    area1 = 0.5D0 * ABS(x*(y_nodes(2)-y_nodes(3)) + &
                        x_nodes(2)*(y_nodes(3)-y) + &
                        x_nodes(3)*(y-y_nodes(2)))

    area2 = 0.5D0 * ABS(x_nodes(1)*(y-y_nodes(3)) + &
                        x*(y_nodes(3)-y_nodes(1)) + &
                        x_nodes(3)*(y_nodes(1)-y))

    area3 = 0.5D0 * ABS(x_nodes(1)*(y_nodes(2)-y) + &
                        x_nodes(2)*(y-y_nodes(1)) + &
                        x*(y_nodes(1)-y_nodes(2)))

    ! Check if the sum of the sub-triangle areas equals the total area
    point_in_triangle = ABS(area - (area1 + area2 + area3)) < 1.0D-10
  END FUNCTION point_in_triangle

  !real function lineintegration(qp_len, x_vec, y_vec, np_len, x_array, y_array ,f)
  ! Given a set of points (x, y) along a line of sight, retruns the line integration of
  ! f. The function f is evaluated IN the closest nodes to points (x, y).
  !  implicit none
  !  real,intent(IN) :: qp_len, np_len
  !  real,dimension(qp_len),intent(IN) :: x_vec,  y_vec
  !  real,dimension(np_len),intent(IN) :: x_array,  y_array
  !  real, intent(IN) :: f
  !  integer :: i, iel, inp, n_ind, f_ind(2)
  !  real*8 :: x, y, dl(qp_len-1), x_qp(qp_len), f_qp(qp_len)

  ! Search closest node and evaluate f
  !  DO i = 1, qp_len
  !     x = x_vec(i)
  !     y = y_vec(i)
  !     n_ind = nodesearch(x, y, np_len, x_array, y_array)
  !     f_ind = findloc(Mesh%T(n_ind))
  !     iel = f_ind(1)
  !     ip = f_ind(2)
  !     x_qp(i) = x_array(n_ind)
  !     f_qp(i) = f((iel-1)*refElPol%N2D + ip)
  !  END DO

  ! Integrate
  !  dl = x_qp(2:qp_len) - x_qp(1:qp_len-1)
  !  lineintegration = sum(0.5*(f_qp(2:qp_len) + f_qp(1:qp_len-1))*dl)

  !end function lineintegration

END MODULE interpolation
