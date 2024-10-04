MODULE GMSH_io_module

  IMPLICIT NONE

CONTAINS


  SUBROUTINE load_gmsh_mesh(gmsh_filename, flip)
    USE globals
    !USE in_out


  !*****************************************************************************80
  !
    INTEGER ( kind = 4 ), ALLOCATABLE            :: element_node(:,:)
    INTEGER ( kind = 4 )                         :: element_num
    INTEGER ( kind = 4 )                         :: element_order
    CHARACTER (*), INTENT(IN)                    :: gmsh_filename
    INTEGER      , INTENT(IN)                    :: flip
    CHARACTER (len = 255)                        :: temp_filename
    INTEGER ( kind = 4 )                         :: m
    INTEGER ( kind = 4 )                         :: node_num
    REAL ( kind = 8 ), ALLOCATABLE               :: node_x(:,:)


    temp_filename      = TRIM(ADJUSTL(gmsh_filename))//'.msh'

    WRITE ( *, '(a)' ) ''
    WRITE ( *, '(a)' ) 'LOADING GMSH FILE'
  !
  !  Get the data size.
  !
    CALL gmsh_size_read ( temp_filename, node_num, m, element_num, &
      element_order )
  !
  !  Print the sizes.
  !
    WRITE ( *, '(a)' ) ''
    WRITE ( *, '(a)' ) '  Node data read from file "' // TRIM ( temp_filename ) // '"'
    WRITE ( *, '(a)' ) ''
    WRITE ( *, '(a,i4)' ) '  Number of nodes = ', node_num
    WRITE ( *, '(a,i4)' ) '  Spatial dimension = ', m
    WRITE ( *, '(a,i4)' ) '  Number of elements = ', element_num
    WRITE ( *, '(a,i4)' ) '  Element order = ', element_order
  !
  !  Allocate memory.
  !

    ALLOCATE ( node_x(1:m,1:node_num) )
    ALLOCATE ( element_node(1:element_order,1:element_num) )
  !
  !  Get the data.
  !

    CALL gmsh_data_read ( temp_filename, m, node_num, node_x, element_order, element_num, element_node, flip )


  !
  !  Clean up.
  !
    DEALLOCATE ( element_node )
    DEALLOCATE ( node_x )
    !stop
    RETURN
  END SUBROUTINE load_gmsh_mesh

  SUBROUTINE gmsh_data_read ( gmsh_filename, node_dim, node_num, node_x, &
    element_order, element_num, element_node, flip )
    USE globals

  !*****************************************************************************80
  !
  !! GMSH_DATA_READ reads data from a GMSH file.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    16 October 2014
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    inputt, character ( len = * ) GMSH_FILENAME, the GMSH filename.
  !
  !    inputt, integer ( kind = 4 ) NODE_DIM, the spatial dimension.
  !
  !    inputt, integer ( kind = 4 ) NODE_NUM, the number of nodes.
  !
  !    Output, real ( kind = 8 ) NODE_X(NODE_DIM,NODE_NUM), the node coordinates.
  !
  !    inputt, integer ( kind = 4 ) ELEMENT_ORDER, the order of the elements.
  !
  !    inputt, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
  !
  !    Output, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM),
  !    the nodes that make up each element.
  !
    USE globals
    USE reference_element, ONLY: generate_fekete_nodes, create_reference_element

    TYPE(Reference_element_type)             :: refEl
    INTEGER, INTENT(in)                      :: flip
    INTEGER ( kind = 4 ) element_num
    INTEGER ( kind = 4 ) element_order
    INTEGER ( kind = 4 ) node_dim
    INTEGER ( kind = 4 ) node_num, n_Tb_IN, n_Tb_OUT, n_Tb_LIM, n_Tb_PUMP, n_Tb_PUFF, n_T_DOM
    INTEGER ( kind = 4 ), DIMENSION(:,:), ALLOCATABLE :: face_info, Tb_Dirichlet, Tb_LEFT, Tb_RIGHT, Tb_UP, Tb_PUMP, Tb_PUFF, Tb_DOWN, Tb_WALL, Tb_LIM, Tb_IN, Tb_OUT, Tb_ULIM, T, Tb
    INTEGER    ( kind = 4 ), DIMENSION(:), ALLOCATABLE   :: temp
    INTEGER   ( kind = 4 ), DIMENSION(:), ALLOCATABLE :: boundaryFlag
    INTEGER ( kind = 4 ) i_Tb_IN, i_Tb_OUT, i_Tb_LIM, i_Tb_PUMP, i_Tb_PUFF,  i_T_DOM
    INTEGER ( kind = 4 ), ALLOCATABLE ::  indices(:)


    INTEGER ( kind = 4 ) element_node(element_order,element_num)
    CHARACTER ( len = * ) gmsh_filename
    INTEGER ( kind = 4 ) i, END
    INTEGER ( kind = 4 ) i4_dummy, i4_dummy_save
    INTEGER ( kind = 4 ) ierror
    INTEGER ( kind = 4 ) indx
    INTEGER ( kind = 4 ) inputt
    INTEGER ( kind = 4 ) input_stat
    INTEGER ( kind = 4 ) j
    INTEGER ( kind = 4 ) k
    INTEGER ( kind = 4 ) length
    INTEGER ( kind = 4 ) level
    REAL    ( kind = 8 ) node_x(node_dim,node_num), x
    LOGICAL flag
    CHARACTER ( len = 255 ) text

    n_Tb_IN = 0
    n_Tb_OUT = 0
    n_Tb_LIM = 0
    n_Tb_PUFF = 0
    n_Tb_PUMP = 0
    n_T_DOM = 0

    i_Tb_IN = 1
    i_Tb_OUT = 1
    i_Tb_LIM = 1
    i_Tb_PUFF = 1
    i_Tb_PUMP = 1
    i_T_DOM = 1


    CALL get_unit ( inputt )

    OPEN ( unit = inputt, file = gmsh_filename, status = 'old', &
      iostat = input_stat )

    IF ( input_stat /= 0 ) THEN
       WRITE ( *, '(a)' ) ''
       WRITE ( *, '(a)' ) 'GMSH_DATA_READ - Fatal error!'
       WRITE ( *, '(a)' ) '  Could not open inputt file "' // &
            TRIM ( gmsh_filename ) // '"'
       STOP 1
    END IF

    level = 0

    DO

       READ ( inputt, '(a)', iostat = input_stat ) text

       IF ( input_stat /= 0 ) THEN
          WRITE ( *, '(a)' ) 'GMSH_DATA_READ - Fatal error!'
          WRITE ( *, '(a)' ) '  Error while seeking node coordinates.'
          STOP 1
       END IF

       IF ( level == 0 ) THEN
          CALL s_begin ( text(1:6), '$Nodes', flag )
          IF (flag)  THEN
          level = 1
          j = 0
          END IF
       ELSE IF ( level == 1 ) THEN
          CALL s_to_i4 ( text, i4_dummy, ierror, length )
        level = 2
       ELSE IF ( level == 2 ) THEN
          CALL s_begin ( text(1:9), '$EndNodes' , flag)
          IF ( flag ) THEN
             EXIT
          ELSE
          j = j + 1
             CALL s_to_i4 ( text, indx, ierror, length )
          text = text(length+1:)
             DO i = 1, node_dim
                CALL s_to_r8 ( text, x, ierror, length )
            text = text(length+1:)
            node_x(i,j) = x
             END DO
          END IF
       END IF

    END DO
  !
  !  Now read element information.
  !
    level = 0

    DO

       READ ( inputt, '(a)', iostat = input_stat ) text

       IF ( input_stat /= 0 ) THEN
          WRITE ( *, '(a)' ) 'GMSH_DATA_READ - Fatal error!'
          WRITE ( *, '(a)' ) '  Error while seeking element connectivity.'
          STOP 1
       END IF

       IF ( level == 0 ) THEN
          CALL s_begin ( text(1:9), '$Elements', flag )
          IF ( flag ) THEN
          level = 1
          j = 0
          END IF
       ELSE IF ( level == 1 ) THEN
          CALL s_to_i4 ( text, i4_dummy, ierror, length )
        level = 2
       ELSE IF ( level == 2 ) THEN
          CALL s_begin ( text(1:12), '$EndElements' , flag)
          IF ( flag) THEN
             EXIT
          ELSE
          j = j + 1
          k = 0
             DO k = 1, 5 ! loop over the first 5 numbers in the row
                CALL s_to_i4 ( text, i4_dummy, ierror, length )

            !! count number of elements in the boundary/domain
                IF (k .EQ. 4) THEN
                   SELECT CASE (i4_dummy)
                   CASE (1)
                  n_Tb_IN = n_Tb_IN + 1
                   CASE (2)
                  n_Tb_OUT = n_Tb_OUT + 1
                   CASE (3)
                  n_Tb_LIM = n_Tb_LIM + 1
                   CASE (4)
                  n_T_DOM = n_T_DOM + 1
                   CASE (5)
                  n_Tb_PUMP = n_Tb_PUMP + 1
                   CASE (6)
                  n_Tb_PUFF = n_Tb_PUFF + 1
                   END SELECT
             ENDIF
             text = text(length+1:)
             END DO

             DO i = 1, element_order
                CALL s_to_i4 ( text, k, ierror, length )
            text = text(length+1:)
            element_node(i,j) = k
             END DO
          END IF
       END IF
    END DO


    rewind(inputt)

    ! element order is equal to p + 1
    ALLOCATE(Tb_IN(n_Tb_IN, element_order))
    ALLOCATE(Tb_OUT(n_Tb_OUT, element_order))
    ALLOCATE(Tb_LIM(n_Tb_LIM, element_order))
    ALLOCATE(Tb_PUFF(n_Tb_PUFF, element_order))
    ALLOCATE(Tb_PUMP(n_Tb_PUMP, element_order))
    ALLOCATE(T(n_T_DOM, (element_order)*(element_order+1)/2))

    !
    !  Now fill Tb_IN, Tb_OUt, Tb_LIM and T
    !
      level = 0

      do

        read ( inputt, '(a)', iostat = input_stat ) text

        if ( input_stat /= 0 ) then
          write ( *, '(a)' ) 'GMSH_DATA_READ - Fatal error!'
          write ( *, '(a)' ) '  Error while seeking element connectivity.'
          stop 1
        end if

        if ( level == 0 ) then
          call s_begin ( text(1:9), '$Elements' , flag)
          if ( flag  ) then
            level = 1
            j = 0
          end if
        else if ( level == 1 ) then
          call s_to_i4 ( text, i4_dummy, ierror, length )
          level = 2
        else if ( level == 2 ) then
          call s_begin ( text(1:12), '$EndElements' , flag)
          if ( flag  ) then
            exit
          else
            j = j + 1
            k = 0
            do k = 1, 5 ! empty loop over the first 4 numbers of the row
              call s_to_i4 ( text, i4_dummy, ierror, length )
              IF(k .eq. 4) THEN
                i4_dummy_save = i4_dummy
              ENDIF
               text = text(length+1:)
            end do

            IF (i4_dummy_save .eq. 4) THEN
            end = (element_order)*(element_order+1)/2
            ELSE
              end = element_order
            ENDIF

            do i = 1, end
              call s_to_i4 ( text, k, ierror, length )

              select case (i4_dummy_save)
               case (1)
                 Tb_IN(i_Tb_IN,i) = k
               case (2)
                 Tb_OUT(i_Tb_OUT,i) = k
                case (3)
                  Tb_LIM(i_Tb_LIM,i) =k
                case (4)
                  T(i_T_DOM,i) = k
                case (5)
                  Tb_PUMP(i_Tb_PUMP,i) =k
                case (6)
                  Tb_PUFF(i_Tb_PUFF,i) =k
              end select
              text = text(length+1:)
            end do

            select case (i4_dummy_save)
             case (1)
               i_Tb_IN = i_Tb_IN + 1
             case (2)
               i_Tb_OUT = i_Tb_OUT + 1
              case (3)
                i_Tb_LIM = i_Tb_LIM + 1
              case (4)
                i_T_DOM = i_T_DOM + 1
              case (5)
                i_Tb_PUMP = i_Tb_PUMP + 1
              case (6)
                i_Tb_PUFF = i_Tb_PUFF + 1
            end select
          end if
        end if
      end do

    close ( unit = inputt )

    IF (element_order-1 .gt. 5) THEN
      allocate(indices(element_order*(element_order+1)/2))
    ENDIF

    ! Reshape Tb_PUFF
    ! shift colomns to the left by one
    IF (n_Tb_PUFF .ne. 0) THEN
      ALLOCATE(temp(n_Tb_PUFF))

      temp = Tb_PUFF(:, 2)
      do j = 1, element_order - 2
        Tb_PUFF(:, j + 1) = Tb_PUFF(:, j + 2)
      enddo
      Tb_PUFF(:,element_order) = temp

      DEALLOCATE(temp)
    ENDIF

    ! Reshape Tb_PUMP
    ! shift colomns to the left by one
    IF (n_Tb_PUMP .ne. 0) THEN
      ALLOCATE(temp(n_Tb_PUMP))

      temp = Tb_PUMP(:, 2)
      do j = 1, element_order - 2
        Tb_PUMP(:, j + 1) = Tb_PUMP(:, j + 2)
      enddo
      Tb_PUMP(:,element_order) = temp

      DEALLOCATE(temp)
    ENDIF

    ! Reshape Tb_OUT
    ! shift last 4 colomns to the left by one
    IF (n_Tb_OUT /= 0) THEN
      ALLOCATE(temp(n_Tb_OUT))

      temp = Tb_OUT(:, 2)
      do j = 1, element_order - 2
        Tb_OUT(:, j + 1) = Tb_OUT(:, j + 2)
      enddo
      Tb_OUT(:,element_order) = temp

      DEALLOCATE(temp)
    ELSE
      WRITE(*,*) "There's something wrong: Tb_OUT doesn't make sense empty."
      stop
    ENDIF

    ! Reshape Tb_LIM
    ! shift colomns to the left by one
    IF (n_Tb_LIM .ne. 0) THEN
      ALLOCATE(temp(n_Tb_LIM))

      temp = Tb_LIM(:, 2)
      do j = 1, element_order - 2
        Tb_LIM(:, j + 1) = Tb_LIM(:, j + 2)
      enddo
      Tb_LIM(:,element_order) = temp

      DEALLOCATE(temp)
    ENDIF

    ! Reshape Tb_IN if n_Tb_IN is not zero
    if ((n_Tb_IN .ne. 0) .and. (flip .eq. 1)) then
      ALLOCATE(temp(n_Tb_IN))
      if ((switch%testcase .lt. 50) .or. (switch%testcase .gt. 59)) then

        ! swap first and second column
        temp = Tb_IN(:, 2)
        Tb_IN(:, 2) =  Tb_IN(:, 1)
        Tb_IN(:, 1) = temp

        ! Swap elements in columns 3 to element order (p+1)
         do j = 1, int((element_order-2)/2)
            temp = Tb_IN(:, j + 2)
            Tb_IN(:, j + 2) = Tb_IN(:, element_order - j + 1)
            Tb_IN(:, element_order - j + 1) = temp
          end do
      endif

      ! shift colomns to the left by one
      temp = Tb_IN(:, 2)
      do j = 1, element_order - 2
        Tb_IN(:, j + 1) = Tb_IN(:, j + 2)
      enddo
      Tb_IN(:,element_order) = temp

      DEALLOCATE(temp)
    end if


    ! Reshape T based on the value of p
    IF(element_order-1 .ge. 5) THEN
      select case (element_order-1)
      case (5)
        indices = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,19,17,21,20,18]
      case (6)
        !indices = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,19,22,23,20,27,28,24,26,25,21]
        indices = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,22,23,20,27,28,24,26,25,21]
      case (7)
        !indices = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,19,22,23,33,34,35,28,32,36,29,31,30,24]
        indices = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,25,26,27,23,33,34,35,28,32,36,29,31,30,24]
      case (8)
        !indices = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,19,22,23,25,28,29,30,31,33,34,35,36,37,38,39,40,41,42,43,44,45]
        indices = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,28,29,30,31,26,39,40,43,41,32,38,45,44,33,37,42,34,36,35,27]
      end select
      T = T(:, indices)
      DEALLOCATE(indices)
    ENDIF

    CALL generate_elemface_info(T,Tb_IN, Tb_LIM, Tb_PUFF, Tb_PUMP, Tb_OUT, element_order, face_info)
    CALL generate_boundary_names(Tb_Dirichlet, Tb_LEFT, Tb_RIGHT, Tb_UP, Tb_DOWN, Tb_WALL, Tb_LIM, Tb_IN, Tb_OUT, Tb_PUFF, Tb_PUMP, Tb_ULIM, Tb, boundaryFlag, element_order)
    CALL load_mesh2global_var(SIZE(node_x,1), SIZE(T,1), SIZE(Tb,1), SIZE(node_x,2), SIZE(T,2), element_order, 0, T, transpose(node_x), Tb, boundaryFlag, face_info)
    CALL create_reference_element(refEl, SIZE(node_x,1), verbose = 0)
    CALL generate_fekete_nodes(Mesh%X,Mesh%T, element_order-1, refEl)

    DEALLOCATE(Tb_IN, Tb_OUT, Tb_LIM, T, Tb_PUFF, Tb_PUMP, boundaryFlag)

    CALL free_reference_element_pol(refEl)

  END


  SUBROUTINE generate_boundary_names(Tb_Dirichlet, Tb_LEFT, Tb_RIGHT, Tb_UP, Tb_DOWN, Tb_WALL, Tb_LIM, Tb_IN, Tb_OUT, Tb_PUFF, Tb_PUMP, Tb_ULIM, Tb, boundaryFlag, element_order)

    INTEGER, DIMENSION(:,:), INTENT(IN), ALLOCATABLE     :: Tb_Dirichlet, Tb_LEFT, Tb_RIGHT, Tb_UP, Tb_DOWN, Tb_WALL, Tb_LIM, Tb_IN, Tb_OUT, Tb_PUFF, Tb_PUMP, Tb_ULIM
    INTEGER, DIMENSION(:), INTENT(OUT), ALLOCATABLE      :: boundaryFlag
    INTEGER, DIMENSION(:,:), INTENT(OUT), ALLOCATABLE    :: Tb
    INTEGER, INTENT(IN)                                  :: element_order
    INTEGER                                              :: n_boundaries
    INTEGER                                              :: s1, s2, index

    n_boundaries = 0
    s1 = 0

    IF(allocated(Tb_Dirichlet) ) THEN
      n_boundaries = n_boundaries + 1
      s1 = s1 + SIZE(Tb_Dirichlet,1)
    ENDIF
    IF(allocated(Tb_LEFT) ) THEN
      n_boundaries = n_boundaries + 1
      s1 = s1 + SIZE(Tb_LEFT,1)
    ENDIF
    IF(allocated(Tb_RIGHT) ) THEN
      n_boundaries = n_boundaries + 1
      s1 = s1 + SIZE(Tb_RIGHT,1)
    ENDIF
    IF(allocated(Tb_UP) ) THEN
      n_boundaries = n_boundaries + 1
      s1 = s1 + SIZE(Tb_UP,1)
    ENDIF
    IF(allocated(Tb_DOWN) ) THEN
      n_boundaries = n_boundaries + 1
      s1 = s1 + SIZE(Tb_DOWN,1)
    ENDIF
    IF(allocated(Tb_WALL) ) THEN
      n_boundaries = n_boundaries + 1
      s1 = s1 + SIZE(Tb_WALL,1)
    ENDIF
    IF(allocated(Tb_IN) ) THEN
      n_boundaries = n_boundaries + 1
      s1 = s1 + SIZE(Tb_IN,1)
    ENDIF
    IF(allocated(Tb_LIM) ) THEN
      n_boundaries = n_boundaries + 1
      s1 = s1 + SIZE(Tb_LIM,1)
    ENDIF
    IF(allocated(Tb_OUT) ) THEN
      n_boundaries = n_boundaries + 1
      s1 = s1 + SIZE(Tb_OUT,1)
    ENDIF
    IF(allocated(Tb_PUFF) ) THEN
      n_boundaries = n_boundaries + 1
      s1 = s1 + SIZE(Tb_PUFF,1)
    ENDIF
    IF(allocated(Tb_PUMP) ) THEN
      n_boundaries = n_boundaries + 1
      s1 = s1 + SIZE(Tb_PUMP,1)
    ENDIF
    IF(allocated(Tb_ULIM) ) THEN
      n_boundaries = n_boundaries + 1
      s1 = s1 + SIZE(Tb_ULIM,1)
    ENDIF

    s2 = element_order

    ALLOCATE(Tb(s1,s2))
    ALLOCATE(boundaryFlag(s1))

    Tb = 0.
    boundaryFlag = 0.


    s1 = 1
    IF(allocated(Tb_Dirichlet) ) THEN
      index = s1+SIZE(Tb_Dirichlet,1)
      Tb(s1:index-1,:) = Tb_Dirichlet
      boundaryFlag(s1:index-1) = 1
      s1 = index
    ENDIF
    IF(allocated(Tb_LEFT) ) THEN
      index = s1+SIZE(Tb_LEFT,1)
      Tb(s1:index-1,:) = Tb_LEFT
      boundaryFlag(s1:index-1) = 2
      s1 = index
    ENDIF
    IF(allocated(Tb_RIGHT) ) THEN
      index = s1+SIZE(Tb_RIGHT,1)
      Tb(s1:index-1,:) = Tb_RIGHT
      boundaryFlag(s1:index-1) = 3
      s1 = index
    ENDIF
    IF(allocated(Tb_UP) ) THEN
      index = s1+SIZE(Tb_UP,1)
      Tb(s1:index-1,:) = Tb_UP
      boundaryFlag(s1:index-1) = 4
      s1 = index
    ENDIF
    IF(allocated(Tb_PUMP) ) THEN
      index = s1+SIZE(Tb_PUMP,1)
      Tb(s1:index-1,:) = Tb_PUMP
      boundaryFlag(s1:index-1) = 5
      s1 = index
    ENDIF
    IF(allocated(Tb_PUFF) ) THEN
      index = s1+SIZE(Tb_PUFF,1)
      Tb(s1:index-1,:) = Tb_PUFF
      boundaryFlag(s1:index-1) = 6
      s1 = index
    ENDIF
    IF(allocated(Tb_IN) ) THEN
      index = s1+SIZE(Tb_IN,1)
      Tb(s1:index-1,:) = Tb_IN
      boundaryFlag(s1:index-1) = 8
      s1 = index
    ENDIF
    IF(allocated(Tb_LIM) ) THEN
      index = s1+SIZE(Tb_LIM,1)
      Tb(s1:index-1,:) = Tb_LIM
      boundaryFlag(s1:index-1) = 7
      s1 = index
    ENDIF
    IF(allocated(Tb_OUT) ) THEN
      index = s1+SIZE(Tb_OUT,1)
      Tb(s1:index-1,:) = Tb_OUT
      boundaryFlag(s1:index-1) = 9
      s1 = index
    ENDIF
    IF(allocated(Tb_ULIM) ) THEN
      index = s1+SIZE(Tb_ULIM,1)
      Tb(s1:index-1,:) = Tb_ULIM
      boundaryFlag(s1:index-1) = 10
      s1 = index
    ENDIF

  ENDSUBROUTINE

  SUBROUTINE generate_elemface_info(T, Tb_IN, Tb_LIM, Tb_PUFF, Tb_PUMP, Tb_OUT, element_order, face_info)
      INTEGER, INTENT(IN)                              :: T(:,:)
      INTEGER, INTENT(IN)                              :: element_order
      INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Tb_IN, Tb_LIM,Tb_PUMP, Tb_PUFF, Tb_OUT
      INTEGER, ALLOCATABLE, INTENT(OUT)                :: face_info(:,:)
      INTEGER, ALLOCATABLE                             :: aux_extfaces(:,:), aux_Tb(:,:)
      INTEGER                                          :: n_Tb_vec(5)
      INTEGER                                          :: nodes(2)
      INTEGER                                          :: n_faces, ifa, iel, n_elements, n_boundaries, el_index, loc_fa, elemFaceInfo_fa
      INTEGER                                          :: i, counter

      n_elements = SIZE(T,1)
      n_Tb_vec = 0
      IF(ALLOCATED(Tb_PUMP)) THEN
        n_Tb_vec(1) = SIZE(Tb_PUMP,1)
      ENDIF
      IF(ALLOCATED(Tb_PUFF)) THEN
        n_Tb_vec(2) = SIZE(Tb_PUFF,1)
      ENDIF
      IF(ALLOCATED(Tb_IN)) THEN
        n_Tb_vec(3) = SIZE(Tb_IN,1)
      ENDIF
      IF(ALLOCATED(Tb_LIM)) THEN
        n_Tb_vec(4) = SIZE(Tb_LIM,1)
      ENDIF
      IF(ALLOCATED(Tb_OUT)) THEN
        n_Tb_vec(5) = SIZE(Tb_OUT,1)
      ENDIF

      n_boundaries = SIZE(n_Tb_vec)

      ALLOCATE(face_info(SUM(n_Tb_vec), 2))

      DO i = 1, n_boundaries
          IF(n_Tb_vec(i) .ne. 0) THEN
              ALLOCATE(aux_extfaces(n_Tb_vec(i), 2))
              ALLOCATE(aux_Tb(n_Tb_vec(i), element_order))
              SELECT CASE(i)
                  CASE (1)
                      aux_Tb = Tb_PUMP
                  CASE (2)
                      aux_Tb = Tb_PUFF
                  CASE (3)
                      aux_Tb = Tb_IN
                  CASE (4)
                      aux_Tb = Tb_LIM
                  CASE (5)
                      aux_Tb = Tb_OUT
              END SELECT
          ELSE
              CYCLE
          ENDIF

          n_faces = n_Tb_vec(i)
          el_index = 0
          loc_fa = 0
          counter = 0

          DO ifa = 1, n_faces
              IF(i .eq. 1) THEN
                  nodes(1) = aux_Tb(ifa, element_order)
                  nodes(2) = aux_Tb(ifa, 1)
              ELSE
                  nodes(1) = aux_Tb(ifa, 1)
                  nodes(2) = aux_Tb(ifa, element_order)
              ENDIF

              DO iel = 1, n_elements
                  counter = COUNT(T(iel, 1:3) .eq. nodes(1))
                  counter = counter + COUNT(T(iel, 1:3) .eq. nodes(2))
                  IF(counter .eq. 2) THEN
                      el_index = iel
                      EXIT
                  ENDIF
                  counter = 0
              ENDDO

              IF(el_index .eq. 0) THEN
                  PRINT *, "Error in generate_elemface_info: element not found. STOP."
                  STOP
              ENDIF

              IF (equality(T(el_index, 1),nodes(1)) .OR. (equality(T(el_index, 1),nodes(2)))) THEN
                  ! Check for loc_fa = 1 or 3
                  IF (equality(T(el_index, 2),nodes(1)) .OR. (equality(T(el_index, 2),nodes(2)))) THEN
                      loc_fa = 1
                  ELSEIF(equality(T(el_index, 3),nodes(1)) .OR. (equality(T(el_index, 3),nodes(2)))) THEN
                      loc_fa = 3
                  END IF
              ELSEIF (equality(T(el_index, 2),nodes(1)) .OR. (equality(T(el_index, 2),nodes(2)))) THEN
                  ! Check for loc_fa = 2 or 1
                  IF (equality(T(el_index, 3),nodes(1)) .OR. (equality(T(el_index, 3),nodes(2)))) THEN
                      loc_fa = 2
                  ELSEIF(equality(T(el_index, 1),nodes(1)) .OR. (equality(T(el_index, 1),nodes(2)))) THEN
                      loc_fa = 1
                  END IF
              ELSEIF (equality(T(el_index, 3),nodes(1)) .OR. (equality(T(el_index, 3),nodes(2)))) THEN
                  ! Check for loc_fa = 3 or elemFaceInfo_fa = 2
                  IF (equality(T(el_index, 1),nodes(1)) .OR. (equality(T(el_index, 1),nodes(2)))) THEN
                      loc_fa = 3
                  ELSEIF(equality(T(el_index, 2),nodes(1)) .OR. (equality(T(el_index, 2),nodes(2)))) THEN
                      elemFaceInfo_fa = 2
                  END IF
              ELSE
                  ! Error condition
                  PRINT *, 'Something is wrong in gmsh_io'
                  STOP
              END IF

              aux_extfaces(ifa,1) = el_index
              aux_extfaces(ifa,2) = loc_fa

          ENDDO

          IF(i .eq. 1) THEN
              face_info(1:n_Tb_vec(i),:) = aux_extfaces(:,:)
          ELSE
              face_info(SUM(n_Tb_vec(1:i-1))+1:SUM(n_Tb_vec(1:i)),:) = aux_extfaces(:,:)
          ENDIF

          DEALLOCATE(aux_extfaces)
          DEALLOCATE(aux_Tb)
      ENDDO

  contains

    function equality (input1,input2) result(flag)
      INTEGER, intent(IN) :: input1
      INTEGER, intent(IN) :: input2
      logical            :: flag

      IF(input1 .eq. input2) THEN
        flag = .true.
      ELSE
        flag = .false.
      ENDIF

      return
    endfunction
  ENDSUBROUTINE

  subroutine convert_gmsh_to_hdf5(h5_filename, Ndim, Nelems, Nextfaces, Nnodes, Nnodesperelem, Nnodesperface, elemType, T, X, Tb, boundaryFlag)
    USE globals
    character(LEN=*), INTENT(IN) :: h5_filename
    INTEGER, INTENT(IN)          :: Ndim, Nelems, Nextfaces, Nnodes, Nnodesperelem, Nnodesperface, elemType
    INTEGER, INTENT(IN)          :: T(:,:), Tb(:,:), boundaryFlag(:)
    REAL, INTENT(IN)             :: X(:,:)
    REAL                         :: temp(SIZE(X,1),SIZE(X,2))
    real                         :: xmin

    temp = X
    !Apply length scale back
    temp = temp*phys%lscale

    xmin = minval(temp(:,1))
    ! Apply shift back if axisymmetric case
    IF ((switch%axisym .and. switch%testcase .ge. 60 .and. switch%testcase .lt. 80)) THEN
      temp(:, 1) = temp(:, 1) - geom%R0
    END IF


    CALL HDF5_save_mesh(h5_filename, Ndim, Nelems, Nextfaces, Nnodes, Nnodesperelem, Nnodesperface, elemType, T, temp, Tb, boundaryFlag)

  ENDSUBROUTINE

  SUBROUTINE HDF5_save_mesh(fname, Ndim, Nelems, Nextfaces, Nnodes, Nnodesperelem, Nnodesperface, elemType, T, X, Tb, boundaryFlag)
    USE HDF5
    USE HDF5_io_module
    USE MPI_OMP
    USE GLOBALS

    character(LEN=*), INTENT(IN) :: fname
    INTEGER                      :: ierr
    INTEGER(HID_T)               :: file_id
    INTEGER, INTENT(IN)          :: Ndim, Nelems, Nextfaces, Nnodes, Nnodesperelem, Nnodesperface, elemType
    INTEGER, INTENT(IN)          :: T(:,:), Tb(:,:), boundaryFlag(:)
    REAL, INTENT(IN)             :: X(:,:)
    CHARACTER(LEN=500)           :: fname_mod, Num
    INTEGER                      :: i

    fname_mod = trim(adjustl(fname))

    i = INDEX(fname_mod, 'P', .TRUE.)
    IF(i .ne. 0) THEN
      WRITE (Num, "(i10)") refElPol%nDeg
      fname_mod(i+1:i+1) = TRIM(ADJUSTL(Num))
    ENDIF

    call HDF5_create(fname_mod, file_id, ierr)
    call HDF5_integer_saving(file_id,Ndim,'Ndim')
    call HDF5_integer_saving(file_id,Nelems,'Nelems')
    call HDF5_integer_saving(file_id,Nextfaces,'Nextfaces')
    call HDF5_integer_saving(file_id,Nnodes,'Nnodes')
    call HDF5_integer_saving(file_id,Nnodesperelem,'Nnodesperelem')
    call HDF5_integer_saving(file_id,Nnodesperface,'Nnodesperface')
    call HDF5_integer_saving(file_id,elemType,'elemType')
    call HDF5_array2D_saving(file_id, real(T), size(T, 1), size(T, 2), 'T')
    call HDF5_array2D_saving(file_id, real(Tb), size(Tb, 1), size(Tb, 2), 'Tb')
    call HDF5_array2D_saving(file_id, X, size(X, 1), size(X, 2), 'X')
    call HDF5_array1D_saving(file_id, real(boundaryFlag), SIZE(boundaryFlag), 'boundaryFlag')
    call HDF5_close(file_id)

    ! Message to confirm succesful creation and filling of file
    IF (MPIvar%glob_id .eq. 0) THEN
      print *, 'Mesh converted from gmsh to hdf5. Output written to file ', trim(adjustl(fname_mod))
      print *, '        '
    END IF
  ENDSUBROUTINE

  !********************************
  ! Loads mesh from an hdf5 file
  ! external file
  !********************************
  SUBROUTINE load_mesh2global_var(Ndim, Nelems, Nextfaces, Nnodes, Nnodesperelem, Nnodesperface, elemType, T, X, Tb, boundaryFlag, face_info)
    USE MPI_OMP
    USE globals
    USE printutils

    INTEGER, INTENT(IN)          :: Ndim, Nelems, Nextfaces, Nnodes, Nnodesperelem, Nnodesperface, elemType
    INTEGER, INTENT(IN)          :: T(:,:), Tb(:,:), boundaryFlag(:),face_info(:,:)
    REAL*8, INTENT(IN)           :: X(:,:)
    character(10)                :: str
    real*8, parameter            :: tol = 1e-6
    real*8                       :: xmin

    IF (utils%printint > 0) THEN
      print *, 'Loading mesh.'
      print *, '        '
    ENDIF

    ALLOCATE (Mesh%T(Nelems, Nnodesperelem))
    ALLOCATE (Mesh%X(Nnodes, Ndim))
    ALLOCATE (Mesh%Tb(Nextfaces, Nnodesperface))
    ALLOCATE (Mesh%boundaryFlag(Nextfaces))
    ALLOCATE (Mesh%face_info(Nextfaces,2))

    Mesh%T = T
    Mesh%Tb = Tb
    Mesh%boundaryFlag = boundaryFlag
    Mesh%X = X
    Mesh%face_info = face_info

    Mesh%Ndim = Ndim
    Mesh%Nnodes = Nnodes
    Mesh%Nelems = Nelems
    Mesh%Nnodesperelem = Nnodesperelem
    Mesh%Nnodesperface = Nnodesperface
    Mesh%elemType = elemType
    Mesh%Nextfaces = Nextfaces

    xmin = minval(Mesh%X(:,1))

    ! Apply shift if axisymmetric case
    IF ((switch%axisym .and. switch%testcase .ge. 60 .and. switch%testcase .lt. 80) .or. (switch%axisym .and. xmin < tol)) THEN
      IF (MPIvar%glob_id .eq. 0) THEN
        WRITE (6, *) "*** Applying translation in axisymmetric case!"
      ENDIF
      Mesh%X(:, 1) = Mesh%X(:, 1) + geom%R0
    END IF

    !Apply length scale
    Mesh%X = Mesh%X/phys%lscale

    Mesh%xmax = maxval(Mesh%X(:, 1))
    Mesh%xmin = minval(Mesh%X(:, 1))
    Mesh%ymax = maxval(Mesh%X(:, 2))
    Mesh%ymin = minval(Mesh%X(:, 2))

    IF (utils%printint > 0) then
      IF (MPIvar%glob_id .eq. 0) THEN
        IF (elemType == 0) then
          WRITE (str, '(A)') 'triangles'
        ELSEIF (elemType == 1) then
          WRITE (str, '(A)') 'quads'
        ELSEIF (elemType == 2) then
          WRITE (str, '(A)') 'thetra'
        ELSEIF (elemType == 3) then
          WRITE (str, '(A)') 'hexa'
        END IF
        WRITE (6, *) '*************************************************'
        WRITE (6, *) '*                    MESH                       *'
        WRITE (6, *) '*************************************************'
        WRITE (6, '(A,I18)') ' Number of dimensions:         ', Ndim
        WRITE (6, '(A,A34)') ' Element type: ', trim(str)
        WRITE (6, '(A,I18)') ' Number of elements:           ', Nelems
        WRITE (6, '(A,I18)') ' Number of nodes:              ', Nnodes
        WRITE (6, '(A,I18)') ' Number of nodes per element:  ', Nnodesperelem
        WRITE (6, '(A,I18)') ' Number of nodes per face:     ', Nnodesperface
        WRITE (6, '(A,I18)') ' Number of exterior faces:     ', Nextfaces
        WRITE (6, *) ' '
        WRITE (6, *) ' '
        IF (utils%printint > 1) THEN
          WRITE (6, *) "Connectivity matrix T:"
          CALL displayMatrixInt(Mesh%T)
          WRITE (6, *) "Boundary connectivity matrix Tb:"
          CALL displayMatrixInt(Mesh%Tb)
        END IF
      ENDIF
    ENDIF

  END SUBROUTINE load_mesh2global_var

  SUBROUTINE gmsh_size_read ( gmsh_filename, node_num, node_dim, element_num, element_order )

    !*****************************************************************************80
    !
    !! GMSH_SIZE_READ reads sizes from a GMSH file.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    16 October 2014
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    inputt, character ( len = * ) GMSH_FILENAME, the GMSH filename.
    !
    !    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes.
    !
    !    Output, integer ( kind = 4 ) NODE_DIM, the spatial dimension.
    !
    !    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
    !
    !    Output, integer ( kind = 4 ) ELEMENT_ORDER, the order of the elements.
    !


    INTEGER ( kind = 4 ) element_num
    INTEGER ( kind = 4 ) element_order
    CHARACTER ( len = * ) gmsh_filename
    INTEGER ( kind = 4 ) ierror
    INTEGER ( kind = 4 ) indx
    INTEGER ( kind = 4 ) inputt
    INTEGER ( kind = 4 ) input_stat
    INTEGER ( kind = 4 ) k
    INTEGER ( kind = 4 ) length
    INTEGER ( kind = 4 ) level
    INTEGER ( kind = 4 ) node_dim
    INTEGER ( kind = 4 ) node_num
    REAL ( kind = 8 ), PARAMETER :: r8_big = 1.0D+15
    LOGICAL flag
    CHARACTER ( len = 255 ) text
    REAL ( kind = 8 ) x
    REAL ( kind = 8 ) x_max
    REAL ( kind = 8 ) x_min
    REAL ( kind = 8 ) y
    REAL ( kind = 8 ) y_max
    REAL ( kind = 8 ) y_min
    REAL ( kind = 8 ) z
    REAL ( kind = 8 ) z_max
    REAL ( kind = 8 ) z_min

    node_num = 0
    node_dim = 0

    x_max = - r8_big
    x_min = + r8_big
    y_max = - r8_big
    y_min = + r8_big
    z_max = - r8_big
    z_min = + r8_big

    CALL get_unit ( inputt )

    OPEN ( unit = inputt, file = gmsh_filename, status = 'old', &
         iostat = input_stat )

    IF ( input_stat /= 0 ) THEN
       WRITE ( *, '(a)' ) ''
       WRITE ( *, '(a)' ) 'GMSH_SIZE_READ - Fatal error!'
       WRITE ( *, '(a)' ) '  Could not open inputt file "' // &
            TRIM ( gmsh_filename ) // '"'
       STOP 1
    END IF

    level = 0

    DO

       READ ( inputt, '(a)', iostat = input_stat ) text

       IF ( level == 0 ) THEN
          CALL s_begin ( text(1:6), '$Nodes' , flag)
          IF ( flag ) THEN
             level = 1
          END IF
       ELSE IF ( level == 1 ) THEN
          CALL s_to_i4 ( text, node_num, ierror, length )
          level = 2
       ELSE IF ( level == 2 ) THEN
          CALL  s_begin ( text(1:9), '$EndNodes' , flag)
          IF ( flag ) THEN
             EXIT
          ELSE
             CALL s_to_i4 ( text, indx, ierror, length )
             text = text(length+1:)
             CALL s_to_r8 ( text, x, ierror, length )
             x_min = MIN ( x_min, x )
             x_max = MAX ( x_max, x )
             text = text(length+1:)
             !
             !  Need to check that we actually were able to read an R8 here.
             !
             CALL s_to_r8 ( text, y, ierror, length )
             y_min = MIN ( y_min, y )
             y_max = MAX ( y_max, y )
             text = text(length+1:)
             CALL s_to_r8 ( text, z, ierror, length )
             text = text(length+1:)
             z_min = MIN ( z_min, z )
             z_max = MAX ( z_max, z )
          END IF
       END IF

    END DO
    !
    !  Make a vHDF5_integer_savinge guess as to the dimensionality of the data.
    !
    node_dim = 3
    IF ( ABS(z_max - z_min) .LE. 1e-12 ) THEN
       node_dim = 2
       IF ( ABS(y_max - y_min) .LE. 1e-12 ) THEN
          node_dim = 1
       END IF
    END IF
    !
    !  Now read element information.
    !
    level = 0

    DO

       READ ( inputt, '(a)', iostat = input_stat ) text

       IF ( level == 0 ) THEN
          CALL  s_begin ( text(1:9), '$Elements', flag )
          IF ( flag ) THEN
             level = 1
          END IF
       ELSE IF ( level == 1 ) THEN
          CALL s_to_i4 ( text, element_num, ierror, length )
          level = 2
       ELSE IF ( level == 2 ) THEN
          CALL s_begin ( text(1:12), '$EndElements' , flag)
          IF (flag ) THEN
             EXIT
          ELSE
             k = 0
             DO
                CALL s_to_i4 ( text, indx, ierror, length )
                text = text(length+1:)

                IF ( ierror /= 0 ) THEN
                   EXIT
                END IF
                k = k + 1
             END DO

             element_order = k - 5
             EXIT
          END IF
       END IF


    END DO

    CLOSE ( unit = inputt )

    RETURN
  END SUBROUTINE gmsh_size_read

  SUBROUTINE gmsh_mesh1d_write ( gmsh_filename, m, node_num, node_x, &
     element_order, element_num, element_node )

    !*****************************************************************************80
    !
    !! GMSH_MESH1D_WRITE writes 1D mesh data as a Gmsh file.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    06 October 2014
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Christophe Geuzaine, Jean-Francois Remacle,
    !    Gmsh: a three-dimensional finite element mesh generator with
    !    built-in pre- and post-processing facilities,
    !    International Journal for Numerical Methods in Engineering,
    !    Volume 79, Number 11, pages 1309-1331, 2009.
    !
    !  Parameters:
    !
    !    inputt, character * ( * ) GMSH_FILENAME, the name of the Gmsh file.
    !
    !    inputt, integer ( kind = 4 ) M, the spatial dimension.
    !
    !    inputt, integer ( kind = 4 ) NODE_NUM, the number of nodes.
    !
    !    inputt, real ( kind = 8 ) NODE_X(M,NODE_NUM), the node coordinates.
    !
    !    inputt, integer ( kind = 4 ) ELEMENT_ORDER, the order of the elements.
    !
    !    inputt, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
    !
    !    inputt, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM),
    !    the nodes that make up each element.
    !


    INTEGER ( kind = 4 ) element_num
    INTEGER ( kind = 4 ) element_order
    INTEGER ( kind = 4 ) m
    INTEGER ( kind = 4 ) node_num

    INTEGER ( kind = 4 ) element
    INTEGER ( kind = 4 ) element_node(element_order,element_num)
    INTEGER ( kind = 4 ) element_type
    CHARACTER * ( * ) gmsh_filename
    INTEGER ( kind = 4 ) gmsh_unit
    INTEGER ( kind = 4 ) node
    REAL ( kind = 8 ) node_x(m,node_num)
    INTEGER ( kind = 4 ) tag_num
    INTEGER ( kind = 4 ) tag1
    !
    !  Enforce 1-based indexing of nodes.
    !
    CALL mesh_base_one ( node_num, element_order, element_num, element_node )
    !
    !  Open the file.
    !
    CALL get_unit ( gmsh_unit )

    OPEN ( unit = gmsh_unit, file = gmsh_filename, status = 'replace' )
    !
    !  Write the data.
    !
    WRITE ( gmsh_unit, '(a)' ) '$MeshFormat'
    WRITE ( gmsh_unit, '(a)' ) '2.2 0 8'
    WRITE ( gmsh_unit, '(a)' ) '$EndMeshFormat'

    WRITE ( gmsh_unit, '(a)' ) '$Nodes'
    WRITE ( gmsh_unit, '(i6)' ) node_num
    DO node = 1, node_num
       WRITE ( gmsh_unit, '(i6,2x,g14.6,a)' ) &
            node, node_x(1:m,node), '  0.0  0.0'
    END DO
    WRITE ( gmsh_unit, '(a)' ) '$EndNodes'

    element_type = 1

    tag_num = 2
    tag1 = 0
    WRITE ( gmsh_unit, '(a)' ) '$Elements'
    WRITE ( gmsh_unit, '(i6)' ) element_num
    DO element = 1, element_num
       WRITE ( gmsh_unit, '(i6,2x,i2,2x,i2,2x,i2,2x,i6,2(2x,i6))' ) &
            element, element_type, tag_num, tag1, element, &
            element_node(1:element_order,element)
    END DO
    WRITE ( gmsh_unit, '(a)' ) '$EndElements'

    CLOSE ( unit = gmsh_unit )

    RETURN
  END SUBROUTINE gmsh_mesh1d_write

  SUBROUTINE gmsh_mesh2d_element_data_example ( element_num, element_order, &
     element_node )

    !*****************************************************************************80
    !
    !! GMSH_MESH2D_ELEMENT_DATA_EXAMPLE returns element data for the 2D example.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    10 October 2014
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    inputt, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
    !
    !    inputt, integer ( kind = 4 ) ELEMENT_ORDER, the order of the elements.
    !
    !    Output, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM),
    !    the indices of the nodes that make up each element.
    !


    INTEGER ( kind = 4 ) element_num
    INTEGER ( kind = 4 ) element_order

    INTEGER ( kind = 4 ) element_node(element_order,element_num)
    INTEGER ( kind = 4 ), DIMENSION ( 3, 24 ) :: element_node_save = &
         RESHAPE ( (/ &
         1,  2,  6, &
         7,  6,  2, &
         2,  3,  7, &
         8,  7,  3, &
         3,  4,  8, &
         9,  8,  4, &
         4,  5,  9, &
         10,  9,  5, &
         6,  7, 11, &
         12, 11,  7, &
         7,  8, 12, &
         13, 12,  8, &
         8,  9, 13, &
         14, 13,  9, &
         9, 10, 14, &
         15, 14, 10, &
         11, 12, 16, &
         17, 16, 12, &
         12, 13, 17, &
         18, 17, 13, &
         16, 17, 19, &
         20, 19, 17, &
         17, 18, 20, &
         21, 20, 18 /), (/ 3, 24 /) )

    CALL i4mat_copy ( element_order, element_num, element_node_save, &
         element_node )

    RETURN
  END SUBROUTINE gmsh_mesh2d_element_data_example

  SUBROUTINE gmsh_mesh2d_element_size_example ( element_num, element_order )

    !*****************************************************************************80
    !
    !! GMSH_MESH2D_ELEMENT_SIZE_EXAMPLE returns element sizes for the 2D example.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    10 October 2014
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
    !
    !    Output, integer ( kind = 4 ) ELEMENT_ORDER, the order of the elements.
    !


    INTEGER ( kind = 4 ) element_num
    INTEGER ( kind = 4 ) element_order

    element_num = 24
    element_order = 3

    RETURN
  END SUBROUTINE gmsh_mesh2d_element_size_example
  SUBROUTINE gmsh_mesh2d_node_data_example ( node_num, node_dim, node_x )

    !*****************************************************************************80
    !
    !! GMSH_MESH2D_NODE_DATA_EXAMPLE returns node data for the 2D example.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    10 October 2014
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    inputt, integer ( kind = 4 ) NODE_NUM, the number of nodes.
    !
    !    inputt, integer ( kind = 4 ) NODE_DIM, the spatial dimension.
    !
    !    Output, real ( kind = 8 ) NODE_X(NODE_DIM,NODE_NUM), the nodal
    !    coordinates.
    !


    INTEGER ( kind = 4 ) node_dim
    INTEGER ( kind = 4 ) node_num

    REAL ( kind = 8 ) node_x(node_dim,node_num)
    REAL ( kind = 8 ), DIMENSION ( 2, 21 ) :: node_x_save = RESHAPE ( (/ &
         0.0, 0.0, &
         1.0, 0.0, &
         2.0, 0.0, &
         3.0, 0.0, &
         4.0, 0.0, &
         0.0, 1.0, &
         1.0, 1.0, &
         2.0, 1.0, &
         3.0, 1.0, &
         4.0, 1.0, &
         0.0, 2.0, &
         1.0, 2.0, &
         2.0, 2.0, &
         3.0, 2.0, &
         4.0, 2.0, &
         0.0, 3.0, &
         1.0, 3.0, &
         2.0, 3.0, &
         0.0, 4.0, &
         1.0, 4.0, &
         2.0, 4.0 /), (/ 2, 21 /) )

    CALL r8mat_copy ( node_dim, node_num, node_x_save, node_x )

    RETURN
  END SUBROUTINE gmsh_mesh2d_node_data_example

  SUBROUTINE gmsh_mesh2d_node_size_example ( node_num, node_dim )

    !*****************************************************************************80
    !
    !! GMSH_MESH2D_NODE_SIZE_EXAMPLE returns node sizes for the 2D example.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    16 October 2014
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes.
    !
    !    Output, integer ( kind = 4 ) NODE_DIM, the spatial dimension.
    !


    INTEGER ( kind = 4 ) node_dim
    INTEGER ( kind = 4 ) node_num

    node_num = 21
    node_dim = 2

    RETURN
  END SUBROUTINE gmsh_mesh2d_node_size_example

  SUBROUTINE gmsh_mesh2d_write ( gmsh_filename, m, node_num, node_x, &
       element_order, element_num, element_node )

    !*****************************************************************************80
    !
    !! GMSH_MESH2D_WRITE writes 2D mesh data as a Gmsh file.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    06 October 2014
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Christophe Geuzaine, Jean-Francois Remacle,
    !    Gmsh: a three-dimensional finite element mesh generator with
    !    built-in pre- and post-processing facilities,
    !    International Journal for Numerical Methods in Engineering,
    !    Volume 79, Number 11, pages 1309-1331, 2009.
    !
    !  Parameters:
    !
    !    inputt, character * ( * ) GMSH_FILENAME, the name of the Gmsh file.
    !
    !    inputt, integer ( kind = 4 ) M, the spatial dimension.
    !
    !    inputt, integer ( kind = 4 ) NODE_NUM, the number of nodes.
    !
    !    inputt, real ( kind = 8 ) NODE_X(M,NODE_NUM), the node coordinates.
    !
    !    inputt, integer ( kind = 4 ) ELEMENT_ORDER, the order of the elements.
    !
    !    inputt, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
    !
    !    inputt, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM),
    !    the nodes that make up each element.
    !


    INTEGER ( kind = 4 ) element_num
    INTEGER ( kind = 4 ) element_order
    INTEGER ( kind = 4 ) m
    INTEGER ( kind = 4 ) node_num

    INTEGER ( kind = 4 ) element
    INTEGER ( kind = 4 ) element_node(element_order,element_num)
    INTEGER ( kind = 4 ) element_type
    CHARACTER * ( * ) gmsh_filename
    INTEGER ( kind = 4 ) gmsh_unit
    INTEGER ( kind = 4 ) node
    REAL ( kind = 8 ) node_x(m,node_num)
    INTEGER ( kind = 4 ) tag_num
    INTEGER ( kind = 4 ) tag1

    !
    !  Enforce 1-based indexing of nodes.
    !
    CALL mesh_base_one ( node_num, element_order, element_num, element_node )
    !
    !  Open the file.
    !
    CALL get_unit ( gmsh_unit )

    OPEN ( unit = gmsh_unit, file = gmsh_filename, status = 'replace' )
    !
    !  Write the data.
    !
    WRITE ( gmsh_unit, '(a)' ) '$MeshFormat'
    WRITE ( gmsh_unit, '(a)' ) '2.2 0 8'
    WRITE ( gmsh_unit, '(a)' ) '$EndMeshFormat'

    WRITE ( gmsh_unit, '(a)' ) '$Nodes'
    WRITE ( gmsh_unit, * ) node_num
    DO node = 1, node_num
       WRITE ( gmsh_unit, '(i6,2x,g14.6,2x,g14.6,a)' ) &
            node, node_x(1:m,node), '  0.0'
    END DO
    WRITE ( gmsh_unit, '(a)' ) '$EndNodes'

    IF ( element_order == 3 ) THEN
       element_type = 2
    ELSE IF ( element_order == 6 ) THEN
       element_type = 9
    ELSE IF ( element_order == 15 ) THEN
       element_type = 2
    END IF

    tag_num = 2
    tag1 = 0
    WRITE ( gmsh_unit, '(a)' ) '$Elements'
    WRITE ( gmsh_unit, * ) element_num
    DO element = 1, element_num
       !write ( gmsh_unit, '(i6,2x,i2,2x,i2,2x,i2,2x,i6,6(2x,i6))' ) &
       WRITE ( gmsh_unit, '(A,1x,i2,1x,i2,1x,i2,1x,i6,15(1x,i6))' ) &
                                  !write(gmsh_unit, '(5(I2, 1x), 4(I3, 1x))') element, element_type, tag_num, 1, 1, element_node(1:element_order, element)
            element, element_type, tag_num, tag1, 1, &
            element_node(1:element_order,element)
    END DO
    WRITE ( gmsh_unit, '(a)' ) '$EndElements'

    CLOSE ( unit = gmsh_unit )

    RETURN
  END SUBROUTINE gmsh_mesh2d_write

  SUBROUTINE gmsh_mesh3d_write ( gmsh_filename, m, node_num, node_x, &
       element_order, element_num, element_node )

    !*****************************************************************************80
    !
    !! GMSH_MESH3D_WRITE writes 3D mesh data as a Gmsh file.
    !
    !  Discussion:
    !
    !    The node ordering for the 20 node element is not standard.
    !
    !    Assuming the vertices are A, B, C and D, Gmsh uses the following ordering:
    !
    !    1:    a
    !    2:        b
    !    3:            c
    !    4:                d
    !    5: (2*a  +b        )/3
    !    6: (  a+2*b        )/3
    !    7: (    2*b+  c    )/3
    !    8: (      b+2*c    )/3
    !    9: (  a    +2*c    )/3
    !   10: (2*a    +  c    )/3
    !   11: (2*a        +  d)/3
    !   12: (  a        +2*d)/3
    !   13: (     b     +2*d)/3
    !   14: (   2*b     +  d)/3
    !   15: (       +  c+2*d)/3
    !   16: (       +2*c+  d)/3
    !   17: (  a+  b+  c    )/3
    !   18: (  a+  b    +  d)/3
    !   19: (      b+  c+  d)/3
    !   20: (  a+      c+  d)/3
    !
    !    Leo Rebholz used the following ordering:
    !
    !    1:    a
    !    2:        b
    !    3:            c
    !    4:                d
    !    5: (2*a  +b        )/3
    !    6: (2*a    +  c    )/3
    !    7: (  a+2*b        )/3
    !    8: (  a    +2*c    )/3
    !    9: (  a+  b+  c    )/3
    !   10: (    2*b+  c    )/3
    !   11: (      b+2*c    )/3
    !   12: (2*a        +  d)/3
    !   13: (   2*b     +  d)/3
    !   14: (       +2*c+  d)/3
    !   15: (  a+  b    +  d)/3
    !   16: (      b+  c+  d)/3
    !   17: (  a+      c+  d)/3
    !   18: (  a        +2*d)/3
    !   19: (     b     +2*d)/3
    !   20: (       +  c+2*d)/3
    !
    !    Since the only 20 node data we have is from Leo, we will assume that
    !    all 20 node inputt data is in Leo's format, and needs to be converted
    !    to the Gmsh convention.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    03 June 2013
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Christophe Geuzaine, Jean-Francois Remacle,
    !    Gmsh: a three-dimensional finite element mesh generator with
    !    built-in pre- and post-processing facilities,
    !    International Journal for Numerical Methods in Engineering,
    !    Volume 79, Number 11, pages 1309-1331, 2009.
    !
    !  Parameters:
    !
    !    inputt, character * ( * ) GMSH_FILENAME, the name of the Gmsh file.
    !
    !    inputt, integer ( kind = 4 ) M, the spatial dimension.
    !
    !    inputt, integer ( kind = 4 ) NODE_NUM, the number of nodes.
    !
    !    inputt, real ( kind = 8 ) NODE_X(M,NODE_NUM), the node coordinates.
    !
    !    inputt, integer ( kind = 4 ) ELEMENT_ORDER, the order of the elements.
    !
    !    inputt, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
    !
    !    inputt, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM),
    !    the nodes that make up each element.
    !


    INTEGER ( kind = 4 ) element_num
    INTEGER ( kind = 4 ) element_order
    INTEGER ( kind = 4 ) m
    INTEGER ( kind = 4 ) node_num

    INTEGER ( kind = 4 ) element
    INTEGER ( kind = 4 ) element_node(element_order,element_num)
    INTEGER ( kind = 4 ) element_type
    CHARACTER * ( * ) gmsh_filename
    INTEGER ( kind = 4 ) gmsh_unit
    INTEGER ( kind = 4 ), DIMENSION ( 20 ) :: leo_to_gmsh = (/ &
         1,  2,  3,  4,  5, &
         7, 10, 11,  8,  6, &
         12, 18, 19, 13, 20, &
         14,  9, 15, 16, 17 /)
    INTEGER ( kind = 4 ) node
    REAL ( kind = 8 ) node_x(3,node_num)
    INTEGER ( kind = 4 ) tag_num
    INTEGER ( kind = 4 ) tag1

    !
    !  Enforce 1-based indexing of nodes.
    !
    CALL mesh_base_one ( node_num, element_order, element_num, element_node )
    !
    !  Open the file.
    !
    CALL get_unit ( gmsh_unit )

    OPEN ( unit = gmsh_unit, file = gmsh_filename, status = 'replace' )
    !
    !  Write the data.
    !
    WRITE ( gmsh_unit, '(a)' ) '$MeshFormat'
    WRITE ( gmsh_unit, '(a)' ) '2.2 0 8'
    WRITE ( gmsh_unit, '(a)' ) '$EndMeshFormat'

    WRITE ( gmsh_unit, '(a)' ) '$Nodes'
    WRITE ( gmsh_unit, '(i6)' ) node_num
    DO node = 1, node_num
       WRITE ( gmsh_unit, '(i6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
            node, node_x(1:m,node)
    END DO
    WRITE ( gmsh_unit, '(a)' ) '$EndNodes'

    IF ( element_order == 4 ) THEN
       element_type = 4
    ELSE IF ( element_order == 10 ) THEN
       element_type = 11
    ELSE IF ( element_order == 20 ) THEN
       element_type = 29
    ELSE IF ( element_order == 15 ) THEN
       element_type = 2
    END IF

    tag_num = 2
    tag1 = 0
    WRITE ( gmsh_unit, '(a)' ) '$Elements'
    WRITE ( gmsh_unit, '(i6)' ) element_num
    DO element = 1, element_num
       IF ( element_order == 20 ) THEN
          WRITE ( gmsh_unit, '(i6,2x,i2,2x,i2,2x,i2,2x,i6,6(2x,i6))' ) &
               element, element_type, tag_num, tag1, element, &
               element_node(leo_to_gmsh(1:element_order),element)
       ELSE
          WRITE ( gmsh_unit, '(i6,2x,i2,2x,i2,2x,i2,2x,i6,6(2x,i6))' ) &
               element, element_type, tag_num, tag1, element, &
               element_node(1:element_order,element)
       END IF
    END DO
    WRITE ( gmsh_unit, '(a)' ) '$EndElements'

    CLOSE ( unit = gmsh_unit )

    RETURN
  END SUBROUTINE gmsh_mesh3d_write

  SUBROUTINE i4mat_copy ( m, n, a1, a2 )

    !*****************************************************************************80
    !
    !! I4MAT_COPY copies an I4MAT.
    !
    !  Discussion:
    !
    !    An I4MAT is a rectangular array of I4's.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    23 October 2010
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    inputt, integer ( kind = 4 ) M, N, the number of rows and columns.
    !
    !    inputt, integer ( kind = 4 ) A1(M,N), the matrix to be copied.
    !
    !    Output, integer ( kind = 4 ) A2(M,N), the copied matrix.
    !


    INTEGER ( kind = 4 ) m
    INTEGER ( kind = 4 ) n

    INTEGER ( kind = 4 ) a1(m,n)
    INTEGER ( kind = 4 ) a2(m,n)

    a2(1:m,1:n) = a1(1:m,1:n)

    RETURN
  END SUBROUTINE i4mat_copy
  SUBROUTINE i4mat_transpose_print ( m, n, a, title )

    !*****************************************************************************80
    !
    !! I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
    !
    !  Discussion:
    !
    !    An I4MAT is a rectangular array of I4's.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    28 December 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    inputt, integer ( kind = 4 ) M, N, the number of rows and columns.
    !
    !    inputt, integer ( kind = 4 ) A(M,N), an M by N matrix to be printed.
    !
    !    inputt, character ( len = * ) TITLE, a title.
    !


    INTEGER ( kind = 4 ) m
    INTEGER ( kind = 4 ) n

    INTEGER ( kind = 4 ) a(m,n)
    CHARACTER ( len = * ) title

    CALL i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

    RETURN
  END SUBROUTINE i4mat_transpose_print
  SUBROUTINE i4mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

    !*****************************************************************************80
    !
    !! I4MAT_TRANSPOSE_PRINT_SOME prints some of the transpose of an I4MAT.
    !
    !  Discussion:
    !
    !    An I4MAT is a rectangular array of I4's.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    10 September 2009
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    inputt, integer ( kind = 4 ) M, N, the number of rows and columns.
    !
    !    inputt, integer ( kind = 4 ) A(M,N), an M by N matrix to be printed.
    !
    !    inputt, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
    !
    !    inputt, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
    !
    !    inputt, character ( len = * ) TITLE, a title.
    !


    INTEGER ( kind = 4 ), PARAMETER :: incx = 10
    INTEGER ( kind = 4 ) m
    INTEGER ( kind = 4 ) n

    INTEGER ( kind = 4 ) a(m,n)
    CHARACTER ( len = 8 ) ctemp(incx)
    INTEGER ( kind = 4 ) i
    INTEGER ( kind = 4 ) i2
    INTEGER ( kind = 4 ) i2hi
    INTEGER ( kind = 4 ) i2lo
    INTEGER ( kind = 4 ) ihi
    INTEGER ( kind = 4 ) ilo
    INTEGER ( kind = 4 ) inc
    INTEGER ( kind = 4 ) j
    INTEGER ( kind = 4 ) j2hi
    INTEGER ( kind = 4 ) j2lo
    INTEGER ( kind = 4 ) jhi
    INTEGER ( kind = 4 ) jlo
    CHARACTER ( len = * ) title

    WRITE ( *, '(a)' ) ' '
    WRITE ( *, '(a)' ) TRIM ( title )

    IF ( m <= 0 .OR. n <= 0 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) '  (None)'
       RETURN
    END IF

    DO i2lo = MAX ( ilo, 1 ), MIN ( ihi, m ), incx

       i2hi = i2lo + incx - 1
       i2hi = MIN ( i2hi, m )
       i2hi = MIN ( i2hi, ihi )

       inc = i2hi + 1 - i2lo

       WRITE ( *, '(a)' ) ' '

       DO i = i2lo, i2hi
          i2 = i + 1 - i2lo
          WRITE ( ctemp(i2), '(i8)' ) i
       END DO

       WRITE ( *, '(''  Row '',10a8)' ) ctemp(1:inc)
       WRITE ( *, '(a)' ) '  Col'
       WRITE ( *, '(a)' ) ' '

       j2lo = MAX ( jlo, 1 )
       j2hi = MIN ( jhi, n )

       DO j = j2lo, j2hi

          DO i2 = 1, inc

             i = i2lo - 1 + i2

             WRITE ( ctemp(i2), '(i8)' ) a(i,j)

          END DO

          WRITE ( *, '(i5,a,10a8)' ) j, ':', ( ctemp(i), i = 1, inc )

       END DO

    END DO

    RETURN
  END SUBROUTINE i4mat_transpose_print_some

  SUBROUTINE mesh_base_one ( node_num, element_order, element_num, element_node )

    !*****************************************************************************80
    !
    !! MESH_BASE_ONE ensures that the element definition is one-based.
    !
    !  Discussion:
    !
    !    The ELEMENT_NODE array contains nodes indices that form elements.
    !    The convention for node indexing might start at 0 or at 1.
    !    Since a FORTRAN90 program will naturally assume a 1-based indexing, it is
    !    necessary to check a given element definition and, if it is actually
    !    0-based, to convert it.
    !
    !    This function attempts to detect 9-based node indexing and correct it.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    18 October 2014
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    inputt, int NODE_NUM, the number of nodes.
    !
    !    inputt, int ELEMENT_ORDER, the order of the elements.
    !
    !    inputt, int ELEMENT_NUM, the number of elements.
    !
    !    inputt/output, int ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM), the element
    !    definitions.
    !


    INTEGER ( kind = 4 ) element_num
    INTEGER ( kind = 4 ) element_order

    INTEGER ( kind = 4 ) element_node(element_order,element_num)
    INTEGER ( kind = 4 ), PARAMETER :: i4_huge = 2147483647
    INTEGER ( kind = 4 ) node_max
    INTEGER ( kind = 4 ) node_min
    INTEGER ( kind = 4 ) node_num

    node_min = + i4_huge
    node_max = - i4_huge

    node_min = MINVAL ( element_node(1:element_order,1:element_num) )
    node_max = MAXVAL ( element_node(1:element_order,1:element_num) )

    IF ( node_min == 0 .AND. node_max == node_num - 1 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' )'MESH_BASE_ONE:'
       WRITE ( *, '(a)' )'  The element indexing appears to be 0-based!'
       WRITE ( *, '(a)' )'  This will be converted to 1-based.'
       element_node(1:element_order,1:element_num) = &
            element_node(1:element_order,1:element_num) + 1
    ELSE IF ( node_min == 1 .AND. node_max == node_num  ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' )'MESH_BASE_ONE:'
       WRITE ( *, '(a)' )'  The element indexing appears to be 1-based!'
       WRITE ( *, '(a)' )'  No conversion is necessary.'
    ELSE
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'MESH_BASE_ONE - Warning!'
       WRITE ( *, '(a)' ) '  The element indexing is not of a recognized type.'
       WRITE ( *, '(a,i8)' ) '  NODE_MIN = ', node_min
       WRITE ( *, '(a,i8)' ) '  NODE_MAX = ', node_max
       WRITE ( *, '(a,i8)' ) '  NODE_NUM = ', node_num
    END IF

    RETURN
  END SUBROUTINE mesh_base_one

  SUBROUTINE r8mat_copy ( m, n, a, b )

    !*****************************************************************************80
    !
    !! R8MAT_COPY copies an R8MAT.
    !
    !  Discussion:
    !
    !    An R8MAT is an MxN array of R8's, stored by (I,J) -> (I+J*M).
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    26 July 2008
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    inputt, integer ( kind = 4 ) M, N, the order of the matrix.
    !
    !    inputt, real ( kind = 8 ) A(M,N), the matrix to be copied.
    !
    !    Output, real ( kind = 8 ) B(M,N), a copy of the matrix.
    !


    INTEGER ( kind = 4 ) m
    INTEGER ( kind = 4 ) n

    REAL ( kind = 8 ) a(m,n)
    REAL ( kind = 8 ) b(m,n)

    b(1:m,1:n) = a(1:m,1:n)

    RETURN
  END SUBROUTINE r8mat_copy

  SUBROUTINE r8mat_transpose_print ( m, n, a, title )

    !*****************************************************************************80
    !
    !! R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
    !
    !  Discussion:
    !
    !    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    14 June 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    inputt, integer ( kind = 4 ) M, N, the number of rows and columns.
    !
    !    inputt, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
    !
    !    inputt, character ( len = * ) TITLE, a title.
    !


    INTEGER ( kind = 4 ) m
    INTEGER ( kind = 4 ) n

    REAL ( kind = 8 ) a(m,n)
    CHARACTER ( len = * ) title

    CALL r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

    RETURN
  END SUBROUTINE r8mat_transpose_print

  SUBROUTINE r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

    !*****************************************************************************80
    !
    !! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
    !
    !  Discussion:
    !
    !    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    10 September 2009
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    inputt, integer ( kind = 4 ) M, N, the number of rows and columns.
    !
    !    inputt, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
    !
    !    inputt, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
    !
    !    inputt, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
    !
    !    inputt, character ( len = * ) TITLE, a title.
    !


    INTEGER ( kind = 4 ), PARAMETER :: incx = 5
    INTEGER ( kind = 4 ) m
    INTEGER ( kind = 4 ) n

    REAL ( kind = 8 ) a(m,n)
    CHARACTER ( len = 14 ) ctemp(incx)
    INTEGER ( kind = 4 ) i
    INTEGER ( kind = 4 ) i2
    INTEGER ( kind = 4 ) i2hi
    INTEGER ( kind = 4 ) i2lo
    INTEGER ( kind = 4 ) ihi
    INTEGER ( kind = 4 ) ilo
    INTEGER ( kind = 4 ) inc
    INTEGER ( kind = 4 ) j
    INTEGER ( kind = 4 ) j2hi
    INTEGER ( kind = 4 ) j2lo
    INTEGER ( kind = 4 ) jhi
    INTEGER ( kind = 4 ) jlo
    CHARACTER ( len = * ) title

    WRITE ( *, '(a)' ) ' '
    WRITE ( *, '(a)' ) TRIM ( title )

    IF ( m <= 0 .OR. n <= 0 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) '  (None)'
       RETURN
    END IF

    DO i2lo = MAX ( ilo, 1 ), MIN ( ihi, m ), incx

       i2hi = i2lo + incx - 1
       i2hi = MIN ( i2hi, m )
       i2hi = MIN ( i2hi, ihi )

       inc = i2hi + 1 - i2lo

       WRITE ( *, '(a)' ) ' '

       DO i = i2lo, i2hi
          i2 = i + 1 - i2lo
          WRITE ( ctemp(i2), '(i8,6x)' ) i
       END DO

       WRITE ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
       WRITE ( *, '(a)' ) '  Col'
       WRITE ( *, '(a)' ) ' '

       j2lo = MAX ( jlo, 1 )
       j2hi = MIN ( jhi, n )

       DO j = j2lo, j2hi

          DO i2 = 1, inc
             i = i2lo - 1 + i2
             WRITE ( ctemp(i2), '(g14.6)' ) a(i,j)
          END DO

          WRITE ( *, '(i5,a,5a14)' ) j, ':', ( ctemp(i), i = 1, inc )

       END DO

    END DO

    RETURN
  END SUBROUTINE r8mat_transpose_print_some

  SUBROUTINE s_begin ( s1, s2, flag )

    !*****************************************************************************80
    !
    !! S_BEGIN is TRUE if one string matches the beginning of the other.
    !
    !  Discussion:
    !
    !    The strings are compared, ignoring blanks, spaces and capitalization.
    !
    !  Example:
    !
    !     S1              S2      S_BEGIN
    !
    !    'Bob'          'BOB'     TRUE
    !    '  B  o b '    ' bo b'   TRUE
    !    'Bob'          'Bobby'   TRUE
    !    'Bobo'         'Bobb'    FALSE
    !    ' '            'Bob'     FALSE    (Do not allow a blank to match
    !                                       anything but another blank string.)
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    20 January 2000
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    inputt, character ( len = * ) S1, S2, the strings to be compared.
    !
    !    Output, logical ( kind = 4 ) S_BEGIN, is TRUE if the strings match up to
    !    the end of the shorter string, ignoring case.
    !


    LOGICAL ( kind = 4 ) flag_ch_eqi
    INTEGER ( kind = 4 ) i1
    INTEGER ( kind = 4 ) i2
    LOGICAL ( kind = 4 ), INTENT(OUT) :: flag
    CHARACTER ( len = * ) s1
    INTEGER ( kind = 4 ) s1_length
    CHARACTER ( len = * ) s2
    INTEGER ( kind = 4 ) s2_length

    s1_length = len_TRIM ( s1 )
    s2_length = len_TRIM ( s2 )
    !
    !  If either string is blank, then both must be blank to match.
    !  Otherwise, a blank string matches anything, which is not
    !  what most people want.
    !
    IF ( s1_length == 0 .OR. s2_length == 0 ) THEN

       IF ( s1_length == 0 .AND. s2_length == 0 ) THEN
          flag = .TRUE.
       ELSE
          flag = .FALSE.
       END IF

       RETURN

    END IF

    i1 = 0
    i2 = 0
    !
    !  Find the next nonblank in S1.
    !
    DO

       DO

          i1 = i1 + 1

          IF ( s1_length < i1 ) THEN
             flag = .TRUE.
             RETURN
          END IF

          IF ( s1(i1:i1) /= ' ' ) THEN
             EXIT
          END IF

       END DO
       !
       !  Find the next nonblank in S2.
       !
       DO

          i2 = i2 + 1

          IF ( s2_length < i2 ) THEN
             flag = .TRUE.
             RETURN
          END IF

          IF ( s2(i2:i2) /= ' ' ) THEN
             EXIT
          END IF

       END DO
       !
       !  If the characters match, get the next pair.
       !
       CALL ch_eqi(s1(i1:i1), s2(i2:i2), flag_ch_eqi)
       IF ( .NOT. flag_ch_eqi) THEN
          EXIT
       END IF

    END DO

    flag = .FALSE.

    RETURN
  END SUBROUTINE s_begin

  SUBROUTINE ch_eqi ( c1, c2, flag )

    !*****************************************************************************80
    !
    !! CH_EQI is a case insensitive comparison of two characters for equality.
    !
    !  Example:
    !
    !    CH_EQI ( 'A', 'a' ) is .TRUE.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    28 July 2000
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    inputt, character C1, C2, the characters to compare.
    !
    !    Output, logical ( kind = 4 ) CH_EQI, the result of the comparison.
    !


    LOGICAL ( kind = 4 ), INTENT(OUT) :: flag
    CHARACTER c1
    CHARACTER c1_cap
    CHARACTER c2
    CHARACTER c2_cap

    c1_cap = c1
    c2_cap = c2

    CALL ch_cap ( c1_cap )
    CALL ch_cap ( c2_cap )

    IF ( c1_cap == c2_cap ) THEN
       flag = .TRUE.
    ELSE
       flag = .FALSE.
    END IF

    RETURN
  END SUBROUTINE ch_eqi

  SUBROUTINE s_to_i4 ( s, ival, ierror, length )

    !*****************************************************************************80
    !
    !! S_TO_I4 reads an I4 from a string.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    28 June 2000
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    inputt, character ( len = * ) S, a string to be examined.
    !
    !    Output, integer ( kind = 4 ) IVAL, the integer value read from the string.
    !    If the string is blank, then IVAL will be returned 0.
    !
    !    Output, integer ( kind = 4 ) IERROR, an error flag.
    !    0, no error.
    !    1, an error occurred.
    !
    !    Output, integer ( kind = 4 ) LENGTH, the number of characters of S
    !    used to make IVAL.
    !


    CHARACTER c
    INTEGER ( kind = 4 ) i
    INTEGER ( kind = 4 ) ierror
    INTEGER ( kind = 4 ) isgn
    INTEGER ( kind = 4 ) istate
    INTEGER ( kind = 4 ) ival
    INTEGER ( kind = 4 ) length
    CHARACTER ( len = * ) s

    ierror = 0
    istate = 0
    isgn = 1
    ival = 0

    DO i = 1, len_TRIM ( s )

       c = s(i:i)
       !
       !  Haven't read anything.
       !
       IF ( istate == 0 ) THEN

          IF ( c == ' ' ) THEN

          ELSE IF ( c == '-' ) THEN
             istate = 1
             isgn = -1
          ELSE IF ( c == '+' ) THEN
             istate = 1
             isgn = + 1
          ELSE IF ( LLE ( '0', c ) .AND. LLE ( c, '9' ) ) THEN
             istate = 2
             ival = ICHAR ( c ) - ICHAR ( '0' )
          ELSE
             ierror = 1
             RETURN
          END IF
          !
          !  Have read the sign, expecting digits.
          !
       ELSE IF ( istate == 1 ) THEN

          IF ( c == ' ' ) THEN

          ELSE IF ( LLE ( '0', c ) .AND. LLE ( c, '9' ) ) THEN
             istate = 2
             ival = ICHAR ( c ) - ICHAR ( '0' )
          ELSE
             ierror = 1
             RETURN
          END IF
          !
          !  Have read at least one digit, expecting more.
          !
       ELSE IF ( istate == 2 ) THEN

          IF ( LLE ( '0', c ) .AND. LLE ( c, '9' ) ) THEN
             ival = 10 * ival + ICHAR ( c ) - ICHAR ( '0' )
          ELSE
             ival = isgn * ival
             length = i - 1
             RETURN
          END IF

       END IF

    END DO
    !
    !  If we read all the characters in the string, see if we're OK.
    !
    IF ( istate == 2 ) THEN
       ival = isgn * ival
       length = len_TRIM ( s )
    ELSE
       ierror = 1
       length = 0
    END IF

    RETURN
  END SUBROUTINE s_to_i4

  SUBROUTINE s_to_r8 ( s, dval, ierror, length )

    !*****************************************************************************80
    !
    !! S_TO_R8 reads an R8 from a string.
    !
    !  Discussion:
    !
    !    The routine will read as many characters as possible until it reaches
    !    the end of the string, or encounters a character which cannot be
    !    part of the number.
    !
    !    Legal inputt is:
    !
    !       1 blanks,
    !       2 '+' or '-' sign,
    !       2.5 blanks
    !       3 integer part,
    !       4 decimal point,
    !       5 fraction part,
    !       6 'E' or 'e' or 'D' or 'd', exponent marker,
    !       7 exponent sign,
    !       8 exponent integer part,
    !       9 exponent decimal point,
    !      10 exponent fraction part,
    !      11 blanks,
    !      12 final comma or semicolon,
    !
    !    with most quantities optional.
    !
    !  Example:
    !
    !    S                 DVAL
    !
    !    '1'               1.0
    !    '     1   '       1.0
    !    '1A'              1.0
    !    '12,34,56'        12.0
    !    '  34 7'          34.0
    !    '-1E2ABCD'        -100.0
    !    '-1X2ABCD'        -1.0
    !    ' 2E-1'           0.2
    !    '23.45'           23.45
    !    '-4.2E+2'         -420.0
    !    '17d2'            1700.0
    !    '-14e-2'         -0.14
    !    'e2'              100.0
    !    '-12.73e-9.23'   -12.73 * 10.0^(-9.23)
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    07 September 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    inputt, character ( len = * ) S, the string containing the
    !    data to be read.  Reading will begin at position 1 and
    !    terminate at the end of the string, or when no more
    !    characters can be read to form a legal real.  Blanks,
    !    commas, or other nonnumeric data will, in particular,
    !    cause the conversion to halt.
    !
    !    Output, real ( kind = 8 ) DVAL, the value read from the string.
    !
    !    Output, integer ( kind = 4 ) IERROR, error flag.
    !    0, no errors occurred.
    !    1, 2, 6 or 7, the inputt number was garbled.  The
    !    value of IERROR is the last type of inputt successfully
    !    read.  For instance, 1 means initial blanks, 2 means
    !    a plus or minus sign, and so on.
    !
    !    Output, integer ( kind = 4 ) LENGTH, the number of characters read
    !    to form the number, including any terminating
    !    characters such as a trailing comma or blanks.
    !

    CHARACTER c
    LOGICAL ( kind = 4 ) flag1, flag2
    REAL ( kind = 8 ) dval
    INTEGER ( kind = 4 ) ierror
    INTEGER ( kind = 4 ) ihave
    INTEGER ( kind = 4 ) isgn
    INTEGER ( kind = 4 ) iterm
    INTEGER ( kind = 4 ) jbot
    INTEGER ( kind = 4 ) jsgn
    INTEGER ( kind = 4 ) jtop
    INTEGER ( kind = 4 ) length
    INTEGER ( kind = 4 ) nchar
    INTEGER ( kind = 4 ) ndig
    REAL ( kind = 8 ) rbot
    REAL ( kind = 8 ) rexp
    REAL ( kind = 8 ) rtop
    CHARACTER ( len = * ) s

    nchar = len_TRIM ( s )

    ierror = 0
    dval = 0.0D+00
    length = -1
    isgn = 1
    rtop = 0
    rbot = 1
    jsgn = 1
    jtop = 0
    jbot = 1
    ihave = 1
    iterm = 0

    DO

       length = length + 1

       IF ( nchar < length+1 ) THEN
          EXIT
       END IF

       c = s(length+1:length+1)
       !
       !  Blank character.
       !
       CALL ch_eqi ( c, 'E' , flag1)
       CALL ch_eqi ( c, 'D' , flag2)
       IF ( c == ' ' ) THEN

          IF ( ihave == 2 ) THEN

          ELSE IF ( ihave == 6 .OR. ihave == 7 ) THEN
             iterm = 1
          ELSE IF ( 1 < ihave ) THEN
             ihave = 11
          END IF
          !
          !  Comma.
          !
       ELSE IF ( c == ',' .OR. c == ';' ) THEN

          IF ( ihave /= 1 ) THEN
             iterm = 1
             ihave = 12
             length = length + 1
          END IF
          !
          !  Minus sign.
          !
       ELSE IF ( c == '-' ) THEN

          IF ( ihave == 1 ) THEN
             ihave = 2
             isgn = -1
          ELSE IF ( ihave == 6 ) THEN
             ihave = 7
             jsgn = -1
          ELSE
             iterm = 1
          END IF
          !
          !  Plus sign.
          !
       ELSE IF ( c == '+' ) THEN

          IF ( ihave == 1 ) THEN
             ihave = 2
          ELSE IF ( ihave == 6 ) THEN
             ihave = 7
          ELSE
             iterm = 1
          END IF
          !
          !  Decimal point.
          !
       ELSE IF ( c == '.' ) THEN

          IF ( ihave < 4 ) THEN
             ihave = 4
          ELSE IF ( 6 <= ihave .AND. ihave <= 8 ) THEN
             ihave = 9
          ELSE
             iterm = 1
          END IF
          !
          !  Scientific notation exponent marker.
       ELSE IF ( flag1 .OR. flag2 ) THEN

          IF ( ihave < 6 ) THEN
             ihave = 6
          ELSE
             iterm = 1
          END IF
          !
          !  Digit.
          !
       ELSE IF (  ihave < 11 .AND. LLE ( '0', c ) .AND. LLE ( c, '9' ) ) THEN

          IF ( ihave <= 2 ) THEN
             ihave = 3
          ELSE IF ( ihave == 4 ) THEN
             ihave = 5
          ELSE IF ( ihave == 6 .OR. ihave == 7 ) THEN
             ihave = 8
          ELSE IF ( ihave == 9 ) THEN
             ihave = 10
          END IF

          CALL ch_to_digit ( c, ndig )

          IF ( ihave == 3 ) THEN
             rtop = 10.0D+00 * rtop + REAL ( ndig, kind = 8 )
          ELSE IF ( ihave == 5 ) THEN
             rtop = 10.0D+00 * rtop + REAL ( ndig, kind = 8 )
             rbot = 10.0D+00 * rbot
          ELSE IF ( ihave == 8 ) THEN
             jtop = 10 * jtop + ndig
          ELSE IF ( ihave == 10 ) THEN
             jtop = 10 * jtop + ndig
             jbot = 10 * jbot
          END IF
          !
          !  Anything else is regarded as a terminator.
          !
       ELSE
          iterm = 1
       END IF
       !
       !  If we haven't seen a terminator, and we haven't examined the
       !  entire string, go get the next character.
       !
       IF ( iterm == 1 ) THEN
          EXIT
       END IF

    END DO
    !
    !  If we haven't seen a terminator, and we have examined the
    !  entire string, then we're done, and LENGTH is equal to NCHAR.
    !
    IF ( iterm /= 1 .AND. length + 1 == nchar ) THEN
       length = nchar
    END IF
    !
    !  Number seems to have terminated.  Have we got a legal number?
    !  Not if we terminated in states 1, 2, 6 or 7!
    !
    IF ( ihave == 1 .OR. ihave == 2 .OR. ihave == 6 .OR. ihave == 7 ) THEN
       ierror = ihave
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'S_TO_R8 - Serious error!'
       WRITE ( *, '(a)' ) '  Illegal or nonnumeric inputt:'
       WRITE ( *, '(a)' ) '    ' // TRIM ( s )
       RETURN
    END IF
    !
    !  Number seems OK.  Form it.
    !
    IF ( jtop == 0 ) THEN
       rexp = 1.0D+00
    ELSE
       IF ( jbot == 1 ) THEN
          rexp = 10.0D+00 ** ( jsgn * jtop )
       ELSE
          rexp = 10.0D+00 ** ( REAL ( jsgn * jtop, kind = 8 ) &
               / REAL ( jbot, kind = 8 ) )
       END IF
    END IF

    dval = REAL ( isgn, kind = 8 ) * rexp * rtop / rbot

    RETURN
  END SUBROUTINE s_to_r8

  SUBROUTINE timestamp ( )

    !*****************************************************************************80
    !
    !! TIMESTAMP prints the current YMDHMS date as a time stamp.
    !
    !  Example:
    !
    !    31 May 2001   9:45:54.872 AM
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    18 May 2013
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    None
    !


    CHARACTER ( len = 8 ) ampm
    INTEGER ( kind = 4 ) d
    INTEGER ( kind = 4 ) h
    INTEGER ( kind = 4 ) m
    INTEGER ( kind = 4 ) mm
    CHARACTER ( len = 9 ), PARAMETER, DIMENSION(12) :: month = (/ &
         'January  ', 'February ', 'March    ', 'April    ', &
         'May      ', 'June     ', 'July     ', 'August   ', &
         'September', 'October  ', 'November ', 'December ' /)
    INTEGER ( kind = 4 ) n
    INTEGER ( kind = 4 ) s
    INTEGER ( kind = 4 ) values(8)
    INTEGER ( kind = 4 ) y

    CALL date_and_TIME ( values = values )

    y = values(1)
    m = values(2)
    d = values(3)
    h = values(5)
    n = values(6)
    s = values(7)
    mm = values(8)

    IF ( h < 12 ) THEN
       ampm = 'AM'
    ELSE IF ( h == 12 ) THEN
       IF ( n == 0 .AND. s == 0 ) THEN
          ampm = 'Noon'
       ELSE
          ampm = 'PM'
       END IF
    ELSE
       h = h - 12
       IF ( h < 12 ) THEN
          ampm = 'PM'
       ELSE IF ( h == 12 ) THEN
          IF ( n == 0 .AND. s == 0 ) THEN
             ampm = 'Midnight'
          ELSE
             ampm = 'AM'
          END IF
       END IF
    END IF

    WRITE ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
         d, TRIM ( month(m) ), y, h, ':', n, ':', s, '.', mm, TRIM ( ampm )

    RETURN
  END SUBROUTINE timestamp

  SUBROUTINE ch_cap ( ch )

    !*****************************************************************************80
    !
    !! CH_CAP capitalizes a single character.
    !
    !  Discussion:
    !
    !    Instead of CHAR and ICHAR, we now use the ACHAR and IACHAR functions,
    !    which guarantee the ASCII collating sequence.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    19 July 1998
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    inputt/output, character CH, the character to capitalize.


    CHARACTER ch
    INTEGER ( kind = 4 ) itemp

    itemp = IACHAR ( ch )

    IF ( 97 <= itemp .AND. itemp <= 122 ) THEN
       ch = ACHAR ( itemp - 32 )
    END IF

    RETURN
  END SUBROUTINE ch_cap

  SUBROUTINE ch_to_digit ( c, digit )

    !*****************************************************************************80
    !
    !! CH_TO_DIGIT returns the integer value of a base 10 digit.
    !
    !  Example:
    !
    !     C   DIGIT
    !    ---  -----
    !    '0'    0
    !    '1'    1
    !    ...  ...
    !    '9'    9
    !    ' '    0
    !    'X'   -1
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    04 August 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    inputt, character C, the decimal digit, '0' through '9' or blank
    !    are legal.
    !
    !    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.
    !    If C was 'illegal', then DIGIT is -1.
    !


    CHARACTER c
    INTEGER ( kind = 4 ) digit

    IF ( LGE ( c, '0' ) .AND. LLE ( c, '9' ) ) THEN

       digit = ICHAR ( c ) - 48

    ELSE IF ( c == ' ' ) THEN

       digit = 0

    ELSE

       digit = -1

    END IF

    RETURN
  END SUBROUTINE ch_to_digit

  SUBROUTINE get_unit ( iunit )

    !*****************************************************************************80
    !
    !! GET_UNIT returns a free FORTRAN unit number.
    !
    !  Discussion:
    !
    !    A "free" FORTRAN unit number is an integer between 1 and 99 which
    !    is not currently associated with an I/O device.  A free FORTRAN unit
    !    number is needed in order to open a file with the OPEN command.
    !
    !    If IUNIT = 0, then no free FORTRAN unit could be found, although
    !    all 99 units were checked (except for units 5, 6 and 9, which
    !    are commonly reserved for console I/O).
    !
    !    Otherwise, IUNIT is an integer between 1 and 99, representing a
    !    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
    !    are special, and will never return those values.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    18 September 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Output, integer ( kind = 4 ) IUNIT, the free unit number.
    !


    INTEGER ( kind = 4 ) i
    INTEGER ( kind = 4 ) ios
    INTEGER ( kind = 4 ) iunit
    LOGICAL lopen

    iunit = 0

    DO i = 1, 99

       IF ( i /= 5 .AND. i /= 6 .AND. i /= 9 ) THEN

          INQUIRE ( unit = i, opened = lopen, iostat = ios )

          IF ( ios == 0 ) THEN
             IF ( .NOT. lopen ) THEN
                iunit = i
                RETURN
             END IF
          END IF

       END IF

    END DO

    RETURN
  END SUBROUTINE get_unit

  SUBROUTINE generate_splines_from_geo_file(filename)
    USE mod_splines
    USE globals
      IMPLICIT NONE
      CHARACTER (*), INTENT(IN)       :: filename

      CHARACTER(len=1000)             :: data_str, line
      INTEGER                         :: i, j, k, ios, n_splines, n_splines_phys, n_points, index1, index2, point
      INTEGER, ALLOCATABLE            :: phys_curve_IN(:), phys_curve_OUT(:), phys_curve_LIM(:), phys_curve_PUMP(:), phys_curve_PUFF(:)
      REAL*8                          :: xcenter, ycenter, x1, x2, y1, y2, angle, RotMat(2,2)
      INTEGER, PARAMETER              :: read_unit = 99


      ! Open the file
      OPEN(unit=read_unit, file=filename, iostat=ios)
      IF ( ios /= 0 ) STOP "Error opening file .geo"

      n_splines = 0
      DO
        READ(read_unit, '(A)', iostat=ios) line
        ! Count the number of splines
        IF(line(:7) .eq. 'BSpline') THEN
          n_splines = n_splines + 1
        ENDIF

        IF (ios /= 0) EXIT
      ENDDO

      ! Go back to top of file
      REWIND (unit = read_unit)

      ALLOCATE(splines(n_splines))

      ! Read the points the each spline contains
      i = 1
      DO
        READ(read_unit, '(A)', iostat=ios) line
        ! Look for the Splines
        IF(line(:7) .eq. 'BSpline') THEN

          ! extract the spline number between the ()
          data_str = line(index(line,'(')+1:index(line,')')-1)

          READ(data_str, *) splines(i)%spline_number

          ! extract the numbers between the {}
          data_str = line(index(line,'{')+1:index(line,'}')-1)

          ! count the number of points
          ! I can actually use the same iterator i cause it's temporary in COUNT, for some reasons i cannot use other iterators other than i
          n_points = count( (/ (data_str(i:i), i=1, len_trim(data_str)) /) .eq. ",") + 1

          splines(i)%n_points = n_points
          ALLOCATE(splines(i)%points_number(n_points))

          ! read the points
          READ(data_str, *) splines(i)%points_number(:)
          i = i + 1
        ENDIF
        IF (ios .ne. 0) EXIT
      ENDDO

      DO i = 1, n_splines
        ALLOCATE(splines(i)%points_coord(splines(i)%n_points, 2))
      ENDDO

      DO i = 1, n_splines
        ! Go back to top of file
        REWIND (unit = read_unit)
        ! loop through file
        DO
          READ(read_unit, '(A)', iostat=ios) line
          ! Look for the points
          IF(line(:5) .eq. 'Point') THEN

            ! extract the numbers between the ()
            data_str = line(index(line,'(')+1:index(line,')')-1)
            ! convert to integer
            READ(data_str, *) point

            ! loop through spline points
            DO j = 1, splines(i)%n_points
              ! if the point belongs to a spline
              IF(splines(i)%points_number(j) .eq. point) THEN
                ! starting string index is given by {
                index1 = index(line,'{') + 1
                index2 = index1

                ! ending string index is given by the second comma
                k = 0
                DO
                  IF(line(index2:index2) .eq. ',') THEN
                    k = k + 1
                  ENDIF
                  IF(k .eq. 2) THEN
                    EXIT
                  ENDIF
                  index2 = index2 + 1
                ENDDO
                READ(line(index1:index2-1), *) splines(i)%points_coord(j,:)
              ENDIF
            ENDDO
          ENDIF
          ! exit when reached end of file
          IF (ios .ne. 0) EXIT
        ENDDO
      ENDDO

      ! Go back to top of file
      REWIND (unit = read_unit)

      ! set boundary type
      DO
        READ(read_unit, '(A)', iostat=ios) line
        ! Count the number of splines
        IF(line(:21) .eq. 'Physical Curve("PUFF"') THEN
          ! extract the numbers between the {}
          data_str = line(index(line,'{')+1:index(line,'}')-1)
          ! count how many splines make up the physical boundary IN
          n_splines_phys = count( (/ (data_str(i:i), i=1, len_trim(data_str)) /) .eq. ",") + 1

          ALLOCATE(phys_curve_PUFF(n_splines_phys))
          ! read the points
          READ(data_str, *) phys_curve_PUFF
        ENDIF
        IF(line(:21) .eq. 'Physical Curve("PUMP"') THEN
          ! extract the numbers between the {}
          data_str = line(index(line,'{')+1:index(line,'}')-1)
          ! count how many splines make up the physical boundary IN
          n_splines_phys = count( (/ (data_str(i:i), i=1, len_trim(data_str)) /) .eq. ",") + 1

          ALLOCATE(phys_curve_PUMP(n_splines_phys))
          ! read the points
          READ(data_str, *) phys_curve_PUMP
        ENDIF
        IF(line(:19) .eq. 'Physical Curve("IN"') THEN
          ! extract the numbers between the {}
          data_str = line(index(line,'{')+1:index(line,'}')-1)
          ! count how many splines make up the physical boundary IN
          n_splines_phys = count( (/ (data_str(i:i), i=1, len_trim(data_str)) /) .eq. ",") + 1

          ALLOCATE(phys_curve_IN(n_splines_phys))
          ! read the points
          READ(data_str, *) phys_curve_IN
        ENDIF
        IF(line(:20) .eq. 'Physical Curve("LIM"') THEN
          ! extract the numbers between the {}
          data_str = line(index(line,'{')+1:index(line,'}')-1)
          ! count how many splines make up the physical boundary IN
          n_splines_phys = count( (/ (data_str(i:i), i=1, len_trim(data_str)) /) .eq. ",") + 1

          ALLOCATE(phys_curve_LIM(n_splines_phys))
          ! read the points
          READ(data_str, *) phys_curve_LIM
        ENDIF

        IF(line(:20) .eq. 'Physical Curve("OUT"') THEN
          ! extract the numbers between the {}
          data_str = line(index(line,'{')+1:index(line,'}')-1)
          ! count how many splines make up the physical boundary IN
          n_splines_phys = count( (/ (data_str(i:i), i=1, len_trim(data_str)) /) .eq. ",") + 1

          ALLOCATE(phys_curve_OUT(n_splines_phys))
          ! read the points
          READ(data_str, *) phys_curve_OUT
        ENDIF

        IF (ios /= 0) EXIT
      ENDDO

      DO i = 1, n_splines
        IF(ALLOCATED(phys_curve_PUMP)) THEN
          IF(ANY(phys_curve_PUMP .eq. splines(i)%spline_number)) THEN
            splines(i)%boundaryFlag = 5
          ENDIF
        ENDIF
        IF(ALLOCATED(phys_curve_PUFF)) THEN
          IF(ANY(phys_curve_PUFF .eq. splines(i)%spline_number)) THEN
            splines(i)%boundaryFlag = 6
          ENDIF
        ENDIF
        IF(ALLOCATED(phys_curve_LIM)) THEN
          IF(ANY(phys_curve_LIM .eq. splines(i)%spline_number)) THEN
            splines(i)%boundaryFlag = 7
          ENDIF
        ENDIF

        IF(ALLOCATED(phys_curve_IN)) THEN
          IF(ANY(phys_curve_IN .eq. splines(i)%spline_number)) THEN
            splines(i)%boundaryFlag = 8
          ENDIF
        ENDIF

        IF(ALLOCATED(phys_curve_OUT)) THEN
          IF(ANY(phys_curve_OUT .eq. splines(i)%spline_number)) THEN
            splines(i)%boundaryFlag = 9
          ENDIF
        ENDIF
      ENDDO

      IF(ALLOCATED(phys_curve_LIM)) DEALLOCATE(phys_curve_LIM)
      IF(ALLOCATED(phys_curve_IN)) DEALLOCATE(phys_curve_IN)
      IF(ALLOCATED(phys_curve_OUT)) DEALLOCATE(phys_curve_OUT)
      IF(ALLOCATED(phys_curve_PUFF)) DEALLOCATE(phys_curve_PUFF)
      IF(ALLOCATED(phys_curve_PUMP)) DEALLOCATE(phys_curve_PUMP)


      IF((switch%testcase .ge. 60) .and. (switch%testcase .le. 80)) THEN
        DO i = 1, n_splines
          splines(i)%points_coord(:,1) = splines(i)%points_coord(:,1) + geom%R0
        ENDDO
      ENDIF


      DO i = 1, n_splines
        splines(i)%xmax = MAXVAL(splines(i)%points_coord(:,1))
        splines(i)%xmin = MINVAL(splines(i)%points_coord(:,1))
        splines(i)%ymax = MAXVAL(splines(i)%points_coord(:,2))
        splines(i)%ymin = MINVAL(splines(i)%points_coord(:,2))
        splines(i)%tot_n_splines = n_splines
      ENDDO

      DO i = 1, n_splines
        ! rotate the points so that the straight boundary aligns to the x axis
        xcenter = SUM(splines(i)%points_coord(:,1))/REAL(splines(i)%n_points)
        ycenter = SUM(splines(i)%points_coord(:,2))/REAL(splines(i)%n_points)

        splines(i)%xcenter = xcenter
        splines(i)%ycenter = ycenter


        splines(i)%points_coord(:,1) = splines(i)%points_coord(:,1) - xcenter
        splines(i)%points_coord(:,2) = splines(i)%points_coord(:,2) - ycenter

        y1 = splines(i)%points_coord(1,2)
        y2 = splines(i)%points_coord(splines(i)%n_points,2)
        x1 = splines(i)%points_coord(1,1)
        x2 = splines(i)%points_coord(splines(i)%n_points,1)

        angle = - ATAN(y2 - y1, x2 - x1)
        splines(i)%angle = angle

        RotMat(1,1) = COS(angle)
        RotMat(1,2) = -SIN(angle)
        RotMat(2,1) = SIN(angle)
        RotMat(2,2) = COS(angle)

        splines(i)%RotMat = RotMat

        splines(i)%points_coord = TRANSPOSE(MATMUL(RotMat,TRANSPOSE(splines(i)%points_coord)))

        splines(i)%points_coord(:,1) = splines(i)%points_coord(:,1) + xcenter
        splines(i)%points_coord(:,2) = splines(i)%points_coord(:,2) + ycenter

        splines(i)%spline = spline(splines(i)%points_coord(:,1), splines(i)%points_coord(:,2))

        RotMat(1,1) = COS(-angle)
        RotMat(1,2) = -SIN(-angle)
        RotMat(2,1) = SIN(-angle)
        RotMat(2,2) = COS(-angle)

        splines(i)%AntiRotMat = RotMat

        splines(i)%points_coord(:,1) = splines(i)%points_coord(:,1) - xcenter
        splines(i)%points_coord(:,2) = splines(i)%points_coord(:,2) - ycenter

        splines(i)%points_coord = TRANSPOSE(MATMUL(RotMat,TRANSPOSE(splines(i)%points_coord)))

        splines(i)%points_coord(:,1) = splines(i)%points_coord(:,1) + xcenter
        splines(i)%points_coord(:,2) = splines(i)%points_coord(:,2) + ycenter

      ENDDO

  end subroutine generate_splines_from_geo_file


END MODULE
