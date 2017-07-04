  !*********************************************************************/
  !*****     Author: D. Dagnino     ************************************/
  !*********************************************************************/
  module m_dag_sort
    
    contains
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    ! Funci�n     : void QuickSort(double *lista, int CantElem, int *indices)
    ! Descripci�n : Esta funci�n ordena una lista (double) con el algoritmo
    !  QuickSort y adem�s guarda los indices antiguoas.
    ! Par�metros  : *lista   -> La lista a ordenar.
    !               CantElem -> La cantidad de elementos que tiene la lista.
    !               *indices -> Es el output de los indices.
    ! Retorna     : La lista ordenada y un array con los indices que tenia cada
    !  elemento antes de ornerar la lista.
    !**********************************************************!
    recursive subroutine quick_sort( size_list, list, order, start )
    ! Quick sort routine from:
    ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
    ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
    ! Modified by Alan Miller to include an associated integer array which gives
    ! the positions of the elements in the original order.
    implicit none
    integer, intent(in)    :: size_list
    real(8), intent(inout) :: list(:)
    integer, intent(inout) :: order(:)
    logical, intent(in)    :: start
    ! Local variable
    integer :: i
    
    if ( start ) then
      do i=1,size_list
        order(i) = i
      end do
    end if

    call quick_sort_1( 1, size_list )
    
    !**********************************************************!
    contains

    !----------------------------------------------------------!
    recursive subroutine quick_sort_1( left_end, right_end )
    ! 
    integer, intent(in) :: left_end, right_end
    ! Local variables
    integer             :: i, j, itemp
    real(8)             :: reference, temp
    integer, parameter  :: max_simple_sort_size = 6

    if (right_end < left_end + max_simple_sort_size) then
      ! Use interchange sort for small lists
      call interchange_sort(left_end, right_end)

    else
      ! Use partition ("quick") sort
      reference = list( (left_end + right_end)/2 )
      i = left_end  - 1
      j = right_end + 1

      do
        ! Scan list from left end until element >= reference is found
        do
          i = i + 1
          if (list(i) >= reference) exit
        end do
        ! Scan list from right end until element <= reference is found
        do
          j = j - 1
          if (list(j) <= reference) exit
        end do


        if (i < j) then
          ! Swap two out-of-order elements
          temp    = list(i)
          list(i) = list(j)
          list(j) = temp

          itemp    = order(i)
          order(i) = order(j)
          order(j) = itemp
        else if (i == j) then
          i = i + 1
          exit
        else
          exit
        end if
      end do

      if ( left_end < j  ) call quick_sort_1(left_end, j)
      if ( i < right_end ) call quick_sort_1(i, right_end)
    end if

    end subroutine quick_sort_1

    !----------------------------------------------------------!
    subroutine interchange_sort( left_end, right_end )

    integer, intent(in) :: left_end, right_end

    ! Local variables
    integer :: i, j, itemp
    real(8)  :: temp

    do i = left_end, right_end - 1
      do j = i+1, right_end
        if (list(i) > list(j)) then
          ! 
          temp    = list(i)
          list(i) = list(j)
          list(j) = temp
          ! 
          itemp    = order(i)
          order(i) = order(j)
          order(j) = itemp
        end if
      end do
    end do

    end subroutine interchange_sort
    
    !**********************************************************!
    end subroutine quick_sort
      
  !***********************************************************************/
  !***********************************************************************/
  !***********************************************************************/
  end module m_dag_sort

