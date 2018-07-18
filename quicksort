    !------------------------------------------------------------------------------
    !> @brief
    !> Indexes an array arr(1:n), i.e., outputs the array indx(1:n) such 
    !> such that arr(indx(j)) is in ascending order for j=1,2...n
    ! 
    !> @author 
    !> D.Gates
    !
    !> @remarks
    !> The pivot selection is done using a median of three approach to avoid a
    !> a worst case scenario in which the array is sorted already. In 
    !> addition, for small partitions, we use insertion sort as it
    !> outperforms all other sort routines for small array sizes. 
    !>
    !> @param[in]  first start index
    !> @param[in]  last  array size
    !> @param[in]  a     input array to be sorted
    !> @param[out] indx  array of sorted indices
    !           
    !---------------------------------------------------------------------------
      recursive subroutine quicksort(first,last,a,indx)
      implicit none
      integer a(*),temp,indx(*)
      integer pivot
      integer i, j, x, t
      integer first, last
      
      if(last-first < 32) then 
          do i = 2+first, last
              x = a(indx(i))
              t = indx(i)
              j = i - 1
              do while (j >= 1)
                  if (a(indx(j)) <= x) exit
                  indx(j+1) = indx(j)
                  j = j - 1
              end do
              indx(j+1) = t
          end do
      else
          pivot = (first + last) / 2
          if (a(indx(pivot)) < a(indx(first))) then
               temp = a(indx(pivot))
               a(indx(pivot)) = a(indx(first))
               a(indx(first)) = temp
          endif
          if (a(indx(last)) < a(indx(first))) then
              temp = a(indx(last))
              a(indx(last)) = a(indx(first))
              a(indx(first)) = temp
          endif
          if (a(indx(pivot)) < a(indx(last))) then 
              temp = a(indx(pivot))
              a(indx(pivot)) = a(indx(last))
              a(indx(last)) = temp
          endif
          pivot = a(indx(last))
          i = first
          j = last
          do
              do while (a(indx(i)) < pivot)
                  i=i+1
              enddo
              do while (pivot < a(indx(j)))
                  j=j-1
              enddo
              if (i >= j) exit
              temp = indx(i);indx(i) = indx(j);indx(j) = temp
              i=i+1
              j=j-1
          enddo
          if (first < i-1) call quicksort(first,i-1,a,indx)
          if (j+1 < last)  call quicksort(j+1,last,a,indx)
      endif
      END subroutine quicksort

