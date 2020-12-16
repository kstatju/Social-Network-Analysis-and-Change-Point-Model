    implicit none
    real(8), allocatable :: data(:), xbar(:,:), s(:,:), dev(:), dev2(:), &
                            fract(:,:,:), window(:)
    real(4), allocatable :: tmax(:,:)
    integer :: nord, nwind, ignore, i, j, k, n, isim, nsim, ij, kl, where, &
               nsurv, last, iglist(5) = (/0, 20, 25, 30, 40 /), ig, ialf, &
               nordlist(5) = (/2, 3, 4, 5, 10/), inord, outer
    real(8) :: one = 1.d0, zero = 0.d0, rat, denom, count, spl, splmax, tsq,&
               tsqmax, getper, alph(5) = (/0.99, 0.995, 0.998, 0.999, 0.9995/)
    logical(1), allocatable :: alive(:)

    logical :: first


    ij = 10
    kl = 11

    nsim = 200000
    nwind = 1000
    nord = 5
    first = .true.

    do outer = 1, 20
      do inord = 5, 1, -1
        nord = nordlist(inord)

        if (.not. first) deallocate (data, xbar, s, dev, dev2, tmax, fract, &
          window, alive)
        first = .false.
        allocate (data(nord), xbar(nord,0:nwind), s(nord,nord), dev(nord), &
          dev2(nord), tmax(nsim,nwind), fract(nwind, 5, 5), window(nsim), &
          alive(nsim))

        iglist(1) = nord + 1

        tmax = -1
        do isim = 1, nsim
          xbar(:,0) = zero
          s = zero
          count = zero
          do n = 1, nwind
            call randn(data, nord, ij, kl)
!           write(*,'(5f8.4)') data
            count = count + one
            dev = data - xbar(:,n-1)
            xbar(:,n) = xbar(:,n-1) + dev / count

            select case (n-iglist(1))

              case(:1)                    ! Build up S
              dev2 = data - xbar(:,n)
              do j = 1, nord
                s(:,j) = s(:,j) + dev * dev2(j)
                end do

              if (n - iglist(1) == 1) then  ! Invert  S
                do j = 1, nord
                  call sweep(s, nord, j, one)
                  end do
                s = -s
                end if

              case(2:)                  ! Production.  Update inverse
              do j = 1, nord
                dev2(j) = dot_product(s(:,j), dev)
                end do
              rat = (count - one) / count
              denom = one + rat * dot_product(dev, dev2)
              do j = 1, nord
                s(:,j) = s(:,j) - rat * dev2 * dev2(j) / denom
                end do

              end select

            if (n > iglist(1)) then
              splmax = 0
              do k = 1, n-1
                dev = xbar(:,k) - xbar(:,n)
                do j = 1, nord
                  dev2(j) = dot_product(s(:,j), dev)
                  end do
                rat = real(n) * real(k) / real(n-k)
                spl = rat * dot_product(dev, dev2)
                tsq = spl / (one - spl) * real(n-2)
!                if (n == nwind) write(*,'(i5, 2f9.3)') k, spl, tsq
                if (spl > splmax) then
                  splmax  = spl
                  where = k
                  end if
                end do

              tsqmax = splmax / (one - splmax) * real(n - 2)
              tmax(isim, n) = tsqmax
!             write(*,'(2i5, 2f12.3)') isim, n, splmax, tsqmax
              if (n < nwind) cycle
              end if
            end do
          end do
!
!        harvest
!
        fract = -99
        do ig = 1, 5
          do ialf = 1, 5
            alive = .true.
            nsurv = nsim
            do i = iglist(ig)+1, nwind
!           write(*,'(10f8.2)') tmax(:,i)
              last = nsurv
              nsurv = 0
              do j = 1, nsim
                if (alive(j)) then
                  nsurv = nsurv + 1
                  window(nsurv) = tmax(j,i)
                  end if
                end do
              if (nsurv > 0) fract(i, ig, ialf) = getper(window, nsurv, alph(ialf))
!             write(*,'(4i8,2f9.6,g15.7)') i, nsurv, ig, ialf, real(nsurv)/real(last), &
!               fract(i, ig, ialf)
              do j = 1, nsim
                alive(j) = alive(j) .and. tmax(j,i) < fract(i, ig, ialf)
                end do
              end do
            end do
          end do

!       open(20, file='fractab.out', action='write', position='append')
        open(20, file='fractab.lng', action='write', position='append')
        do ialf = 1, 5
          write(*,'('' p='', i5, '' alpha = '',f8.6/'' n     ignore'',i3,10i14)') &
            nord, alph(ialf), iglist
          do i = iglist(1) + 1, nwind
            write(*, '(i5,3x,5g14.5)') i, fract(i, :, ialf)
            write(20, '(3i5,5g14.5)') nord, ialf, i, fract(i, :, ialf)
            end do
          end do
        close(20)
        end do
      end do
    end program

