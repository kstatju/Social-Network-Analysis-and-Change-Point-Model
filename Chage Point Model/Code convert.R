# implicit none
# real(8), allocatable :: data(:), xbar(:,:), s(:,:), dev(:), dev2(:), &
#   fract(:,:,:), window(:)
# real(4), allocatable :: tmax(:,:)
# integer :: nord, nwind, ignore, i, j, k, n, isim, nsim, ij, kl, where, &
#   nsurv, last, iglist(5) = (/0, 20, 25, 30, 40 /), ig, ialf, &
#   nordlist(5) = (/2, 3, 4, 5, 10/), inord, outer
# real(8) :: one = 1.d0, zero = 0.d0, rat, denom, count, spl, splmax, tsq,&
#   tsqmax, getper, alph(5) = (/0.99, 0.995, 0.998, 0.999, 0.9995/)
# logical(1), allocatable :: alive(:)
# 
# logical :: first

iglist = c(0, 20, 25, 30, 40)
nordlist = c(2, 3, 4, 5, 10)
alph = c(0.99, 0.995, 0.998, 0.999, 0.9995)
1-alph
ij = 10
kl = 11
zero = 0.0
one = 1.0
nsim = 200000
nwind = 1000
nord = 5
first = TRUE

for (outer in 1:20){
  for (inord in 5:1){
    nord = nordlist[inord]
    
    # if (!first){deallocate (data, xbar, s, dev, dev2, tmax, fract, &
    #                                window, alive)}
    
    first = FALSE
    # allocate (data(nord), xbar(nord,0:nwind), s(nord,nord), dev(nord), &
    #             dev2(nord), tmax(nsim,nwind), fract(nwind, 5, 5), window(nsim), &
    #             alive(nsim))
    
    iglist[1] = nord + 1
    xbar = array(data = NA, dim = c(nord,nwind+1))
    s = array(data = NA, dim = c(nord,nord))
    tmax = array(data = -1, dim = c(nsim,nwind))
    alive = array(data = TRUE, dim = nsim)
    window = array(data = NA, dim = nsim)
    fract = array(data = NA, dim = c(nwind, 5, 5))
    dev2 = array(data = NA, dim = nord)
    for (isim in 1:nsim){
      xbar[,1] = zero
      s[,] = zero
      count = zero
      for (n in 2:(nwind+1)){
        #call randn(data, nord, ij, kl) ********************************
        data = rnorm(nord, ij, kl)
        # !           write(*,'(5f8.4)') data
        count = count + one
        dev = data - xbar[,n-1]
        xbar[,n] = xbar[,n-1] + dev / count
        

        
         if(n-iglist[1] <= 1) {                   #! Build up S
            dev2[] = data - xbar[,n]
            for (j in 1:nord){
              s[,j] = s[,j] + dev * dev2[j]
            }
            
            if (n - iglist[1] == 1) {  #! Invert  S
              for (j in 1:nord){
                call sweep(s, nord, j, one) ****************************************
              }
              s = -s
            }
          }
          
          if(n-iglist[1] >= 2) {                 # ! Production.  Update inverse
            for (j in 1:nord){
              dev2[j] = s[,j] %*% dev
            }
            rat = (count - one) / count
            denom = one + rat * (dev %*% dev2)
            for (j in 1:nord){
              s[,j] = s[,j] - rat * dev2 * dev2[j] / denom
            }
        
          }
        
        if (n > iglist[1]) {
          splmax = 0
          for (k in 1:(n-1)){
            dev = xbar[,k] - xbar[,n]
            for (j in 1:nord){
              dev2[j] = s[,j] %*% dev
            }
            rat = n * k / (n-k)
            spl = rat * (dev %*% dev2)
            tsq = spl / (one - spl) * (n-2)
            # !                if (n == nwind) write(*,'(i5, 2f9.3)') k, spl, tsq
            if (spl > splmax) {
              splmax  = spl
              where = k
            }
          }
          
          tsqmax = splmax / (one - splmax) * (n - 2)
          tmax[isim, n] = tsqmax
          # !             write(*,'(2i5, 2f12.3)') isim, n, splmax, tsqmax
          if (n < nwind){ 
            cycle 
          }
        }
      }
    }
    
    
    
    
    
    
    
    
    # !
    #   !        harvest
    # !
      fract[,,] = -99
    for (ig in 1:5){
      for (ialf in 1:5){
        alive[] = TRUE
        nsurv = nsim
        for (i in (iglist[ig]+1):nwind){
          # !           write(*,'(10f8.2)') tmax(:,i)
          last = nsurv
          nsurv = 0
          for (j in 1:nsim){
            if (alive[j]) {
              nsurv = nsurv + 1
              window[nsurv] = tmax[j,i]
            }
          }
          if (nsurv > 0){
            fract[i, ig, ialf] = getper(window, nsurv, alph[ialf]) ******************************
          }
          # !             write(*,'(4i8,2f9.6,g15.7)') i, nsurv, ig, ialf, real(nsurv)/real(last), &
          #   !               fract(i, ig, ialf)
          for (j in 1:nsim){
            alive[j] = alive[j] * tmax[j,i] < fract[i, ig, ialf]
          }
        }
      }
    }
    
    # !       open(20, file='fractab.out', action='write', position='append')
    #open(20, file='fractab.lng', action='write', position='append')
    for (ialf in 1:5){
      line=c(nord, alph[ialf], iglist)
      write(line,file="C:/Users/ka746940/Desktop/UCF/STA 6908 - Edgard Maboudou/R code/myfile.txt",append=TRUE)
      #write(*,'('' p='', i5, '' alpha = '',f8.6/'' n     ignore'',i3,10i14)') &
      #  nord, alph(ialf), iglist
      
      for (i in (iglist[1] + 1):nwind){
        
        line1=c(i, fract[i, , ialf])
        write(line1,file="C:/Users/ka746940/Desktop/UCF/STA 6908 - Edgard Maboudou/R code/myfile1.txt",append=TRUE)
        
        line2=c(nord, ialf, i, fract[i, , ialf])
        write(line2,file="C:/Users/ka746940/Desktop/UCF/STA 6908 - Edgard Maboudou/R code/myfile2.txt",append=TRUE)
        
        #write(*, '(i5,3x,5g14.5)') i, fract(i, :, ialf)
        #write(20, '(3i5,5g14.5)') nord, ialf, i, fract(i, :, ialf)
      }
    }
  }
}
