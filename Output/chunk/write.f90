program abc
implicit none
integer * 1 :: i
integer * 1 :: j
integer * 1 :: k
integer * 1 :: l

READ(*,*) i, j, k, l
open(1, file='test,txt')
write(1,*) i, j, k, l

end program abc
