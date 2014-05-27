! -*- f90 -*-
! FILE: fortran_routines.f90
!  f2py -c -m fortran_routines fortran_routines.f90
!  This file will have a collection of fortran routines that can be interfaced to python
!  To be used with PyAnalog
! ======================

Subroutine rank_analog(trainField,fcstField,trainnum,iNum,jNum,&
     startLat,endLat,startLon,endLon,allLats,allLons,window,outRanks)
INTEGER, INTENT(IN) :: iNum,jNum,startLat,endLat,startLon,endLon,window,trainnum
real, intent(IN), dimension(iNum) :: allLats
real, intent(IN), dimension(jNum) :: allLons
real, INTENT(IN), DIMENSION(iNum,jNum) :: fcstField
real, INTENT(IN), DIMENSION(trainnum,iNum,jNum) :: trainField
real, intent(out), dimension(trainnum,iNum,jNum) :: outRanks

! --- Now, some other variables we'll need
integer :: i,j,np,startLatIdx,startLonIdx,endLonIdx,endLatIdx,n_grdpts
real, DIMENSION(:,:,:),allocatable :: ranks
real, DIMENSION(:,:,:),allocatable :: allData
integer, DIMENSION(:),allocatable :: dummyIdx,dummy_finalRanks,dummy_ranks
real, dimension(:,:),allocatable :: trainData
real, dimension(:), allocatable :: dummy_array,rankDiffs

!f2py intent(in) iNum,jNum,trainnum,startLat,endLat,startLon,endLon,window,allLats,allLons,trainField,fcstField
!f2py depend(iNum) allLats
!f2py depend(jNum) allLons
!f2py depend(iNum,jNum) fcstField
!f2py depend(trainnum,iNum,jNum) trainField
!f2py intent(out) outRanks
!f2py depend(trainnum,iNum,jNum) outRanks

outRanks(:,:,:) = -9999.9

! --- Get latitude indices to iterate over
do i = 1,iNum
   if (allLats(i) .eq. startLat) then
      startLatIdx = i
   elseif (allLats(i) .eq. endLat) then
      endLatIdx = i
   end if
end do

! --- Get longitude indices to iterate over
do j = 1,jNum
   if (allLons(j) .eq. startLon) then
      startLonIdx = j
   elseif (allLons(j) .eq. endLon) then
      endLonIdx = j
   end if
end do

! --- Start big coefficient loop
n_grdpts = ((window*2)+1)*((window*2)+1)

! --- Let's allocate some arrays
allocate(allData(trainnum+1,iNum,jNum)) ! --- train data + forecast data
allocate(ranks(trainnum+1,iNum,jNum)) ! --- Array of rankings
allocate(dummy_ranks(trainnum+1)) ! --- temporary ranking array
allocate(dummy_array(trainnum+1)) ! --- temporary data array
allocate(trainData(trainnum+1,((window*2)+1)*((window*2)+1))) ! --- 2-D training data
allocate(rankDiffs(trainnum)) ! --- mean absolute difference of ranks
allocate(dummy_finalRanks(trainnum)) ! --- rankings of rankDiffs array
allocate(dummyIdx(trainnum)) ! --- indices of training dates

! --- We need to put the training and fcst data in one array
allData(:trainnum,:,:) = trainField(:,:,:) ! --- Training data
allData(trainnum+1,:,:) = fcstField(:,:) ! --- Forecast data

! --- Ranking each grid point in this loop
do j = 1,jNum
   do i = 1,iNum
      dummy_array(:) = allData(:,i,j)
      dummy_ranks(:)=0
      call real_rank(dummy_array,size(dummy_array),dummy_ranks)
      call ties(dummy_ranks,dummy_array,size(dummy_array),ranks(:,i,j))
   end do
end do


! --- Now time to find ranking differences at each grid point given a window of grid points
do j = startLonIdx,endLonIdx
   do i = startLatIdx,endLatIdx

      rankDiffs(:) = 999999
      ! --- First, extract data
      trainData(:,:) = reshape(ranks(:,i-window:i+window,j-window:j+window),(/trainnum+1,n_grdpts/))

      ! --- Find sum of absolute value of rank differences
      do np = 1,trainnum
         rankDiffs(np) = sum(abs(trainData(np,:)-trainData(trainnum+1,:)),dim=1)
      end do

      ! --- Ranking the ranked diffs to make it easier to get probabilities
      call real_rank(rankDiffs,size(rankDiffs),dummy_finalRanks)
      do np = 1,size(dummy_finalRanks)
         outRanks(dummy_finalRanks(np),i,j) = np-1
      end do


   end do
end do

! --- Time to deallocate those arrays
deallocate(allData) ! --- train data + forecast data

deallocate(ranks) ! --- Array of rankings

deallocate(dummy_ranks) ! --- temporary ranking array

deallocate(trainData) ! --- 2-D training data

deallocate(rankDiffs) ! --- mean absolute difference of ranks

deallocate(dummy_finalRanks) ! --- rankings of rankDiffs array

deallocate(dummyIdx) ! --- indices of training dates

deallocate(dummy_array) ! --- temporary data array

return

end Subroutine rank_analog

Subroutine rmse_analog(trainField,fcstField,trainnum,iNum,jNum,&
     startLat,endLat,startLon,endLon,allLats,allLons,window,outRanks)
INTEGER, INTENT(IN) :: iNum,jNum,startLat,endLat,startLon,endLon,window,trainnum
real, intent(IN), dimension(iNum) :: allLats
real, intent(IN), dimension(jNum) :: allLons
real, INTENT(IN), DIMENSION(iNum,jNum) :: fcstField
real, INTENT(IN), DIMENSION(trainnum,iNum,jNum) :: trainField
real, intent(out), dimension(trainnum,iNum,jNum) :: outRanks

! --- Now, some other variables we'll need
integer :: i,j,np,startLatIdx,startLonIdx,endLonIdx,endLatIdx,n_grdpts
real, DIMENSION(:,:,:),allocatable :: allData
integer, DIMENSION(:),allocatable :: dummy_finalRanks
real, dimension(:,:),allocatable :: trainData
real, dimension(:), allocatable :: rankDiffs

!f2py intent(in) iNum,jNum,trainnum,startLat,endLat,startLon,endLon,window,allLats,allLons,trainField,fcstField
!f2py depend(iNum) allLats
!f2py depend(jNum) allLons
!f2py depend(iNum,jNum) fcstField
!f2py depend(trainnum,iNum,jNum) trainField
!f2py intent(out) outRanks
!f2py depend(trainnum,iNum,jNum) outRanks

outRanks(:,:,:) = -9999.9

! --- Get latitude indices to iterate over
do i = 1,iNum
   if (allLats(i) .eq. startLat) then
      startLatIdx = i
   elseif (allLats(i) .eq. endLat) then
      endLatIdx = i
   end if
end do

! --- Get longitude indices to iterate over
do j = 1,jNum
   if (allLons(j) .eq. startLon) then
      startLonIdx = j
   elseif (allLons(j) .eq. endLon) then
      endLonIdx = j
   end if
end do

! --- Start big coefficient loop
n_grdpts = ((window*2)+1)*((window*2)+1)

! --- Let's allocate some arrays
allocate(allData(trainnum+1,iNum,jNum)) ! --- train data + forecast data
allocate(trainData(trainnum+1,((window*2)+1)*((window*2)+1))) ! --- 2-D training data
allocate(rankDiffs(trainnum)) ! --- mean absolute difference of ranks
allocate(dummy_finalRanks(trainnum)) ! --- rankings of rankDiffs array

! --- We need to put the training and fcst data in one array
allData(:trainnum,:,:) = trainField(:,:,:) ! --- Training data
allData(trainnum+1,:,:) = fcstField(:,:) ! --- Forecast data

! --- Now time to find ranking differences at each grid point given a window of grid points
do j = startLonIdx,endLonIdx
   do i = startLatIdx,endLatIdx

      rankDiffs(:) = 999999
      ! --- First, extract data
      trainData(:,:) = reshape(allData(:,i-window:i+window,j-window:j+window),(/trainnum+1,n_grdpts/))

      ! --- Find sum of absolute value of rank differences
      do np = 1,trainnum
         call rmse(trainData(np,:),trainData(trainnum+1,:),n_grdpts,rankDiffs(np))
      end do

      ! --- Ranking the ranked diffs to make it easier to get probabilities
      call real_rank(rankDiffs,size(rankDiffs),dummy_finalRanks)
      do np = 1,size(dummy_finalRanks)
         outRanks(dummy_finalRanks(np),i,j) = np-1
      end do

   end do
end do

! --- Time to deallocate those arrays
deallocate(allData) ! --- train data + forecast data

deallocate(trainData) ! --- 2-D training data

deallocate(rankDiffs) ! --- mean absolute difference of ranks

deallocate(dummy_finalRanks) ! --- rankings of rankDiffs array

return

end Subroutine rmse_analog


Subroutine mae_analog(trainField,fcstField,trainnum,iNum,jNum,&
     startLat,endLat,startLon,endLon,allLats,allLons,window,outRanks)
INTEGER, INTENT(IN) :: iNum,jNum,startLat,endLat,startLon,endLon,window,trainnum
real, intent(IN), dimension(iNum) :: allLats
real, intent(IN), dimension(jNum) :: allLons
real, INTENT(IN), DIMENSION(iNum,jNum) :: fcstField
real, INTENT(IN), DIMENSION(trainnum,iNum,jNum) :: trainField
real, intent(out), dimension(trainnum,iNum,jNum) :: outRanks

! --- Now, some other variables we'll need
integer :: i,j,np,startLatIdx,startLonIdx,endLonIdx,endLatIdx,n_grdpts
real, DIMENSION(:,:,:),allocatable :: allData
integer, DIMENSION(:),allocatable :: dummy_finalRanks
real, dimension(:,:),allocatable :: trainData
real, dimension(:), allocatable :: rankDiffs

!f2py intent(in) iNum,jNum,trainnum,startLat,endLat,startLon,endLon,window,allLats,allLons,trainField,fcstField
!f2py depend(iNum) allLats
!f2py depend(jNum) allLons
!f2py depend(iNum,jNum) fcstField
!f2py depend(trainnum,iNum,jNum) trainField
!f2py intent(out) outRanks
!f2py depend(trainnum,iNum,jNum) outRanks

outRanks(:,:,:) = -9999.9

! --- Get latitude indices to iterate over
do i = 1,iNum
   if (allLats(i) .eq. startLat) then
      startLatIdx = i
   elseif (allLats(i) .eq. endLat) then
      endLatIdx = i
   end if
end do

! --- Get longitude indices to iterate over
do j = 1,jNum
   if (allLons(j) .eq. startLon) then
      startLonIdx = j
   elseif (allLons(j) .eq. endLon) then
      endLonIdx = j
   end if
end do

! --- Start big coefficient loop
n_grdpts = ((window*2)+1)*((window*2)+1)

! --- Let's allocate some arrays
allocate(allData(trainnum+1,iNum,jNum)) ! --- train data + forecast data
allocate(trainData(trainnum+1,((window*2)+1)*((window*2)+1))) ! --- 2-D training data
allocate(rankDiffs(trainnum)) ! --- mean absolute difference of ranks
allocate(dummy_finalRanks(trainnum)) ! --- rankings of rankDiffs array

! --- We need to put the training and fcst data in one array
allData(:trainnum,:,:) = trainField(:,:,:) ! --- Training data
allData(trainnum+1,:,:) = fcstField(:,:) ! --- Forecast data

! --- Now time to find ranking differences at each grid point given a window of grid points
do j = startLonIdx,endLonIdx
   do i = startLatIdx,endLatIdx

      rankDiffs(:) = 999999
      ! --- First, extract data
      trainData(:,:) = reshape(allData(:,i-window:i+window,j-window:j+window),(/trainnum+1,n_grdpts/))

      ! --- Find sum of absolute value of rank differences
      do np = 1,trainnum
         call mae(trainData(np,:),trainData(trainnum+1,:),n_grdpts,rankDiffs(np))
      end do

      ! --- Ranking the ranked diffs to make it easier to get probabilities
      call real_rank(rankDiffs,size(rankDiffs),dummy_finalRanks)
      do np = 1,size(dummy_finalRanks)
         outRanks(dummy_finalRanks(np),i,j) = np-1
      end do

   end do
end do

! --- Time to deallocate those arrays
deallocate(allData) ! --- train data + forecast data

deallocate(trainData) ! --- 2-D training data

deallocate(rankDiffs) ! --- mean absolute difference of ranks

deallocate(dummy_finalRanks) ! --- rankings of rankDiffs array

return

end Subroutine mae_analog

Subroutine corr_analog(trainField,fcstField,trainnum,iNum,jNum,&
     startLat,endLat,startLon,endLon,allLats,allLons,window,outRanks)
INTEGER, INTENT(IN) :: iNum,jNum,startLat,endLat,startLon,endLon,window,trainnum
real, intent(IN), dimension(iNum) :: allLats
real, intent(IN), dimension(jNum) :: allLons
real, INTENT(IN), DIMENSION(iNum,jNum) :: fcstField
real, INTENT(IN), DIMENSION(trainnum,iNum,jNum) :: trainField
real, intent(out), dimension(trainnum,iNum,jNum) :: outRanks

! --- Now, some other variables we'll need
integer :: i,j,np,startLatIdx,startLonIdx,endLonIdx,endLatIdx,n_grdpts
real, DIMENSION(:,:,:),allocatable :: allData
integer, DIMENSION(:),allocatable :: dummy_finalRanks
real, dimension(:,:),allocatable :: trainData
real, dimension(:), allocatable :: rankDiffs

!f2py intent(in) iNum,jNum,trainnum,startLat,endLat,startLon,endLon,window,allLats,allLons,trainField,fcstField
!f2py depend(iNum) allLats
!f2py depend(jNum) allLons
!f2py depend(iNum,jNum) fcstField
!f2py depend(trainnum,iNum,jNum) trainField
!f2py intent(out) outRanks
!f2py depend(trainnum,iNum,jNum) outRanks

outRanks(:,:,:) = -9999.9

! --- Get latitude indices to iterate over
do i = 1,iNum
   if (allLats(i) .eq. startLat) then
      startLatIdx = i
   elseif (allLats(i) .eq. endLat) then
      endLatIdx = i
   end if
end do

! --- Get longitude indices to iterate over
do j = 1,jNum
   if (allLons(j) .eq. startLon) then
      startLonIdx = j
   elseif (allLons(j) .eq. endLon) then
      endLonIdx = j
   end if
end do

! --- Start big coefficient loop
n_grdpts = ((window*2)+1)*((window*2)+1)

! --- Let's allocate some arrays
allocate(allData(trainnum+1,iNum,jNum)) ! --- train data + forecast data
allocate(trainData(trainnum+1,((window*2)+1)*((window*2)+1))) ! --- 2-D training data
allocate(rankDiffs(trainnum)) ! --- mean absolute difference of ranks
allocate(dummy_finalRanks(trainnum)) ! --- rankings of rankDiffs array

! --- We need to put the training and fcst data in one array
allData(:trainnum,:,:) = trainField(:,:,:) ! --- Training data
allData(trainnum+1,:,:) = fcstField(:,:) ! --- Forecast data

! --- Now time to find ranking differences at each grid point given a window of grid points
do j = startLonIdx,endLonIdx
   do i = startLatIdx,endLatIdx

      rankDiffs(:) = 999999
      ! --- First, extract data
      trainData(:,:) = reshape(allData(:,i-window:i+window,j-window:j+window),(/trainnum+1,n_grdpts/))

      ! --- Find sum of absolute value of rank differences
      do np = 1,trainnum
         call correlation(trainData(np,:),trainData(trainnum+1,:),n_grdpts,rankDiffs(np))
      end do

      ! --- Ranking the ranked diffs to make it easier to get probabilities
      call real_rank(rankDiffs,size(rankDiffs),dummy_finalRanks)
      do np = 1,size(dummy_finalRanks)
         outRanks(dummy_finalRanks(np),i,j) = np-1
      end do

   end do
end do

! --- Time to deallocate those arrays
deallocate(allData) ! --- train data + forecast data

deallocate(trainData) ! --- 2-D training data

deallocate(rankDiffs) ! --- mean absolute difference of ranks

deallocate(dummy_finalRanks) ! --- rankings of rankDiffs array

return

end Subroutine corr_analog

Subroutine rank_analog_point(trainField,fcstField,trainnum,iNum,jNum,&
     closeLat,closeLon,allLats,allLons,window,outRanks)
INTEGER, INTENT(IN) :: iNum,jNum,closeLat,closeLon,window,trainnum
real, intent(IN), dimension(iNum) :: allLats
real, intent(IN), dimension(jNum) :: allLons
real, INTENT(IN), DIMENSION(iNum,jNum) :: fcstField
real, INTENT(IN), DIMENSION(trainnum,iNum,jNum) :: trainField
real, intent(out), dimension(trainnum) :: outRanks

! --- Now, some other variables we'll need
integer :: i,j,np,startLatIdx,startLonIdx,n_grdpts
real, DIMENSION(:,:,:),allocatable :: ranks
real, DIMENSION(:,:,:),allocatable :: allData
integer, DIMENSION(:),allocatable :: dummyIdx,dummy_finalRanks,dummy_ranks
real, dimension(:,:),allocatable :: trainData
real, dimension(:), allocatable :: dummy_array,rankDiffs

!f2py intent(in) iNum,jNum,trainnum,closeLat,closeLon,window,allLats,allLons,trainField,fcstField
!f2py depend(iNum) allLats
!f2py depend(jNum) allLons
!f2py depend(iNum,jNum) fcstField
!f2py depend(trainnum,iNum,jNum) trainField
!f2py intent(out) outRanks
!f2py depend(trainnum) outRanks

outRanks(:) = -9999.9

! --- Get latitude indices to iterate over
do i = 1,iNum
   if (allLats(i) .eq. closeLat) then
      startLatIdx = i
   end if
end do

! --- Get longitude indices to iterate over
do j = 1,jNum
   if (allLons(j) .eq. closeLon) then
      startLonIdx = j
   end if
end do

! --- Start big coefficient loop
n_grdpts = ((window*2)+2)*((window*2)+2)

! --- Let's allocate some arrays
allocate(allData(trainnum+1,iNum,jNum)) ! --- train data + forecast data
allocate(ranks(trainnum+1,iNum,jNum)) ! --- Array of rankings
allocate(dummy_ranks(trainnum+1)) ! --- temporary ranking array
allocate(dummy_array(trainnum+1)) ! --- temporary data array
allocate(trainData(trainnum+1,((window*2)+1)*((window*2)+1))) ! --- 2-D training data
allocate(rankDiffs(trainnum)) ! --- mean absolute difference of ranks
allocate(dummy_finalRanks(trainnum)) ! --- rankings of rankDiffs array
allocate(dummyIdx(trainnum)) ! --- indices of training dates

! --- We need to put the training and fcst data in one array
allData(:trainnum,:,:) = trainField(:,:,:) ! --- Training data
allData(trainnum+1,:,:) = fcstField(:,:) ! --- Forecast data

! --- Ranking each grid point in this loop
do j = 1,jNum
   do i = 1,iNum
      dummy_array(:) = allData(:,i,j)
      dummy_ranks(:)=0
      call real_rank(dummy_array,size(dummy_array),dummy_ranks)
      call ties(dummy_ranks,dummy_array,size(dummy_array),ranks(:,i,j))
   end do
end do

rankDiffs(:) = 999999
! --- First, extract data
! --- We assume the closest Lat/Lon combo is the NW corner of the box
trainData(:,:) = reshape(ranks(:,i-window:i+window+1,j-window:j+window+1),(/trainnum+1,n_grdpts/))

! --- Find sum of absolute value of rank differences
do np = 1,trainnum
   call rank_corr(trainData(np,:),trainData(trainnum+1,:),n_grdpts,rankDiffs(np))
end do

! --- Ranking the ranked diffs to make it easier to get probabilities
call real_rank(rankDiffs,size(rankDiffs),dummy_finalRanks)
do np = 1,size(dummy_finalRanks)
   outRanks(dummy_finalRanks(np)) = np-1
end do


! --- Time to deallocate those arrays
deallocate(allData) ! --- train data + forecast data

deallocate(ranks) ! --- Array of rankings

deallocate(dummy_ranks) ! --- temporary ranking array

deallocate(trainData) ! --- 2-D training data

deallocate(rankDiffs) ! --- mean absolute difference of ranks

deallocate(dummy_finalRanks) ! --- rankings of rankDiffs array

deallocate(dummyIdx) ! --- indices of training dates

deallocate(dummy_array) ! --- temporary data array

return

end Subroutine rank_analog_point

Subroutine rmse_analog_point(trainField,fcstField,trainnum,iNum,jNum,&
     closeLat,closeLon,allLats,allLons,window,outRanks)
INTEGER, INTENT(IN) :: iNum,jNum,closeLat,closeLon,window,trainnum
real, intent(IN), dimension(iNum) :: allLats
real, intent(IN), dimension(jNum) :: allLons
real, INTENT(IN), DIMENSION(iNum,jNum) :: fcstField
real, INTENT(IN), DIMENSION(trainnum,iNum,jNum) :: trainField
real, intent(out), dimension(trainnum) :: outRanks

! --- Now, some other variables we'll need
integer :: i,j,np,startLatIdx,startLonIdx,endLonIdx,endLatIdx,n_grdpts
real, DIMENSION(:,:,:),allocatable :: allData
integer, DIMENSION(:),allocatable :: dummy_finalRanks
real, dimension(:,:),allocatable :: trainData
real, dimension(:), allocatable :: rankDiffs

!f2py intent(in) iNum,jNum,trainnum,closeLat,closeLon,window,allLats,allLons,trainField,fcstField
!f2py depend(iNum) allLats
!f2py depend(jNum) allLons
!f2py depend(iNum,jNum) fcstField
!f2py depend(trainnum,iNum,jNum) trainField
!f2py intent(out) outRanks
!f2py depend(trainnum) outRanks

outRanks(:) = -9999.9

! --- Get latitude indices to iterate over
do i = 1,iNum
   if (allLats(i) .eq. closeLat) then
      startLatIdx = i
   end if
end do

! --- Get longitude indices to iterate over
do j = 1,jNum
   if (allLons(j) .eq. closeLon) then
      startLonIdx = j
   end if
end do

! --- Start big coefficient loop
n_grdpts = ((window*2)+2)*((window*2)+2)

! --- Let's allocate some arrays
allocate(allData(trainnum+1,iNum,jNum)) ! --- train data + forecast data
allocate(trainData(trainnum+1,((window*2)+1)*((window*2)+1))) ! --- 2-D training data
allocate(rankDiffs(trainnum)) ! --- mean absolute difference of ranks
allocate(dummy_finalRanks(trainnum)) ! --- rankings of rankDiffs array

! --- We need to put the training and fcst data in one array
allData(:trainnum,:,:) = trainField(:,:,:) ! --- Training data
allData(trainnum+1,:,:) = fcstField(:,:) ! --- Forecast data

! --- Now time to find ranking differences at each grid point given a window of grid points

rankDiffs(:) = 999999
! --- First, extract data
trainData(:,:) = reshape(allData(:,i-window:i+window,j-window:j+window),(/trainnum+1,n_grdpts/))

! --- Find sum of absolute value of rank differences
do np = 1,trainnum
   call rmse(trainData(np,:),trainData(trainnum+1,:),n_grdpts,rankDiffs(np))
end do

! --- Ranking the ranked diffs to make it easier to get probabilities
call real_rank(rankDiffs,size(rankDiffs),dummy_finalRanks)
do np = 1,size(dummy_finalRanks)
   outRanks(dummy_finalRanks(np)) = np-1
end do


! --- Time to deallocate those arrays
deallocate(allData) ! --- train data + forecast data

deallocate(trainData) ! --- 2-D training data

deallocate(rankDiffs) ! --- mean absolute difference of ranks

deallocate(dummy_finalRanks) ! --- rankings of rankDiffs array

return

end Subroutine rmse_analog_point

Subroutine mae_analog_point(trainField,fcstField,trainnum,iNum,jNum,&
     closeLat,closeLon,allLats,allLons,window,outRanks)
INTEGER, INTENT(IN) :: iNum,jNum,closeLat,closeLon,window,trainnum
real, intent(IN), dimension(iNum) :: allLats
real, intent(IN), dimension(jNum) :: allLons
real, INTENT(IN), DIMENSION(iNum,jNum) :: fcstField
real, INTENT(IN), DIMENSION(trainnum,iNum,jNum) :: trainField
real, intent(out), dimension(trainnum) :: outRanks

! --- Now, some other variables we'll need
integer :: i,j,np,startLatIdx,startLonIdx,endLonIdx,endLatIdx,n_grdpts
real, DIMENSION(:,:,:),allocatable :: allData
integer, DIMENSION(:),allocatable :: dummy_finalRanks
real, dimension(:,:),allocatable :: trainData
real, dimension(:), allocatable :: rankDiffs

!f2py intent(in) iNum,jNum,trainnum,closeLat,closeLon,window,allLats,allLons,trainField,fcstField
!f2py depend(iNum) allLats
!f2py depend(jNum) allLons
!f2py depend(iNum,jNum) fcstField
!f2py depend(trainnum,iNum,jNum) trainField
!f2py intent(out) outRanks
!f2py depend(trainnum) outRanks

outRanks(:) = -9999.9

! --- Get latitude indices to iterate over
do i = 1,iNum
   if (allLats(i) .eq. closeLat) then
      startLatIdx = i
   end if
end do

! --- Get longitude indices to iterate over
do j = 1,jNum
   if (allLons(j) .eq. closeLon) then
      startLonIdx = j
   end if
end do

! --- Start big coefficient loop
n_grdpts = ((window*2)+2)*((window*2)+2)

! --- Let's allocate some arrays
allocate(allData(trainnum+1,iNum,jNum)) ! --- train data + forecast data
allocate(trainData(trainnum+1,((window*2)+1)*((window*2)+1))) ! --- 2-D training data
allocate(rankDiffs(trainnum)) ! --- mean absolute difference of ranks
allocate(dummy_finalRanks(trainnum)) ! --- rankings of rankDiffs array

! --- We need to put the training and fcst data in one array
allData(:trainnum,:,:) = trainField(:,:,:) ! --- Training data
allData(trainnum+1,:,:) = fcstField(:,:) ! --- Forecast data

! --- Now time to find ranking differences at each grid point given a window of grid points

rankDiffs(:) = 999999
! --- First, extract data
trainData(:,:) = reshape(allData(:,i-window:i+window,j-window:j+window),(/trainnum+1,n_grdpts/))

! --- Find sum of absolute value of rank differences
do np = 1,trainnum
   call mae(trainData(np,:),trainData(trainnum+1,:),n_grdpts,rankDiffs(np))
end do

! --- Ranking the ranked diffs to make it easier to get probabilities
call real_rank(rankDiffs,size(rankDiffs),dummy_finalRanks)
do np = 1,size(dummy_finalRanks)
   outRanks(dummy_finalRanks(np)) = np-1
end do


! --- Time to deallocate those arrays
deallocate(allData) ! --- train data + forecast data

deallocate(trainData) ! --- 2-D training data

deallocate(rankDiffs) ! --- mean absolute difference of ranks

deallocate(dummy_finalRanks) ! --- rankings of rankDiffs array

return

end Subroutine mae_analog_point

Subroutine corr_analog_point(trainField,fcstField,trainnum,iNum,jNum,&
     closeLat,closeLon,allLats,allLons,window,outRanks)
INTEGER, INTENT(IN) :: iNum,jNum,closeLat,closeLon,window,trainnum
real, intent(IN), dimension(iNum) :: allLats
real, intent(IN), dimension(jNum) :: allLons
real, INTENT(IN), DIMENSION(iNum,jNum) :: fcstField
real, INTENT(IN), DIMENSION(trainnum,iNum,jNum) :: trainField
real, intent(out), dimension(trainnum) :: outRanks

! --- Now, some other variables we'll need
integer :: i,j,np,startLatIdx,startLonIdx,endLonIdx,endLatIdx,n_grdpts
real, DIMENSION(:,:,:),allocatable :: allData
integer, DIMENSION(:),allocatable :: dummy_finalRanks
real, dimension(:,:),allocatable :: trainData
real, dimension(:), allocatable :: rankDiffs

!f2py intent(in) iNum,jNum,trainnum,closeLat,closeLon,window,allLats,allLons,trainField,fcstField
!f2py depend(iNum) allLats
!f2py depend(jNum) allLons
!f2py depend(iNum,jNum) fcstField
!f2py depend(trainnum,iNum,jNum) trainField
!f2py intent(out) outRanks
!f2py depend(trainnum) outRanks

outRanks(:) = -9999.9

! --- Get latitude indices to iterate over
do i = 1,iNum
   if (allLats(i) .eq. closeLat) then
      startLatIdx = i
   end if
end do

! --- Get longitude indices to iterate over
do j = 1,jNum
   if (allLons(j) .eq. closeLon) then
      startLonIdx = j
   end if
end do

! --- Start big coefficient loop
n_grdpts = ((window*2)+2)*((window*2)+2)

! --- Let's allocate some arrays
allocate(allData(trainnum+1,iNum,jNum)) ! --- train data + forecast data
allocate(trainData(trainnum+1,((window*2)+1)*((window*2)+1))) ! --- 2-D training data
allocate(rankDiffs(trainnum)) ! --- mean absolute difference of ranks
allocate(dummy_finalRanks(trainnum)) ! --- rankings of rankDiffs array

! --- We need to put the training and fcst data in one array
allData(:trainnum,:,:) = trainField(:,:,:) ! --- Training data
allData(trainnum+1,:,:) = fcstField(:,:) ! --- Forecast data

! --- Now time to find ranking differences at each grid point given a window of grid points

rankDiffs(:) = 999999
! --- First, extract data
trainData(:,:) = reshape(allData(:,i-window:i+window,j-window:j+window),(/trainnum+1,n_grdpts/))

! --- Find sum of absolute value of rank differences
do np = 1,trainnum
   call correlation(trainData(np,:),trainData(trainnum+1,:),n_grdpts,rankDiffs(np))
end do

! --- Ranking the ranked diffs to make it easier to get probabilities
call real_rank(rankDiffs,size(rankDiffs),dummy_finalRanks)
do np = 1,size(dummy_finalRanks)
   outRanks(dummy_finalRanks(np)) = np-1
end do


! --- Time to deallocate those arrays
deallocate(allData) ! --- train data + forecast data

deallocate(trainData) ! --- 2-D training data

deallocate(rankDiffs) ! --- mean absolute difference of ranks

deallocate(dummy_finalRanks) ! --- rankings of rankDiffs array

return

end Subroutine corr_analog_point


Subroutine real_rank (XDONT, NOBS, IRNGT)
!
! From the ORDERPACK 2.0 code: http://www.fortran-2000.com/rank/
!
! __________________________________________________________
!   MRGRNK = Merge-sort ranking of an array
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! _________________________________________________________
      integer, intent(in) :: NOBS
      Real, Dimension (NOBS), Intent (In) :: XDONT
      Integer, Dimension (NOBS), Intent (Out) :: IRNGT
! __________________________________________________________
      Real :: XVALA, XVALB
!
      Integer, Dimension (NOBS) :: JWRKT
      Integer :: LMTNA, LMTNC, IRNG1, IRNG2
      Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
     ! print*, "Here!"
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XDONT(IIND-1) <= XDONT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo(NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      Do
         If (NVAL <= 2) Exit
!
!   Loop on merges of A and B into C
!
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
!
!   1 2 3
!
               If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
!
!   1 3 2
!
               If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
!
!   1 2 3 4
!
            If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
!
!   1 3 x x
!
            If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
!
!   3 x x x
!
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         Exit
      End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B. This line may be activated when the
!   initial set is already close to sorted.
!
!          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
!
!  One steps in the C subset, that we build in the final rank array
!
!  Make a copy of the rank array for the merge iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
!
            XVALA = XDONT (JWRKT(IINDA))
            XVALB = XDONT (IRNGT(IINDB))
!
            Do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (XVALA > XVALB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XVALB = XDONT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XVALA = XDONT (JWRKT(IINDA))
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
      Return
!
End Subroutine real_rank


subroutine ties(ranks,vals,n_ranks,outranks)
! --- Subroutine designed to account for ties in rankings, something
! --- this other ranking program I have doesn't do (ugh).
! --- Not very optimized, but... meh...
! --- Written by F. Alvarez, v. 1: 5/13

integer, intent(in) :: n_ranks
integer, intent(in) :: ranks(n_ranks)
real, intent(in) :: vals(n_ranks)
real, intent(out) :: outranks(n_ranks)
integer :: idx,begin_idx,idx2,intcount,intcount2
real :: curr_val,counter,rank_sum

counter=0.
intcount = 0
curr_val = 0.
rank_sum=0.
begin_idx=1
outranks(:)=0.
do idx = 1,n_ranks
   if (vals(ranks(idx)) .eq. curr_val) then
      counter=counter+1.
      intcount = intcount+1
      rank_sum=rank_sum+real(idx)
   else if (vals(ranks(idx)) .ne. curr_val) then
      if (counter .gt.0) then
         intcount2=1
         idx2 = (idx-intcount)
         do while (intcount2 .le. intcount)
            outranks(ranks(idx2)) = (rank_sum/counter)
            idx2=idx2+1
            intcount2=intcount2+1
         end do
      else if (idx.eq.1) then
         counter=0.
         intcount = 0
         curr_val = vals(ranks(idx))
         rank_sum=0.
         begin_idx=idx
      else if (counter .eq. 0 .and. (idx.ne.n_ranks) .and. (idx.ne.1)) then
         outranks(ranks(idx-1)) = idx-1
      else if ((idx.eq.n_ranks) .and. (counter .eq. 0)) then
         outranks(ranks(idx-1)) = idx-1
         outranks(ranks(idx)) = idx
      else if ((idx.eq.n_ranks) .and. (counter .gt. 0)) then
         rank_sum = rank_sum+real(idx)
         intcount2=1
         idx2 = (idx-intcount)
         do while (intcount2 .le. intcount)
            outranks(ranks(idx2)) = (rank_sum/counter)
            idx2=idx2+1
            intcount2=intcount2+1
         end do
         outranks(ranks(idx)) = idx

      end if
      counter=0.
      intcount = 0
      curr_val = vals(ranks(idx))
      rank_sum=0.
      begin_idx=idx
   end if
end do

return

end subroutine ties

subroutine rmse(train,fcst,shape,outData)
! --- Subroutine to calculate root mean squared error
integer,intent(in) :: shape
real,intent(in) :: train(shape),fcst(shape)
real,intent(out) :: outData

outdata = sqrt(sum((train(:)-fcst(:))**2,dim=1)/shape)

return

end subroutine rmse

subroutine mae(train,fcst,shape,outData)
! --- Subroutine to calculate mean absolute error
integer,intent(in) :: shape
real,intent(in) :: train(shape),fcst(shape)
real,intent(out) :: outData

outdata = sum(abs(train(:)-fcst(:)),dim=1)/shape

return

end subroutine mae

Subroutine correlation(trainData,fcstData,shape,corr,sig)

! --- Subroutine to find the correlation (and significance) of two arrays
! --- In variables:
! --- Shape - integer, shape of 1D arrays
! --- fcstData(shape) - real, array with data to correlate against
! --- trainData(shape) - real, array with data to find correlation with
!
! --- Out variables:
! --- Corr - real, value of correlation between two arrays

integer, intent(in) :: shape
real,dimension(shape),intent(in) :: trainData,fcstData
real,intent(out) :: corr

! --- Other variables...
real, dimension(shape):: xmn, ymn ! --- mean of data
real, dimension(shape):: xdev, ydev ! --- demeaned data
real, dimension(shape) :: xdevydev, xdevxdev, ydevydev !--- multiplication of demeaned data
real :: Sx,Sy,CovXY ! --- Other assorted variables we'll need

! --- Find mean of each array
xmn = SUM( fcstData, DIM = 1 ) / REAL(shape)
ymn = SUM( trainData, DIM = 1 ) / REAL(shape)

! --- Remove mean from data
xdev(:) = fcstData(:) - xmn
ydev(:) = trainData(:) - ymn

! --- Multiply de-meaned data together
xdevydev(:) = xdev(:) * ydev(:)

! --- Find covariance
COVxy = SUM( xdevydev, DIM = 1 ) / REAL(shape)

! --- Standard deviation stuff
xdevxdev(:) = xdev(:) * xdev(:)
Sx = SQRT( SUM( xdevxdev, DIM = 1 ) / REAL(shape) )

ydevydev(:) = ydev(:) * ydev(:)
Sy = SQRT( SUM( ydevydev, DIM = 1 ) / REAL(shape) )

! --- Correlation time!
corr = COVxy / ( Sx * Sy )


! --- Not sure if we really need significance, but whatever
if ( corr .ne. 1.0 ) then
   sig = corr * SQRT( (REAL(shape) - 2) / ( 1 - corr*corr) )
else
   sig = 1.0
end if

return

end Subroutine correlation

Subroutine rank_corr(trainData,fcstData,shape,corr,sig)

! --- Subroutine to find the correlation (and significance) of two arrays
! --- In variables:
! --- Shape - integer, shape of 1D arrays
! --- fcstData(shape) - real, array with data to correlate against
! --- trainData(shape) - real, array with data to find correlation with
!
! --- Out variables:
! --- Corr - real, value of correlation between two arrays

integer, intent(in) :: shape
real,dimension(shape),intent(in) :: trainData,fcstData
real,intent(out) :: corr

! --- Other variables...
real, dimension(shape):: xmn, ymn ! --- mean of data
real, dimension(shape):: xdev, ydev ! --- demeaned data
real, dimension(shape) :: xdevydev, xdevxdev, ydevydev !--- multiplication of demeaned data
real :: Sx,Sy,CovXY ! --- Other assorted variables we'll need

! --- Find mean of each array
xmn = SUM( fcstData, DIM = 1 ) / REAL(shape)
ymn = SUM( trainData, DIM = 1 ) / REAL(shape)

! --- Remove mean from data
xdev(:) = fcstData(:) - xmn
ydev(:) = trainData(:) - ymn

! --- Multiply de-meaned data together
xdevydev(:) = xdev(:) * ydev(:)

! --- Find covariance
COVxy = SUM( xdevydev, DIM = 1 ) / REAL(shape)

! --- Standard deviation stuff
xdevxdev(:) = xdev(:) * xdev(:)
Sx = SQRT( SUM( xdevxdev, DIM = 1 ) / REAL(shape) )

ydevydev(:) = ydev(:) * ydev(:)
Sy = SQRT( SUM( ydevydev, DIM = 1 ) / REAL(shape) )

! --- Correlation time!
corr = COVxy / ( Sx * Sy )


! --- Not sure if we really need significance, but whatever
if ( corr .ne. 1.0 ) then
   sig = corr * SQRT( (REAL(shape) - 2) / ( 1 - corr*corr) ) 
else
   sig = 1.0
end if

return

end Subroutine rank_corr
